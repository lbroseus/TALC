//
//  Jellyfish.cpp
//  LReader
//
//  Created by Lucile Broseus on 16/10/2017.
//  Copyright © 2017 LuB. All rights reserved.
//

/****************************************************************************
 * Version TALC_1.01
 ****************************************************************************
 ****************************************************************************
 * Updates:
 * 27/01/2019: canonical mode removed: kmers must be directional
 * 27/01/2019: coloured Kmers are canonical (read and stored in both directions)
 * 27/01/2019: add a threshold to colour a kmer
 ****************************************************************************/
 
#include "Jellyfish.hpp"

#include <stdio.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "utils.hpp"
#include "Settings.hpp"

using namespace std;
using namespace seqan;

typedef Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;
//---------------------------------------------------------------------//

extern unsigned int gp_MIN_COUNT;


extern std::string queryMode;
extern std::string io_pathToJF;
extern bool gp_useJunctions;

extern std::string logFile;

extern bool DEBUG;


static std::vector<char > Dict = {'A', 'C', 'G', 'T'};
static unsigned int colouredCountThr = 10000;
//---------------------------------------------------------------------//

//----------------------------------------------------------------------/

 std::vector<TSeq> getKmers(TSeq& Seq, unsigned int& kmerSize)
{
    unsigned readLength = length(Seq);
    std::vector<TSeq > kmers;
    TSeq tmp;

    for(unsigned startPosition(0); startPosition < readLength-kmerSize+1; startPosition++)
    {	
	 tmp = infix(Seq, startPosition, startPosition + kmerSize);
     kmers.push_back( tmp );
    }
if( length(Seq)-kmerSize+1 != kmers.size() ) std::cout << "PB when chunking kmers " << std::endl; 
    return kmers;
}

//----------------------------------------------------------------------/



/*************************
 * appendSuccessors()
 * output all the four right neighbours of a given kmer
 *************************/
 
std::string appendSuccessors(TSeq& kmer)
{
    std::string result;

    TSeq rest(kmer);
    std::string tmp;

    erase(rest,0);
    for(unsigned int i(1); i <length(kmer); ++i)
    {
        tmp+= *toCString(rest);
        erase(rest,0);
    }
    
    if(length(tmp)+1 != length(kmer)) std::cout << "Pb with tmp: " << tmp;
    
    for(unsigned int b(0); b<Dict.size(); b++)
    {
        result += tmp + Dict[b] + " ";
    }
    return result;
}

std::vector<TSeq> getSuccessors(TSeq& kmer, Direction direction)
{
    std::vector<TSeq> result;

    for(unsigned int b(0); b<Dict.size(); b++)
    {
        result.push_back( formNextKmer(kmer, Dict[b], direction) );
    }

    return result;
}

/*************************
 * appendPredecessors()
 * output all the four left neighbours of a given kmer
 *************************/
 
std::string appendPredecessors(TSeq& kmer)
{
    std::string result;
    
    TSeq rest(kmer);
	std::string tmp;
	
    for(unsigned int i(0); i <length(kmer)-1; ++i)
    {
        tmp += *toCString(rest);
        erase(rest,0);
    }  
    for(unsigned int b(0); b<Dict.size(); b++)
    {
        result += Dict[b] + tmp + " ";
    }

    return result;
}

std::string appendKmers(std::vector<TSeq>& kmers)
{
    std::string result;
    unsigned int size(0);

    for(unsigned int i(0); i<kmers.size(); ++i)
    {
        TSeq tmp(kmers[i]);
        size = length(tmp);

        for(unsigned int j(0); j<size;++j)
        {
          result += *toCString(tmp);
          erase(tmp,0);
        }
        result += " ";
    }
    return result;
}

//----------------------------------------------------------------------/

std::map<TSeq, colouredCount > buildDBGfromDump(std::string& countsTable, std::string& junctionCountsTable)
{
    std::map<TSeq, colouredCount > dBG;
    
	std::fstream in;
	in.open(countsTable);
	
    int onlineCounter(0);
	int addedjcounter(0); 
	
	std::string line;
	while(getline(in, line))
	{
	 std::istringstream iss(line); // Create an input stream to read values
	 std::string kmer0;
	 std::string count;
	    
	    if(iss >> kmer0 >> count) //Read the input delimited by whitespace chars
	     {
		  TSeq kmer = kmer0;
		  if(std::stoi(count)>=gp_MIN_COUNT) dBG.insert( std::pair<TSeq, colouredCount>(kmer,std::make_pair(std::stoi(count),0)) );
		  onlineCounter++;
		 }
		else std::cout << "Issue with count table when building dBG..." << std::endl;
	 if((onlineCounter>=15000000) & (onlineCounter%15000000==0)) std::cout << onlineCounter << " keys inserted." << std::endl;	
    }   
    in.close();
  
  if( gp_useJunctions ){ 
	std::cout << "Junction mode activated." << std::endl;  
	
    in.open(junctionCountsTable);
    
	while(getline(in, line))
	{
		std::istringstream iss2(line); // Create an input stream to read your values
		std::string jmer0;
	    std::string jcount;
	    
	    if(iss2 >> jmer0 >> jcount) //Read the input delimted by whitespace chars
	    {
		  if( dBG.count(jmer0)>0 ) 
		  {
			  TSeq jmer = jmer0;
			  std::get<1>(dBG[jmer]) = std::stoi(jcount);
		  }
		  //else 
		  //{
		  //  TSeq jmer = jmer0;
		  //  dBG.insert( std::pair<TSeq, colouredCount>(jmer,std::make_pair(0, std::stoi(jcount))) );
		  //  addedjcounter++;
		  //}
		}
		else std::cout << "Error when building dBG..." << std::endl;	
    }
    std::cout << "Nbr of junction kmers absent from SR data: " << addedjcounter << std::endl;
   }//NED-IF useJunctions
   
	return dBG;
}

std::map<TSeq, colouredCount > buildCDBG(int mode, std::string& countsTable, std::string& junctionCountsTable)
{
    std::map<TSeq, colouredCount > dBG;
    std::string jmer0;
	std::string jcount;
	
    std::fstream in;
    
    int onlineCounter(0);
	int actualCounter(0);
	 
	std::string line;
	
	in.open(countsTable);

	while(getline(in, line))
	{
		std::istringstream iss(line); // Create an input stream to read values
		std::string kmer0;
	    std::string count;
	    
	   if(iss >> kmer0 >> count) //Read the input delimited by whitespace chars
	    {
		  TSeq kmer = kmer0;
		  if(std::stoi(count)>=gp_MIN_COUNT) 
		   {
			 dBG.insert( std::pair<TSeq, colouredCount>(kmer,std::make_pair(std::stoi(count),0)) );
		     actualCounter++;
		   }   
		   onlineCounter++;
		 if( (onlineCounter>=15000000) & (onlineCounter%15000000==0)) std::cout << onlineCounter << " keys evaluated." << std::endl;	
		}
	  else std::cout << "Error when building dBG..." << std::endl;
    }
    
    in.close();
  
  if( gp_useJunctions )
   { 
	std::cout << "Junction mode activated." << std::endl;  
    in.open(junctionCountsTable);
    
	while(getline(in, line))
	 {
		std::istringstream iss2(line); // Create an input stream to read your values
	    if(iss2 >> jmer0 >> jcount) //Read the input delimted by whitespace chars
	     {
		  TSeq jmer(jmer0);
		  if( (std::stoi(jcount)<colouredCountThr) & (dBG.count(jmer)>0)  ) std::get<1>(dBG[jmer]) = std::stoi(jcount);
		  reverseComplement(jmer);
		  if( (std::stoi(jcount)<colouredCountThr) & (dBG.count(jmer)>0) ) std::get<1>(dBG[jmer]) = std::stoi(jcount);
		 }
	    else std::cout << "Error when building dBG..." << std::endl;	
	 }
   }//END-IF useJunctions
std::cout << "There were " << onlineCounter << " k-mers retrieved from database." << std::endl;
std::cout << "In the whole, we have kept " << actualCounter << "k-mers, whose counts were over the specified threshold." << std::endl;

	return dBG;
}

//----------------------------------------------------------------------/

std::vector<colouredCount > getNextCounts(TSeq& kmer, Direction direction, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable)
{
    std::vector<colouredCount > colouredCounts;
    if( queryMode=="jellyfish2-only" || queryMode=="jellyfish2-colour" ) colouredCounts = getNextCountsFromJF(kmer, direction, countsTable, junctionCountsTable);
	else colouredCounts = getNextCountsFromDBG(kmer, direction, dBG);
		
    return colouredCounts;
}

std::vector<colouredCount > getNextCountsFromDBG(TSeq& kmer, Direction direction, colouredDBG& dBG)
{
   std::vector<colouredCount > colouredCounts;
   std::vector<TSeq> nextKmers;
    
   nextKmers = getSuccessors(kmer, direction); 
  
    for(unsigned int k(0); k < nextKmers.size(); k++)
    {		
	 if(dBG.count(nextKmers[k])>0) colouredCounts.push_back(dBG.at(nextKmers[k]));
	 else colouredCounts.push_back(std::make_pair(0,0));
	}
    return colouredCounts;
}

std::vector<colouredCount > getNextCountsFromJF(TSeq& kmer, Direction direction, std::string& countsTable, std::string& junctionCountsTable)
{
	std::vector<colouredCount > colouredCounts;
	std::string nextKmers;

     if(direction == LEFT) nextKmers = appendPredecessors(kmer);
     else                  nextKmers = appendSuccessors(kmer); 
     
    std::vector<unsigned int> counts;
    std::vector<unsigned int> jcounts;
    FILE * stream;
    const int max_buffer = 1024;
    char buffer[max_buffer];
    std::string cmd;

    cmd = io_pathToJF + "/" + "jellyfish" + " " + "query" + " " + countsTable + " " + nextKmers;
    cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
                std::string stri(buffer);
                counts.push_back(atoi(stri.substr(length(kmer), stri.length()).c_str()));        
           }
        pclose(stream);
    }
    else std::cout << "DEBUG: Pb with stream " << std::endl;
    
    if( gp_useJunctions )
    {
	 cmd = io_pathToJF + "/" + "jellyfish" + " " + "query" + " " + junctionCountsTable + " " + nextKmers;
     cmd.append(" 2>&1");
     
     stream = popen(cmd.c_str(), "r");
     if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
                std::string stri(buffer);
                jcounts.push_back(atoi(stri.substr(length(kmer), stri.length()).c_str()));        
            }
        pclose(stream);
        //buffer.clear();
     }	
	}
    
    for(unsigned int b(0); b<counts.size(); b++)
    {
	  if( gp_useJunctions ) colouredCounts.push_back(std::make_pair(counts[b],jcounts[b]));
	  else colouredCounts.push_back(std::make_pair(counts[b],0));
	}
if( colouredCounts.size() != 4) std::cout << "PB in nextCounts " << colouredCounts.size() << std::endl; 

    return colouredCounts;
}

//----------------------------------------------------------------------/

int getOutDegree(TSeq& kmer, Direction direction, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable)
{
   std::vector<colouredCount > nextCounts = getNextCounts(kmer,direction, dBG, countsTable, junctionCountsTable);
   int outDegree(0);
   for(unsigned int i(0); i<nextCounts.size(); i++) 
    { 
	  if(std::get<0>(nextCounts[i])>=gp_MIN_COUNT) ++outDegree;
	}
if(outDegree>4) std::cout << "PB in outdegree: " << nextCounts.size() << std::endl;
   return outDegree;
}

//----------------------------------------------------------------------/

colouredCount getCount(TSeq& kmer, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable)
{
	colouredCount count;

	if( queryMode=="jellyfish2-only" || queryMode=="jellyfish2-colour" ) count = getCountFromJF(kmer, countsTable, junctionCountsTable);
	else count = getCountFromDBG(kmer, dBG);
	
	return count;
}

colouredCount getCountFromDBG(TSeq& kmer, colouredDBG& dBG)
{
	colouredCount counts;
    if((dBG.count(kmer))>0) counts = dBG.at(kmer);
    else counts = std::make_pair(0,0);
    return counts;
}

colouredCount getCountFromJF(TSeq& kmer, std::string& countsTable, std::string& junctionCountsTable)
{
    unsigned int count(0);
    unsigned int jfcolour(0);
    TSeq tmp = kmer;
    std::string s;

    for(unsigned int i(0); i < length(kmer); ++i)
    {
        s+= *toCString(tmp);
        erase(tmp,0);
    }

    FILE * stream;
    const int max_buffer = 1024;
    char buffer[max_buffer];

    std::string cmd;
    cmd = io_pathToJF + "/" + "jellyfish" + " " + "query" + " " + countsTable + " " + s;
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream){
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
                std::string stri(buffer);
                count = atoi(stri.substr(length(kmer), stri.length()).c_str());
            }
        pclose(stream);
        //buffer.clear();
    }
    
    if( gp_useJunctions )
    {
      cmd = io_pathToJF + "/" + "jellyfish" +" "+ "query" + " " + junctionCountsTable + " " + s;
      cmd.append(" 2>&1");
      
      stream = popen(cmd.c_str(), "r");
      if (stream){
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
                std::string stri(buffer);
                jfcolour = atoi(stri.substr(length(kmer), stri.length()).c_str());
            }
        pclose(stream);
        //buffer.clear();
      }	
	}//END-IF useJunctions

    return std::make_pair(count, jfcolour);
}

//----------------------------------------------------------------------/

std::vector<colouredCount> getLRCountsInSR(TSeq& Seq, unsigned int& kmerSize, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable)
{
//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION GETLRCOUNTS SR" << std::endl;
//--------------------------------------------------------------------->DEBUG		
	std::vector<colouredCount > counts;
	bool junctions = (queryMode=="jellyfish2-only");
	
	if( queryMode=="jellyfish2-only" || queryMode=="jellyfish2-colour" ) counts = getLRCountsInSRFromJF(Seq, kmerSize, countsTable, junctionCountsTable); 								
	else counts = getLRCountsInSRFromDBG(Seq, kmerSize, dBG);
	
	return counts;
}

std::vector<colouredCount> getLRCountsInSRFromDBG(TSeq& Seq, unsigned int& kmerSize, colouredDBG& dBG)
{
    std::vector<colouredCount > counts;
    std::vector<TSeq> kmers = getKmers(Seq, kmerSize);
       
    for(unsigned int position(0); position < kmers.size(); position++)
    {	 
      if( dBG.count(kmers[position])>0 ) counts.push_back(dBG.at(kmers[position]));
      else counts.push_back(std::make_pair(0,0));
	}
    return counts;
}

std::vector<colouredCount> getLRCountsInSRFromJF(TSeq& Seq,unsigned int& kmerSize, std::string& countsTable, std::string& junctionCountsTable)
{
//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION GETLRCOUNTFROMJF" << std::endl;
//--------------------------------------------------------------------->DEBUG	
	
	std::vector<colouredCount > counts;
	std::vector<TSeq> kmers = getKmers(Seq, kmerSize);

    std::vector<unsigned int > jfcounts;
    std::vector<unsigned int > jfcolours;
    
    FILE * stream;
    const int max_buffer = 1024;
    char buffer[max_buffer];
    std::string cmd;
    cmd = io_pathToJF + "/" + "jellyfish" + " " + "query" + " " + countsTable + " " + (appendKmers(kmers));
    cmd.append(" 2>&1");
        
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
             std::string stri(buffer);
              jfcounts.push_back(atoi(stri.substr(kmerSize, stri.length()).c_str()));
            }
        pclose(stream);      
    }
    
    if( gp_useJunctions )
    {
      cmd = io_pathToJF + "/" + "jellyfish" + " " + "query" + " " + junctionCountsTable + " " + (appendKmers(kmers));
      cmd.append(" 2>&1");

      stream = popen(cmd.c_str(), "r");
      if (stream){
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
            {
              std::string stri(buffer);
              jfcolours.push_back(atoi(stri.substr(kmerSize, stri.length()).c_str()));
            }
        pclose(stream);
      }	
	}//END-IF useJunctions
    
    for(unsigned int pos(0); pos<jfcounts.size(); pos++)
    {
	  if( gp_useJunctions ) counts.push_back( std::make_pair(jfcounts[pos], jfcolours[pos] ) );
	  else counts.push_back( std::make_pair(jfcounts[pos], 0) );
	}
    std::cout << "# of recorded counts: " << counts.size() <<  " " << jfcounts.size() << std::endl;
    return counts;
}




