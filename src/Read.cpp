//
//  Read.cpp

//  Created by Lucile Broseus in April 2018. The sky was cloudy.
//  Copyright © 2018&2019 LuB. All rights reserved.
//

/****************************************************************************
 * Version TALC_1.01
 ****************************************************************************
 ****************************************************************************
 * Updates:
 * 27/01/2019: removed messages used for test
 * 20/08/2019: tried to solve a bug in findInRegion().
 ****************************************************************************/

#include "Read.hpp"

#include <iostream>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <string>
#include <cmath>
#include <time.h>
#include <map>
#include <set>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "Settings.hpp"
#include "Explorer.hpp"

#include "Jellyfish.hpp"
#include "utils.hpp"
#include "io.hpp"

using namespace seqan;

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;

//---------------------------------------------------------------------//

extern unsigned int gp_MIN_COUNT;
extern double gp_SR_ERROR_RATE;
extern unsigned int K;
extern double gp_ALPHA;

extern std::string queryMode;
extern colouredDBG SR_DBG;

extern std::string io_pathToJF;
extern std::string io_JFcountsTable;
extern std::string io_JFjunctcountsTable;

extern std::string io_logFile;

extern bool DEBUG_USER;
extern bool DEBUG_TEST;

static bool DEBUG_READ(true);

Read::~Read()
{
}

//Constructors

Read::Read() : m_id(make_empty_CharString()), m_sequence(make_empty_Dna5String()),
											  m_scorr(make_empty_CharString()),
                                              m_coverage(make_empty_vecOfColouredPairs()),
                                              m_priorLambda_noise(gp_MIN_COUNT), 
                                              m_LRcoverage(make_empty_vector()),
                                              m_InKmersPositions(make_empty_vectOfKmPos()),
                                              m_newInnerStructure(make_empty_inner()),
                                              m_head(make_empty_border()),
                                              m_tail(make_empty_border()),
                                              m_nbInKmers(0), m_nbSolidKmers(0)
{

}

Read::Read(CharString& id) : m_id(id), m_sequence(make_empty_Dna5String()),
                                       m_scorr(make_empty_CharString()),
                                       m_correction(make_empty_Dna5String()),
                                       m_coverage(make_empty_vecOfColouredPairs()),
                                       m_priorLambda_noise(gp_MIN_COUNT), 
                                       m_LRcoverage(make_empty_vector()),
                                       m_InKmersPositions(make_empty_vectOfKmPos()),
                                       m_newInnerStructure(make_empty_inner()),
                                       m_head(make_empty_border()),
                                       m_tail(make_empty_border()),
                                       m_nbInKmers(0), m_nbSolidKmers(0)
{

}

Read::Read(CharString& id, TSeq& sequence) : m_id(id), 
                                             m_sequence(sequence), 
                                             m_scorr(make_empty_CharString()),
                                             m_correction(make_empty_Dna5String()),
                                             m_coverage(make_empty_vecOfColouredPairs()), 
                                             m_priorLambda_noise(gp_MIN_COUNT),
                                             m_LRcoverage(make_empty_vector()),
                                             m_InKmersPositions(make_empty_vectOfKmPos()),
                                             m_newInnerStructure(make_empty_inner()),
                                             m_head(make_empty_border()),
                                             m_tail(make_empty_border()),                                             
                                             m_nbInKmers(0), m_nbSolidKmers(0)
{
	
}

//Methods

seqan::CharString Read::getName()
{
    return m_id;
}

TSeq Read::getSeq()
{
    return m_sequence;
}

TSeq Read::getCorrSeq()
{
    return m_correction;
}

void Read::setCorrSeq(TSeq& corrSeq)
{
    m_correction = corrSeq;
}

CharString Read::getScorr()
{
    return m_scorr;
}

void Read::setScorr(CharString& scorr)
{
    m_scorr = scorr;
}

void Read::updateScorr(CharString scorrNewPart)
{
	append(m_scorr, scorrNewPart);
}


int Read::getLength()
{
	return length(m_sequence);
}

int Read::getCorrLength()
{
	return length(m_correction);
}


//--------------------------------------------------------------------->DEBUG
// FUNCTION reCoverage()
// query database for the count of each k-mer in the raw sequence
//--------------------------------------------------------------------->DEBUG

bool Read::reCoverage()
{
//--------------------------------------------------------------------->DEBUG
//std::cout << "FUNCTION RECOVERAGE" << std::endl;
//--------------------------------------------------------------------->DEBUG
  m_nbInKmers = 0;
    	
  m_coverage = getLRCountsInSR(m_sequence, K, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable); 
  
  if( DEBUG_READ )
   {
	std::cout << "[DEBUG_READ]: raw coverage" <<  std::endl;  
    for(unsigned int i(0); i<m_coverage.size();i++) std::cout << std::get<0>(m_coverage[i]) << " ";
    std::cout << std::endl;
   }

  if(!empty(m_coverage)) for(unsigned int p(0); p<m_coverage.size(); p++) if(std::get<0>(m_coverage[p])>gp_MIN_COUNT) m_nbInKmers++; 

  //if( DEBUG_USER ) std::cout << "[DEBUG READ]: nb solid kmers: " << m_nbInKmers << std::endl;  

    return (m_nbInKmers>0);     
}

std::vector<colouredCount > Read::getCoverage()
{
	return m_coverage;
}

void Read::setPriorNoise(double lambda)
{
	m_priorLambda_noise = lambda;
}


void Read::displaySegmentation()
{
	std::cout << "[DEBUG_READ] SEGMENTATION: " << m_InKmersPositions.size() << " IN-regions" << std::endl;
    for(int re(0); re<m_InKmersPositions.size(); re++) std::cout << "IN-Region #" << re << " " << std::get<0>(m_InKmersPositions[re]) << "-" << std::get<1>(m_InKmersPositions[re]) << std::endl; 
}

bool Read::setInitialStructure()
{
 TSeq temp;
 unsigned int len(0);
 m_head = std::make_pair(temp, ABSENT);
 m_tail = std::make_pair(temp, ABSENT);
 
 if( !empty(m_InKmersPositions) )
  {
   if( std::get<0>(m_InKmersPositions[0])>0 )
    {
	 //std::cout << "WEAK HEAD DETECTED" << std::endl;
	 temp = extractWeakBorderSequence(m_sequence, std::get<0>(m_InKmersPositions[0]), K, HEAD);
	 m_head = std::make_pair(temp, UNCORRECTED); 
	 len += length(temp);
    }
   }
   if( std::get<1>(m_InKmersPositions.back())+1<m_coverage.size() )
    {
	 //std::cout << "WEAK TAIL DETECTED: " << std::get<1>(m_InKmersPositions.back()) << std::endl;
	 temp = extractWeakBorderSequence(m_sequence, std::get<1>(m_InKmersPositions.back()), K, TAIL);
	 m_tail = std::make_pair(temp, UNCORRECTED); 
	 len += length(temp);
    }
  for(unsigned int i(0); i+1<m_InKmersPositions.size(); i++)
   {
    //std::cout << "INNER #" << i << std::endl;
	temp = extractSolidSequence(m_sequence, std::get<0>(m_InKmersPositions[i]), std::get<1>(m_InKmersPositions[i]),K) ;
	m_newInnerStructure.push_back(std::make_pair(temp, std::get<2>(m_InKmersPositions[i]))); 
	len += length(temp);
	if(std::get<0>(m_InKmersPositions[i+1])>std::get<1>(m_InKmersPositions[i])+K) temp = extractWeakSequence(m_sequence, std::get<1>(m_InKmersPositions[i]), std::get<0>(m_InKmersPositions[i+1]), K) ;
	else clear(temp);
	
	m_newInnerStructure.push_back(std::make_pair(temp, UNCORRECTED)); 
	len += length(temp);

   }
  temp = extractSolidSequence(m_sequence, std::get<0>(m_InKmersPositions.back()), std::get<1>(m_InKmersPositions.back()), K) ;
  m_newInnerStructure.push_back(std::make_pair(temp, std::get<2>(m_InKmersPositions.back()))); 
  len += length(temp);

if( len != length(m_sequence) ) std::cout << "[WARNING] Read: PB in Read: initStructure: " << len << " vs " << length(m_sequence) << std::endl;

  return (len == length(m_sequence));
}

bool Read::defineStructure2()
{
  bool checok(true);
  double thr(gp_MIN_COUNT);

 checok = findINRegions(m_InKmersPositions, m_coverage, K);
 if( DEBUG_READ ) std::cout << "[DEBUG_READ]: NB IN-REG FOUND: " << m_InKmersPositions.size() << std::endl;
 
 thr = computeSeqErrorThreshold(m_coverage);
 this->setPriorNoise(thr);

 analyzeINRegions(m_InKmersPositions, m_sequence, m_coverage, thr, K);
 
 if( !empty(m_InKmersPositions) ) checok &= this->setInitialStructure();
 	
 return checok;
}

void Read::display()
{		
	std::cout << "[DEBUG_MODE]: Read " << m_id << std::endl;
	std::cout << "[DEBUG_MODE]: Length: " << this->getLength() << std::endl;
}

/***********************************************************************
 * Correct sequence according to settings  
 ***********************************************************************/
 
std::pair<TSeq, Status > Read::getSolidRegion(std::tuple<unsigned int, unsigned int, Status > coordinates)
{
  TSeq seq =  extractSolidSequence(m_sequence, std::get<0>(coordinates), std::get<1>(coordinates), K);
  return std::make_pair(seq, std::get<2>(coordinates));
}

void Read::updateINNER(Explorer& myExplorer, int& reg)
{
    m_InKmersPositions[reg] = myExplorer.getLEFTHandPositions();
	m_InKmersPositions[reg+1] = myExplorer.getRIGHTHandPositions();
	
    m_newInnerStructure[2*reg+1] = myExplorer.getWeakSeq();

    m_newInnerStructure[2*reg] = getSolidRegion(myExplorer.getLEFTHandPositions());
    m_newInnerStructure[2*(reg+1)] = getSolidRegion(myExplorer.getRIGHTHandPositions()); 
}

void Read::updateHEAD(Explorer& myExplorer)
{
//if( DEBUG_USER ) std::cout << "UPDATE HEAD: " << std::get<0>(myExplorer.getRIGHTHandPositions()) << "-" << std::get<1>(myExplorer.getRIGHTHandPositions()) << std::endl;	
  m_InKmersPositions[0] = myExplorer.getRIGHTHandPositions();
  m_newInnerStructure[0] = getSolidRegion(myExplorer.getRIGHTHandPositions());
  m_head = myExplorer.getWeakSeq();
}

void Read::updateTAIL(Explorer& myExplorer)
{
  m_InKmersPositions.back() = myExplorer.getLEFTHandPositions();
  m_newInnerStructure.back() = getSolidRegion(myExplorer.getLEFTHandPositions());
  m_tail = myExplorer.getWeakSeq();
}

void Read::updateCorrSeq()
{
	clear(m_correction);
	append(m_correction, std::get<0>(m_head));
	for(int reg(0); reg<m_newInnerStructure.size();reg++) append(m_correction, std::get<0>(m_newInnerStructure[reg]));
	append(m_correction, std::get<0>(m_tail));
} 

TSeq Read::concatenateNewStructure()
{
 TSeq newSequence;
 for(unsigned int i(0); i<m_newInnerStructure.size(); ++i) append(newSequence, std::get<0>(m_newInnerStructure[i]));
	
 return newSequence;
}

void Read::correct2()
{
  Explorer myExplorer(this->getSeq(), this->getCoverage(), m_priorLambda_noise);
  bool success;
  
  if( m_InKmersPositions.size()>0 )
  {	
	try
	 {		 
	   //create explorer to correct the current region
	  for(int reg(0); reg<m_InKmersPositions.size()-1; reg++)
       {
//std::cout << "INNER REGION #" << reg << std::endl; 		  
		 myExplorer.initializeINNER(m_InKmersPositions[reg], m_InKmersPositions[reg+1], RIGHT);	 
	     success = myExplorer.searchBridge();
	     if(!success) 
	      {
//std::cout << "TRYING LEFT FOR REGION #" << reg << std::endl; 		  			  
		  myExplorer.initializeINNER(m_InKmersPositions[reg], m_InKmersPositions[reg+1], LEFT);	 
	      myExplorer.searchBridge();
//std::cout << "INNER REGION #" << reg << " correction achieved." << std::endl; 		  	       
	      }
         this->updateINNER(myExplorer, reg);
//std::cout << "END INNER REGION #" << reg << std::endl;	             
	    }
	  if( (std::get<1>(m_head)!=ABSENT) & (length(std::get<0>(m_head))<=500) ) 
	    {
//std::cout << "WEAK HEAD" << std::endl; 		   
		 myExplorer.initializeHEAD( m_InKmersPositions[0] ); 	 
		 success = myExplorer.searchEdge();
	     if(success ) this->updateHEAD(myExplorer);    
	    }
	   if( (std::get<1>(m_tail)!=ABSENT) & (length(std::get<0>(m_tail))<=500) )
	    {
//std::cout << "WEAK TAIL" << std::endl; 		   
		 myExplorer.initializeTAIL( m_InKmersPositions.back() );
	     success = myExplorer.searchEdge();
	     if( success ) this->updateTAIL(myExplorer);
	    }
      }//END-TRY
     catch(std::string const& chaine)
      {
	   this->display(); 	 
	   std::cout << chaine << std::endl;
	   throwToLog(this->getName(), chaine, io_logFile);
	  } 
    }
  else std::cout << "Read " << m_id << " has no defined solid region. No correction could be performed." << std::endl;

  this->updateCorrSeq();
}

/***********************************************************************
 * Compute statistics  
 ***********************************************************************/

//read_id raw_length head_length tail_length nbInKmers nbSolidKmers nbSolidReg corr_length (nbInKmers2?)

void setBasicReadStatsHeader(std::string& outFileName)
{
   std::ofstream outputFile;
    
   outputFile.open(outFileName, std::ios_base::trunc);
     outputFile << "read_name" << "\t"
                << "raw_length" << "\t"
                << "whead_length" << "\t"
                << "wtail_length" << "\t"
                << "nbInKmers" << "\t"
                << "nbSolidKmers" << "\t"
                << "nbSolidReg" << "\t"
                << "nbInWeakReg" << "\t"
                << "nbInCorrReg" << "\t"
                << "CorrHead?" << "\t"
                << "CorrHeadLen" << "\t"
                << "CorrTail?" << "\t"
                << "CorrTailLen" << "\t"
                << "Corrlength" << "\t"	
                << "nbInKmers2" << "\n";
   outputFile.close();            
}


void Read::outputBasicReadStats(std::string& outFileName)
{
   int nbSReg(0);
   unsigned int nbInKmersBefore(0);
 
   if( !empty(m_InKmersPositions) ) for(unsigned int i(0); i<m_InKmersPositions.size(); i++)  nbInKmersBefore += std::get<1>(m_InKmersPositions[i])-std::get<0>(m_InKmersPositions[i])+1;
   nbSReg = m_InKmersPositions.size();
   std::ofstream outputFile;
   outputFile.open(outFileName, std::ios_base::app);      
   outputFile << "\n" << m_id << "\t"                                 //read name
               << length(m_sequence) << "\t"                          //raw length
               << nbInKmersBefore << "\t"                                   //span of prior solid regions
               << nbSReg << "\t"                                      //initial " of solid regions
               << length(m_correction);                               //length after correction
    outputFile.close();
}

//Function findINRegions()
//Screens the kmer count vector and records start+end kmer position of regions of the LR that are IN the SR-dBG
//Fills StartEndKmers with the position of the 1st/last kmer of each IN region
//Returns whether there is at least one IN region

bool findINRegions(std::vector<std::tuple<unsigned int, unsigned int, Status > >& StartEndKmers, 
                   std::vector<colouredCount >& counts,
                   unsigned int& kmerSize)
{
  unsigned int pos(0);
  unsigned int current_startKmer(0);
  unsigned int current_count(0);
  
  bool state(false);                 //records whether the kmer #pos is IN or OUTSIDE the SR-dBG
  
  if( DEBUG_READ ) std::cout << "[DEBUG_READ] Searching for Solid regions...width:" << counts.size() <<  std::endl;
  
  if( counts.size() > 1 )
  {
	   try
	    { 
           while( pos<counts.size() )
			{
				current_count = std::get<0>( counts[pos] );
				//If we enter the SR-dBG    
				if( (current_count>=gp_MIN_COUNT) & (state==false) )
				{
				current_startKmer = pos;
				state = true;
				}
				//else if we were in the SR-dBG and get out
				else if( (current_count<gp_MIN_COUNT) & (state==true) )
				{
				StartEndKmers.push_back( std::make_tuple(current_startKmer, pos-1, UNEXPECTED) );
				state = false;
				}
				++pos;
			}//END-WHILE
			
			//Adressing the last kmer if still in the SR-dBG
			if( state==true ) StartEndKmers.push_back( std::make_tuple(current_startKmer, pos-1, UNEXPECTED) );
		}
		catch(std::string const& chaine)
        { 
		  std::cout << "[ERROR] in findINRegions" << std::endl;
		  std::cout << chaine << std::endl;
		}	
  }else{
    std::cout << "[WARNING] in findInRegions: count vector seems unexpectedly empty..." << std::endl;
  }
  
  //if( DEBUG_USER ) std::cout << "[DEBUG USER] " << StartEndKmers.size() << " IN-Regions found." << std::endl;
   
   return !empty(StartEndKmers);	
}

//Extract counts of in kmers
//Compute a robust mean of these counts and a thresold for SR sequencing error counts
double computeSeqErrorThreshold( std::vector<colouredCount >& counts)
{
	std::vector<unsigned int > INcounts;
	double robMean( gp_MIN_COUNT );
	
	unsigned int first;
	unsigned int last;
	
	for(unsigned int pos(0); pos<counts.size(); pos++) 
	 {
	  if( std::get<0>(counts[pos])>=gp_MIN_COUNT ) INcounts.push_back( std::get<0>(counts[pos]) );
	 }
	sort(INcounts.begin(), INcounts.end());
	
	INcounts.size()>10 ? first = (unsigned int)(0.15*(double)INcounts.size()) : first = 0; 
    INcounts.size()>10 ? last = (unsigned int)(0.90*(double)INcounts.size()) : last = INcounts.size(); 
	
	//Take the mean of the counts removing the first 10% and last 5%
	for(unsigned int i(first); i<last; i++) robMean += INcounts[i];
	robMean /= last-first;  
	
if( DEBUG_USER || DEBUG_TEST ) std::cout << "ROB MEAN " << robMean << " THR: " << robMean*gp_SR_ERROR_RATE << std::endl;
        
	return robMean*gp_SR_ERROR_RATE;
	
}

//Checks whether each IN-Region satisfies the following criteria:
// can be continued both on the right and left side (outdegree of edge kmer >=1, 2 if not primary first kmer)
// space between any two consecutive IN-Regions must be >=0
//replaces StartEndKmers with the new positions
void analyzeINRegions(std::vector<std::tuple<unsigned int, unsigned int, Status > >& StartEndKmers,
                      TSeq& refSequence,
                      std::vector<colouredCount >& counts, 
                      double& solidityThr, unsigned int& kmerSize)
{	
 std::vector<std::tuple<unsigned int, unsigned int, Status > > newStartEndKmers;
 TSeq kmer;
 bool OK(true);
 unsigned int new_start_pos(0);
 unsigned int new_end_pos(0);
 Status regionStatus(UNEXPECTED);
 unsigned int c(0);
 int span(0);
 
   for( unsigned int reg(0); reg<StartEndKmers.size(); reg++)
    {
	 span = 0;
	 //Adressing IN-Region #reg
	 //Are there SR-kmers continuing the region both on the left and right hand sides?
	 OK = true;
	 //LEFT SIDE
	 new_start_pos = std::get<0>(StartEndKmers[reg]);
	 kmer = getKmerAt(refSequence, new_start_pos, K);
	 //If not, reduce the region by one kmer and check again
	 if( !(StartEndKmers.size()==reg+1 & empty(newStartEndKmers)) & (getOutDegree(kmer, LEFT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable)==0) & (new_start_pos!=0) )
	  {
	   OK = false;
	   while( (new_start_pos<std::get<1>(StartEndKmers[reg])) & !OK )
	    {
		  ++new_start_pos;
		  kmer = getKmerAt(refSequence, new_start_pos, K);
		  if( getOutDegree(kmer, LEFT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable)>1 ) OK = true;
		}//END-WHILE
	  } 
	 new_end_pos = std::get<1>(StartEndKmers[reg]); 
	 if( OK & !(StartEndKmers.size()==reg+1 & empty(newStartEndKmers)) )
	  {
	   //RIGHT SIDE
	   kmer =  getKmerAt(refSequence, new_end_pos, K);
	   if( (getOutDegree(kmer, RIGHT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable)==0) & (new_end_pos!=counts.size()-1) )
	    {
		 OK = false;
	     while( (new_end_pos>std::get<0>(StartEndKmers[reg])) & !OK )
	     {
		  --new_end_pos;
		  kmer = getKmerAt(refSequence, new_end_pos, K);
		  if( getOutDegree(kmer, RIGHT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable)>1 ) OK = true;
		 }//END-WHILE
	    }	
	   }
	 //if OK region can be extended further on both side in the SR-dBG  
	 //Check whether the distance with the next in region in long enough
	 if( OK ) 
	  {
	    if( (reg+1<StartEndKmers.size()) ) span =  ( (int)std::get<0>(StartEndKmers[reg+1])-(int)(new_end_pos + kmerSize) );
	    //If overlap with following region, skip the current region and extend the following one
	    if( span<0 ) 
	     { 
		  if( (int)std::get<1>(StartEndKmers[reg+1])+span>=(int)std::get<0>(StartEndKmers[reg+1]) ) std::get<0>(StartEndKmers[reg+1]) -=  span; 
		  else //merge current region with next region, skip current region
		   {
			std::get<0>(StartEndKmers[reg+1]) = new_start_pos;
			OK = false;
		   } 
	     }
	    if( OK )
	     {
		  c = 0;
		  for( unsigned int i(new_start_pos); i<=new_end_pos;i++ ) c < std::get<0>(counts[i]) ? c = std::get<0>(counts[i]) : c += 0;
		  !isExpectedbyMyModel((unsigned int)c, (unsigned int)solidityThr, gp_ALPHA, UNEXPECTED) ? regionStatus = EXPECTED : regionStatus = LOWCOUNT;
		  if(regionStatus == EXPECTED) newStartEndKmers.push_back( std::make_tuple(new_start_pos, new_end_pos, regionStatus) ); 
		 }
	    }//END-IF OK
    }//END-FOR each IN-Region
    
  if( !empty(newStartEndKmers) ) StartEndKmers = newStartEndKmers;	 
}

Read create_empty_Read()
{
	Read r;
	return r;
}
