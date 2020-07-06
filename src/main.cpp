//
//  main.cpp
//  TALC
//
//  Created by Lucile Broseus on the 06/10/2017.
//  Copyright © 2017&2018&2019 LuB. All rights reserved!!!!!
//

/****************************************************************************
 * Date: 01/01/2019
 * Version TALC 1.01
 ****************************************************************************
 * Object:
 * First version to be shared
 ****************************************************************************
 * Updates:
 * Added DEBUG_USER mode
 * Changed options for query mode: memory as default mode, coloured iff file specified
 ****************************************************************************/

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

#include <omp.h>

#include "utils.hpp"
#include "io.hpp"
#include "Jellyfish.hpp"  
#include "Settings.hpp"
#include "Read.hpp"

using namespace seqan;

typedef Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;

extern unsigned int gp_MIN_COUNT;
extern unsigned int K;
extern unsigned int gp_MAX_NB_COMPETING_PATHS;

extern unsigned int gp_minTailWidth;

extern std::string queryMode;
extern bool gp_useJunctions;
extern colouredDBG SR_DBG;
extern bool gp_reverse;

extern CharString io_seqFile;
extern CharString io_outFile;
extern std::string io_pathToJF;
extern std::string io_JFcountsTable;
extern std::string io_JFjunctcountsTable;

extern std::string io_logFile;
extern std::string io_statFile;

extern int g_NTHREADS;

extern bool DEBUG_USER;
extern bool DEBUG_MODE;
extern bool DEBUG_TEST;


/***********************************************************************
 * MAIN
 ***********************************************************************/

int main(int argc, const char * argv[]) {
	
	std::cout << "******************************************************" << std::endl;
	std::cout << "* TALC : Transcriptome-Aware Long Read Correction    *" << std::endl;
	std::cout << "*----------------------------------------------------*" << std::endl;
  std::cout << "*                                                    *" << std::endl;
  std::cout << "* Kmers are assumed directional                      *" << std::endl;
  std::cout << "******************************************************" << std::endl;

//parse arguments

	 ArgumentParser parser("TALC: Transcriptome-Aware Long Read Correction");
	 setShortDescription(parser, "Hybrid Long Read Correction using Short Read coverage");
	 setVersion(parser, "1.01");
	 setDate(parser, "September 2019");
    
    //INPUT sequences

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "Input fastQ/A file containing Long Reads"));
    
    //OUTPUTs
    
    addOption(parser, seqan::ArgParseOption("o", "output", "Prefix to be used for output files",
                                            seqan::ArgParseArgument::STRING, "TEXT"));
    setDefaultValue(parser, "output", "out");
   
    //----------------------------------------------------------------------------------------//
    // Parameters for the De Bruijn Graph

    addOption(parser, seqan::ArgParseOption("k", "kmerSize", "length k of k-mers.",
                                            seqan::ArgParseArgument::INTEGER, "UNSIGNED INT"));
    setRequired(parser, "k", true);
    setMaxValue(parser, "k", "30");
    setMinValue(parser, "k", "18");

    
    // Parameters for query mode
    
    addOption(parser, seqan::ArgParseOption("qm", "query-mode", "Mode that should be used to query kmers (advised: memory). ", seqan::ArgParseArgument::STRING, "TEXT"));
    setDefaultValue(parser, "query-mode", "memory");
    setValidValues(parser, "query-mode", "memory jellyfish2");
   
    addOption(parser, seqan::ArgParseOption("SR", "SRCounts", "Short reads kmer counts in format .dump file (memory mode)  or .jf file (jf2 mode), obtained from Jellyfish2",
                                            seqan::ArgParseArgument::STRING, "TEXT"));
    setRequired(parser, "SRCounts", true);  
                                  
    addOption(parser, seqan::ArgParseOption("j", "junctions", "k-mers flanking junctions and their counts.",
                                            seqan::ArgParseArgument::STRING, "TEXT"));      
    setRequired(parser, "junctions", false);                                  

    addOption(parser, seqan::ArgParseOption("jf2", "pathToJF2", "Specifies where the Jellyfish2 program should be found.",
                                            seqan::ArgParseArgument::STRING, "TEXT"));
    setDefaultValue(parser, "jf2", "");  
          
    //----------------------------------------------------------------------------------------//
    //Parameters and thresholds for correction
    
    addOption(parser, seqan::ArgParseOption("MIN_INNER_SCORE", "MIN_INNER_SCORE", "Minimum %ID score required for the best-correction-path-between-any-two-solid-regions to be kept.",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(parser, "MIN_INNER_SCORE", "0.7");
    setMaxValue(parser, "MIN_INNER_SCORE", "0.9");
    setMinValue(parser, "MIN_INNER_SCORE", "0.3");
    
    addOption(parser, seqan::ArgParseOption("MIN_BORDER_SCORE", "MIN_BORDER_SCORE", "Minimum %ID score required for the best-correction-path-on-left-or-right-borders to be kept.",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(parser, "MIN_BORDER_SCORE", "0.7");
    setMaxValue(parser, "MIN_BORDER_SCORE", "0.9");
    setMinValue(parser, "MIN_BORDER_SCORE", "0.5");
      
    //----------------------------------------------------------------------------------------//
    //Parameters and thresholds for correction
    
    addOption(parser, seqan::ArgParseOption("MIN_COUNT", "MIN_COUNT", "Minimal count for a k-mer to be kept in the SR-de Bruijn Graph.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "MIN_COUNT", "2"); 
    setMinValue(parser, "MIN_COUNT", "2"); 

    addOption(parser, seqan::ArgParseOption("SR_ERROR_RATE", "SR_ERROR_RATE", "Prior estimate of the error rate in short reads.",
                                            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "SR_ERROR_RATE", "0.025");
    setMaxValue(parser, "SR_ERROR_RATE", "0.1");
    setMinValue(parser, "SR_ERROR_RATE", "0.01");
      
    addOption(parser, seqan::ArgParseOption("WINDOW_SIZE", "WINDOW_SIZE", "[ADVANCED] Size of the window. The larger, the more fresh air.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "WINDOW_SIZE", "9"); 
    setMinValue(parser, "WINDOW_SIZE", "6"); 
    
    addOption(parser, seqan::ArgParseOption("MAX_NB_BRANCHES", "MAX_NB_BRANCHES", "[ADVANCED] Maximal number of competing branches.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "MAX_NB_BRANCHES", "7"); 
    setMinValue(parser, "MAX_NB_BRANCHES", "5"); 
    
    //----------------------------------------------------------------------------------------//
    //EXPECTATION INTERVALS:  
    
    //ALPHA COEFFICIENT
    addOption(parser, seqan::ArgParseOption("ALPHA_FOR_PRED", "ALPHA_FOR_PRED", "[ADVANCED] Coefficient used to build count confidence interval: [count +/- ALPHA*sqrt(count)].",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(parser, "ALPHA_FOR_PRED", "2.57"); 
    setMinValue(parser, "ALPHA_FOR_PRED", "0.67");  
                             
    //----------------------------------------------------------------------------------------//
    //Other
    
    addOption(parser, seqan::ArgParseOption("t", "num_threads", "number of threads",
                                            seqan::ArgParseArgument::INTEGER, "UNSIGNED INT"));
    setDefaultValue(parser, "t", 1);
    setMinValue(parser, "num_threads", "1");
    
    addOption(parser, seqan::ArgParseOption("DEBUG_MODE", "DEBUG_MODE", "Will activate output for debugging purposes",
                                            seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("rev", "reverse", "If set, long reads will be reverse complemented before correction by short reads."));                                       
                                    
std::cout << "[TALC]: Parsing arguments" << std::endl; 
seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
if (res != seqan::ArgumentParser::PARSE_OK) return res == seqan::ArgumentParser::PARSE_ERROR;
else
 {
  bool verbose(true);
  setSettings(parser, verbose); 
  setBasicReadStatsHeader(io_statFile);
  
  if( g_NTHREADS>1 ) DEBUG_USER = false; 
 
  //All about (long read) sequences to be corrected
  StringSet<CharString> ids;
  StringSet<TSeq > mySeqs;
  StringSet<CharString > myScorrs;
        
  clock_t mapstart, mapend, corrstart, corrend;
  CharString scorr;
         
 std::cout << "[TALC]: Attempting to load sequences." << std::endl;
 corrstart = clock();
  
 if( loadSeqData(ids, mySeqs, io_seqFile)==0 )
  {
	std::cout << "[TALC]: Hmm...it seems the sequence file is OK." << std::endl;
    std::cout << "[TALC]: " << length(ids) << " long read(s) loaded" << std::endl;
       
  if( queryMode=="memory" )
   {   	   
    if( (queryMode == "memory") &  gp_useJunctions ) std::cout << "[TALC]: Building the SR-cdBG from count files: " << io_JFcountsTable << " and " << io_JFjunctcountsTable << std::endl;
    else  std::cout << "[TALC]: Building the SR-dBG from count file: " << io_JFcountsTable << std::endl;
    
    mapstart = clock();
    
    SR_DBG =  buildCDBG(1, io_JFcountsTable, io_JFjunctcountsTable);
    decolourRepeatsFromDBG(SR_DBG, K);
    
    mapend = clock();
    
    std::cout << "[TALC]: SR-dBG created. Time needed: " << ((mapend - mapstart)/CLOCKS_PER_SEC)/60 << " min." <<std::endl;
    std::cout << "[TALC]: SR-dBG contains " << SR_DBG.size() << " nodes." << std::endl; 
   } 
   	
  if( queryMode=="jellyfish2" || !empty(SR_DBG))
   {  
    omp_set_num_threads( g_NTHREADS );
    
    std::cout << "[TALC]: " << "Good news, there are nodes in the de Bruijn Graph." << std::endl;
    std::cout << "[TALC]: " << "Maybe we can try and correct some long reads, then?" << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(long unsigned int r=0; r<length(ids); r++)
    {
     clock_t start, end;                     //records time taken for correction
          
     //Cut an store Left and Right A-tails  
     if( gp_reverse ) reverseComplement(mySeqs[r]);
      
     //if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: Cutting A or T tails" << std::endl;
     //mySeqs[r] = trimAT(mySeqs[r], gp_minTailWidth);
 
     Read myLRead(ids[r], mySeqs[r]);
     if( DEBUG_USER ) myLRead.display();
         
     start = clock(); 
     if( myLRead.getLength()>K )
      { 
	   try
	    {     		  	 
         //Call k-mer count database to get the coverage of the LR in the SR-dBG
         if( myLRead.reCoverage() )
          {

           if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: Coverage recovered" << std::endl;

           if( myLRead.defineStructure2() )
           { 
            if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: STRUCTURE DEFINED" << std::endl;
            if( DEBUG_USER ) myLRead.displaySegmentation();

            myLRead.correct2();

            if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: CORRECTION PERFORMED" << std::endl;

           //Once LR has been corrected, append the left and right tails to its new sequence 
           //append(leftTail, myLRead.getCorrSeq());
           //append(leftTail, rightTail)	;             	                  
           // mySeqs[r] = leftTail;
           mySeqs[r] = myLRead.getCorrSeq();
          if( gp_reverse ) reverseComplement(mySeqs[r]);
               
            if( DEBUG_USER ) std::cout << mySeqs[r] << std::endl; 
         
           }else throwToLog(myLRead.getName(), "Unable to define convenient structure.", io_logFile);  
         
         }else{
          if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: No solid kmer found." << std::endl;
          throwToLog(myLRead.getName(), "No solid kmer could be found.", io_logFile); 
          if( DEBUG_USER ) std::cout << "[DEBUG_MODE]: Saved in log" << std::endl;
        }

	    }catch(std::string const& chaine)
        {
          myLRead.display(); 
          std::cout << chaine << std::endl;
          throwToLog(myLRead.getName(), chaine, io_logFile);
		    }	
	  //Append some stats to the *.stat_basics.txt file	
      //myLRead.outputBasicReadStats(io_statFile);                                                        
       }       
     end = clock();      
    } //END-FOR all imported reads
   //Output results as a fasta file 
   outputSeqData(ids, mySeqs, io_outFile);  
   corrend = clock();
   std::cout << "[TALC]: " << "Looks like we are done now." << std::endl;
   std::cout << "[TALC]: The correction needed: " << (corrend - corrstart)/(CLOCKS_PER_SEC*60) << " min to be completed." << std::endl;     

    return 0;  
   } //END if dBF construction ok (when required)
    else
    {
	  std::cout << "[TALC]: The de Bruijn Graph is empty...Correction aborted." << std::endl;
	  return 1;
	} 
   }  //END-IF data can be loaded
  else std::cout << "[TALC]: ISSUE WITH INPUT FILES" << std::endl; 
  }//END-IF COMMAND LINE PARSING IS OK
}
