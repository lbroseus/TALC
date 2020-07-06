//
//  Settings.cpp
//

#include "Settings.hpp"
#include "utils.hpp"

#include <stdio.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <string>
#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>


using namespace seqan;

typedef Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;

/***********************************************************************
 * Global variables
 ***********************************************************************/
 
unsigned int K;

unsigned int gp_MIN_COUNT;

double gp_ALPHA;
unsigned int gp_WINDOW_SIZE;
double gp_SR_ERROR_RATE(0.025);
double gp_LR_ERROR_RATE;
double gp_MIN_INNER_SCORE;
double gp_MIN_BORDER_SCORE;

unsigned int gp_MAX_NB_COMPETING_PATHS(7);

unsigned int gp_minTailWidth(3);

std::string queryMode;
bool gp_useJunctions;
colouredDBG SR_DBG;
bool gp_reverse;

CharString io_seqFile;
CharString io_outFile;
std::string io_pathToJF;
std::string io_JFcountsTable;
std::string io_JFjunctcountsTable;
std::string io_LRcountsTable;

std::string io_logFile;
std::string io_statFile;

int g_NTHREADS;

bool DEBUG_MODE;
bool DEBUG_TEST(true);
bool DEBUG_USER(true);

bool DEBUG;

/***********************************************************************/


void setSettings(seqan::ArgumentParser& parser, bool& verbose){
    
    std::string outPrefix;
    
    std::string mode;
	//INPUT Long read sequences
    getArgumentValue(io_seqFile, parser, 0);  
    //kmer counts in short reads:
    getOptionValue(io_JFcountsTable, parser, "SRCounts");
    //junction kmer counts in short reads:
    getOptionValue(io_JFjunctcountsTable, parser, "junctions");
    //kmer sizes
    getOptionValue(K, parser, "k");
    //query mode
    getOptionValue(queryMode, parser, "qm");
    
                             // Optional arguments //
    //OUTPUT sequences

    getOptionValue(outPrefix, parser, "o");
    io_outFile = outPrefix + ".fa";
    io_statFile = outPrefix + ".stats_basics.txt";
    io_logFile = outPrefix + ".log";

    //Jellyfish2 query program
    getOptionValue(io_pathToJF, parser, "jf2");
    
    //score thresholds for correction:
    getOptionValue(gp_SR_ERROR_RATE, parser, "SR_ERROR_RATE");
    getOptionValue(gp_MIN_INNER_SCORE, parser, "MIN_INNER_SCORE");
    getOptionValue(gp_MIN_BORDER_SCORE, parser, "MIN_BORDER_SCORE");

    gp_useJunctions = isSet(parser, "junctions");
    gp_reverse = isSet(parser, "reverse");

    
    //Parameters for dBG exploration
    getOptionValue(gp_MIN_COUNT, parser, "MIN_COUNT");
    getOptionValue(gp_MAX_NB_COMPETING_PATHS, parser, "MAX_NB_BRANCHES");
    getOptionValue(gp_ALPHA, parser, "ALPHA_FOR_PRED");
    getOptionValue(gp_WINDOW_SIZE, parser, "WINDOW_SIZE");

    //OTHER
    getOptionValue(g_NTHREADS, parser, "t");
    
    DEBUG_MODE = isSet(parser, "DEBUG_MODE");

    if( verbose ) displaySettings();
    outputConfig(outPrefix);
}

void displaySettings()
{
	std::cout << "******************************" << std::endl;
	std::cout << "*           SETTINGS         *" << std::endl;
	std::cout << "******************************" << std::endl;
    std::cout << "INPUTS SUMMARY" << std::endl;
    std::cout << "INPUT sequences: " << io_seqFile << std::endl;
    std::cout << "Kmer size: " << K << std::endl;
    std::cout << "Query mode:" << queryMode << std::endl;
    std::cout << "Use junctions? " << gp_useJunctions << std::endl;
    std::cout << "Specified table with SR kmer counts: " << io_JFcountsTable << std::endl;
    
    if( gp_useJunctions )
     {
	  std::cout << "Junction mode activated " << std::endl;
	  std::cout << "Specified table with junction kmer counts: " << io_JFcountsTable << std::endl;
	 }
    std::cout << "---------------------------" << std::endl;
    std::cout << "MIN_INNER_SCORE: " << gp_MIN_INNER_SCORE << std::endl;
    std::cout << "MIN_BORDER_SCORE: " << gp_MIN_BORDER_SCORE << std::endl;  
    std::cout << "MAX_COMPETING_PATH: " << gp_MAX_NB_COMPETING_PATHS << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "OUTPUT: " << io_outFile<< std::endl;
    std::cout << "STATS: " << io_statFile << std::endl;
    std::cout << "****************************" << std::endl;
    std::cout << "EXPLORATION PARAMETERS" << std::endl;
    std::cout << "ALPHA: " << gp_ALPHA << std::endl;
    std::cout << "MIN_COUNT: " << gp_MIN_COUNT << std::endl;
    std::cout << "WINDOW_SIZE: " << gp_WINDOW_SIZE << std::endl;
    std::cout << "MAX_NB_BRANCHES: " << gp_MAX_NB_COMPETING_PATHS << std::endl;
    std::cout << "NTHREADS: " << g_NTHREADS << std::endl;
    std::cout << "****************************" << std::endl;
}

// Add starting time
void outputConfig(std::string& outPrefix)
{  
	std::string out;
    out = outPrefix + ".config.txt";
	
	std::ofstream outputFile;
    outputFile.open(out, std::ios_base::trunc);      
    outputFile << "TALC: Parameters used for sample: " << outPrefix << "\n"
               << "****************************" << "\n"
               << "INPUT=" << io_seqFile << "\n"
               << "OUTPUT=" << outPrefix << "\n"
               << "STATS=" << io_statFile << "\n"
               << "****************************" << "\n"
               << "KmerSize=" << K << "\n"
               << "Junction mode activated? " << gp_useJunctions  << "\n"
               << "queryMode=" << queryMode << "\n"
               << "****************************" << "\n"
               << "MIN_INNER_SCORE=" << gp_MIN_INNER_SCORE << "\n"
               << "MIN_BORDER_SCORE=" << gp_MIN_BORDER_SCORE << "\n"  
               << "MAX_NB_BRANCHES=" << gp_MAX_NB_COMPETING_PATHS << "\n"
               << "ALPHA=" << gp_ALPHA << "\n"
               << "MIN_SR_COUNT=" << gp_MIN_COUNT << "\n"
               << "WINDOW_SIZE=" << gp_WINDOW_SIZE << "\n"
               << "****************************" << std::endl;
    outputFile.close();
}
