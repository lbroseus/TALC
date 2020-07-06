//
//  Settings.hpp
//
//  -> Class that store chosen parameters and paths to needed files along the correction step
//  -> Called by Explorer objects

#ifndef Settings_hpp
#define Settings_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/arg_parse.h>

#include "Settings.hpp"

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;

/***********************************************************************/
   
void setSettings(seqan::ArgumentParser& parser, bool& verbose);
    
void displaySettings();
    
void outputConfig(std::string& outPrefix);
    
 
#endif /* Settings_hpp */
