//
//  Jellyfish.hpp
//  LReader
//
//  Created by Lucile Broseus on 16/10/2017.
//  Copyright © 2017 LuB. All rights reserved.
//

#ifndef Jellyfish_hpp
#define Jellyfish_hpp

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

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;
typedef std::pair<seqan::Dna5, unsigned int > Base; 

//----------------------------------------------------------------------/                                                                                         

std::map<TSeq, colouredCount > buildDBGfromDump(std::string& countsTable, std::string& junctionCountsTable);

std::map<TSeq, colouredCount > buildCDBG(int mode, std::string& countsTable, std::string& junctionCountsTable);								
									            
//----------------------------------------------------------------------/

std::vector<colouredCount > getNextCounts(TSeq& kmer,
                                          Direction direction,
                                          colouredDBG& dBG, 
                                          std::string& countsTable, std::string& junctionCountsTable);
                                          
std::vector<colouredCount > getNextCountsFromDBG(TSeq& kmer, Direction direction, colouredDBG& dBG);                                          

std::vector<colouredCount > getNextCountsFromJF(TSeq& kmer, Direction direction,std::string& countsTable, std::string& junctionCountsTable);

int getOutDegree(TSeq& kmer, Direction direction, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable); 
                                                      
//----------------------------------------------------------------------/
                                                       
colouredCount getCount(TSeq& kmer, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable);                                                      
                                                       
colouredCount getCountFromDBG(TSeq& kmer, colouredDBG& dBG);                                                  
									            															    
colouredCount getCountFromJF(TSeq& kmer, std::string& countsTable, std::string& junctionCountsTable);

//----------------------------------------------------------------------/

std::vector<colouredCount> getLRCountsInSR(TSeq& Seq, unsigned int& kmerSize, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable);

std::vector<colouredCount> getLRCountsInSRFromDBG(TSeq& Seq, unsigned int& kmerSize, colouredDBG& dBG);
                                             
std::vector<colouredCount> getLRCountsInSRFromJF(TSeq& Seq, unsigned int& kmerSize, std::string& countsTable, std::string& junctionCountsTable);  
                                                                                                                                                   
//----------------------------------------------------------------------/

std::vector<TSeq> getSuccessors(TSeq& kmer, Direction direction);


#endif /* Jellyfish_hpp */
