//
//  utils.hpp
//  LReader
//
//  Created by Lucile Broseus on 10/10/2017.
//  Copyright © 2017 LuB. All rights reserved.
//

/****************************************************************************
 * Date: 10/10/2017
 *
 ****************************************************************************
 * Object: some functions to:
 * - query KMC (k-mer counts table);
 * - output statistics;
 * - maybe more?
 ****************************************************************************
 * Updates:
 *
 ****************************************************************************/

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/find.h>

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;
typedef std::pair<seqan::Dna5, unsigned int > Base;
typedef std::tuple<TSeq, unsigned int, unsigned int > anchorTuple;
typedef std::pair<TSeq, colouredCount > tipNode; 


typedef seqan::StringSet<TSeq> TStringSet;                    // container for strings
typedef seqan::StringSet<TSeq, seqan::Dependent<> > TDepStringSet;   // dependent string set
typedef seqan::Graph<seqan::Alignment<TDepStringSet> > TAlignGraph;

enum Confidence { SURE, UNSURE };
enum Complexity { LOW, HIGH };
enum Outcome { NO, DEADEND, SUCCESS, ABORTION , BRANCHING, CYCLE};
enum Status { EXPECTED, UNEXPECTED, LOWCOUNT, SUPPORTED, CORRECTED, UNCORRECTED, ABSENT, BREAKPOINT};
enum Bifurcation { MAJOR, MINOR };
enum Direction { LEFT, RIGHT };
enum Location { HEAD, INNER, TAIL, UNKNOWN };

std::string getStdoutFromCommand(std::string cmd);

unsigned int computeQuantile(std::vector<unsigned int > const& counts, unsigned int const& whichQuantile);

int outputSequence(seqan::Dna5String mySeq, seqan::CharString Id, seqan::CharString outFileName);

seqan::CharString make_empty_CharString();

TSeq make_empty_Dna5String();

tipNode make_empty_tipNode();

std::vector<unsigned int > make_vector(unsigned int x);

std::vector<unsigned int > make_empty_vector();

std::vector< std::tuple<int,Status,int> > make_empty_vect_of_NodeInfo();

std::vector<bool > make_empty_bool_vector();

std::vector<std::pair<int, int> > make_empty_vecOfIntPairs();

std::vector<anchorTuple > make_empty_vecOfAnchors();

std::vector<colouredCount > make_vecOfColouredPairs(colouredCount counts);

std::vector<colouredCount > make_empty_vecOfColouredPairs();

colouredCount make_empty_colouredCount();

std::vector<std::pair<seqan::Dna5, colouredCount > > make_vecOfColouredBases(seqan::Dna5 base, colouredCount counts);

std::vector<std::pair<seqan::Dna5, colouredCount > > make_empty_vecOfColouredBases();

std::vector<TSeq > make_empty_vectOfTSeq();

std::deque<int > make_empty_IntDeque();

std::deque<Status > make_empty_StatusDeque();

DBG make_empty_DBG();

std::pair<TSeq, Status> make_empty_border();
std::vector<std::pair<TSeq, Status> > make_empty_inner();

std::tuple<unsigned int, unsigned int, Status > make_empty_KmPos();

std::vector<std::tuple<unsigned int, unsigned int, Status > > make_empty_vectOfKmPos();

bool descComparePair(std::pair<seqan::Dna5, unsigned int> pair1,
                     std::pair<seqan::Dna5, unsigned int> pair2);
                     
bool descCompareColouredPair(std::pair<seqan::Dna5, colouredCount > pair1,
                     std::pair<seqan::Dna5, colouredCount > pair2);  
                     
TSeq cutATail(TSeq& Seq, unsigned int& minWidth);

TSeq trimAT(TSeq& Seq, unsigned int& minWidth);
                     
double predIntLow(double count, double coeff);

double predIntUp(double count, double coeff);

bool isExpected(unsigned int& next_count,
                 unsigned int& current_count,
                 double coeff);  
                 
TSeq formNextKmer(TSeq kmer, seqan::Dna5 new_base, Direction direction);

double overlapAlign(TSeq reference, TSeq candidate, int mode);  

std::pair<double, double> computeOverlapScore(TSeq& reference, TSeq& sequence, Direction direction);       

std::tuple<double, int, int> InnerEditAlignment(TSeq& reference, TSeq& pathSeq);  
std::tuple<double, int, int> BorderEditAlignment(TSeq& reference, TSeq& pathSeq, Direction& direction);

// SEPT 18: new functions

TSeq getKmerAt(TSeq& refSequence, unsigned int position, unsigned int& kmerSize);
 
TSeq extractSolidSequence(TSeq& refSequence, unsigned int start_kmPosition, unsigned int end_kmPosition, unsigned int& kmerSize);
TSeq extractWeakSequence(TSeq& refSequence, unsigned int end_kmPosition, unsigned int start_kmPosition, unsigned int& kmerSize);
TSeq extractWeakBorderSequence(TSeq& refSequence, unsigned int kmPosition, unsigned int& kmerSize, Location m_location);      

double absdistance(TSeq& seq1, TSeq& seq2);

void decolourRepeatsFromDBG(colouredDBG& dBG, unsigned int& kmerSize);

#endif /* utils_hpp */
