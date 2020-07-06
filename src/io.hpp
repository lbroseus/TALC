//
//  io.hpp
//
//  Functions related to the importation and exportation of data...
//

#ifndef io_hpp
#define io_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

typedef seqan::Dna5String TSeq;

int loadSeqData(seqan::StringSet<seqan::CharString>& ids,
                seqan::StringSet<TSeq>& Seqs,
                seqan::CharString& seqFile);
                
int outputSeqData(seqan::StringSet<seqan::CharString>& ids,
                  seqan::StringSet<TSeq>& Seqs,
                  seqan::CharString& outFile);               

int outputSequalData(seqan::StringSet<seqan::CharString>& ids,
                  seqan::StringSet<TSeq>& Seqs,
                  seqan::StringSet<seqan::CharString>& Quals,
                  seqan::CharString& outFile);
                  
void throwToLog(seqan::CharString seqName, std::string chaine, std::string outFile);

std::vector<int > make_empty_IntVector();

std::vector<double > make_empty_DoubVector();

std::vector<float > make_empty_FloatVector();

std::string concatIntVectToString(std::vector<int > const& vector);

#endif /* io_hpp */
