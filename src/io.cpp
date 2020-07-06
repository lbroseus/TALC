//
//  io.cpp
//
//  Functions related to the importation and exportation of data
//
//  Created by Lucile Broseus on 26/10/2017.
//  Copyright © 2017 LuB. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "io.hpp"

using namespace seqan;

typedef Dna5String TSeq;

int loadSeqData(StringSet<CharString>& ids,
                StringSet<TSeq>& Seqs,
                CharString& seqFile)
{
    SeqFileIn seqFileIn;

    if (!open(seqFileIn, toCString(seqFile)))
    {
        std::cerr << "ERROR: Could not open file " << seqFile << "\n";
        return 1;
    }
    try
    {
        readRecords(ids, Seqs, seqFileIn);
        return 0;
    }
    catch (Exception const& e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

}

int outputSeqData(StringSet<CharString>& ids,
                  StringSet<TSeq>& Seqs,
                  CharString& outFile)
{
    std::cout << "Specified output file name: " << outFile << std::endl;

SeqFileOut seqFileOut;

    if (!open(seqFileOut, toCString(outFile)))
    {
        std::cerr << "ERROR: Could not open the file " << outFile << "\n";
        return 1;
    }

    try
    {
        writeRecords(seqFileOut, ids, Seqs);
    }
    catch (Exception const& e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

int outputSequalData(StringSet<CharString>& ids,
                  StringSet<TSeq>& Seqs,
                  seqan::StringSet<seqan::CharString>& Quals,
                  CharString& outFile)
{
    std::cout << "Specified output file name: " << outFile << std::endl;

SeqFileOut seqFileOut;

    if (!open(seqFileOut, toCString(outFile)))
    {
        std::cerr << "ERROR: Could not open the file " << outFile << "\n";
        return 1;
    }

    try
    {
        writeRecords(seqFileOut, ids, Seqs, Quals);
    }
    catch (Exception const& e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

void throwToLog(seqan::CharString seqName, std::string chaine, std::string outFile)
{ 
	std::ofstream outputFile;
    outputFile.open(outFile, std::ios_base::app);      
    outputFile << "[Read: " << seqName << " ]: " << chaine << std::endl;
    outputFile.close();
}

/////////////


template<class Type > std::vector<Type > make_empty_vector(Type x)
{
	std::vector<Type > t;
	return t;
}

std::vector<int > make_empty_IntVector()
{
    std::vector<int > v;

    return v;
}

std::vector<double > make_empty_DoubVector()
{
    std::vector<double > v;

    return v;
}

std::vector<float > make_empty_FloatVector()
{
    std::vector<float > v;

    return v;
}



