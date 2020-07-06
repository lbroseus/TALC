//
//  utils.cpp
//  LReader
//
//  Created by Lucile Broseus on 10/10/2017.
//  Copyright © 2017&2018&2019 LuB. All rights reserved.
//
/****************************************************************************
 * Version TALC 1.01
 ****************************************************************************
 ****************************************************************************
 * Updates:
 * 22/08/2019: fixed cutATail() function
 ****************************************************************************/

#include "utils.hpp"

#include "Settings.hpp"

#include <iostream>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/find.h>
#include <seqan/align.h>

using namespace seqan;

typedef Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<Dna5, colouredCount > colouredBase; 
typedef std::map<TSeq, unsigned int > DBG;
typedef std::pair<seqan::Dna5, unsigned int > Base; 
typedef std::tuple<TSeq, unsigned int, unsigned int > anchorTuple;
typedef std::pair<TSeq, colouredCount > tipNode; 


typedef Align<TSeq, ArrayGaps> TAlign;     // align type
typedef Row<TAlign>::Type TRow;                 // gapped sequence type
	    typedef Iterator<TRow>::Type TRowIterator;
	    
//---------------------------------------------------------------------//

extern unsigned int gp_MIN_COUNT;
extern unsigned int K;

extern std::string queryMode;
extern colouredDBG SR_DBG;
extern DBG LR_DBG;

extern std::string pathToJF;
extern std::string JFcountsTable;
extern std::string JFjunctcountsTable;
extern std::string LRcountsTable;

extern std::string logFile;


extern std::string logFile;

static std::vector<char > Dict = {'A', 'C', 'G', 'T'};

//---------------------------------------------------------------------//	    
	    
/*****************************************************
 * TOADD
 *    GRAPH
 * -> fonction pour l'initialisation du graphe de correction
 * -> fonction sortie graph.dot
 *
 *    CORRECTION
 * -> setSource(): k-1 mer d'une région solide, d'où partir
 * -> setTarget(): k-mer de la région solide suivante, à viser
 * ->
 *
 *    QUERY
 * -> query KMC or squeakr count tab for successors
 *****************************************************/



/*****************************************************
 *
 *****************************************************/

std::string getStdoutFromCommand(std::string cmd) {
    
    std::string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}

unsigned int computeQuantile(std::vector<unsigned int > const& counts, unsigned int const& whichQuantile)
{
    //remove zero counts
    std::vector<unsigned int > posCounts;
    for(unsigned int i(0); i<counts.size(); i++)
    {
        if(counts[i]>0) posCounts.push_back(counts[i]>0);
    }
    
    sort(posCounts);
    
    return posCounts[posCounts[(posCounts.size()/4)*whichQuantile-1]];
}

//Screen beginning/end of a sequence for AAA or TTT tails
//Returns the AAA tail if exists
//Cut the corresponding tail from Seq (seq is then modified)

// A/T 5' prime Border trimming
// minWidth: minimum tail length to consider it a tail
TSeq cutATail(TSeq& Seq, unsigned int& minWidth)
{
        unsigned int pos;
        CharString trimmedSeq;
        char firstBase;
        bool goesOn;
        

        pos = 0;
        firstBase = Seq[pos];
        goesOn =( (firstBase == (char)'A') ||  (firstBase == (char)'T') );
          
        while(goesOn & (pos<length(Seq)) ) 
        {
          pos++;
          goesOn = (Seq[pos]==firstBase);
        }
           
        if( pos>(minWidth-1) ) trimmedSeq = suffix(Seq, pos);
        else trimmedSeq = Seq;

        return trimmedSeq;    
}

TSeq trimAT(TSeq& Seq, unsigned int& minWidth)
{
    TSeq newSeq(Seq);

    reverse(newSeq);
    newSeq = cutATail(newSeq, minWidth);
    reverse(newSeq);
    newSeq = cutATail(newSeq, minWidth);

    return( newSeq );
}

int outputSequence(Dna5String mySeq, CharString Id, CharString outFileName)
{
    std::cout << "Specified output file name: " << outFileName << std::endl;

    SeqFileOut seqFileOut;
    
    if (!open(seqFileOut, toCString(outFileName)))
    {
        std::cerr << "ERROR: Could not open the file " << outFileName << "\n";
        return 1;
    }
    
    try
    {
        writeRecord(seqFileOut, Id, mySeq);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

CharString make_empty_CharString()
{
    CharString v;

    return v;
}


TSeq make_empty_Dna5String()
{
    TSeq v;
    return v;
}

tipNode make_empty_tipNode()
{
 tipNode t;
 return t;
}
	


std::vector<unsigned int > make_vector(unsigned int x)
{
    std::vector<unsigned int > v = {x};
    
    return v;
}

std::vector<unsigned int > make_empty_vector()
{
    std::vector<unsigned int > v;
    
    return v;
}

std::vector< std::tuple<int,Status,int> > make_empty_vect_of_NodeInfo()
{
	std::vector< std::tuple<int,Status,int> > vector;
	return vector;
}

std::vector<bool > make_empty_bool_vector()
{
    std::vector<bool > v;
    
    return v;
}

std::vector<std::pair<int, int> > make_empty_vecOfIntPairs()
{
    std::vector<std::pair<int ,int > > v;
    
    return v;
}

std::vector<TSeq > make_empty_vectOfTSeq()
{
    std::vector<TSeq > v;
    
    return v;
}

std::vector<anchorTuple > make_empty_vecOfAnchors()
{
    std::vector<anchorTuple > v;
    
    return v;
}

std::vector<colouredCount > make_vecOfColouredPairs(colouredCount counts)
{
    std::vector<colouredCount> v = {counts};
    return v;
}

std::vector<colouredCount > make_empty_vecOfColouredPairs()
{
    std::vector<colouredCount > v;
    return v;
}

colouredCount make_empty_colouredCount()
{
    colouredCount v;
    return v;
}

std::vector<std::pair<seqan::Dna5, colouredCount > > make_vecOfColouredBases(seqan::Dna5 base, colouredCount counts)
{
    std::vector<std::pair<seqan::Dna5, colouredCount > > v = {std::make_pair(base,counts)};
    return v;
}

std::vector<std::pair<seqan::Dna5, colouredCount > > make_empty_vecOfColouredBases()
{
    std::vector<std::pair<seqan::Dna5, colouredCount > > v;
    return v;
}
 
std::deque<int > make_empty_IntDeque()
{
    std::deque<int > deq;
    return deq;
}

std::deque<Status > make_empty_StatusDeque()
{
    std::deque<Status > deq;
    return deq;
}

std::tuple<unsigned int, unsigned int, Status > make_empty_KmPos()
{
  std::tuple<unsigned int, unsigned int, Status > pos;
  return pos;	
}

std::vector<std::tuple<unsigned int, unsigned int, Status > > make_empty_vectOfKmPos()
{
  std::vector<std::tuple<unsigned int, unsigned int, Status > > pos;
  return pos;	
}

DBG make_empty_DBG()
{
	DBG dBG;
	return dBG;
}

std::pair<TSeq, Status> make_empty_border()
{
	std::pair<TSeq, Status> border;
	return border;
}

std::vector<std::pair<TSeq, Status> > make_empty_inner()
{
	std::vector<std::pair<TSeq, Status> > inners;
	return inners;
}

bool descComparePair(std::pair<Dna5, unsigned int> pair1,
                     std::pair<Dna5, unsigned int> pair2)
{
    return ( std::get<1>(pair1)>std::get<1>(pair2) );
}

bool descCompareColouredPair(std::pair<Dna5, colouredCount> pair1,
                     std::pair<Dna5, colouredCount> pair2)
{
    return ( std::get<1>(pair1)>std::get<1>(pair2) );
}

/* Computation of prediction interval lower/upper bounds */

double predIntLow(double count, double coeff)
{
    return (count - coeff*sqrt(count));
}

double predIntUp(double count, double coeff)
{
    return (count + coeff*sqrt(count));
}

/* Check whether a given count falls into the predicted interval */

bool isExpected(unsigned int& next_count,
                 unsigned int& current_count,
                 double coeff)

{
  return ( (predIntLow((double)current_count, coeff)<=(double)next_count)
        & ((double)next_count <= predIntUp((double)current_count,coeff)) );
}

/* Return next kmer from a new base */

TSeq formNextKmer(TSeq kmer, Dna5 new_base, Direction direction)
{
	TSeq newKmer;
	
	if(direction == RIGHT)
	 {
      newKmer = suffix(kmer, 1);
      appendValue(newKmer, new_base);
     }
    else
     {
	  TSeq tmp(new_base);
	  append(tmp,kmer);
	  newKmer = prefix(tmp, length(kmer));
	 } 

    return newKmer;
}

double overlapAlign(TSeq reference, TSeq candidate, int mode)
{
	
    TStringSet sequences;
    double score(0);
    
    //Scoring scheme:
    int match = 1;
    int mismatch = -1;
    int indel = -1;
    
    if( length(reference)>=length(candidate) )
    {
      appendValue(sequences, reference);
      appendValue(sequences, candidate);
    }
    else
    {
       appendValue(sequences, reference);
       appendValue(sequences, candidate);
    }
    
    TAlignGraph alignG(sequences);
     
    switch( mode )
    {
     //Global alignment, penalize both left and right gaps 
     case 0:
     score = globalAlignment(alignG, Score<int, Simple>(match, mismatch, indel),
                                AlignConfig<false, false, false, false>(), LinearGaps());  
     break;
    
     case 1:
     //Correction from left to right, do not penalize gap on the right                            
     score = globalAlignment(alignG, Score<int, Simple>(match, mismatch, indel),
                                AlignConfig<false, false, true, true>(), LinearGaps());           
     break;
      //Correction from right to left, do not penalize gap on the left                            
     case 2:
     score = globalAlignment(alignG, Score<int, Simple>(match, mismatch, indel),
                                AlignConfig<true, true, false, false>(), LinearGaps());
     break; 
      //Global alignment, de-penalize both left and right gaps 
     default :
     score = globalAlignment(alignG, Score<int, Simple>(match, mismatch, indel), AlignConfig<true, true, true, true>(), LinearGaps()); 
    }    
    
    return score;
}

std::pair<double, double> computeOverlapScore(TSeq& reference, TSeq& sequence, Direction direction)
{
   double score(0);
   double regscore(0);
   int mode;
   
   TSeq truncRef;
   if(length(reference)>length(sequence)*1.1) truncRef = prefix(reference, length(sequence)*1.2);
   else truncRef = reference;
   switch( direction )
	{
	 case LEFT : 
	   mode = 2;
	 break;
	 case RIGHT : 
	   mode = 1;
	 break;
   } 
  	 
  score = overlapAlign(truncRef, sequence, mode);
  regscore = score/length(truncRef);
  
  return std::pair<double, double>(score, regscore);
}

// Compute edit distance + longest ins + longuest del;
std::tuple<double, int, int> InnerEditAlignment(TSeq& reference, TSeq& pathSeq)
{
  double score;
  int longuestINS(0);
  int longuestDEL(0);
  int lengthCount(0);
  bool stop(true);
  
  int len1 = (int)length(reference);
  int len2 = (int)length(pathSeq);
  
  TAlign align;
  resize(rows(align), 2);
  
  if( len1>len2 )
   {
    assignSource(row(align, 0), reference);
    assignSource(row(align, 1), pathSeq);
   }
  else 
   {
    assignSource(row(align, 0), pathSeq);
    assignSource(row(align, 1), reference);
   }
  
  score = (double)globalAlignment(align, Score<int, Simple>(0, -1, -1));

  
  TRow & row1 = row(align, 0);
  TRow & row2 = row(align, 1);
  TRowIterator it = begin(row1), itEnd = end(row1);

//LONGUEST INSERTION	
	stop = false;
    it = begin(row1), itEnd = end(row1);
    while(it != itEnd)
    {
      if(isGap(it)) ++lengthCount;
      else
       {
		if(longuestINS<lengthCount) longuestINS = lengthCount;
		++lengthCount = 0;
	   }
      ++it;
     }
   
//LONGUEST DELETION
	stop = false;
	lengthCount = 0;
    it = begin(row2), itEnd = end(row2);
    while(it != itEnd)
    {
      if(isGap(it)) ++lengthCount;
      else
       {
		if(longuestDEL<lengthCount) longuestDEL = lengthCount; 
		lengthCount = 0;
	   }
      ++it;
     }
// std::cout << "DEL LENGTH: " << longuestDEL << std::endl;    
    
   if( len1>len2 ) return std::make_tuple(score, longuestINS, longuestDEL);
   else return std::make_tuple(score, longuestDEL, longuestINS);
}

//OUTPUT edit distance nbIns nbDel in pathSeq
std::tuple<double, int, int> BorderEditAlignment(TSeq& reference, TSeq& pathSeq, Direction& direction)
{
//std::cout << "BORDER EDIT ALIGNMENT" << std::endl;	
  double score(-50000);
  int borderINS(0);
  int borderDEL(0);
  int gapCount(0);
  
  TAlign align;
  resize(rows(align), 2);
  
  int len1 = (int)length(reference);
  int len2 = (int)length(pathSeq);
  bool stop(false);
  
  if( (len1>0) & (len2>0) )
  {
    if( len1>len2 )
   {
    assignSource(row(align, 0), reference);
    assignSource(row(align, 1), pathSeq);
   }
  else 
   {
    assignSource(row(align, 0), pathSeq);
    assignSource(row(align, 1), reference);
   }
  
  score = (double)globalAlignment(align, Score<int, Simple>(0, -1, -1), LinearGaps());

  
  //Insertion/Deletion length computation
  TRow & row1 = row(align, 0);
  TRow & row2 = row(align, 1);
  
  TRowIterator it = begin(row1), itEnd = end(row1);
  if( direction == RIGHT )
   {
    while( (itEnd!=it) & !stop)
    {
      if(isGap(itEnd)) ++gapCount;
      else stop = true;
      --itEnd;             
    }
   }
  else
   {
	while( (it!=itEnd) & !stop)
    {
      if(isGap(it)) ++gapCount;
      else stop = true;
      ++it;             
    }
   }
  borderINS = gapCount;
  
  it = begin(row2), itEnd = end(row2);
  if( direction == RIGHT )
   {
    while( (itEnd!=it) & !stop)
    {
      if(isGap(itEnd)) ++gapCount;
      else stop = true;
      --itEnd;             
    }
   }
  else
   {
	while( (it!=itEnd) & !stop)
    {
      if(isGap(it)) ++gapCount;
      else stop = true;
      ++it;             
    }
   }
  borderDEL = gapCount;

 if( len1 > len2 ) return std::make_tuple(score, borderINS, borderDEL);
 else return std::make_tuple(score, borderDEL, borderINS);
}
 else
  {
   std::cout << "PB in BORDER EDIT ALIGNMENT" << std::endl;
   return std::make_tuple(score, 100, 100);
  }
	
}

/***********************************************************************
 * SEPTEMBRE 
 * Fonctions pour améliorer l'étape de définition des régions SOLID/WEAK
 ***********************************************************************/

TSeq getKmerAt(TSeq& refSequence, unsigned int position, unsigned int& kmerSize)
{
	TSeq answer = infix(refSequence, position, position + kmerSize);
	return answer;
}

//Return the part of the reference sequence that contains its kmers #start_kmPosition to #end_kmPosition (both included) 
TSeq extractSolidSequence(TSeq& refSequence, unsigned int start_kmPosition, unsigned int end_kmPosition, unsigned int& kmerSize)
{
	TSeq answer = infix(refSequence, start_kmPosition, end_kmPosition + kmerSize);
	return answer;
}

TSeq extractWeakSequence(TSeq& refSequence, unsigned int end_kmPosition, unsigned int start_kmPosition, unsigned int& kmerSize)
{
	TSeq answer = infix(refSequence, end_kmPosition+kmerSize, start_kmPosition);
	return answer;
}

TSeq extractWeakBorderSequence(TSeq& refSequence, unsigned int kmPosition, unsigned int& kmerSize, Location m_location)
{
	TSeq answer;
	if(m_location == HEAD) answer = prefix(refSequence, kmPosition);
	else answer = suffix(refSequence, kmPosition+kmerSize);
	return answer;
}

double absdistance(TSeq& seq1, TSeq& seq2)
{
  return (double)(abs((int)length(seq1)-(int)length(seq2)));
}

//Uncolour repeat kmers
void decolourRepeatsFromDBG(colouredDBG& dBG, unsigned int& kmerSize)
{
	TSeq kmerToRemove;	
	for(unsigned int b(0); b<Dict.size(); b++)
	 {
	  Dna5 base(Dict[b]);
	  clear(kmerToRemove);
	  for(unsigned int n(0); n<kmerSize; n++) appendValue(kmerToRemove, base);
	  std::cout << kmerToRemove << " " <<  dBG.count(kmerToRemove) << std::endl;
	  if( dBG.count(kmerToRemove)>0 ) std::get<1>(dBG[kmerToRemove]) = 0;
	}   
}
