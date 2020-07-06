//
//  Trail.cpp

//  Created by Lucile Broseus.
//  Copyright © 2019 LuB. All rights reserved.
//
/****************************************************************************
 * Version TALC 1.01
 ****************************************************************************
 ****************************************************************************
 * Updates:
 * 20/08/2019: remove Repeat warning
 ****************************************************************************/

#include "Trail.hpp"

#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>

#include "utils.hpp"
#include "Jellyfish.hpp"

using namespace seqan;
  
typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<Dna5, colouredCount > colouredBase; 
typedef std::pair<TSeq, colouredCount > tipNode; 

extern bool DEBUG_TEST;
extern bool DEBUG_USER;
extern unsigned int K;

Trail::~Trail()
{
	
}


Trail::Trail() : m_sequence(make_empty_Dna5String()), 
                 m_lastStep(make_empty_tipNode()),    
                 m_lastScore(0), m_nbFailuresInARow(0),
                 m_nbBreakpoints(0), m_distance(0), m_leftAnchor(-1), m_rightAnchor(-1)
{
    
}

//Constructs inner nodes
Trail::Trail(TSeq kmer, colouredCount& ccount) : m_sequence(kmer),
                                          m_lastStep(std::make_pair(kmer, ccount)), m_lastScore(0), 
                                          m_nbFailuresInARow(0),
                                          m_nbBreakpoints(0), m_distance(0),
                                          m_leftAnchor(-1), m_rightAnchor(-1)
                                          
{
    
}

Trail::Trail(Trail path, TSeq seq) : m_sequence(seq),
                                     m_lastStep(make_empty_tipNode()),  m_lastScore(path.getLastScore()),
                                     m_nbFailuresInARow(path.getNbFailuresInARow()),
                                     m_nbBreakpoints(path.getNbBreakpoints()), m_distance(path.getDistance()),
                                     m_leftAnchor(path.getLeftAnchor()), m_rightAnchor(path.getRightAnchor())
{
    
}

Trail::Trail(Trail path, Dna5 newBase, Direction& direction, unsigned int count) : m_sequence(addNewBase(path.getSeq(), newBase, direction)),
                                         m_lastStep(makeTipNode(path.getLastKmer(), newBase, direction, count)),  
                                         m_lastScore(path.getLastScore()),
                                         m_nbFailuresInARow(path.getNbFailuresInARow()),
                                         m_nbBreakpoints(path.getNbBreakpoints()), m_distance(path.getDistance()),
                                         m_leftAnchor(path.getLeftAnchor()), m_rightAnchor(path.getRightAnchor())
{
    
}

//Methods

TSeq Trail::getSeq()
{
 return m_sequence;
}

void Trail::setSeq(TSeq seq)
{
 m_sequence = seq;
}

unsigned int Trail::getLength()
{
 return length(m_sequence);
}

tipNode Trail::getLastStep()
{
 return m_lastStep;
}

void Trail::setLastStep(tipNode lastStep)
{
 m_lastStep = lastStep;
}

unsigned int Trail::getLastCount()
{
 return std::get<0>(std::get<1>(m_lastStep)); 
}

TSeq Trail::getLastKmer()
{
 return std::get<0>(m_lastStep); 
}
    
double Trail::getLastScore()
{
 return m_lastScore;
}

void Trail::setLastScore(double& score)
{
 m_lastScore = score;
}
//Compute standard alignment score of path against reference sequence.
// length(reference) should be larger than length(path)
void Trail::DPscore(TSeq& reference)
{
	Align<TSeq > align;
    resize(rows(align), 2);
    
    assignSource(row(align, 0), reference);
    assignSource(row(align, 1), this->getSeq());

	m_lastScore = localAlignment(align, Score<int, Simple>(1, 0, 0), LinearGaps());;
}

void Trail::Overlapscore(TSeq& reference, Direction& direction)
{
    TStringSet sequences;
    TSeq candidate(this->getSeq());

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
     
   switch( direction )
    {
     case LEFT:
     //Correction from left to right, do not penalize gap on the right                            
     m_lastScore = globalAlignment(alignG, Score<int, Simple>(4, -3, -2),
                                AlignConfig<false, false, true, true>(), LinearGaps());           
     break;
      //Correction from right to left, do not penalize gap on the left                            
     case RIGHT:
     m_lastScore = globalAlignment(alignG, Score<int, Simple>(4, -3, -2),
                                AlignConfig<true, true, false, false>(), LinearGaps());
     }
}

unsigned int Trail::getNbFailuresInARow()
{
 return m_nbFailuresInARow;
}

void Trail::recordFailure()
{
 m_nbFailuresInARow++;
}

void Trail::eraseFailures()
{
 m_nbFailuresInARow = 0;
}


//reference is assumed to be anchored
bool Trail::seedAndExtend(TSeq& reference, Direction& direction, int xdrop, unsigned int MAX_FAILURES)
{
 TSeq candidate(this->getSeq());
 bool ok(true);
 
 std::tuple<TSeq, TSeq, int, double, bool > extension_res;
 	   
 extension_res = getSeedAndExtension(reference, candidate, xdrop, direction, K) ;
 
 //check whether total # of errors allowed is overcome
 ok = ( length(std::get<1>(extension_res))==length(candidate) );
 if(!ok) this->recordFailure();
 else this->eraseFailures();
 
 //#of errors compared to the reference, score computed as ..?
 this->setLastScore( std::get<3>( extension_res ) );
 
 if( direction == RIGHT ) this->setRightAnchor( std::get<2>( extension_res ) );
 else this->setLeftAnchor( std::get<2>( extension_res ) );
 
 ok = (this->getNbFailuresInARow()<=MAX_FAILURES);
 ok &= !std::get<4>( extension_res );
 return ok;
}

int Trail::getNbBreakpoints()
{
	return m_nbBreakpoints;
}

void Trail::setNbBreakpoints(int nb)
{
	m_nbBreakpoints = nb;
}

void Trail::recordBreakpoint()
{
	m_nbBreakpoints++;
}

double Trail::getDistance()
{
 return m_distance; 
}

void Trail::recordDistance(double distance)
{
 m_distance += distance;
}

int Trail::getRightAnchor()
{
  return m_rightAnchor;
}

void Trail::setRightAnchor(int position)
{
  m_rightAnchor = position;
}

int Trail::getLeftAnchor()
{
  return m_leftAnchor;
}

void Trail::setLeftAnchor(int position)
{
  m_leftAnchor = position;
}

void Trail::setReachedAim(int aim, Direction& direction)
{
  if(direction == RIGHT) m_rightAnchor = aim;
  else m_leftAnchor = aim;
}

//INPUT: list of the kmers chosen as targets
//Check whether last kmer matches one of the aim-kmers
//If so, records the hit kmer into the Trail object
//OUTPUT: returns whether an aim has been reached
bool Trail::checkAims(std::vector<anchorTuple > aims, Direction& direction)
{
 bool isEqual(false);
 unsigned int i(0);
 while( (!isEqual) & (i<aims.size()) )
  {
   isEqual=(std::get<0>(m_lastStep)==std::get<0>(aims[i]));
   ++i;
  }
 if(isEqual) this->setReachedAim( std::get<1>(aims[i-1]), direction );

 return isEqual;
}



bool Trail::ThinkIveAlreadyGotThere(TSeq history)
{
 int alreadyOccurred = -1;
 
 if(length(history)>length(std::get<0>(m_lastStep)))
  {
   Finder<TSeq> finder(history);
   Pattern<CharString, Horspool> pattern(std::get<0>(m_lastStep));
        
   while (find(finder, pattern) & (alreadyOccurred<0)) alreadyOccurred = beginPosition(finder);
  }
  //if(alreadyOccurred>0) std::cout << "REPEATED KMER: " << alreadyOccurred << std::endl;  
    return alreadyOccurred>0;
}


std::vector<colouredCount > Trail::whatsNext(Direction& direction, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable)
{
//std::cout << "[NEXT STEP]: " << std::get<0>(m_lastStep) << std::endl;
 return getNextCounts(std::get<0>(m_lastStep), direction, dBG, countsTable, junctionCountsTable);
}

//---------------------------------------------------------------------//

TSeq addNewBase(TSeq sequence, Dna5 newBase, Direction& direction)
{
 TSeq tmp(newBase);
  
 switch( direction )
  {
   case LEFT: 
	   append(tmp, sequence);
	   sequence = tmp;
   break;   
   case RIGHT:
       tmp = sequence;
	   appendValue(sequence, newBase);
   break;
  }  
//std::cout << "[ADD NEW BASE] " << newBase << " -> " << sequence << std::endl;
 return sequence;
}

tipNode makeTipNode(TSeq kmer, Dna5 base, Direction& direction, unsigned int count)
{
	TSeq next_kmer = formNextKmer(kmer, base, direction);
	tipNode result = std::make_pair(next_kmer, std::make_pair(count,0));
//std::cout << "[TIP NODE] " << next_kmer << std::endl;
	return result;
}


std::tuple<TSeq, TSeq, int, double, bool> getSeedAndExtension(TSeq& reference,TSeq& candidate, int& xdrop, Direction& direction, unsigned int& seedSize)
{
  TSeq seq1;
  TSeq seq2;
  
  TSeq refExtension;
  TSeq histExtension;
  int posOnRef(-1);                //should be positive when +
  bool state(true);                //if history is shorter than reference or not
  Align<TSeq > align;
  Score<int, Simple > scoringScheme(0, -1, -1);
  bool stopThere(false);

  double score(0);
  
  int limit(0);
//std::cout << "[DEBUG] Ref:  " << length(reference) << " candidate: " <<  length(candidate) << " xdrop: " << xdrop << "\n";
  if( length(reference) < length(candidate) )
   {
	seq1 = candidate;
	seq2 = reference;
	state = false;
   }
  else
   {
	seq1 = reference;
	seq2 = candidate;  
   }
 
  if( direction == RIGHT )
   {
    Seed<Simple> seedR(0, 0, seedSize-1, seedSize-1);
    extendSeed(seedR, seq1, seq2, EXTEND_RIGHT, scoringScheme, xdrop, GappedXDrop());
    if(endPositionV(seedR)>length(seq2)) std::cout << "BUG rightext: endPos "  << endPositionV(seedR) << " seq2: " << length(seq2) << std::endl;
    if( state ) 
     {
		 histExtension = prefix(candidate, endPositionV(seedR)); 
		 refExtension = prefix(reference, endPositionH(seedR)); 
		 posOnRef = endPositionH(seedR);		 
	 }	 
    else 
     {
      histExtension = prefix(candidate, endPositionH(seedR));
      refExtension = prefix(reference, endPositionV(seedR)); 
      posOnRef = endPositionV(seedR);	       
     }
   }
  else
   {
    Seed<Simple> seedL(length(seq1)-seedSize, length(seq2)-seedSize, length(seq1)-1, length(seq2)-1);
    extendSeed(seedL, seq1, seq2, EXTEND_LEFT, scoringScheme, xdrop, GappedXDrop());
    if(beginPositionV(seedL)>=length(seq2)) std::cout << "BUG leftext" << std::endl;
    if (state ) 
     {
	  histExtension = suffix(candidate, beginPositionV(seedL));
	  refExtension = suffix(reference, beginPositionH(seedL));
	  posOnRef = beginPositionH(seedL);
     }
    else 
     {
	  histExtension = suffix(candidate, beginPositionH(seedL));
	  refExtension = suffix(reference, beginPositionV(seedL));
	  posOnRef = beginPositionV(seedL);
     }
   }
   
//--------------------------------------------------------------------->DEBUG
 if( std::max(length(refExtension), length(histExtension))>=seedSize )
   {
    resize(rows(align), 2);
    if( length(histExtension) <= length(refExtension) )
     {
	  assignSource(row(align, 0), refExtension);
      assignSource(row(align, 1), histExtension);   
     }
    else
     {
	  assignSource(row(align, 0), histExtension);
      assignSource(row(align, 1), refExtension);     
     }  
//std::cout << "[DEBUG] Ref:  " << length(refExtension) << " candidate: " <<  length(histExtension) << "\n";      
   score = globalAlignment(align,  scoringScheme);     

   //std::cout << "[DEBUG] Alignment score: " << score << "\n";
    //std::cout << "xdrop: " << xdrop << "\n";
    //std::cout << "history's length: " << length(history) << "\n";
   //}
  } 
 else
  {
   std::cout << "[DEBUG] BUG IN ALIGNMENT " << length(reference) << " candidate: " <<  length(candidate) << "\n";  
   score = (-1)*xdrop;  
   stopThere = true;
  }
//--------------------------------------------------------------------->DEBUG 
  return std::make_tuple(refExtension, histExtension, posOnRef, score, stopThere);
}

