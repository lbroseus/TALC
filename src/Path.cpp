//
//  Path.cpp

//  Created by Lucile Broseus.
//  Copyright © 2018 LuB. All rights reserved.
//

#include "Path.hpp"

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

using namespace seqan;
  
typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<Dna5, colouredCount > colouredBase; 

typedef StringSet<TSeq> TStringSet;                    // container for strings
typedef StringSet<TSeq, Dependent<> > TDepStringSet;   // dependent string set
typedef Graph<Alignment<TDepStringSet> > TAlignGraph;       // alignment graph

extern bool DEBUG_align;
extern bool DEBUG_info;
extern bool DEBUG_correction;

Path::~Path()
{
}

//Constructors

Path::Path() : m_sequence(make_empty_Dna5String()), 
               m_lastStep(make_empty_Dna5String()),    
               m_lastScore(0), m_lastHitLandmark(0),       
               m_ExpectedWays(make_empty_vecOfColouredBases()), 
               m_Breakpoints(make_empty_vecOfColouredBases()),
               m_UnexpectedWays(make_empty_vecOfColouredBases()), 
               m_nbVisits(0), m_nbFailures(0), m_nbBreakpoints(0), m_nbUnpredBranch(0), m_start(-1), m_reachedAim(-1)
{
    
}

//Constructs first nodes
Path::Path(int whichStart) : m_sequence(make_empty_Dna5String()),
                                       m_lastStep(make_empty_Dna5String()), m_lastScore(0), m_lastHitLandmark(0),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), m_nbFailures(0), m_xdropenalty(2), 
                                       m_nbBreakpoints(0), m_nbUnpredBranch(0),
                                       m_start(whichStart), m_reachedAim(-1),
                                       m_status(EXPECTED)
{
    
}

//Constructs inner nodes
Path::Path(TSeq seq, int whichStart) : m_sequence(seq),
                                       m_lastStep(make_empty_Dna5String()), m_lastScore(0), m_lastHitLandmark(0),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), 
                                       m_nbFailures(0),  m_nbFailuresInARow(0), m_xdropenalty(2),  
                                       m_nbBreakpoints(0), m_nbUnpredBranch(0),
                                       m_start(whichStart), m_reachedAim(-1),
                                       m_status(EXPECTED)
{
    
}

Path::Path(TSeq seq, int whichStart, Status status) : m_sequence(seq),
                                       m_lastStep(make_empty_Dna5String()), m_lastScore(0), m_lastHitLandmark(0),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), 
                                       m_nbFailures(0), m_nbFailuresInARow(0), m_xdropenalty(2), 
                                       m_nbBreakpoints(0), m_nbUnpredBranch(0),
                                       m_start(whichStart), m_reachedAim(-1),
                                       m_status(status)
{
    
}

//Constructs last node when an aim has been reached
Path::Path(TSeq seq, int whichStart, int whichAim) : m_sequence(seq),
                                       m_lastStep(make_empty_Dna5String()), m_lastScore(0), m_lastHitLandmark(0),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), 
                                       m_nbFailures(0), m_nbFailuresInARow(0), m_xdropenalty(2),  
                                       m_nbBreakpoints(0), m_nbUnpredBranch(0),
                                       m_start(whichStart), m_reachedAim(whichAim),
                                       m_status(EXPECTED)
{
    
}

//Constructs new node heriting from previous path
Path::Path(Path path, TSeq seq) : m_sequence(seq),
                                       m_lastStep(make_empty_Dna5String()),  m_lastScore(path.getLastScore()), m_lastHitLandmark(path.getLastHitLandmark()),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), m_nbFailures(path.getNbFailures()), m_nbFailuresInARow(path.getNbFailuresInARow()),
                                       m_xdropenalty(path.getXdropenalty()), 
                                       m_nbBreakpoints(path.getNbBreakpoints()), m_nbUnpredBranch(path.getNbUnpred()),
                                       m_start(path.getStart()), m_reachedAim(-1),
                                       m_status(path.getStatus())
{
    
}

Path::Path(Path path, TSeq seq, Status newStatus) : m_sequence(seq),
                                       m_lastStep(make_empty_Dna5String()), m_lastScore(path.getLastScore()), m_lastHitLandmark(path.getLastHitLandmark()),
                                       m_ExpectedWays(make_empty_vecOfColouredBases()), 
                                       m_Breakpoints(make_empty_vecOfColouredBases()),
                                       m_UnexpectedWays(make_empty_vecOfColouredBases()), 
                                       m_nbVisits(0), m_nbFailures(path.getNbFailures()), m_nbFailuresInARow(path.getNbFailuresInARow()), 
                                       m_xdropenalty(path.getXdropenalty()), 
                                       m_nbBreakpoints(path.getNbBreakpoints()), m_nbUnpredBranch(path.getNbUnpred()),
                                       m_start(path.getStart()), m_reachedAim(-1),
                                       m_status(newStatus)
{
    
}


//Methods

TSeq Path::getSeq()
{
    return m_sequence;
}

void Path::setSeq(TSeq seq)
{
    m_sequence = seq;
}

int Path::getLength()
{
	return length(m_sequence);
}

void Path::addBase(seqan::Dna5 base)
{
    appendValue(m_sequence, base);
}

TSeq Path::getLastStep()
{
    return m_lastStep;
}

void Path::setLastStep(TSeq step)
{
    m_lastStep = step;
}

double Path::getLastScore()
{
    return m_lastScore;
}

void Path::setLastScore(double score)
{
    m_lastScore = score;
}

unsigned int Path::getLastHitLandmark()
{
  return m_lastHitLandmark;
}

void Path::setLastHitLandmark(unsigned int nb)
{
	m_lastHitLandmark = nb;
}

std::vector<colouredBase > Path::getExpectedWays()
{
	return m_ExpectedWays;
}
std::vector<colouredBase > Path::getBreakpoints()
{
	return m_Breakpoints;
}
std::vector<colouredBase > Path::getUnexpectedWays()
{
	return m_UnexpectedWays;
}
    
void Path::setExpectedWays(std::vector<colouredBase > ExpectedWays)
{
	m_ExpectedWays = ExpectedWays;
}
void Path::setBreakpoints(std::vector<colouredBase > Breakpoints)
{
	m_Breakpoints = Breakpoints;
}
void Path::setUnpredBases(std::vector<colouredBase > UnexpectedWays)
{
	m_UnexpectedWays = UnexpectedWays;
}

void Path::pruneUnpredBases()
{
	clear( m_UnexpectedWays );
}
           
unsigned int Path::getNbVisits()
{
    return m_nbVisits;
}

void Path::setNbVisits(unsigned int nbVisits)
{
    m_nbVisits = nbVisits;
}

void Path::recordVisit()
{
    m_nbVisits++;
}

void Path::update(seqan::Dna5 next_base, Status status, Direction direction)
{
  TSeq tmp(next_base);
  
	switch( direction )
	 {
	   case LEFT: 
	     append(tmp, m_sequence);
	     m_sequence = tmp;
	   break;   
	   case RIGHT:
	      append(m_sequence, tmp);
	   break;
    }   
       m_status = status;
}

Status Path::getStatus()
{
    return m_status;
}

void Path::setStatus(Status status)
{
    m_status = status;
}

int Path::getStart()
{
    return m_start;
}

void Path::setStart(int whichStart)
{
    m_start = whichStart;
}

int Path::getNbFailures()
{
	return m_nbFailures;
}

void Path::recordFailure()
{
	m_nbFailures++;
	m_nbFailuresInARow++;
}

void Path::eraseFailure(int mode)
{
	if(mode == 0) m_nbFailures = std::max(0,m_nbFailures-1);
	else if(mode = 1) m_nbFailuresInARow = 0;
	else
	{
	 m_nbFailures = std::max(0,m_nbFailures-1);
	 m_nbFailuresInARow = 0;
	}
}

int Path::getNbFailuresInARow()
{
	return m_nbFailuresInARow;
}

int Path::getXdropenalty()
{
	return m_xdropenalty;
}

void Path::setXdrop(int& xdrop)
{
	m_xdropenalty = xdrop;
}

void Path::updateXdropenalty(int pen)
{
	m_xdropenalty += pen;
}

int Path::getNbBreakpoints()
{
	return m_nbBreakpoints;
}

void Path::recordBreakpoint()
{
	m_nbBreakpoints++;
}

int Path::getNbUnpred()
{
	return m_nbUnpredBranch;
}

void Path::recordUnpred()
{
	m_nbUnpredBranch++;
}

int Path::getReachedAim()
{
    return m_reachedAim;
}

void Path::setReachedAim(int aim)
{
    m_reachedAim = aim;
}

/***********************************************************************
 * seemsOK()
 * Takes ANCHORED history sequence and ANCHORED weak sequence
 * ANCHOR is the seed
 * Seed and extend both sequence from the anchor
 * xdrop parameter is function of history's seq length
 * 
 ***********************************************************************/
 
//Here seq2 should be shorter than seq1
//Perform seed&enxtend between seq1 and seq2 with allowed #errors=xdrop>=0
std::tuple<TSeq, TSeq, int > performExtension(TSeq& seq1,TSeq& seq2, int xdrop, Direction& direction, unsigned int& seedSize)
{

  TSeq refExtension;
  TSeq histExtension;
  int posOnRef(-1);                //should be positive when +
  bool state(true);                //if history is shorter than reference or not
  Score<int, Simple > scoringScheme1(0, -1, -1);

  int limit(0);
 
  if( direction == RIGHT )
   {
    Seed<Simple> seedR(0, 0, seedSize-1, seedSize-1);
    extendSeed(seedR, seq1, seq2, EXTEND_RIGHT, scoringScheme1, xdrop, GappedXDrop());
    if(endPositionV(seedR)>length(seq2)) std::cout << "BUG rightext: endPos "  << endPositionV(seedR) << " seq2: " << length(seq2) << std::endl;
    histExtension = prefix(seq2, endPositionV(seedR)); 
    refExtension = prefix(seq1, endPositionH(seedR)); 
	posOnRef = endPositionH(seedR);		 ;
   }
  else
   {
    Seed<Simple> seedL(length(seq1)-seedSize, length(seq2)-seedSize, length(seq1)-1, length(seq2)-1);
    extendSeed(seedL, seq1, seq2, EXTEND_LEFT, scoringScheme1, xdrop, GappedXDrop());
    if(beginPositionV(seedL)>=length(seq2)) std::cout << "BUG leftext: " << beginPositionV(seedL) << " vs. " << length(seq2) << std::endl;
	histExtension = suffix(seq2, beginPositionV(seedL));
	refExtension = suffix(seq1, beginPositionH(seedL));
	posOnRef = beginPositionH(seedL);   
   }
 
 return std::make_tuple(refExtension, histExtension, posOnRef);

}

 double performGlobalAlignment(TSeq seq1, TSeq seq2)
{ 
 Score<int, Simple > scoringScheme(0, -1, -1);
 Align<TSeq > align;
 double score;
 
 resize(rows(align), 2);
    if( length(seq2) <= length(seq1) )
     {
	  assignSource(row(align, 0), seq1);
      assignSource(row(align, 1), seq2);   
     }
    else
     {
	  assignSource(row(align, 0), seq2);
      assignSource(row(align, 1), seq1); 
     }      
  score = globalAlignment(align, scoringScheme);
if(DEBUG_align) std::cout << "DEBUG| SHORTER BORDER ALIGNMENT" << align << std::endl;
  return score;
}

double performSemiGlobalAlignment(TSeq seq1, TSeq seq2, Direction direction)
{ 
 Score<int, Simple > scoringScheme(0, -1, -1);
 double score;
 TStringSet sequences;
 
 appendValue(sequences, seq1);
 appendValue(sequences, seq2);
 TAlignGraph alignG(sequences);

 if(direction == RIGHT) score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>());
 else  if(direction == LEFT) score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<true, true, false, false>());

//std::cout << "DEBUG| SHORTER PATH SEMI-GLOBAL ALIGNMENT" << score << std::endl;
  return score;
}

double performLocalAlignment(TSeq seq1, TSeq seq2)
{ 
 Score<int, Simple > scoringScheme(0, -1, -1);
 Align<TSeq > align;
 double score;
 
 resize(rows(align), 2);
    if( length(seq2) <= length(seq1) )
     {
	  assignSource(row(align, 0), seq1);
      assignSource(row(align, 1), seq2);   
     }
    else
     {
	  assignSource(row(align, 0), seq2);
      assignSource(row(align, 1), seq1); 
     }      
  score = localAlignment(align, scoringScheme);
//if(DEBUG_align) std::cout << "DEBUG| SHORTER INNER ALIGNMENT" << align << std::endl;
  return score;
}
 
double computePercentID(TSeq seq1, TSeq seq2)
{
 Score<int, Simple > scoringScheme(1, 0, 0);
 Align<TSeq > align;
 double score;
 double len;
 
 resize(rows(align), 2);
    if( length(seq2) <= length(seq1) )
     {
	  assignSource(row(align, 0), seq1);
      assignSource(row(align, 1), seq2); 
      len = (double)length(seq1);  
     }
    else
     {
	  assignSource(row(align, 0), seq2);
      assignSource(row(align, 1), seq1);
      len = (double)length(seq2);   
     }      
  score = localAlignment(align, scoringScheme)/len;
  return score;
}

std::tuple<TSeq, TSeq, int > findWhereItEnds(TSeq& extGapSeq, TSeq& shortDEADENDpath, int xdrop1, Direction& direction, unsigned int& kmerSize)
{
  int j(0);
  bool goFurther(true);
  double score(0);
  std::tuple<TSeq, TSeq, int > extension_results;
  std::tuple<TSeq, TSeq, int > new_extension_results;
  
  new_extension_results = performExtension(extGapSeq, shortDEADENDpath, xdrop1, direction, kmerSize);

  do
   {
	 ++j;
	 extension_results = new_extension_results;
	 new_extension_results = performExtension(extGapSeq, shortDEADENDpath, xdrop1-j, direction, kmerSize);
	 if( length(std::get<1>(new_extension_results))<length(std::get<1>(extension_results)) ) goFurther = false;
//std::cout << "DEBUG| new pos " << std::get<2>(extension_results)  << std::endl;
   }while( goFurther & (xdrop1>=j) );

std::cout << "DEBUG| INIT xdrop " << xdrop1 << " NOW: " << xdrop1-j+1 << std::endl;
  return  std::make_tuple( std::get<0>(extension_results), std::get<1>(extension_results), std::get<2>(extension_results));	    
}

 //extGapSeq: weak region + anchor
std::tuple<TSeq, TSeq, int, double> getExtension(TSeq& extGapSeq,TSeq& history, int& xdrop, Direction& direction, unsigned int seedSize, int mode)
{
  TSeq seq1;
  TSeq seq2;
  
  TSeq refExtension;
  TSeq histExtension;
  int posOnRef(-1);                //should be positive when +
  bool state(true);                //if history is shorter than reference or not
  Align<TSeq > align;
  Score<int, Simple > scoringScheme(0, -1, -1);

  Seed<Simple> seed;
  double score(-200000);
  
  int limit(0);
 
  if( length(extGapSeq) < length(history) )
   {
	seq1 = history;
	seq2 = extGapSeq;
	state = false;
   }
  else
   {
	seq1 = extGapSeq;
	seq2 = history;  
   }
 
  if( direction == RIGHT )
   {
    Seed<Simple> seedR(0, 0, seedSize-1, seedSize-1);
    extendSeed(seedR, seq1, seq2, EXTEND_RIGHT, scoringScheme, xdrop, GappedXDrop());
    if(endPositionV(seedR)>length(seq2)) std::cout << "BUG rightext: endPos "  << endPositionV(seedR) << " seq2: " << length(seq2) << std::endl;
    if( state ) 
     {
		 histExtension = prefix(history, endPositionV(seedR)); 
		 refExtension = prefix(extGapSeq, endPositionH(seedR)); 
		 posOnRef = endPositionH(seedR);		 
	 }	 
    else 
     {
      histExtension = prefix(history, endPositionH(seedR));
      refExtension = prefix(extGapSeq, endPositionV(seedR)); 
      posOnRef = endPositionV(seedR);	       
     }
    seed = seedR;
   }
  else
   {
    Seed<Simple> seedL(length(seq1)-seedSize, length(seq2)-seedSize, length(seq1)-1, length(seq2)-1);
    extendSeed(seedL, seq1, seq2, EXTEND_LEFT, scoringScheme, xdrop, GappedXDrop());
    if(beginPositionV(seedL)>=length(seq2)) std::cout << "BUG leftext" << std::endl;
    if (state ) 
     {
	  histExtension = suffix(history, beginPositionV(seedL));
	  refExtension = suffix(extGapSeq, beginPositionH(seedL));
	  posOnRef = beginPositionH(seedL);
     }
    else 
     {
	  histExtension = suffix(history, beginPositionH(seedL));
	  refExtension = suffix(extGapSeq, beginPositionV(seedL));
	  posOnRef = beginPositionV(seedL);
     }
    seed = seedL;
   }
   
//--------------------------------------------------------------------->DEBUG
	   
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
   score = globalAlignment(align,  scoringScheme);     
    //std::cout << "Resulting extendSeed alignment\n" << align << "\n";
    //std::cout << "Alignment score: " << score << "\n";
    //std::cout << "xdrop: " << xdrop << "\n";
    //std::cout << "history's length: " << length(history) << "\n";
   //} 
//--------------------------------------------------------------------->DEBUG 
  return std::make_tuple(refExtension, histExtension, posOnRef, score);
}

bool Path::seemsOK(TSeq& weakSeq, TSeq& history, TSeq& currentAnchor, Direction& direction, unsigned int& kappa, double& errorRate)
{ 
//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION seemsOK 1" << std::endl;
//--------------------------------------------------------------------->DEBUG	
    bool ok(true);
    int trimLength(0);
    std::tuple<TSeq, TSeq, int, double > extension_res;
    TSeq extendGap;
    Align<TSeq > align;
    Score<int, Simple > scoringScheme(0, -1, -1);
    Seed<Simple> seed;
    int xdrop(0);
    
    xdrop = std::max(1, (this->getXdropenalty())) + 2;
//std::cout << "[PATH] :" << "XDROP: " << xdrop << std::endl;   

    
    if(direction == RIGHT)
     {
      append(extendGap, currentAnchor);
      append(extendGap, weakSeq);

     }
    else
     {
	  append(extendGap, weakSeq); 
	  append(extendGap, currentAnchor);
	 } 
	   
    extension_res = getExtension(extendGap, history, xdrop, direction, length(currentAnchor), 0) ;
    ok = ( length(std::get<1>(extension_res))==length(history) );
    this->setLastScore( std::get<3>( extension_res ) );
    //Record the number of errors between history and its extension on extendGap
    xdrop = (int)((-1)*performSemiGlobalAlignment(extendGap, history, direction));
    this->setXdrop( xdrop );

//std::cout << "[PATH] :" << "New score: " << std::get<3>( extension_res ) << std::endl;   
//std::cout << "[PATH] :" << "New xdrop: " << xdrop << std::endl;   

//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION seemsOK 2" << std::endl;
// if(!ok) std::cout << "Diagnostic: " << ok << std::endl;
//--------------------------------------------------------------------->DEBUG	
    return (ok);
}

bool Path::checkInnerPath(unsigned int kmerSize)
{
	bool diagnostic;
	diagnostic = (m_start>=0);
	if(!diagnostic) std::cout << "Inner path lacks a specified start kmer!" << std::endl;
	else 
	{ 
	  diagnostic &= (!empty(m_ExpectedWays) || !empty(m_Breakpoints) || !empty(m_UnexpectedWays));
	  if(!diagnostic) std::cout << "Inner path has no specified continuation!" << std::endl;
	  else 
	  {
		diagnostic &= (length(m_lastStep)==kmerSize);
	    if(!diagnostic) std::cout << "Last step is not a kmer! -> length " << length(m_lastStep) << std::endl;   
	  }
    }
  return diagnostic;
}

bool Path::checkInitPath(unsigned int kmerSize)
{
	bool diagnostic;
	diagnostic = (m_start>=0);
	if(!diagnostic) std::cout << "Inner path lacks a specified start kmer!" << std::endl;
	else 
	{ 
	  diagnostic &= (length(m_lastStep)==kmerSize);
	  if(!diagnostic) std::cout << "Last step is not a kmer! -> length " << length(m_lastStep) << std::endl;	    
    }
  return diagnostic;
}

bool Path::checkNewBranch()
{
	bool diagnostic;
	diagnostic = (m_start>=0);
	if(!diagnostic) std::cout << "Inner path lacks a specified start kmer!" << std::endl;
	else 
	{ 
	  diagnostic &= (length(m_sequence)==1);
	  if(!diagnostic) std::cout << "Pb in seq. initiation, its length is " << length(m_sequence) << std::endl;	    
    }
  return diagnostic;
}

void Path::displayPath()
{
	std::cout << "******************* Path info ************************" << "\n"
	          << "Whole sequence: " << m_sequence << "\n"
	          << "Starting from position " << m_start << "\n"
	          << "Reached aim (?) " << m_reachedAim << "\n"
	          << "Last step (?) " << m_lastStep << "\n"
	          << "# of possible way-outs " << m_ExpectedWays.size()+ m_Breakpoints.size() + m_UnexpectedWays.size() << "\n"
	          << "Types:  " << m_status << "\n"
	          << "Nb times visited: " << m_nbVisits << "\n"
	          << "******************************************************" << "\n"
	          << std::endl;
}

TSeq catPaths(std::vector<Path >& Nodes, Direction& direction)
{
  TSeq fullPath;

  if(Nodes.size()>0)
   {
    Path tmp;
    
    if(direction == RIGHT)
    {
     for(int i(0); i < Nodes.size(); i++)
     {
        tmp = Nodes[i];
        append(fullPath, tmp.getSeq());
     }
    }
   else 
    {
     for(int i(Nodes.size()-1); i>=0; i--)
      {
        tmp = Nodes[i];
        append(fullPath, tmp.getSeq());
      }
    }
  }
 return fullPath;
}

