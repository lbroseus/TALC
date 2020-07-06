//
//  Explorer.cpp
//
//  Created by Lucile Broseus in April 2018. 
//  Copyright © 2018 LuB. All rights reserved.
//
//--------------------------------------------------------------------//
//                        CLASS DESCRIPTION
//
// Defined from two consecutive solid regions + the weak one
// in-between: SOURCE - GAP - TARGET
//
// Comes with functions that implement an exploration of the de Bruijn 
// graph from SOURCE to TARGET
//--------------------------------------------------------------------//

/****************************************************************************
 * Version TALC 1.01
 ****************************************************************************
 ****************************************************************************
 * Updates:
 * 27/01/2019: removed messages used for test
 * SomeDay/02/2019: Borders: stricter ID score when in case of pruning
 * 28/01/2019: Borders: relax thr in simple cases (few paths)
 ****************************************************************************/

#include "Explorer.hpp"

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <seqan/seeds.h>
#include <seqan/find.h>
#include <seqan/align.h>

#include "utils.hpp"
#include "Jellyfish.hpp"
#include "Settings.hpp"
#include "Trail.hpp"
#include "Trajectory.hpp"

using namespace seqan;

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<Dna5, colouredCount > colouredBase; 
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;
typedef std::tuple<unsigned int, unsigned int, Status> kmerStretch; 

//---------------------------------------------------------------------//

extern unsigned int gp_MIN_COUNT;
extern unsigned int K;

extern std::string queryMode;

extern colouredDBG SR_DBG;

extern std::string io_pathToJF;
extern std::string io_JFcountsTable;
extern std::string io_JFjunctcountsTable;

extern double gp_ALPHA;
extern unsigned int gp_WINDOW_SIZE;
extern unsigned int gp_MAX_NB_COMPETING_PATHS;
extern double gp_MIN_INNER_SCORE;
extern double gp_MIN_BORDER_SCORE;

extern double gp_SR_ERROR_RATE;

extern std::string io_logFile;

extern bool DEBUG;

extern bool DEBUG_TEST;
extern bool DEBUG_USER;

static unsigned int p_MIN_START_ANCHORS(3);
static unsigned int p_MAX_START_ANCHORS(5);
static unsigned int p_NB_LOW_COUNT_ALLOWED(2);
static unsigned int p_MAX_IN_COUNT(100000);

unsigned int p_MAX_NB_OF_BORDER_PATHS(75);
unsigned int p_MAX_NB_OF_INNER_PATHS(50);

unsigned int p_CHECK_INTERVAL(6);
double p_ALLOWED_FAILURE_RATE(0.3);

int p_MAX_NB_INNER_FAILURES(3);
int p_MAX_NB_BORDER_FAILURES(3);

int p_MAX_EXP_INS_LENGTH(5);
int p_MAX_EXP_DEL_LENGTH(5);

unsigned int p_MAX_REGION_LENGTH(1000);

std::vector<seqan::Dna5 > Bases = {'A', 'C', 'G', 'T'};

static unsigned int maxNbCompetingPaths(0);

Explorer::~Explorer()
{
}

Explorer::Explorer() :                m_sequence(make_empty_Dna5String()),
                                      m_weakSequence(make_empty_border()),
                                      m_coverage(make_empty_vecOfColouredPairs()),
                                      m_priorLambda_noise(gp_MIN_COUNT),
                                      m_LEFT_KMpositions(make_empty_KmPos()),
                                      m_RIGHT_KMpositions(make_empty_KmPos()),
                                      m_LEFT_anchors(make_empty_vecOfAnchors()),
                                      m_RIGHT_anchors(make_empty_vecOfAnchors()),
                                      m_location(UNKNOWN),
                                      m_direction(LEFT),
                                      m_nbGoodPathsFound(0), m_nbVeryGoodPathsFound(0),
                                      m_complexRegion(false),
                                      m_fullPaths(make_empty_vectOfTraj()), 
                                      m_longPaths(make_empty_vectOfTraj()),
                                      m_shortPaths(make_empty_vectOfTraj()),
                                      m_scorr(make_empty_CharString())
{
	
} 

Explorer::Explorer(TSeq refSequence, std::vector<colouredCount > coverage, double lambda) :      
                                      m_sequence(refSequence),
                                      m_weakSequence(make_empty_border()),
                                      m_coverage( coverage ), 
                                      m_priorLambda_noise(lambda),
                                      m_LEFT_KMpositions(make_empty_KmPos()),
                                      m_RIGHT_KMpositions(make_empty_KmPos()),
                                      m_LEFT_anchors(make_empty_vecOfAnchors()),
                                      m_RIGHT_anchors(make_empty_vecOfAnchors()),
                                      m_location(UNKNOWN),
                                      m_direction(RIGHT),
                                      m_nbGoodPathsFound(0), m_nbVeryGoodPathsFound(0),
                                      m_complexRegion(false),
                                      m_fullPaths(make_empty_vectOfTraj()), 
                                      m_longPaths(make_empty_vectOfTraj()),
                                      m_shortPaths(make_empty_vectOfTraj()),
                                      m_scorr(make_empty_CharString())
{
	
} 

//Methods

void Explorer::reset()
{
	m_weakSequence = make_empty_border();
	m_LEFT_KMpositions = make_empty_KmPos();
	m_RIGHT_KMpositions = make_empty_KmPos();
	
	m_LEFT_anchors = make_empty_vecOfAnchors();
	m_RIGHT_anchors = make_empty_vecOfAnchors();
		
	m_location = UNKNOWN;
	m_direction = RIGHT;
	m_nbGoodPathsFound = 0;
    m_nbVeryGoodPathsFound = 0;
    
    m_fullPaths = make_empty_vectOfTraj();
    m_longPaths = make_empty_vectOfTraj();
    m_shortPaths = make_empty_vectOfTraj();
    m_scorr = make_empty_CharString();
    
}

void Explorer::setCoverage(std::vector<colouredCount > counts)
{
	m_coverage = counts;
}

void Explorer::setSequence(TSeq sequence)
{
   m_sequence = sequence;
}

kmerStretch Explorer::getLEFTHandPositions()
{
	return m_LEFT_KMpositions;
}
kmerStretch Explorer::getRIGHTHandPositions()
{
	return m_RIGHT_KMpositions;
}

void Explorer::setLEFTHandPositions(kmerStretch positions)
{
	m_LEFT_KMpositions = positions;
}

void Explorer::setRIGHTHandPositions(kmerStretch  positions)
{
	m_RIGHT_KMpositions = positions;
}



std::pair<TSeq, Status > Explorer::getWeakSeq()
{
  return m_weakSequence;
}

unsigned int Explorer::getWeakLength()
{
  return length(std::get<0>(m_weakSequence));
}


void Explorer::setWeakSequence()
{ 
 if(m_location == INNER) std::get<0>(m_weakSequence) = extractWeakSequence(m_sequence, std::get<1>(m_LEFT_KMpositions), std::get<0>(m_RIGHT_KMpositions), K);
 else if(m_location == HEAD) std::get<0>(m_weakSequence) = extractWeakBorderSequence(m_sequence, std::get<0>(m_RIGHT_KMpositions), K, m_location);
 else std::get<0>(m_weakSequence) = extractWeakBorderSequence(m_sequence, std::get<1>(m_LEFT_KMpositions), K, m_location);
 
 std::get<1>(m_weakSequence) = UNCORRECTED;
 	
}

void Explorer::initializeINNER(kmerStretch startKmPos, kmerStretch endKmPos, Direction direction)
{
    this->reset();
    
    this->setLocation( INNER );
    this->setDirection( direction );
    
    this->setLEFTHandPositions( startKmPos );
    this->setRIGHTHandPositions( endKmPos );
    
	this->setWeakSequence();
	
	this->anchorLEFTHandSide();
	this->anchorRIGHTHandSide();
    		
}

void Explorer::initializeHEAD(kmerStretch kmPos)
{
    this->reset();
    
    this->setLocation( HEAD );
    this->setDirection( LEFT );
    this->setRIGHTHandPositions( kmPos );
    
	this->setWeakSequence();
	
	this->anchorRIGHTHandSide();
    		
}

void Explorer::initializeTAIL(kmerStretch kmPos)
{
    this->reset();
    
    this->setLocation( TAIL );
    this->setDirection( RIGHT );
    
    this->setLEFTHandPositions( kmPos );
	this->setWeakSequence();
	
	this->anchorLEFTHandSide();
    		
}


Direction Explorer::getDirection()
{
	return m_direction;
}

void Explorer::setLocation(Location location)
{
	m_location = location;
}

void Explorer::setDirection(Direction direction)
{
	m_direction = direction;
}

CharString Explorer::getScorr()
{
	return m_scorr;
}

void Explorer::setScorr(CharString scorr)
{
	m_scorr = scorr;
}

Trajectory Explorer::getBestBridge()
{
	int bestOne(0);
	if( empty(m_fullPaths) ) return Trajectory();
	else
	 {
	   bestOne = findBestBridge(m_fullPaths);
	   return m_fullPaths[ findBestBridge(m_fullPaths) ];
	 } 
}

Trajectory Explorer::sortOutBestBorder()
{
  std::pair<bool, unsigned int > result;
  unsigned int best(0);
  bool isConvenient(false);
	
  result = findBestBORDER(m_longPaths, gp_MIN_BORDER_SCORE);
  if( std::get<0>(result) ) return m_longPaths[std::get<1>(result)];
  else
   {
	//std::cout << "NO LONG BORDER PATH" << std::endl;
	result = findBestBORDER(m_shortPaths, gp_MIN_BORDER_SCORE);
	if( std::get<0>(result) ) return m_shortPaths[std::get<1>(result)];
    else
     {
	  // std::cout << "NO CONVENIENT SHORT PATH EITHER..." << std::endl;
	  return Trajectory();
     }
    } 	
}


/************************************************************************
 * Behaviour when getting back to previous branches
 * getBack():
 * 
 * re-initiate exploration on last recorded possible ways
 * when complexity is HIGH, ignore unpredicted ways
 * otherwise explore any remaining way
 * Return TRUE if exploration can go on, FALSE otherwise
 * if TRUE the new way to be followed should have been initiated
 ***********************************************************************/ 


std::vector<unsigned int >  findNearest(unsigned int& cc, std::vector<colouredCount >& nextCounts)
{
  std::vector<unsigned int > result;
  unsigned int index(0);
  while( (std::get<0>(nextCounts[index])<gp_MIN_COUNT) & index<nextCounts.size() ) ++index;
  for(unsigned int i(index+1); i<nextCounts.size();i++) if( (std::get<0>(nextCounts[i])>=gp_MIN_COUNT) & (abs((int)cc-(int)std::get<0>(nextCounts[index]))>abs((int)cc-(int)std::get<0>(nextCounts[i])))) index=i;
  result.push_back(index);
  for(unsigned int i(index+1); i<nextCounts.size();i++) if( (std::get<0>(nextCounts[i])>=gp_MIN_COUNT) & (abs((int)cc-(int)std::get<0>(nextCounts[index]))==abs((int)cc-(int)std::get<0>(nextCounts[i])))) result.push_back(i);

  return( result );
}

std::vector<unsigned > findMax(unsigned int& cc, std::vector<colouredCount >& nextCounts)
{
  std::vector<unsigned int > result;
  unsigned int index(0);
  while( (std::get<0>(nextCounts[index])<gp_MIN_COUNT) & index<nextCounts.size() ) ++index;
  for(unsigned int i(index+1); i<nextCounts.size();i++) if( std::get<0>(nextCounts[index])<std::get<0>(nextCounts[i]) ) index=i;
  result.push_back(index);
  for(unsigned int i(index+1); i<nextCounts.size();i++) if( std::get<0>(nextCounts[index])==std::get<0>(nextCounts[i]) ) result.push_back(i);

  return( result );	
}

bool isIn(std::vector<unsigned int >& set, unsigned int& element)
{
	bool isIn(false);
	for(unsigned int i(0); i<set.size();i++) if(element == set[i]) isIn = true;
	return isIn;
}


/* Order pairs <nucleotide, count> according to count */

int distance(const colouredCount& cc, const std::pair<Dna5, colouredCount >& cnextc )
{
	return abs((int)std::get<0>(cc)-(int)std::get<0>(std::get<1>(cnextc)));
}

void sortByNearest(std::vector<std::pair<Dna5, colouredCount > >& ways, colouredCount& cc)
{
 std::sort(ways.begin(), ways.end(),
          [cc](const std::pair<Dna5, colouredCount >& lhs, const std::pair<Dna5, colouredCount >& rhs){ return distance(cc, lhs) < distance(cc, rhs); });               
}

/************************************************************************
 * Behaviour at branching points
 * recordBranch():
 * 
 * if the total # of ways to be explored is over the p_MAX_PATH threshold exploration has high complexity
 *    -> in that case if there is no GPred base to follow, STOP EXPLORATION
 *    -> otherwise record all the ways that are available
 *    -> first follow the LPred ones if exist, then the GPred ones
 *       then possibly the (most frequent first) unpred ones
 * Evaluate the local branching rate; if is over threshold, prune UnexpectedWays of last nodes
 ***********************************************************************/                                  


int distance_anchors(const double cc, const anchorTuple& anc )
{
	return abs((int)cc-(int)std::get<2>(anc));
}

void sortAnchorsByNearest(double cc, std::vector<std::tuple<TSeq, unsigned int, unsigned int > >& anchors)
{
 std::sort(anchors.begin(), anchors.end(),
          [cc](const anchorTuple& lhs, const anchorTuple& rhs){ return distance_anchors(cc, lhs) < distance_anchors(cc, rhs); });               
}

void Explorer::anchorLEFTHandSide()
{
	bool goFurther(true);
    unsigned int nbKmers = std::get<1>(m_LEFT_KMpositions)-std::get<0>(m_LEFT_KMpositions) + 1;
    unsigned int len = length(m_sequence);
    unsigned int pivot = std::get<1>(m_LEFT_KMpositions);
    unsigned int limit = std::get<0>(m_LEFT_KMpositions);
    unsigned int lowCounter(1);
    int degree(0);
    
    //LEFT ANCHORS (on the right hand side of the LEFT-hand source region
    TSeq anchor;
    std::vector<unsigned int > anchorPos;
    double current_count = (double)std::get<0>(m_coverage[pivot]);
    double next_count(0);
    unsigned int j(pivot);
    //First screening of the stretch of solid kmers to detect coverage variations: at least one anchor by coverage breakpoint   
    anchorPos.push_back(pivot);
    while( goFurther & (j>=limit+1) )
    {
	 next_count = std::get<0>(m_coverage[j-1]);	
	 if(next_count>=gp_MIN_COUNT & next_count<p_MAX_IN_COUNT) goFurther = isExpectedbyMyLastNode(next_count, current_count, gp_ALPHA);   
     else goFurther = false;
     
     if(!goFurther & (current_count>=gp_MIN_COUNT) & (next_count>=gp_MIN_COUNT) & (next_count<p_MAX_IN_COUNT) ) //stopped by lack of count homogeneity, add another anchor
    {
     anchorPos.push_back(j-1);
     goFurther = true;
     current_count = next_count;
    } 
     --j;
    }//END-WHILE 
   //-------------------------------------------------->START anchor sequence extraction       
    if(empty(anchorPos)) std::cout << "Pb during right anchors definition - first count: " << std::get<0>(m_coverage[0]) << std::endl;
    else //Anchors' sequence extraction
     {
       for(int anc(0); anc<anchorPos.size(); anc++)
        {
	     anchor = getKmerAt(m_sequence, anchorPos[anc], K);
	     if(length(anchor) != K) std::cout << "PB in right anchors " << std::endl;
	     degree = getOutDegree(anchor, RIGHT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);
	     if((anc==0) || (anc!=0 & degree>1)) m_LEFT_anchors.push_back(std::make_tuple(anchor, anchorPos[anc], std::get<0>(m_coverage[anc])));
	    }   			
      }
    //if not enough anchors, add some  
    if( empty(m_RIGHT_anchors) ) anchorPos.push_back(pivot);
    if( m_LEFT_anchors.size()<std::min(p_MIN_START_ANCHORS, nbKmers) )
     {
	  j = pivot;
	  goFurther = true;
	  while( (j>=limit+1) &  (m_LEFT_anchors.size()<std::min(p_MIN_START_ANCHORS, nbKmers)) )
	  {
		for(unsigned int i(0); i<m_LEFT_anchors.size(); --i) goFurther &= (std::get<1>(m_LEFT_anchors[i]) != (j-1) );
		if(  goFurther )
		 {
		  anchor = getKmerAt(m_sequence, j-1, K);		  
	      if(length(anchor) != K) std::cout << "PB in right anchors " << std::endl;
	      degree = getOutDegree(anchor, RIGHT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);
	      if( degree>1 ) m_LEFT_anchors.push_back(std::make_tuple(anchor, j-1, std::get<0>(m_coverage[j-1])));
		 }
		--j;
	  }//END-WHILE
	 } 
  sortAnchorsByNearest(m_priorLambda_noise/gp_SR_ERROR_RATE, m_LEFT_anchors);	 

}      

void Explorer::anchorRIGHTHandSide()
{
	bool goFurther(true);
    unsigned int nbKmers = std::get<1>(m_RIGHT_KMpositions)-std::get<0>(m_RIGHT_KMpositions) + 1;
    unsigned int len = length(m_sequence);
    unsigned int pivot = std::get<0>(m_RIGHT_KMpositions);
    unsigned int limit = std::get<1>(m_RIGHT_KMpositions);
    unsigned int lowCounter(1);
    
    //RIGHT ANCHORS (on the right hand side of the RIGHT-hand source region
    TSeq anchor;
    std::vector<unsigned int > anchorPos;
    double current_count = (double)std::get<0>(m_coverage[pivot]);
    double next_count(0);
    unsigned int j(pivot);
    //First screening of the stretch of solid kmers to detect coverage variations: at least one anchor by coverage breakpoint   
    anchorPos.push_back(pivot);
    while( goFurther & ((j+1)<=limit) )
    {
	 next_count = std::get<0>(m_coverage[j+1]);	
	 if(next_count>=gp_MIN_COUNT & next_count<p_MAX_IN_COUNT) goFurther = isExpectedbyMyLastNode(next_count, current_count, gp_ALPHA);
     else  goFurther = false;            
     
     if(!goFurther & (current_count>=gp_MIN_COUNT) & (next_count>=gp_MIN_COUNT) & (next_count<p_MAX_IN_COUNT) ) //stopped by lack of count homogeneity, add another anchor
    {
     anchorPos.push_back(j+1);
     goFurther = true;
     current_count = next_count;
    } 
     ++j;
    }//END-WHILE 
    
    if(empty(anchorPos)) std::cout << "Pb during right anchors definition - first count: " << std::get<0>(m_coverage[0]) << std::endl;
    else //Anchors' sequence extraction
     {
       for(int anc(0); anc<anchorPos.size(); anc++)
        {
	     anchor = getKmerAt(m_sequence, anchorPos[anc], K);
	     if(length(anchor) != K) std::cout << "PB in right anchors " << std::endl;
	     int degree = getOutDegree(anchor, LEFT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);
	     if((anc==0) || (anc!=0 & degree>1)) m_RIGHT_anchors.push_back(std::make_tuple(anchor, anchorPos[anc], std::get<0>(m_coverage[anc])));
	    }   			
      }
    //if not enough anchors, add some  
    if( empty(m_RIGHT_anchors) ) anchorPos.push_back(pivot);
    if( m_RIGHT_anchors.size()<std::min(p_MIN_START_ANCHORS, nbKmers) )
     {
	  j = pivot;
	  goFurther = true;
	  while( ((j+1)<=limit) &  (m_RIGHT_anchors.size()<std::min(p_MIN_START_ANCHORS, nbKmers)) )
	  {
		for(unsigned int i(0); i<m_RIGHT_anchors.size(); --i) goFurther &= (std::get<1>(m_RIGHT_anchors[i]) != (j+1) );
		if(  goFurther )
		 {
		  anchor = getKmerAt(m_sequence, j+1, K);		  
	      if(length(anchor) != K) std::cout << "PB in right anchors " << std::endl;
	      int degree = getOutDegree(anchor, LEFT, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);
	      if( degree>1 ) m_RIGHT_anchors.push_back(std::make_tuple(anchor, j+1, std::get<0>(m_coverage[j+1])));
		 }
		--j;
	  }//END-WHILE
	 }	  	 
  sortAnchorsByNearest(m_priorLambda_noise/gp_SR_ERROR_RATE, m_RIGHT_anchors);	 
}  

//Inner walk
 void Explorer::oneMoreStep(TSeq& reference, std::vector<Trail>& competingPaths, unsigned int& stepCounter, Direction& direction, 
                           unsigned int& PATH_MAXLENGTH, colouredDBG& SR_DBG, std::string& io_JFcountsTable, std::string& io_JFjunctcountsTable)
{
 std::vector<Trail > newCompetingPaths;
 
 std::vector<colouredCount > nextCounts;                                // where the 4 following counts (A, C, G, T) are stored
 std::vector<std::pair<Status, double > > nodeTags;
 std::vector<unsigned int > indexOfKeptPaths;
 bool cycle(false);
 bool aimReached(false);
 bool complex(false);
 //For each Trail in competingPath:
 // query next-kmers 
 // filter bases
 //create corresponding Path+1 -> newCompetingPath
 for(unsigned int t(0); t<competingPaths.size(); t++)
  {
	complex = (competingPaths.size()>gp_MAX_NB_COMPETING_PATHS);
    nextCounts = competingPaths[t].whatsNext(m_direction, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable); 
    tagNextNodes(nodeTags, nextCounts, competingPaths[t].getLastCount(), gp_ALPHA, m_priorLambda_noise, complex);

    //creating corresponding paths if needed (no dead-end)
    for(unsigned int i(0); i<nodeTags.size();i++)
     {
	   if(std::get<0>(nodeTags[i])!=UNEXPECTED)
	    {
		 newCompetingPaths.push_back( Trail(competingPaths[t], Bases[i], m_direction, std::get<0>(nextCounts[i])));		 
		 if(std::get<0>(nodeTags[i])==BREAKPOINT) newCompetingPaths.back().recordBreakpoint();
         newCompetingPaths.back().recordDistance(std::get<1>(nodeTags[i]));
         if(m_direction == RIGHT) aimReached = newCompetingPaths.back().checkAims(m_RIGHT_anchors, m_direction);
         else aimReached = newCompetingPaths.back().checkAims(m_LEFT_anchors, m_direction);
         if(aimReached){
			 this->recordBridge(newCompetingPaths.back());
			 if( newCompetingPaths.back().getLength()>length(reference) )
			 {
			  newCompetingPaths.pop_back();
		     }
		  }
		 else
		  { 
           cycle = newCompetingPaths.back().ThinkIveAlreadyGotThere(competingPaths[t].getSeq());
           if(cycle) newCompetingPaths.pop_back();
	      }
	    }   
     }
   }
//if( empty(newCompetingPaths) ) std::cout << "[BRIDGE] STOP at step " << stepCounter << " " << cycle << std::endl;

   	complex = (newCompetingPaths.size()>gp_MAX_NB_COMPETING_PATHS);
   ++stepCounter;
    
  unsigned int len = newCompetingPaths.size();
  maxNbCompetingPaths = std::max(maxNbCompetingPaths, len);
  
  //Scoring and ranking
  if( complex & (stepCounter%p_CHECK_INTERVAL==0) )
   {
    scoreBridges(newCompetingPaths, stepCounter, reference, m_direction);
    m_complexRegion |= doABitOfGardening(indexOfKeptPaths, newCompetingPaths);
    clear(competingPaths);
    for(unsigned int t(0); t<indexOfKeptPaths.size();t++) competingPaths.push_back(newCompetingPaths[indexOfKeptPaths[t]]);
   }
  else competingPaths = newCompetingPaths;

//if( empty(competingPaths) ) std::cout << "[BRIDGE] END-STOP at step " << stepCounter << " " << cycle << std::endl;
 
}

//Border walk
void Explorer::oneMoreStepInTheDark(int& xdrop, TSeq& reference, std::vector<Trail>& competingPaths, unsigned int& stepCounter, 
                                    Direction& direction, unsigned int& PATH_MAXLENGTH, 
                                    colouredDBG& SR_DBG, std::string& io_JFcountsTable, std::string& io_JFjunctcountsTable)
{
 std::vector<Trail > newCompetingPaths;
 
 std::vector<colouredCount > nextCounts;                                // where the 4 following counts (A, C, G, T) are stored
 std::vector<std::pair<Status, double > > nodeTags;
 std::vector<unsigned int > indexOfKeptPaths;
 bool complex(false);
 bool cycle(false);
 unsigned int counter(0);
 //For each Trail in competingPath:
 // query next-kmers 
 // filter bases
 //create corresponding Path+1 -> newCompetingPath
 for(unsigned int t(0); t<competingPaths.size(); t++)
  {
	counter = 0;
	complex = (competingPaths.size()>7);
    nextCounts = competingPaths[t].whatsNext(m_direction, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable); 
    tagNextNodes(nodeTags, nextCounts, competingPaths[t].getLastCount(), gp_ALPHA, m_priorLambda_noise, complex);

    //creating corresponding paths if needed (no dead-end)
    for(unsigned int i(0); i<nodeTags.size();i++)
     {
	   if(std::get<0>(nodeTags[i])!=UNEXPECTED)
	    {
		 ++counter;
		 newCompetingPaths.push_back( Trail(competingPaths[t], Bases[i], m_direction, std::get<0>(nextCounts[i])));
		 if(std::get<0>(nodeTags[i])==BREAKPOINT) newCompetingPaths.back().recordBreakpoint();
		 newCompetingPaths.back().recordDistance(std::get<1>(nodeTags[i]));
		 cycle = newCompetingPaths.back().ThinkIveAlreadyGotThere(competingPaths[t].getSeq());
         if( cycle || (stepCounter+1>PATH_MAXLENGTH) )
          {
if(DEBUG_TEST) std::cout << "[EDGE]: path stop at: " << stepCounter << " | cycle? " << cycle << std::endl;	  
		   newCompetingPaths.back().seedAndExtend(reference, m_direction, xdrop, p_MAX_NB_BORDER_FAILURES);
		   this->recordEdge(newCompetingPaths.back(), reference);
		   newCompetingPaths.pop_back();
		  }
	    }  
     }
    if(counter==0) //no continuation-> dead-end path
	 {
	  if(DEBUG_TEST) std::cout << "[DEAD-END]: reached by path: " << t+1 << " out of " << competingPaths.size() << std::endl;
	  competingPaths[t].seedAndExtend(reference, m_direction, xdrop, p_MAX_NB_BORDER_FAILURES);
	  this->recordEdge(competingPaths[t], reference);
	 }  
   }
   ++stepCounter;
   
if( DEBUG_TEST & empty(newCompetingPaths) ) std::cout << "TOTAL NB OF COMPETING PATHS BEFORE:" <<  competingPaths.size() << std::endl; 
if( DEBUG_TEST & empty(newCompetingPaths) ) std::cout << "STOP at step " << stepCounter << " " << cycle << std::endl;
  
  unsigned int len = newCompetingPaths.size();
  maxNbCompetingPaths = std::max(maxNbCompetingPaths, len);
  
  //Scoring and ranking
  if( (stepCounter%p_CHECK_INTERVAL==0) || (newCompetingPaths.size()>=p_MAX_NB_OF_BORDER_PATHS) )
   {
    scoreEdges(xdrop, newCompetingPaths, stepCounter, reference, m_direction);
    
    if( newCompetingPaths.size()>5 )
     {
	  m_complexRegion |= doABitOfGardening(indexOfKeptPaths, newCompetingPaths);
      clear(competingPaths);
      for(unsigned int t(0); t<indexOfKeptPaths.size();t++) competingPaths.push_back(newCompetingPaths[indexOfKeptPaths[t]]);
     }
    else competingPaths = newCompetingPaths;
    
   }
  else competingPaths = newCompetingPaths; 
}

void scoreBridges(std::vector<Trail>& newCompetingPaths, unsigned int& stepCounter, TSeq& reference, Direction& direction)
{
 TSeq truncatedReference;
 int bound(0);
 if( direction==RIGHT )
  {
   bound = K+stepCounter+gp_WINDOW_SIZE;
   if(bound>=length(reference)) truncatedReference = reference;
   else truncatedReference = prefix(reference, bound);
  } 
 else 
  {
   bound = (int)length(reference)-(int)K-(int)stepCounter-(int)gp_WINDOW_SIZE;
   if(bound<0) truncatedReference = reference;
   else truncatedReference = suffix(reference, bound);
  }
 for(unsigned int(j); j<newCompetingPaths.size(); j++) newCompetingPaths[j].Overlapscore(truncatedReference, direction);
}


void Explorer::scoreEdges(int& xdrop, std::vector<Trail>& newCompetingPaths, unsigned int& stepCounter, TSeq& reference, Direction& direction)
{
 bool boolean(true);
 std::vector<Trail> newSelectedPaths;
 std::vector<Trail> trashPaths;
 int new_xdrop(0);
 int current_xdrop(0);
 if( !empty(newCompetingPaths) )
  {
	xdrop += 2;
	for(unsigned int t(0); t<newCompetingPaths.size(); t++)
	 {
	  boolean = newCompetingPaths[t].seedAndExtend(reference, direction, xdrop, p_MAX_NB_BORDER_FAILURES);
	  if(!boolean) trashPaths.push_back(newCompetingPaths[t]);
	  else
	   {
		newSelectedPaths.push_back(newCompetingPaths[t]);
	    current_xdrop = newCompetingPaths[t].getLastScore() * (-1);
	    if( (new_xdrop>current_xdrop) || (new_xdrop==0)) new_xdrop = current_xdrop;
	   }
     }

  xdrop = new_xdrop;
  newCompetingPaths = newSelectedPaths;  

if( empty(newCompetingPaths) ) //all paths have been stopped and where the last best candidates
 {
  boolean = (stepCounter+K>length(reference));
  for(unsigned int t(0); t<trashPaths.size(); t++) this->recordEdge(trashPaths[t], reference);
  }
 }
}

void sortByScore(std::vector<std::tuple<unsigned int, double, double> >& index)
{
 std::sort(index.begin(), index.end(),
          [](const std::tuple<unsigned int, double, double>  & lhs, const std::tuple<unsigned int, double, double>& rhs){ return std::get<1>(lhs) > std::get<1>(rhs); });               
}

void sortByLikelihood(std::vector<std::tuple<unsigned int, double, double> >& index)
{
 std::sort(index.begin(), index.end(),
          [](const std::tuple<unsigned int, double, double>  & lhs, const std::tuple<unsigned int, double, double>& rhs){ return std::get<2>(lhs) < std::get<2>(rhs); });               
}

void sortByMinRank(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> >& ranking)
{
 std::sort(ranking.begin(), ranking.end(),
          [](const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>   & lhs, const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int > & rhs){ return std::get<3>(lhs) < std::get<3>(rhs); });               
}

void sortByMaxScore(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> >& ranking)
{
 std::sort(ranking.begin(), ranking.end(),
          [](const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>   & lhs, const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int > & rhs){ return std::get<1>(lhs) < std::get<1>(rhs); });               
}

void sortByMinDist(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> >& ranking)
{
 std::sort(ranking.begin(), ranking.end(),
          [](const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>   & lhs, const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> & rhs){ return std::get<2>(lhs) < std::get<2>(rhs); });               
}


bool doABitOfGardening(std::vector<unsigned int >&  indexOfKeptPaths, std::vector<Trail>& newCompetingPaths)
{
 clear(indexOfKeptPaths);
 
 std::vector<std::tuple<unsigned int, double, double> > indexOfTrails1;
  std::vector<std::tuple<unsigned int, double, double> > indexOfTrails2;

  std::vector<unsigned int > rankWithTies1;
 std::vector<unsigned int > rankWithTies2;
 //                     //index - score rank - break rank - min rank - sum ranks
 std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> > rankings;
 std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> > newrankings;

 bool ties(true);
 bool isComplex(false); 
 unsigned int nb;
 unsigned int s(0);

 nb = newCompetingPaths.size();
 nb = std::min(nb, gp_MAX_NB_COMPETING_PATHS);
 for(unsigned t(0); t<newCompetingPaths.size();t++)
  { 
   rankings.push_back(std::make_tuple(t,0, 0, 0));
   rankWithTies1.push_back(t);
   rankWithTies2.push_back(t);
   indexOfTrails1.push_back(std::make_tuple(t, newCompetingPaths[t].getLastScore(), newCompetingPaths[t].getDistance()));
   indexOfTrails2.push_back(std::make_tuple(t, newCompetingPaths[t].getLastScore(), newCompetingPaths[t].getDistance()));
  } 
 
 sortByScore(indexOfTrails1);
 std::get<1>(rankings[std::get<0>(indexOfTrails1[0])]) = rankWithTies1[0];
 for(unsigned int t(1); t<rankWithTies1.size(); t++)
  {
   if(std::get<1>(indexOfTrails1[t])==std::get<1>(indexOfTrails1[t-1])) rankWithTies1[t]=rankWithTies1[t-1];
   else rankWithTies1[t]=rankWithTies1[t-1]+1;
   std::get<1>(rankings[std::get<0>(indexOfTrails1[t])]) = rankWithTies1[t];
  } 

 sortByLikelihood(indexOfTrails2);
 std::get<2>(rankings[std::get<0>(indexOfTrails2[0])]) = rankWithTies2[0];
 for(unsigned int t(1); t<rankWithTies2.size(); t++)
  {
   if(std::get<2>(indexOfTrails2[t])==std::get<2>(indexOfTrails2[t-1])) rankWithTies2[t]=rankWithTies2[t-1];
   else rankWithTies2[t]=rankWithTies2[t-1]+1;
   std::get<2>(rankings[std::get<0>(indexOfTrails2[t])]) = rankWithTies2[t];
  } 

 for(unsigned int t(0); t<rankings.size(); t++)
  {
   std::get<3>(rankings[t]) = std::get<1>(rankings[t])+std::get<2>(rankings[t]);
   if( (std::get<3>(rankings[t])==0) || (newCompetingPaths.size()<=gp_MAX_NB_COMPETING_PATHS) ) newrankings.push_back(rankings[t]);
  }
     //else, order path by similarity score
 if( empty( newrankings ) ) //best path according to score do not match best paths according to likelihood
  {
	sortByMaxScore( rankings );
	
    s = 0;
    ties = false;
    do
    {
     if( (s<=nb) || ties )  newrankings.push_back(rankings[s]);
     if(s<rankings.size()-1) ties = (std::get<1>(rankings[s+1])==std::get<1>(rankings[s]));
     ++s;
    }while( ((s<=nb) || ties) & (s<rankings.size()) );
      
 if( newrankings.size()>gp_MAX_NB_COMPETING_PATHS )
  {
	if(std::get<1>(newrankings[0])!=std::get<1>(newrankings[gp_MAX_NB_COMPETING_PATHS]))
	 {
	   newrankings.pop_back();
	   ties=true;
	   while( ( newrankings.size()>=gp_MAX_NB_COMPETING_PATHS) & ties)
	    {
		 ties = (std::get<1>(newrankings.back())==std::get<1>(newrankings[newrankings.size()-2]));
		 ties |= ( newrankings.size()>=gp_MAX_NB_COMPETING_PATHS);
		 if( ties )  newrankings.pop_back();
		}
     }
    if( newrankings.size()>gp_MAX_NB_COMPETING_PATHS & std::get<1>(newrankings[0])==std::get<1>(newrankings[gp_MAX_NB_COMPETING_PATHS]) )
     {
	  isComplex = true;
	  sortByMinDist( newrankings );
	  for(unsigned int t(0); t<gp_MAX_NB_COMPETING_PATHS; t++) indexOfKeptPaths.push_back(std::get<0>(newrankings[t]));
     }
  }
if( empty(newrankings) ) std::cout << "[ISSUE WHEN GARDENING]: no path kept !" << std::endl;
 for(unsigned int t(0); t<newrankings.size();  t++) indexOfKeptPaths.push_back(std::get<0>(newrankings[t]));
 }
 else for(unsigned int t(0); t<newrankings.size();  t++) indexOfKeptPaths.push_back(std::get<0>(newrankings[t]));

 return isComplex;
}


bool Explorer::searchBridge()
{

//--------------------------------------------------------------------//
 bool pathHasBeenFound(false);
 unsigned int PATH_MAXLENGTH(0); 
 unsigned int totalNbOfBranches(0);
 unsigned int limit; 
 unsigned int stepCounter(0);
 double minScore(0);
 double diff(0);
    
 std::vector<anchorTuple > aims;
 std::vector<anchorTuple > anchors;
 
 colouredCount currentCount;
 TSeq currentGap;
 TSeq currentAnchor;
 unsigned int whichStart(0);
 TSeq currentTarget;
 TSeq currentRefSeq;
 TSeq bestPath;
 unsigned int bestOne;
 std::pair<TSeq, bool>  result;

 std::vector<Trail> competingPaths;

 //--------------------------------------------------------------------//

  anchors = (m_direction == LEFT) ? m_RIGHT_anchors : m_LEFT_anchors;      
  aims = (m_direction == LEFT) ? m_LEFT_anchors : m_RIGHT_anchors; 
   
	   	
  limit = anchors.size();
  limit = std::min(limit, p_MAX_START_ANCHORS);

 for(int s(0); s < limit ; s++)  
  {
   if( !pathHasBeenFound )
    {
     clear(competingPaths);
     clear(m_fullPaths);
	 stepCounter = 0;
	 maxNbCompetingPaths = 0;
     whichStart = std::get<1>(anchors[s]);
     currentTarget = (m_direction == RIGHT) ? extractSolidSequence(m_sequence, std::get<0>(m_RIGHT_KMpositions), std::get<1>(m_RIGHT_KMpositions), K): extractSolidSequence(m_sequence, std::get<0>(m_LEFT_KMpositions), std::get<1>(m_LEFT_KMpositions), K);
     
     //function initiateGap
     clear(currentGap);     //should remain empty if following conditions are not met
     if( (m_direction == RIGHT) & (whichStart+K<std::get<0>(m_RIGHT_KMpositions))) currentGap = extractWeakSequence(m_sequence, whichStart, std::get<0>(m_RIGHT_KMpositions), K);
     else if( (m_direction == LEFT) & (std::get<1>(m_LEFT_KMpositions)+K<whichStart)) currentGap = extractWeakSequence(m_sequence, std::get<1>(m_LEFT_KMpositions), whichStart, K);
     PATH_MAXLENGTH = (int)(1.2*length(currentGap)+3*K);
     
     currentAnchor = std::get<0>(anchors[s]);
     currentCount = m_coverage[whichStart];
     competingPaths.push_back(Trail(currentAnchor, currentCount));

     if( m_direction == RIGHT )
      {
		currentRefSeq = currentAnchor;
		append(currentRefSeq, currentGap);
		append(currentRefSeq, currentTarget);
		competingPaths[0].setLeftAnchor(whichStart);
	  }
	 else
	  {
	   currentRefSeq = currentTarget;
	   append(currentRefSeq, currentGap);
	   append(currentRefSeq, currentAnchor);
	   competingPaths[0].setRightAnchor(whichStart);
	  }
    
  while( (!empty(competingPaths)) & (competingPaths.size()<=p_MAX_NB_OF_INNER_PATHS) & (stepCounter<PATH_MAXLENGTH))
   { 
     this->oneMoreStep(currentRefSeq, competingPaths, stepCounter, m_direction , PATH_MAXLENGTH, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);                                  
   }

 if(!empty(m_fullPaths) & !pathHasBeenFound)
  {
	std::vector<unsigned int > index;
	unsigned int limit = std::get<1>(m_RIGHT_KMpositions);
	for(unsigned int t(0); t<m_fullPaths.size(); t++){
		m_fullPaths[t].scoreSequence(currentRefSeq);
		
		if(m_fullPaths[t].cutAnchors(m_location, limit, K)) index.push_back(t);
	}
	
	if(m_fullPaths.size() != index.size()) 
	 {
	  m_shortPaths = m_fullPaths;
	  clear(m_fullPaths);
	  for(unsigned int t(0); t<index.size(); t++) m_fullPaths.push_back(m_shortPaths[t]);
     }
	
  if(!empty(m_fullPaths) )
  {
    //select best trajectory and polish
    bestOne = findBestBridge(m_fullPaths);
    bestPath = m_fullPaths[bestOne].getSeq();

    diff = (double)this->getWeakLength()-(double)length(bestPath);

  if( (this->getWeakLength()>=200) || m_complexRegion ) minScore = std::max(0.75, gp_MIN_INNER_SCORE);
  else minScore = gp_MIN_INNER_SCORE;
 
  if(  ( diff<this->getWeakLength()*0.05 || ((this->getWeakLength()<6) & (length(bestPath)<6)) ) & (m_fullPaths[bestOne].getIDScore()>=gp_MIN_INNER_SCORE) )
  {
     std::get<1>(m_LEFT_KMpositions) = m_fullPaths[bestOne].getLeftAnchor();
	 std::get<0>(m_RIGHT_KMpositions) = m_fullPaths[bestOne].getRightAnchor();	
     std::get<0>(m_weakSequence) = m_fullPaths[bestOne].getSeq();  
	 std::get<1>(m_weakSequence) = CORRECTED; 
     pathHasBeenFound = true;
   }
  }
 }
   
 }//END-IF no bridge has been found yet
  
 }//END-FOR on selected start points
 
  return pathHasBeenFound; 
}


bool Explorer::searchEdge()
{	
 bool pathHasBeenFound(false);
 unsigned int PATH_MAXLENGTH(0); 
 unsigned int totalNbBranches(0);
 unsigned int limit; 
 unsigned int stepCounter(0);
 double diff(0);
 double minScore(0);
 
 std::vector<anchorTuple > anchors;
 colouredCount currentCount;
 TSeq currentGap;
 unsigned int whichStart(0);
 TSeq currentAnchor;
 TSeq currentRefSeq;
 Trajectory winner;
 int xdrop;
  
 anchors = (m_direction == LEFT) ? m_RIGHT_anchors : m_LEFT_anchors; 
 limit = anchors.size();
 limit = std::min(limit,p_MAX_START_ANCHORS);
  
 std::vector<Trail> competingPaths;

 //--------------------------------------------------------------------//
 
  anchors = (m_direction == LEFT) ? m_RIGHT_anchors : m_LEFT_anchors;      
	   	
  limit = anchors.size();
  limit = std::min(limit, p_MAX_START_ANCHORS);

 for(int s(0); s < limit ; s++)  
  {
if( DEBUG_TEST ) std::cout << "[EDGE] Trying anchor: " << s+1 << " kmer #" << std::get<1>(anchors[s]) << " out of " << limit << std::endl;    		
	 stepCounter = 0;
	 maxNbCompetingPaths = 0;
	 clear(competingPaths);
	 
	 xdrop = (int)p_CHECK_INTERVAL*p_ALLOWED_FAILURE_RATE + 1;

     whichStart = std::get<1>(anchors[s]);
     
     currentGap = extractWeakBorderSequence(m_sequence, whichStart, K, m_location);
     PATH_MAXLENGTH = (int)(1.2*length(currentGap)+2*K);
     
     currentAnchor = std::get<0>(anchors[s]);
     currentCount = m_coverage[whichStart];
     competingPaths.push_back(Trail(currentAnchor, currentCount));

     if( m_direction == RIGHT )
      {
		currentRefSeq = currentAnchor;
		append(currentRefSeq, currentGap);
		competingPaths[0].setLeftAnchor(whichStart);
	  }
	 else
	  {
	   currentRefSeq = currentGap;
	   append(currentRefSeq, currentAnchor);
	   competingPaths[0].setRightAnchor(whichStart);
	  }
     
  while( (!empty(competingPaths)) & (competingPaths.size()<=p_MAX_NB_OF_INNER_PATHS) & (stepCounter<PATH_MAXLENGTH))
   {	 
     this->oneMoreStepInTheDark(xdrop, currentRefSeq, competingPaths, stepCounter, m_direction , PATH_MAXLENGTH, SR_DBG, io_JFcountsTable, io_JFjunctcountsTable);                                  
   }
 }//END-FOR on selected start points
   
 if( !empty(m_shortPaths) || !empty(m_longPaths) )
   {
     winner = this->sortOutBestBorder();
     diff = (double)this->getWeakLength()-(double)winner.getLength();
   
  if( (this->getWeakLength()>=300) || m_complexRegion ) minScore = std::max(0.75, gp_MIN_BORDER_SCORE);
  else minScore = gp_MIN_BORDER_SCORE;
  
  if( (diff<this->getWeakLength()*0.05 || ((this->getWeakLength()<6) & (winner.getLength()<6))) & (winner.getIDScore()>=minScore) )
   {   
	pathHasBeenFound = true;
    std::get<0>(m_weakSequence) = winner.getSeq();  
	std::get<1>(m_weakSequence) = CORRECTED; 
    if( m_location == TAIL) std::get<1>(m_LEFT_KMpositions) = winner.getLeftAnchor();
	else std::get<0>(m_RIGHT_KMpositions) = winner.getRightAnchor();
   }
  //else std::cout << "[EDGE] Paths found were not convenient: " << winner.getIDScore() << " - " << diff <<  std::endl;
 } 
 
  return pathHasBeenFound; 
}
 

/************************************************************************
 * Behaviour during the exploration of borders
 * followBorderPath()
 * 
 * The function iterates as long as there is only one possible base to follow and the path is not too long.
 * otherwise until it is faced to a DEAD-END (0) or a BRANCHING POINT (>=2)
 *    -> in that case if there is no GPred base to follow, STOP EXPLORATION
 *    -> otherwise record all the ways that are available
 *    -> first follow the LPred ones if exist, then the GPred ones
 *       then possibly the (most frequent first) unpred ones
 ***********************************************************************/
 
 
void Explorer::recordBridge(Trail& trail)
{
 Trajectory myBridge = Trajectory(trail);
 m_fullPaths.push_back(myBridge);
}

void Explorer::recordEdge(Trail& trail, TSeq& reference)
{
 bool shorter(true);
 Trajectory myTip = Trajectory(trail);

 myTip.trim(K, p_CHECK_INTERVAL, trail.getNbFailuresInARow(), m_direction);

 shorter = (myTip.getLength()<=length(reference));
 myTip.reshape(reference, K, m_direction, shorter);
 if( myTip.cutAnchors(m_location, 0, K) )
  {
   if( shorter ) m_shortPaths.push_back(myTip);
   else m_longPaths.push_back(myTip);
  }
  else if(DEBUG_USER) std::cout << "[DEBUG EDGE] Trajectory ignored due to overlap" << std::endl;
}

/************************************************************************/

void trim(TSeq& history, int nbBases, Direction direction)
{
	TSeq trimmedSeq;
	if(direction == RIGHT) trimmedSeq = prefix(history, length(history)-nbBases-1);
	else trimmedSeq = suffix(history, nbBases);
	
	history = trimmedSeq;
}


TSeq extendGap(TSeq& gap, TSeq& anchor, Location& location)
{
  TSeq extendedGap;
  if(location == HEAD) 
   {
	 extendedGap = gap;
	 append(extendedGap, anchor);   
   } 
  else 
   {
	 extendedGap = anchor;
	 append(extendedGap, gap);
   }
 return extendedGap;
}

TSeq extendINNERGap(TSeq& gap, TSeq& anchor, TSeq& target, Direction& direction)
{
	TSeq extendedGap;
	if( direction == RIGHT )
	 {
	   extendedGap = anchor;
	   append(extendedGap, gap);
	   append(extendedGap, target); 
	 }
	else
	 {
	  extendedGap = target; 
	  append(extendedGap, gap);
	  append(extendedGap, anchor); 
	 } 
  return extendedGap;
}


void Explorer::display()
{
	std::cout << "******************* Explorer info ************************" << "\n"
		      << "Sequence length: " << length(m_sequence) << "\n"
	          << "Coverage size: " << m_coverage.size() << "\n"
	          << "Location: " << m_location << "\n"
	          << "Direction: " << m_direction << "\n"
	          << "LEFT positions: " << std::get<0>(m_LEFT_KMpositions) << "-" << std::get<1>(m_LEFT_KMpositions) << " -> " << std::get<2>(m_LEFT_KMpositions) << "\n"
	          << "RIGHT positions: " << std::get<0>(m_RIGHT_KMpositions) << "-" << std::get<1>(m_RIGHT_KMpositions) << " -> " << std::get<2>(m_RIGHT_KMpositions) << "\n"
	          << "Nb LEFT anchors: " << m_LEFT_anchors.size() << "\n"
	          << "Nb RIGHT anchors: " << m_RIGHT_anchors.size() << "\n"
	          << "Nb recorded bridges " << m_fullPaths.size() << "\n"	          
	          << "********************************************************" << "\n"
	          << std::endl;       	          
}


//returns true if nextc falls withing the area defined as "EXPECTED" or "NOISE" by myself
bool isExpectedbyMyModel(unsigned int nextc, unsigned int cc, double& broad_coeff, Status classe)
{
  //if count <=3 : Poisson CI bounds computed following WALD CC method
  // [Z_1-alpha/2-sqrt(x_-0.02) ;  ( Z_1-alpha/2+Z_+sqrt(x_0.96)]^2)
  // SMALLISH COUNTS
  if( (cc<=3) & (classe==UNEXPECTED) ) return ((double)nextc<=((double)(cc+0.5)+gp_ALPHA*sqrt((double)(cc+0.5))));
  else if(  (cc<=3) & (classe==EXPECTED) ) return ((double)nextc>=((double)(cc-0.5)+(1-gp_ALPHA)*sqrt((double)(cc-0.5))));
  
  //Begaud's method
  //[(sqrt(x_0.2)+Z_alpha1/2)² ; (sqrt(x_0.96)+Z_alpha2/2)²
  else if( (cc>3) & (classe==UNEXPECTED) ) return ((double)nextc<=std::pow((gp_ALPHA/2+sqrt((double)(cc+0.96))),2));
  else return ((double)nextc>=std::pow((gp_ALPHA/2-sqrt((double)(cc+0.02))),2)); 
  //else if( (cc>3) & (classe==EXPECTED) ) return ((double)nextc>=std::pow((gp_ALPHA/2-sqrt((double)(cc+0.02))),2)); 
}

bool isExpectedbyMyLastNode(unsigned int nextc, unsigned int cc, double& broad_coeff)
{
  bool isExpected(true);	
  if( (cc<=3) ) 
   {
	isExpected &= ((double)nextc<=((double)(cc+0.5)+gp_ALPHA*sqrt((double)(cc+0.5))));
	isExpected &= ((double)nextc>=((double)(cc-0.5)+(1-gp_ALPHA)*sqrt((double)(cc-0.5))));
   }
  
  //Begaud's method
  //[(sqrt(x_0.2)+Z_alpha1/2)² ; (sqrt(x_0.96)+Z_alpha2/2)²
  if( (cc>3) )
   {
	isExpected &= ((double)nextc<=std::pow((gp_ALPHA/2+sqrt((double)(cc+0.96))),2));
	isExpected &= ((double)nextc>=std::pow((gp_ALPHA/2-sqrt((double)(cc+0.02))),2)); 
	} 
 return isExpected;
}


Explorer make_empty_Explorer()
{
	Explorer explorer;
	return explorer;
}

void tagNextNodes(std::vector<std::pair<Status, double > >& nodeTags, std::vector<colouredCount >& nextCounts, unsigned int count, double& broad_coeff, double& priorNoiseLevel, bool& complex)
{
  //Counting # of non zero neighbour k-mers
  clear(nodeTags);
  int counter(0);
  double dist(0);
  unsigned int nextc(0);
  unsigned int lambda_noise;
  unsigned int nbExpected(0);
  unsigned int nbBreakpoints(0);
  unsigned int nbUnexpected(0);
    
    for(unsigned int i(0); i < nextCounts.size(); i++) if( (int)std::get<0>(nextCounts[i])>=gp_MIN_COUNT) counter++;
	//If there are counts>gp_MIN_COUNT
    if(counter>0)
    {		
	 lambda_noise = (int)((double)count*gp_SR_ERROR_RATE);   

	for(unsigned int b(0); b < nextCounts.size(); b++)
      {  
	   nextc = std::get<0>(nextCounts[b]);
	   dist = abs((double)count-(double)nextc)/sqrt((double)count);
	   if( nextc>=gp_MIN_COUNT )
	    {
		 //if predicted or is the only possibility
         if( isExpectedbyMyModel(nextc, count, broad_coeff, EXPECTED) || (counter == 1) ) 
         {
		  nodeTags.push_back(std::make_pair(EXPECTED, dist));
		  ++nbExpected;
	     }
         else if( (lambda_noise>=gp_MIN_COUNT))
          {
		   if( !isExpectedbyMyModel(nextc, lambda_noise, broad_coeff, UNEXPECTED) || (std::get<1>(nextCounts[b])>0))
		    {
			 nodeTags.push_back(std::make_pair(BREAKPOINT, dist));
			 ++nbBreakpoints;
		    }
		   else
		    {
			 nodeTags.push_back(std::make_pair(UNEXPECTED, dist));
			 ++nbUnexpected; 
		    }
          }
         else
          {
		   nodeTags.push_back(std::make_pair(BREAKPOINT, dist));
		   ++nbBreakpoints;
	      }
	     }//END-IF next count in SR-DBG
	    else nodeTags.push_back(std::make_pair(UNEXPECTED, dist));
       }//END-FOR all non zero bases
    }//END-IF AT LEAST ONE BASE WITH COUNT>bound 
if( (nbExpected==0) & (nbBreakpoints==1))
 {
  for(unsigned int t(0); t<nodeTags.size();t++) if(std::get<0>(nodeTags[t])==BREAKPOINT) std::get<0>(nodeTags[t])=EXPECTED;
 }   	 
if( (nbExpected==1) & (nbUnexpected>0) & !complex) //check sum of unexpected nodes is below noise threshold
 {
   counter=0;
   unsigned int index(0);
   
   for(unsigned int i(0); i<nextCounts.size(); i++)
    {
	 if(std::get<0>(nodeTags[i])==UNEXPECTED)
	  {
	   if(counter==0) index = i;
	   counter += std::get<0>(nextCounts[i]);
	   if( std::get<0>(nextCounts[index]) < std::get<0>(nextCounts[i]) ) index=i;
	  }
    }
   if( !isExpectedbyMyModel(counter, lambda_noise, broad_coeff, UNEXPECTED) ) std::get<0>(nodeTags[index]) = BREAKPOINT;
  }
}
