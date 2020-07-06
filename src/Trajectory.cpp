//
//  Trajectory.cpp

//  Created by Lucile Broseus.
//  Copyright © 2018 LuB. All rights reserved.
//

#include "Trajectory.hpp"

#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "utils.hpp"
#include "Trail.hpp"

using namespace seqan;
  
typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;

extern bool DEBUG;
extern bool DEBUG_TEST;

Trajectory::~Trajectory()
{
}

//Constructors
	
Trajectory::Trajectory() : m_sequence(make_empty_Dna5String()), 
                           m_leftAnchor(0), m_rightAnchor(0), m_outcome(NO), m_direction(RIGHT), m_lastScore(0), m_score(-200000), m_IDscore(-1), m_nbBreakpoints(0), m_distance(0)
{
	
}
      
Trajectory::Trajectory(Trail& trail):
                       m_sequence(trail.getSeq()), 
                       m_leftAnchor(trail.getLeftAnchor()), m_rightAnchor(trail.getRightAnchor()) , m_outcome(NO), m_direction(RIGHT), m_lastScore(trail.getLastScore()), m_score(-200000), m_IDscore(0), m_nbBreakpoints(trail.getNbBreakpoints()), m_distance(trail.getDistance()/(trail.getLength()+0.01))
{
	
}

                       
                    
//Methods

TSeq Trajectory::getSeq()
{
	return m_sequence;
}
    
void Trajectory::setSeq(TSeq seq)
{
	m_sequence = seq;
}

unsigned int Trajectory::getLength()
{
	return length(m_sequence);
}
           
unsigned int Trajectory::getLeftAnchor()
{
	return m_leftAnchor;
}
    
void Trajectory::setLeftAnchor(unsigned int lanchor)
{
	m_leftAnchor = lanchor;
}
    
unsigned int Trajectory::getRightAnchor()
{
	return m_rightAnchor;
}
    
void Trajectory::setRightAnchor(unsigned int ranchor)
{
	m_rightAnchor = ranchor;
}

void Trajectory::trim(unsigned int& minSize, unsigned int& intervalLength, unsigned int nbFailuresInARow, Direction& direction)
{
  unsigned int nbBases(0);
  int position(0);
  TSeq newEdge;
  nbBases = nbFailuresInARow * intervalLength;
  
  if(this->getLength()>=nbBases+minSize)
   {
	if( direction == RIGHT )
	 {
	  position = this->getLength()-nbBases;
	  newEdge = prefix(this->getSeq(), position);
     }
    else
     {
	  position = nbBases;
	  newEdge = suffix(this->getSeq(), position);
	 }
   }
  else newEdge = this->getSeq();
  
  m_sequence = newEdge;
}

void Trajectory::reshape(TSeq& reference, unsigned int& kmerSize, Direction& direction, bool& shorter)
{
 std::tuple<TSeq, TSeq, int, double > extension_results;
 TSeq newSeq;
 TSeq tmp;
 double id_score;
 int xdrop1(0);
 
 //decrease xdrop1 in extension until just before extGapSeq reduction
 //0 -> SR-path
 //1 -> reference sequence
 tmp = this->getSeq();
 xdrop1 = (int)this->getLastScore()*(-1);

 if( !shorter ) //extend weak reference on path sequence
  {
   extension_results = findStopPosition(tmp, reference, xdrop1, direction, kmerSize);
   if( direction == LEFT ) newSeq = suffix(tmp, std::get<2>(extension_results));
   else newSeq = prefix(tmp, std::get<2>(extension_results)); 
  }
 else
  {
   extension_results = findStopPosition(reference, tmp, xdrop1, direction, kmerSize);
   if( direction == LEFT )
    {
	 newSeq = prefix(reference, std::get<2>(extension_results));
	 append(newSeq, tmp);
    }
   else
    {
	  tmp = suffix(reference, std::get<2>(extension_results));
	  newSeq = this->getSeq();
	  append(newSeq, tmp);
	} 	  
  }
 
 id_score = computePercentID(std::get<0>(extension_results), std::get<1>(extension_results));
 this->setIDScore(id_score);
 this->setScore(std::get<3>(extension_results));
 this->setSeq(newSeq);
 
}

bool Trajectory::cutAnchors(Location location, unsigned int limit, unsigned int kmerSize)
{
  TSeq truncSeq;
  TSeq oldSeq = this->getSeq();
  bool isOK(true);
  unsigned int len = this->getLength();

  try
   {
        switch( location )
         {
                 case HEAD:  //remove the last kmerSize bases
                  if(len>kmerSize) truncSeq = prefix(oldSeq, len-kmerSize);
if(length(truncSeq)!=len-kmerSize) std::cout << "Could not unanchor HEAD" << std::endl;           
                 break; //remove the first kmerSize bases                         
                 case TAIL:
                  if(len>kmerSize) truncSeq = suffix(oldSeq, kmerSize);
if( (len>kmerSize) & (length(truncSeq)!=len-kmerSize) ) std::cout << "Error when unanchoring TAIL" << std::endl;        
                 break;                           
                 case INNER:
                  if(len>=2*kmerSize)
                   {
                        truncSeq = infix(oldSeq, kmerSize, len-kmerSize);
if(length(truncSeq)!=(len-2*kmerSize)) std::cout << "Could not unanchor INNER" << std::endl;                            
               }
              else if( (len<2*kmerSize) & len>kmerSize )
               {
                 //std::cout << "ANCHOR OVERLAP " << len << " " << this->getRightAnchor() << std::endl; 
                 if(this->getRightAnchor()+2*kmerSize-len<=limit)
                  { 
                   this->setRightAnchor(this->getRightAnchor()+2*kmerSize-len);
                   //std::cout << "New right position " << this->getRightAnchor() << std::endl;  
			      }
			     else
			      {
				   // std::cout << "[TRAJ DISCARDED] ANCHOR OVERLAP " << len << " " << this->getRightAnchor() << " " << limit << std::endl; 
				   isOK = false; 
			      }   
               }
               else isOK = false;
              break;    
        }
   }
  catch(std::string const& chaine)
   {
        std::cout << chaine << std::endl;
        std::cout << "Input History: " << this->getLength() << std::endl;
        std::cout << "After truncation: " << truncSeq << std::endl;
    truncSeq = "N";
    isOK = false;
   }    
 this->setSeq(truncSeq);  
 
 return isOK;                  
}


Direction Trajectory::getDirection()
{
	return m_direction;
}

double Trajectory::getLastScore()
{
	return m_lastScore;
}

double Trajectory::getScore()
{
	return m_score;
}
    
void Trajectory::setScore(double score)
{
	m_score = score;
}

void Trajectory::setIDScore(double score)
{
	m_IDscore = score;
}

void Trajectory::scoreSequence(TSeq& reference)
{
  m_score = computeEditDistance(reference, m_sequence);
  m_IDscore = computeIDScore(reference, m_sequence)/std::max(length(reference), length(m_sequence));
}

//++: Add anchors before scoring

void Trajectory::IDscore(TSeq& reference)
{
    double score;
	int len = std::max(length(reference), length(m_sequence));
	
	score = computeIDScore(reference, m_sequence);
	m_IDscore = score/len; 	    
}

double Trajectory::getIDScore()
{
	return m_IDscore;
}
    
void Trajectory::overlapScore(TSeq refGap, int mode)
{
	if(length(refGap)>0) m_score = overlapAlign(refGap, m_sequence, mode)/length(refGap);
	else m_score = -1;
}

int Trajectory::getNbBreakpoints()
{
 return m_nbBreakpoints;
}

double Trajectory::getMeanDistance()
{
 return m_distance;
}
      
//Other

//Check if there are ex aequo path, then decide upon EditDistance?
//Use total number of breakpoint path chunks
//Then path are somehow ordered by likelihood (min abs distance count)
unsigned int findBestBridge(std::vector<Trajectory >& trajectories)
{
  unsigned int index(0);
  std::vector<unsigned int > exAequo1;
  for(unsigned int i(1); i<trajectories.size(); i++ )
   {
    if(trajectories[i].getScore()>trajectories[index].getScore()) index=i;
   } 
//If there are several traj. with same minimal edit distance, classify by desc. # of hits    
  for(unsigned int i(index+1); i<trajectories.size(); i++ )
   {
    if(trajectories[i].getScore()==trajectories[index].getScore()) exAequo1.push_back(i);
   }  
  if(!empty(exAequo1))
   {
	for(unsigned int i(0); i<exAequo1.size(); i++ )
     {
	  if(trajectories[exAequo1[i]].getMeanDistance()>trajectories[index].getMeanDistance() ) index=exAequo1[i];
     } 
   }   
    return index;	
}

//Returns whether there is a long deadend and if so, the index of the best one in the input vector.
std::pair<bool, unsigned int> findBestBORDER(std::vector<Trajectory >& trajectories, double MIN_SCORE)
{
//std::cout << "[FIND BEST EDGE]: " << trajectories.size() << std::endl;
  unsigned int index(0);
  std::vector<unsigned int > exAequo1;
  bool isConvenient = true;
  
  if(empty(trajectories)) isConvenient = false;
  else{
  for(unsigned int i(1); i<trajectories.size(); i++ )
   {
    if(trajectories[i].getScore()>trajectories[index].getScore()) index=i;
   } 
  
//If there are several traj. with same minimal edit distance, classify by desc. # of hits    
  for(unsigned int i(index+1); i<trajectories.size(); i++ )
   {
    if(trajectories[i].getScore()==trajectories[index].getScore()) exAequo1.push_back(i);
   }  
  if(!empty(exAequo1))
   {
	for(unsigned int i(0); i<exAequo1.size(); i++ )
     {
	  if(trajectories[exAequo1[i]].getMeanDistance()>trajectories[index].getMeanDistance()) index=exAequo1[i];
     }
    }  
 }
   return std::make_pair(isConvenient, index);
}

//Returns whether there is a short deadend and if so, the index of the best one in the input vector.
double computeIDScore(TSeq& gap, TSeq& history)
{
//--------------------------------------------------------------------->DEBUG
//std::cout << "FUNCTION COMPUTEIDSCORE" << std::endl;
//--------------------------------------------------------------------->DEBUG	

   double score(-1);
   Align<TSeq > align;
   resize(rows(align), 2);
    
   if(length(gap)==0) std::cout << "PB3" << std::endl;
   if(length(history)==0) std::cout << "PB32" << std::endl;
    
   unsigned int len1 = length(gap);
   unsigned int len2 = length(history);
   unsigned int max_len = std::max(len1,len2);
  
   if( (len1>0) & (len2>0) )
   {
    if( len1>len2 )
     {
      assignSource(row(align, 0), gap);
      assignSource(row(align, 1), history);
     }
    else 
     {
      assignSource(row(align, 0), history);
      assignSource(row(align, 1), gap);
    }     
    try
      {
        score = localAlignment(align, Score<int, Simple>(1, 0, 0), LinearGaps());
      }  
     catch(std::string const& chaine)
     {
	  std::cout << chaine << std::endl; 
	  score = -1;
	  std::cout << "Attributed score: " << score << std::endl;                          
     }                        
   }
    else
    {
	 score = -1;
	 std::cout << "Attributed score: " << score << std::endl; 
	}

     return score;
}

double computeEditDistance(TSeq& reference, TSeq& history)
{
//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION COMPUTEEDITDISTANCE" << std::endl;
//--------------------------------------------------------------------->DEBUG	
	
   TStringSet sequences;
   double score(-100000);
    
   if(length(reference)==0) std::cout << "PB41" << std::endl;
   if(length(history)==0) std::cout << "PB42" << std::endl;
    
   if( (length(history)>0) & (length(reference)>0) )
    {  
     if(length(reference)>=length(history))
     {
      appendValue(sequences, reference);
      appendValue(sequences, history);
     }
     else
     {
      appendValue(sequences, history);
      appendValue(sequences, reference);
     }
     TAlignGraph alignG(sequences);
     try
      {
       score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1));                       
      }                             
     catch(std::string const& chaine)
      {
	   std::cout << chaine << std::endl; 
	   std::cout << "Attributed score: " << score << std::endl;                          
      }
	 }
	else
    {
	 score = (-100000);
	 std::cout << "Attributed score: " << score << std::endl; 
	}

    return score;
}

double computeEditDist(TSeq& reference, TSeq& history, int mode)
{
//--------------------------------------------------------------------->DEBUG
// std::cout << "FUNCTION COMPUTEEDITDISTANCE" << std::endl;
//--------------------------------------------------------------------->DEBUG	
	
   TStringSet sequences;
   double score(-100000);
    
   if(length(reference)==0) std::cout << "PB41" << std::endl;
   if(length(history)==0) std::cout << "PB42" << std::endl;
    
   if( (length(history)>0) & (length(reference)>0) )
    {  
     if(length(reference)>=length(history))
     {
      appendValue(sequences, reference);
      appendValue(sequences, history);
     }
     else
     {
      appendValue(sequences, history);
      appendValue(sequences, reference);
     }
     TAlignGraph alignG(sequences);
     try
      {
        if(mode == 1) score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1)); 
        else score = globalAlignment(alignG, Score<int, Simple>(1, -1, -1));                        
      }                             
     catch(std::string const& chaine)
      {
	   std::cout << chaine << std::endl; 
	   std::cout << "Attributed score: " << score << std::endl;                          
      }
	 }
	else
    {
	 score = (-100000);
	 std::cout << "Attributed score: " << score << std::endl; 
	}

    return score;
}

std::vector<Trajectory > make_empty_vectOfTraj()
{
	std::vector<Trajectory > traj;
	return traj;
}


std::tuple<TSeq, TSeq, int, double > findStopPosition(TSeq& reference, TSeq& shorterPath, int xdrop, Direction& direction, unsigned int& kmerSize)
{
  int j(0);
  int xdrop1(xdrop);
  bool goFurther(true);
  double score(0);
  std::tuple<TSeq, TSeq, int, double, bool> extension_results;
  std::tuple<TSeq, TSeq, int, double, bool> new_extension_results;
  
  new_extension_results = getSeedAndExtension(reference, shorterPath, xdrop1, direction, kmerSize);

  do
   {
	 ++j;
	 --xdrop1;
	 extension_results = new_extension_results;
	 new_extension_results = getSeedAndExtension(reference, shorterPath, xdrop1, direction, kmerSize);
	 if( length(std::get<1>(new_extension_results))<length(std::get<1>(extension_results)) ) goFurther = false;
   }while( goFurther & (xdrop1>0) );

  return std::make_tuple(std::get<0>(extension_results), std::get<1>(extension_results), std::get<2>(extension_results), std::get<3>(extension_results));	    
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
