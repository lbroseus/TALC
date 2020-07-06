//
//  Trajectory.hpp
//
//--------------------------------------------------------------------//
//                        CLASS DESCRIPTION
//
// Stores information needed to define a complete path
// -> ie each successful grpah exploration results in a BRIDGE between SOURCE and TARGET
//    or in a candidate EDGE whose characteristics are stored as a TRAJECTORY object.
// Used by Explorer objects.
//--------------------------------------------------------------------//

#ifndef Trajectory_hpp
#define Trajectory_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "utils.hpp"
#include "Trail.hpp"

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;

class Trajectory
{
public:

    ~Trajectory();

	//Constructors
	
	Trajectory();
      	    
    Trajectory(Trail& trail);
 
    //Methods

    TSeq getSeq();   
    void setSeq(TSeq seq);
    
    TSeq getExtension();
    void setExtension(TSeq extension);
    
    unsigned int getLength();
          
    unsigned int getLeftAnchor();
    
    void setLeftAnchor(unsigned int lanchor);
    
    unsigned int getRightAnchor();
    
    void setRightAnchor(unsigned int ranchor);
    
    void trim(unsigned int& minSize, unsigned int& intervalLength, unsigned int nbFailuresInARow, Direction& direction);
    void reshape(TSeq& reference, unsigned int& kmerSize, Direction& direction, bool& shorter);

    bool cutAnchors(Location location, unsigned int limit, unsigned int kmerSize);
    
    Outcome getOutcome();
    
    Direction getDirection();
    
    double getLastScore();

    double getScore();   
    void setScore(double score);
    void IDscore(TSeq& reference);
    double getIDScore();
    void setIDScore(double score);
    void scoreSequence(TSeq& reference);

    void overlapScore(TSeq refGap, int mode);
      
    int getNbBreakpoints();
     
    double getMeanDistance();

private:

    TSeq m_sequence;          
  
    unsigned int m_leftAnchor;         //index of the leftAnchor
    unsigned int m_rightAnchor;        //index of the rightAnchor
    
    double m_score;           //Score used to compare this trajectory to the others; the way scores are computed is a parameter
    double m_IDscore;         //ID score of m_sequence to the raw sequence (that was to be corrected)
    
    int m_nbBreakpoints;       //Number of nodes in the trajectory detected as change point in the model
    double m_lastScore;
    double m_distance;
    
    Outcome m_outcome;           //whether the trajectory is a DEADEND, has been ABORTED or is a bridge (SUCCESS)
    Direction m_direction;
        
};

std::tuple<TSeq, TSeq, int, double > findStopPosition(TSeq& reference, TSeq& shorterPath, int xdrop1, Direction& direction, unsigned int& kmerSize);

//returns -1 when no successful Trajectories could be found...
unsigned int findBestBridge(std::vector<Trajectory >& trajectories);
std::pair<bool, unsigned int> findBestBORDER(std::vector<Trajectory >& BORDERpaths, double MIN_SCORE);

double computeIDScore(TSeq& gap, TSeq& history);  

double computeEditDistance(TSeq& reference, TSeq& history);
double computeEditDist(TSeq& reference, TSeq& history, int mode);
double computePercentID(TSeq seq1, TSeq seq2);


std::vector<Trajectory > make_empty_vectOfTraj();


#endif /* Trajectory_hpp */
