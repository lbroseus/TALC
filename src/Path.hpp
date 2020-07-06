//
//  Path.hpp

//  Created by Lucile Broseus on a rather sunny afternoon
//  Copyright © 2018 LuB. All rights reserved.
//

#ifndef Path_hpp
#define Path_hpp

#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <utility>
#include <algorithm>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seeds.h>

#include "utils.hpp"
//#include "Trajectory.hpp"

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<seqan::Dna5, colouredCount > colouredBase; 

class Path
{
public:

    ~Path();
    
    //Constructors
    
    Path();
    
    Path(int whichStart); 
    
    Path(TSeq seq, int whichStart);
    
    Path(TSeq seq, int whichStart, Status status);
       
    Path(TSeq seq, int whichStart, int whichAim);
    
    Path(Path path, TSeq seq);
    
    Path(Path path, TSeq seq, Status newStatus);
    
    //Path(Trajectory trajectory);
           
    //Methods
    TSeq getSeq(); 
    void setSeq(TSeq seq);
    
    int getLength();
    
    void addBase(seqan::Dna5 base);
    
    TSeq getLastStep();   
    void setLastStep(TSeq lastStep);
    
    double getLastScore();
    void setLastScore(double score);
    
    unsigned int getLastHitLandmark();
    void setLastHitLandmark(unsigned int nb);
    
    std::vector<std::pair<seqan::Dna5, colouredCount > > getWhatsNext();  
    void setWhatsNext(std::vector<std::pair<seqan::Dna5, colouredCount > > outPaths);
    
    std::vector<colouredBase > getExpectedWays();
    std::vector<colouredBase > getBreakpoints();
    std::vector<colouredBase > getUnexpectedWays();
    
    void setExpectedWays(std::vector<colouredBase > ExpectedWays);
    void setBreakpoints(std::vector<colouredBase > Breakpoints);
    void setUnpredBases(std::vector<colouredBase > UnexpectedWays);
    
    void pruneUnpredBases();
           
    unsigned int getNbVisits();    
    void setNbVisits(unsigned int nbVisits);    
    void recordVisit();
    
    Status getStatus();
    void setStatus(Status status);
    void update(seqan::Dna5 next_base, Status status, Direction direction);

    int getNbFailures();
    void recordFailure();
    void eraseFailure(int mode);
    
    int getNbFailuresInARow();
    
    int getXdropenalty();
    void setXdrop(int& xdrop);
    void updateXdropenalty(int pen);
    
    int getNbBreakpoints();
    void recordBreakpoint();
    
    int getNbUnpred();
    void recordUnpred();
     
    int getStart();
    void setStart(int start);

    int getReachedAim();
    void setReachedAim(int aim);
    
    bool seemsOK(TSeq& weakSeq, TSeq& history, TSeq& currentAnchor,Direction& direction, unsigned int& kappa, double& failRate);
    
    bool checkInnerPath(unsigned int kmerSize);    
    bool checkInitPath(unsigned int kmerSize);
    bool checkNewBranch();
    
    void displayPath();
            
private:
    
    TSeq m_sequence;       
    
    TSeq m_lastStep;
    double m_lastScore;
    unsigned int m_lastHitLandmark;
                           
    int m_nbFailures;
    int m_nbFailuresInARow;
    int m_xdropenalty;
    
    int m_nbBreakpoints;
    int m_nbUnpredBranch;   
    
    std::vector<colouredBase > m_ExpectedWays;
    std::vector<colouredBase > m_Breakpoints;
    std::vector<colouredBase > m_UnexpectedWays;
    
    Status m_status;
    
    unsigned int m_nbVisits;
  
    int m_start;  
    int m_reachedAim;
    
};

TSeq catPaths(std::vector<Path >& Nodes, Direction& direction);

double performGlobalAlignment(TSeq seq1, TSeq seq2);
double performSemiGlobalAlignment(TSeq seq1, TSeq seq2, Direction direction);
double performLocalAlignment(TSeq seq1, TSeq seq2);
double computePercentID(TSeq seq1, TSeq seq2);
std::tuple<TSeq, TSeq, int > findWhereItEnds(TSeq& extGapSeq, TSeq& shortDEADENDpath, int xdrop1, Direction& direction, unsigned int& kmerSize);
//Outputs: TSeq: extension and int position where extension stops on extGapSeq 
std::tuple<TSeq, TSeq, int, double> getExtension(TSeq& extGapSeq,TSeq& history, int& xdrop, Direction& direction, unsigned int seedSize, int mode);


#endif /* Path_hpp */

