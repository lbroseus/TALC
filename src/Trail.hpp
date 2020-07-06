//
//  Trail.hpp

//  Created by Lucile Broseus on a chilly but sunny morning
//  Copyright © 2018 LuB. All rights reserved.


#ifndef Trail_hpp
#define Trail_hpp

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

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<seqan::Dna5, colouredCount > colouredBase; 
typedef std::pair<TSeq, colouredCount > tipNode; 

class Trail
{
public:

    ~Trail();
    
    //Constructors
    
    Trail();
        
    Trail(TSeq kmer,  colouredCount& ccount);
    
    Trail(Trail path, TSeq seq);
    
    Trail(Trail path, seqan::Dna5 newBase, Direction& direction, unsigned int count);
    
                   
    //Methods
    
    TSeq getSeq(); 
    void setSeq(TSeq seq);
    
    unsigned int getLength();
    
    void addBase(seqan::Dna5 base);
    
    tipNode getLastStep();   
    void setLastStep(tipNode lastStep);
    
    unsigned int getLastCount();
    TSeq getLastKmer();
    
    std::vector<colouredCount > whatsNext(Direction& direction, colouredDBG& dBG, std::string& countsTable, std::string& junctionCountsTable);

    double getLastScore();
    void setLastScore(double& score);
    void DPscore(TSeq& reference);
    void Overlapscore(TSeq& reference, Direction& direction);

    unsigned int getNbFailuresInARow();
    void recordFailure();
    void eraseFailures();
    
    bool seedAndExtend(TSeq& reference, Direction& direction, int xdrop, unsigned int MAX_FAILURES);

    int getNbBreakpoints();
    void setNbBreakpoints(int nb);
    void recordBreakpoint();
    
    double getDistance();
    void recordDistance(double dist);
     
    int getLeftAnchor();
    void setLeftAnchor(int position);

    int getRightAnchor();
    void setRightAnchor(int position);
    
    void setReachedAim(int aim, Direction& direction);
    
    
    bool ThinkIveAlreadyGotThere(TSeq history);
    bool checkAims(std::vector<anchorTuple > aims, Direction& direction);
        
    void displayPath();
            
private:
    
    TSeq m_sequence;       
    
    tipNode m_lastStep;
    double m_lastScore;
    unsigned int m_nbFailuresInARow;
                           
    int m_nbBreakpoints;
    double m_distance;
   
    int m_leftAnchor;  
    int m_rightAnchor;
    
};

TSeq addNewBase(TSeq sequence, seqan::Dna5 newBase, Direction& direction);

tipNode makeTipNode(TSeq kmer, seqan::Dna5 base, Direction& direction, unsigned int count);

std::tuple<TSeq, TSeq, int, double, bool> getSeedAndExtension(TSeq& reference,TSeq& candidate, int& xdrop, Direction& direction, unsigned int& seedSize);


#endif /* Trail_hpp */

