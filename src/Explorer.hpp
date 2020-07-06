//
//  Explorer.hpp
//
//--------------------------------------------------------------------//
//                        CLASS DESCRIPTION
//
// Class that manages the correction of weak regions
//
// Defined from two consecutive solid regions + the weak one
// in-between: SOURCE - GAP - TARGET
// Or for HEAD-TAIL : GAP - SOURCE / SOURCE - GAP
// Comes with functions that implement an exploration of the de Bruijn 
// graph from SOURCE to TARGET or from SOURCE to WE-DONT-KNOW-WHERE
//--------------------------------------------------------------------//

#ifndef Explorer_hpp
#define Explorer_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <queue>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "utils.hpp"
#include "Jellyfish.hpp"

#include "Trajectory.hpp"
#include "Settings.hpp"
#include "Trail.hpp"


typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::pair<seqan::Dna5, colouredCount > colouredBase; 
typedef std::map<TSeq, colouredCount > colouredDBG;
typedef std::map<TSeq, unsigned int > DBG;
typedef std::tuple<unsigned int, unsigned int, Status> kmerStretch; 
typedef std::tuple<TSeq, unsigned int, unsigned int > anchorTuple;


class Explorer
{
 public:
 
    ~Explorer();
    
    //Constructors
    
    Explorer();
    Explorer(TSeq refSequence, std::vector<colouredCount > coverage, double lambda);

    //Methods
    
    void reset();
    
    void setCoverage(std::vector<colouredCount > counts);
    void setSequence(TSeq sequence);
    
    std::pair<TSeq, Status > getWeakSeq();
	void setWeakSequence();
	unsigned int getWeakLength();
	
	kmerStretch getLEFTHandPositions();
	kmerStretch getRIGHTHandPositions();
	
	void initializeINNER(kmerStretch startKmPos,kmerStretch endKmPos, Direction direction); 
	void initializeHEAD(kmerStretch kmPos);
	void initializeTAIL(kmerStretch kmPos);
	
	
	void setLocation(Location location);

	Direction getDirection();
	void setDirection(Direction direction);
		
	seqan::CharString getScorr();
	void setScorr(seqan::CharString scorr);
	
//--------------------------------------------------------------------//
//               FUNCTIONS FOR DE BRUIJN GRAPH EXPLORATION
//
// FUNCTION findEdge():
//  wraps the whole graph exploration process for HEAD/TAIL correction
// FUNCTION findBridge():
//  wraps graph exploration for the correction of inner regions
//  call iterations of simple paths explorations with: follow*Path()
//
// FUNCTION followBorderPath():
// performs exploration of border paths whose nodes are not bifurcations
// FUNCTION followInnerPath():
// performs exploration of inner paths whose nodes are not bifurcations

// FUNCTION backToLastBranch():
//  applied when a path finished, to go back to the previous bifurcation
//  mode and possibly try another continuation
// FUNCTION recordBranch()
//  applied when a bifurcation node is met. 
//  Stores and classifies possible ways (bases on which to continue)
//--------------------------------------------------------------------//
    void setLEFTHandPositions(kmerStretch  positions);
    void setRIGHTHandPositions(kmerStretch positions);
    
    void anchorLEFTHandSide();
	void anchorRIGHTHandSide();
		
	void oneMoreStep(TSeq& reference, std::vector<Trail>& competingPaths, unsigned int& stepCounter, 
	                 Direction& direction, unsigned int& PATH_MAXLENGTH, 
                     colouredDBG& SR_DBG, std::string& io_JFcountsTable, std::string& io_JFjunctcountsTable);
   void oneMoreStepInTheDark(int& xdrop, TSeq& reference, std::vector<Trail>& competingPaths, unsigned int& stepCounter, Direction& direction, 
                           unsigned int& PATH_MAXLENGTH, colouredDBG& SR_DBG, std::string& io_JFcountsTable, std::string& io_JFjunctcountsTable);                  
                     
    bool searchBridge();
    bool searchEdge();
    
    void scoreEdges(int& xdrop, std::vector<Trail>& newCompetingPaths, unsigned int& stepCounter, TSeq& reference, Direction& direction);

    Trajectory sortOutBestBorder();

//--------------------------------------------------------------------//
//               FUNCTIONS FOR PATH EVALUATION
//
//--------------------------------------------------------------------//
	
    
    void recordBridge(Trail& trail);
    void recordEdge(Trail& trail, TSeq& reference);

	Trajectory getBestBridge();

//Functions intended to spot possible programming errors of Explorer objects
	void display();
         
 private:
 
    TSeq m_sequence;
    std::pair<TSeq, Status> m_weakSequence;
    std::vector<colouredCount > m_coverage; 
    double m_priorLambda_noise;
      
    std::tuple<unsigned int, unsigned int, Status> m_LEFT_KMpositions;
    std::tuple<unsigned int, unsigned int, Status> m_RIGHT_KMpositions;
    
    std::vector<anchorTuple > m_LEFT_anchors;
    std::vector<anchorTuple > m_RIGHT_anchors;
            
    Location m_location;                                              //Type of the GAP: either HEAD, TAIL or INNER
    Direction m_direction;                                            //Direction chosen for the exploration (LEFT: from right to left, RIGHT: from left to right)
        
    int m_nbGoodPathsFound;
    int m_nbVeryGoodPathsFound;
    bool m_complexRegion;
    
        
    std::vector<Trajectory > m_fullPaths;                             //Stores bridge paths extracted from the SR-dBG 
    std::vector<Trajectory > m_longPaths;                             //Stores long paths extracted from the SR-dBG
    std::vector<Trajectory > m_shortPaths;                            //Stores shorter extracted from the SR-dBG
    
    seqan::CharString m_scorr; 
};


                  
void tagNextNodes(std::vector<std::pair<Status, double > >& nodeTags, std::vector<colouredCount >& nextCounts, unsigned int count, double& broad_coeff, double& priorNoiseLevel, bool& complex);
                     
bool doABitOfGardening(std::vector<unsigned int >& indexOfKeptPaths, std::vector<Trail>& newCompetingPaths);

void scoreBridges(std::vector<Trail>& newCompetingPaths, unsigned int& stepCounter, TSeq& reference, Direction& direction);
                                      
std::pair<seqan::Dna5, colouredCount> getNextBase(int whichOne,
                                 std::vector<colouredBase >  ExpectedWays,
                                 std::vector<colouredBase >  Breakpoints,
                                 std::vector<colouredBase >  UnexpectedWays);
                                 
bool isExpectedbyMyModel(unsigned int nextc, unsigned int cc, double& broad_coeff, Status type);
bool isExpectedbyMyLastNode(unsigned int nextc, unsigned int cc, double& broad_coeff);
                             
                
//--------------------------------------------------------------------//
//               FUNCTIONS FOR BORDER PATH CONTROL
//--------------------------------------------------------------------//                


TSeq extendINNERGap(TSeq& gap, TSeq& anchor, TSeq& target, Direction& direction);

Explorer make_empty_Explorer();





//--------------------------------------------------------------------//
//               TO BE IMPLEMENTED
//--------------------------------------------------------------------// 

        
#endif /* Explorer_hpp */ 
