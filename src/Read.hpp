//
//  Read.hpp
//
//  
//

#ifndef Read_hpp
#define Read_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include "utils.hpp"

#include "Settings.hpp"
#include "Explorer.hpp"
#include "Trajectory.hpp"

typedef seqan::Dna5String TSeq;
typedef std::pair<unsigned int, unsigned int > colouredCount;
typedef std::map<TSeq, colouredCount > colouredDBG;

class Read
{

public:

   ~Read();

    //Constructors
    
    Read();
    
    Read(seqan::CharString& id);

    Read(seqan::CharString& id, TSeq& sequence);
    
    //Methods
    
    seqan::CharString getName();
    TSeq getSeq();
    
    TSeq getCorrSeq();
    void setCorrSeq(TSeq& corrSeq);
    
    seqan::CharString getScorr();
    void setScorr(seqan::CharString& scorr);
    void updateScorr(seqan::CharString scorrNewPart);
    
    int getLength();
    int getCorrLength();
    
    bool reCoverage();
    std::vector<colouredCount > getCoverage();
    
    void setPriorNoise(double lambda);
        
    void display();
    
    void displaySegmentation();
       
    void updateCorrSeq();
    
    void correct2();
          
    void outputBasicReadStats(std::string& outFileName);
    
    TSeq concatenateNewStructure();
    bool setInitialStructure();
    bool defineStructure2();
    //void update2(Explorer myExplorer, int& reg, Location location);
    std::pair<TSeq, Status > getSolidRegion(std::tuple<unsigned int, unsigned int, Status > coordinates);  
    
    void updateINNER(Explorer& myExplorer, int& reg);
    void updateHEAD(Explorer& myExplorer);
    void updateTAIL(Explorer& myExplorer);
    void setMoreReadStatsHeader(std::string& outFileName);
    //void outputMoreReadStats(std::string& outFileName);
    
private:

    seqan::CharString m_id;                   //read name
    
	TSeq m_sequence;                          //raw read sequence
	TSeq m_correction;                        //corrected sequence
	seqan::CharString m_scorr;                //ID scores of the corrected sequence to the raw sequence
	
	std::vector<colouredCount > m_coverage;   //SR kmer count of the raw sequence
	
	double m_priorLambda_noise;
		
	std::vector<unsigned int > m_LRcoverage; //counts in the LR-dBG
	
	std::vector< std::tuple<unsigned int, unsigned int, Status > > m_InKmersPositions; 
		
	std::vector<std::pair<TSeq, Status> > m_newInnerStructure;
	std::pair<TSeq, Status> m_head;
	std::pair<TSeq, Status> m_tail;
	//Statistics
	
	int m_nbInKmers;
	int m_nbSolidKmers;
	
};

void setBasicReadStatsHeader(std::string& outFileName);

Read create_empty_Read();

//---------------------------------------------------------------------//
// SEPT 18: new functions
bool findINRegions(std::vector<std::tuple<unsigned int, unsigned int, Status > >& StartEndKmers, 
                   std::vector<colouredCount >& counts,
                   unsigned int& kmerSize);

double computeSeqErrorThreshold( std::vector<colouredCount >& counts);

void analyzeINRegions(std::vector<std::tuple<unsigned int, unsigned int, Status > >& StartEndKmers,
                      TSeq& refSequence,
                      std::vector<colouredCount >& counts, 
                      double& solidityThr, unsigned int& kmerSize);
                    
//---------------------------------------------------------------------//                                        
                       



#endif /* Read_hpp */
