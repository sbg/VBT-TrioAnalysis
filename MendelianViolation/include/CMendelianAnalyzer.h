//
//  CMendelianAnalyzer.h
//  VCFComparison
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef C_MENDELIAN_ANALYZER_H_
#define C_MENDELIAN_ANALYZER_H_


#include "SConfig.h"
#include "CMendelianVariantProvider.h"
#include "CMendelianResultLog.h"
#include "CMendelianTrioMerger.h"
#include "SChrIdTriplet.h"
#include "CMendelianDecider.h"
#include "ENoCallMode.h"
#include <thread>

namespace mendelian
{

class CMendelianAnalyzer
{
    
public:
    
    //Default Constructor
    CMendelianAnalyzer();
    
    //Executes the mendelian violation analyzer
    int run(int argc, char** argv);
    
private:
    
    //Read the parameters if the execution mode is mendelian. If all mandatory parameters are set, return true.
    bool ReadParameters(int argc, char** argv);
    
        
    //A Function to process mendelian violation pipeline for given chromosome id
    void ProcessChromosome(const std::vector<SChrIdTriplet>& a_rChromosomeIds);
        
    //Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
    int AssignJobsToThreads(int a_nThreadCount);
    
    //Prints the help menu at console
    void PrintHelp() const;
    
    //Config file for father-child comparison which is passed to m_fatherChildProvider
    SConfig m_fatherChildConfig;
    
    //Config file for mother-child comparison which is passed to m_motherChildProvider
    SConfig m_motherChildConfig;
    
    //No call mode of the variant
    ENoCallMode m_noCallMode;
    
    //Variant provider for parent-child comparison
    CMendelianVariantProvider m_provider;
    
    //Mendelian decider to assess mendelian decision of variants
    CMendelianDecider m_mendelianDecider;

    //Trio merger for outputing
    CMendelianTrioMerger m_trioWriter;
    
    //Result log for mendelian comparison
    CMendelianResultLog m_resultLog;
    
    //Best Paths written by each thread for each unique chromosome exists [Between father and child]
    std::vector<core::CPath> m_aBestPathsFatherChildGT;
    std::vector<core::CPath> m_aBestPathsFatherChildAM;
    
    //Best Paths written by each thread for each unique chromosome exists [Between mother and child]
    std::vector<core::CPath> m_aBestPathsMotherChildGT;
    std::vector<core::CPath> m_aBestPathsMotherChildAM;

};

}

#endif //C_MENDELIAN_ANALYZER_H_
