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
#include "CPath.h"
#include <thread>
#include "CResultLog.h"
#include "SMendelianThreadParam.h"

class CMendelianAnalyzer
{
    
public:
    
    //Executes the mendelian violation analyzer
    void run(int argc, char** argv);
    
    //Read the parameters if the execution mode is mendelian. If all mandatory parameters are set, return true.
    bool ReadParameters(int argc, char** argv);
    
    //A thread function to process best path algorithm in allele match with the given parameters
    //void ThreadFunc(const std::vector<SMendelianThreadParam>& a_rThreadParameters);
    
    //A function that perform merge operations after best path algorith m for mother-child and father-child is called
    //void MergeFunc(int a_nChromosomeId);

    void MergeFunc(int a_nChromosomeId);

    //A Function to process mendelian violation pipeline for given chromosome id
    void ProcessChromosome(const std::vector<int>& a_rChromosomeIds);
    
private:
    
    //Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
    int AssignJobsToThreads(int a_nThreadCount);
        
    //Config file for father-child comparison which is passed to m_fatherChildProvider
    SConfig m_fatherChildConfig;
    
    //Config file for mother-child comparison which is passed to m_motherChildProvider
    SConfig m_motherChildConfig;
    
    //Variant provider for parent-child comparison
    CMendelianVariantProvider m_provider;

    //Result log for mendelian comparison
    CResultLog m_resultLog;
    
    //Best Paths written by each thread for each unique chromosome exists [Between father and child]
    CPath m_aBestPathsFatherChildGT[CHROMOSOME_COUNT];
    CPath m_aBestPathsFatherChildAM[CHROMOSOME_COUNT];

    
    //Best Paths written by each thread for each unique chromosome exists [Between mother and child]
    CPath m_aBestPathsMotherChildGT[CHROMOSOME_COUNT];
    CPath m_aBestPathsMotherChildAM[CHROMOSOME_COUNT];

    
    //Thread pool we have for multitasking by per chromosome
    std::thread *m_pThreadPool;
};


#endif //C_MENDELIAN_ANALYZER_H_
