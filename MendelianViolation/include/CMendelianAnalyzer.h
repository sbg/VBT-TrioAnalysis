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
#include "CMendelianTrioMerger.h"
#include "EMendelianDecision.h"
#include "ENoCallMode.h"

class CMendelianAnalyzer
{
    
public:
    
    //Executes the mendelian violation analyzer
    void run(int argc, char** argv);
    
private:
    
    //Read the parameters if the execution mode is mendelian. If all mandatory parameters are set, return true.
    bool ReadParameters(int argc, char** argv);
    
    //A function that perform merge operations after best path algorith m for mother-child and father-child is called
    void MergeFunc(int a_nChromosomeId);
    
    //Return the syncpointlist for given comparison. Writes to the last parameter
    void GetSyncPointList(int a_nChrId, bool a_bIsFatherChild, std::vector<CSyncPoint>& a_rSyncPointList, bool a_bIsGT = false);
    
    //Check sync points which child excluded contains 0 Allele variant. If that 0 allele is playable for the parent, we mark variants as compliant, violation otherwise
    void CheckFor0Path(int a_nChrId,
                       bool a_bIsFatherChild,
                       std::vector<const CVariant*>& a_rOvarList,
                       std::vector<const CVariant*>& a_rViolationList,
                       std::vector<const CVariant*>& a_rCompliantList,
                       bool a_bIsUpdateDecisionList = true);
    
    //Check each unique vars and seek for 0 path at the requested side and fills the a_rSideDecisionList
    void CheckUniqueVars(EMendelianVcfName a_checkSide, int a_nChrId, const std::vector<const CVariant*>& a_rVariantList, std::vector<bool>& a_rSideDecisions);
    
    //Check sync points which child excluded contains 0 Allele variant. If that 0 allele is playable for the parent, we mark variants as compliant, violation otherwise
    void CheckFor0PathFor00(int a_nChrId,
                            bool a_bIsFatherChild,
                            std::vector<const CVariant*>& a_rOvarList,
                            std::vector<const CVariant*>& a_rViolationList,
                            std::vector<const CVariant*>& a_rCompliantList);
    
    
    //Check sync points which child excluded variant is 0/0. If that is playable for the both parent, we mark variants as compliant, violation otherwise
    void CheckFor00Child(int a_nChrId,
                         std::vector<const CVariant*>& a_rOvarList,
                         std::vector<const CVariant*>& a_rViolationList,
                         std::vector<const CVariant*>& a_rCompliantList,
                         bool a_bIsGTMatch);
        
    //A Function to process mendelian violation pipeline for given chromosome id
    void ProcessChromosome(const std::vector<int>& a_rChromosomeIds);
    
    //Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
    int AssignJobsToThreads(int a_nThreadCount);
        
    //Config file for father-child comparison which is passed to m_fatherChildProvider
    SConfig m_fatherChildConfig;
    
    //Config file for mother-child comparison which is passed to m_motherChildProvider
    SConfig m_motherChildConfig;
    
    //No call mode of the variant
    ENoCallMode m_noCallMode;
    
    //Variant provider for parent-child comparison
    CMendelianVariantProvider m_provider;

    //Trio merger for outputing
    CMendelianTrioMerger m_trioWriter;
    
    //Result log for mendelian comparison
    CResultLog m_resultLog;
    
    //Best Paths written by each thread for each unique chromosome exists [Between father and child]
    CPath m_aBestPathsFatherChildGT[CHROMOSOME_COUNT];
    CPath m_aBestPathsFatherChildAM[CHROMOSOME_COUNT];

    //Best Paths written by each thread for each unique chromosome exists [Between mother and child]
    CPath m_aBestPathsMotherChildGT[CHROMOSOME_COUNT];
    CPath m_aBestPathsMotherChildAM[CHROMOSOME_COUNT];
    
    //Menndelian compliant/violation decision of each child variant
    std::vector<EMendelianDecision> m_aChildDecisions[CHROMOSOME_COUNT];
    std::vector<EMendelianDecision> m_aFatherDecisions[CHROMOSOME_COUNT];
    std::vector<EMendelianDecision> m_aMotherDecisions[CHROMOSOME_COUNT];
    
    
    //Thread pool we have for multitasking by per chromosome
    std::thread *m_pThreadPool;
};


#endif //C_MENDELIAN_ANALYZER_H_
