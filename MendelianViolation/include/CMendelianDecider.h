//
//  CMendelianDecider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/24/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_MENDELIAN_DECIDER_H_
#define _C_MENDELIAN_DECIDER_H_

#include "CPath.h"
#include "EMendelianVcfName.h"
#include "EMendelianDecision.h"
#include "ENoCallMode.h"
#include "CSyncPoint.h"

namespace mendelian
{
    
    class CMendelianVariantProvider;
    class CMendelianResultLog;
    struct SChrIdTriplet;
    
class CMendelianDecider
{
    
public:
    
    CMendelianDecider(std::vector<core::CPath>& a_aBestPathsFatherChildGT,
                      std::vector<core::CPath>& a_aBestPathsFatherChildAM,
                      std::vector<core::CPath>& a_aBestPathsMotherChildGT,
                      std::vector<core::CPath>& a_aBestPathsMotherChildAM,
                      CMendelianVariantProvider& a_rProvider,
                      CMendelianResultLog& a_rResultLog);
    
    //Sets the nocall mode selected
    void SetNocallMode(ENoCallMode a_nMode);
    
    //A function that perform merge operations after best path algorith for mother-child and father-child are called
    void MergeFunc(SChrIdTriplet& a_rTriplet,
                   std::vector<EMendelianDecision>& a_rMotherDecisions,
                   std::vector<EMendelianDecision>& a_rFatherDecisions,
                   std::vector<EMendelianDecision>& a_rChildDecisions);
    
    //Return the decision lists of requested sample and chromosome
    std::vector<EMendelianDecision> GetDecisionList(EMendelianVcfName a_checkSide, const SChrIdTriplet& a_rTriplet) const;
    
    
private:
    
    
    //Return the syncpointlist for given comparison. Writes to the last parameter
    void GetSyncPointList(SChrIdTriplet& a_rTriplet,
                          bool a_bIsFatherChild,
                          std::vector<core::CSyncPoint>& a_rSyncPointList,
                          bool a_bIsGT = false);
    
    //Check sync points which child excluded contains 0 Allele variant. If that 0 allele is playable for the parent, we mark variants as compliant, violation otherwise
    void CheckFor0Path(SChrIdTriplet& a_rTriplet,
                       bool a_bIsFatherChild,
                       std::vector<const CVariant*>& a_rOvarList,
                       std::vector<const CVariant*>& a_rViolationList,
                       std::vector<const CVariant*>& a_rCompliantList,
                       std::vector<EMendelianDecision>& a_rParentDecisions,
                       bool a_bIsUpdateDecisionList = true);
    
    //Check each unique vars and seek for 0 path at the requested side and fills the a_rSideDecisionList
    void CheckUniqueVars(EMendelianVcfName a_checkSide,
                         SChrIdTriplet& a_rTriplet,
                         const std::vector<const CVariant*>& a_rVariantList,
                         std::vector<bool>& a_rSideDecisions,
                         const std::vector<EMendelianDecision>& a_rParentSelfDecisions,
                         const std::vector<EMendelianDecision>& a_rParentOtherDecisions,
                         const std::vector<EMendelianDecision>& a_rChildDecisions);
    
    //Check each parent variant and assign all unassigned parent variants as violation or consistent
    void AssignDecisionToParentVars(EMendelianVcfName a_checkSide,
                                    SChrIdTriplet& a_rTriplet, std::vector<EMendelianDecision>& a_rParentDecisions,
                                    std::vector<EMendelianDecision>& a_rChildDecisions);
    
    //Check sync points which child excluded contains 0 Allele variant. If that 0 allele is playable for the parent, we mark variants as compliant, violation otherwise
    void CheckFor0PathFor00(SChrIdTriplet& a_rTriplet,
                            bool a_bIsFatherChild,
                            std::vector<const CVariant*>& a_rOvarList,
                            std::vector<const CVariant*>& a_rViolationList,
                            std::vector<const CVariant*>& a_rCompliantList,
                            std::vector<EMendelianDecision>& a_rMotherDecisions,
                            std::vector<EMendelianDecision>& a_rFatherDecisions);
    
    //Check sync points which child excluded variant is 0/0. If that is playable for the both parent, we mark variants as compliant, violation otherwise
    void CheckFor00Child(SChrIdTriplet& a_rTriplet,
                         std::vector<const CVariant*>& a_rOvarList,
                         std::vector<const CVariant*>& a_rViolationList,
                         std::vector<const CVariant*>& a_rCompliantList,
                         bool a_bIsGTMatch,
                         std::vector<EMendelianDecision>& a_rMotherDecisions,
                         std::vector<EMendelianDecision>& a_rFatherDecisions);

    //Report short output table (Non 0/0 child variants only)
    void ReportChildChromosomeData(SChrIdTriplet& a_rTriplet,
                                   std::vector<const CVariant*>& a_rCompliants,
                                   std::vector<const CVariant*>& a_rViolations);
    
    //No-call mode selected by user (Default is explicit)
    ENoCallMode m_nocallMode;
    
    //Best Paths written by each thread for each unique chromosome exists [Between father and child]
    std::vector<core::CPath>& m_aBestPathsFatherChildGT;
    std::vector<core::CPath>& m_aBestPathsFatherChildAM;
    
    //Best Paths written by each thread for each unique chromosome exists [Between mother and child]
    std::vector<core::CPath>& m_aBestPathsMotherChildGT;
    std::vector<core::CPath>& m_aBestPathsMotherChildAM;
    
    //Variant provider for parent-child comparison
    CMendelianVariantProvider& m_provider;
    
    //Result log for mendelian comparison
    CMendelianResultLog& m_resultLog;

    
};
    
}

#endif //_C_MENDELIAN_DECIDER_H_
