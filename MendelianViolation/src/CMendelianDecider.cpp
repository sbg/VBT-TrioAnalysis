//
//  CMendelianDecider.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/24/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianDecider.h"
#include "CMendelianVariantProvider.h"
#include "SChrIdTriplet.h"
#include "CVariantIterator.h"
#include "ENoCallMode.h"
#include <iostream>
#include "CMendelianResultLog.h"
#include <algorithm>

using namespace mendelian;

CMendelianDecider::CMendelianDecider(std::vector<core::CPath>& a_aBestPathsFatherChildGT,
                                    std::vector<core::CPath>& a_aBestPathsFatherChildAM,
                                     std::vector<core::CPath>& a_aBestPathsMotherChildGT,
                                     std::vector<core::CPath>& a_aBestPathsMotherChildAM,
                                     CMendelianVariantProvider& a_rProvider,
                                     CMendelianResultLog& a_rResultLog) :
m_aBestPathsFatherChildGT(a_aBestPathsFatherChildGT),
m_aBestPathsFatherChildAM(a_aBestPathsFatherChildAM),
m_aBestPathsMotherChildGT(a_aBestPathsMotherChildGT),
m_aBestPathsMotherChildAM(a_aBestPathsMotherChildAM),
m_provider(a_rProvider),
m_resultLog(a_rResultLog)
{
}

void CMendelianDecider::SetNocallMode(ENoCallMode a_nMode)
{
    m_nocallMode = a_nMode;
}


//Checks if the given two range is overlapping
extern bool isOverlap(int left1, int right1, int left2, int right2);


//Compare variants according to id for sort operation
bool variantCompare(const CVariant* v1, const CVariant* v2)
{
    return v1->m_nId < v2->m_nId;
}

void CMendelianDecider::CheckUniqueVars(EMendelianVcfName a_checkSide, SChrIdTriplet& a_rTriplet,
                                        const std::vector<const CVariant*>& a_rVariantList,
                                        std::vector<bool>& a_rSideDecisions,
                                        const std::vector<EMendelianDecision>& a_rParentSelfDecisions,
                                        const std::vector<EMendelianDecision>& a_rChildDecisions)
{
    std::vector<const CVariant*> varListToCheckChild = m_provider.GetVariantList(eCHILD, a_rTriplet.m_nCid);
    std::vector<const CVariant*> varListToCheckParent = a_checkSide == eMOTHER ? m_provider.GetVariantList(eFATHER, a_rTriplet.m_nFid) : m_provider.GetVariantList(eMOTHER, a_rTriplet.m_nMid);
    std::vector<const CVariant*> varListToCheckSelf = a_checkSide == eFATHER ? m_provider.GetVariantList(eFATHER, a_rTriplet.m_nFid) : m_provider.GetVariantList(eMOTHER, a_rTriplet.m_nMid);
    
    
    unsigned int varItrChild = 0;
    unsigned int varItrSelf = 0;
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
    {
        
        //Check if the variant is already marked
        if(a_rParentSelfDecisions[a_rVariantList[k]->m_nId] != eUnknown)
        {
            a_rSideDecisions[k] = a_rParentSelfDecisions[a_rVariantList[k]->m_nId] == eCompliant ? true : false;
            continue;
        }
        
        //If variant does not have 0 allele
        else if(a_rVariantList[k]->m_genotype[0] != 0 && a_rVariantList[k]->m_genotype[1] != 0)
        {
            a_rSideDecisions[k] = false;
            continue;
        }
        
        //Skip child variants until end position of child variant is more than our current Parent variant start position
        while(varItrChild <  varListToCheckChild.size() && varListToCheckChild[varItrChild]->m_nEndPos <= a_rVariantList[k]->m_nStartPos)
            varItrChild++;
        
        //Skip self variants until end position of child variant is more than our current Parent variant start position
        while(varItrSelf <  varListToCheckSelf.size() && varListToCheckSelf[varItrSelf]->m_nEndPos <= a_rVariantList[k]->m_nStartPos)
            varItrSelf++;
        
        unsigned int counterChild = 0;
        unsigned int counterSelf = 0;
        
        std::vector<const CVariant*> childSideMatches;
        std::vector<const CVariant*> selfSideMatches;
        
        //Get All overlapping variants in child side
        while(varItrChild + counterChild < varListToCheckChild.size() && varListToCheckChild[varItrChild+counterChild]->m_nStartPos < a_rVariantList[k]->m_nEndPos)
        {
            if(isOverlap(a_rVariantList[k]->m_nStartPos, a_rVariantList[k]->m_nEndPos, varListToCheckChild[varItrChild+counterChild]->m_nStartPos, varListToCheckChild[varItrChild+counterChild]->m_nEndPos))
                childSideMatches.push_back(varListToCheckChild[varItrChild+counterChild]);
            counterChild++;
        }
        
        //Get All overlapping variants in self side
        while(varItrSelf + counterSelf < varListToCheckSelf.size() && varListToCheckSelf[varItrSelf+counterSelf]->m_nStartPos < a_rVariantList[k]->m_nEndPos)
        {
            if(isOverlap(a_rVariantList[k]->m_nStartPos, a_rVariantList[k]->m_nEndPos, varListToCheckSelf[varItrSelf+counterSelf]->m_nStartPos, varListToCheckSelf[varItrSelf+counterSelf]->m_nEndPos))
                selfSideMatches.push_back(varListToCheckSelf[varItrSelf+counterSelf]);
            counterSelf++;
        }
        
        bool selfCheck = true;
        bool childCheck = true;
        
        //=====STEP 1: Check self side=====
        //check if all variants is 0/x form
        for(unsigned int m=0; m < selfSideMatches.size(); m++)
        {
            if(selfSideMatches[m]->m_genotype[0] != 0 && selfSideMatches[m]->m_genotype[1] != 0)
            {
                selfCheck = false;
                break;
            }
        }
        
        //===== STEP 2: Check child Side =========
        for(unsigned int m=0; m < childSideMatches.size(); m++)
        {
            if(a_rChildDecisions[varListToCheckChild[varItrChild+m]->m_nId] == eViolation)
            {
                childCheck = false;
                break;
            }
            
            else if(childSideMatches[m]->m_genotype[0] != 0 && childSideMatches[m]->m_genotype[1] != 0)
            {
                childCheck = false;
                break;
            }
        }
        
        a_rSideDecisions[k] = selfCheck && childCheck;
        
    }
    
}

void CMendelianDecider::GetSyncPointList(SChrIdTriplet& a_rTriplet, bool a_bIsFatherChild, std::vector<core::CSyncPoint>& a_rSyncPointList, bool a_bIsGT)
{
    std::vector<const core::COrientedVariant*> pBaseIncluded;
    std::vector<const core::COrientedVariant*> pCalledIncluded;
    
    std::vector<const CVariant*> pBaseExcluded;
    std::vector<const CVariant*> pCalledExcluded;
    
    core::CPath *pPath;
    core::CPath *pPathSync;
    
    if(a_bIsGT == false)
    {
        pPath = a_bIsFatherChild ? &m_aBestPathsFatherChildAM[a_rTriplet.m_nTripleIndex] : &m_aBestPathsMotherChildAM[a_rTriplet.m_nTripleIndex];
        core::CPath *pPathGT = a_bIsFatherChild ? &m_aBestPathsFatherChildGT[a_rTriplet.m_nTripleIndex] : &m_aBestPathsMotherChildGT[a_rTriplet.m_nTripleIndex];
        pPathSync = pPathGT;
        
        std::vector<const CVariant*> excludedVarsBase = m_provider.GetVariantList(a_bIsFatherChild ? eFATHER : eMOTHER, a_bIsFatherChild ? a_rTriplet.m_nFid : a_rTriplet.m_nMid, pPathGT->m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsCalled = m_provider.GetVariantList(eCHILD, a_rTriplet.m_nCid, pPathGT->m_calledSemiPath.GetExcluded());
        
        pBaseIncluded = pPath->m_baseSemiPath.GetIncludedVariants();
        pCalledIncluded = pPath->m_calledSemiPath.GetIncludedVariants();
        
        pBaseExcluded = m_provider.GetVariantList(excludedVarsBase, pPath->m_baseSemiPath.GetExcluded());
        pCalledExcluded = m_provider.GetVariantList(excludedVarsCalled, pPath->m_calledSemiPath.GetExcluded());
    }
    else
    {
        pPath = a_bIsFatherChild ? &m_aBestPathsFatherChildGT[a_rTriplet.m_nTripleIndex] : &m_aBestPathsMotherChildGT[a_rTriplet.m_nTripleIndex];
        pPathSync = pPath;
        
        pBaseIncluded = pPath->m_baseSemiPath.GetIncludedVariants();
        pCalledIncluded = pPath->m_calledSemiPath.GetIncludedVariants();
        
        pBaseExcluded = m_provider.GetVariantList(a_bIsFatherChild ? eFATHER : eMOTHER, a_bIsFatherChild ? a_rTriplet.m_nFid : a_rTriplet.m_nMid, pPath->m_baseSemiPath.GetExcluded());
        pCalledExcluded = m_provider.GetVariantList(eCHILD, a_rTriplet.m_nCid, pPath->m_calledSemiPath.GetExcluded());
    }
    
    
    unsigned int baseIncludedItr = 0;
    unsigned int baseExcludedItr = 0;
    unsigned int calledIncludedItr = 0;
    unsigned int calledExcludedItr = 0;
    
    for(unsigned int k = 0; k < pPathSync->m_aSyncPointList.size(); k++)
    {
        core::CSyncPoint ssPoint;
        ssPoint.m_nStartPosition = k > 0 ? pPathSync->m_aSyncPointList[k-1] : 0;
        ssPoint.m_nEndPosition = pPathSync->m_aSyncPointList[k];
        ssPoint.m_nIndex = (int)k;
        int bound = ssPoint.m_nStartPosition == ssPoint.m_nEndPosition ? 1 : 0;
        
        while(baseIncludedItr < pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() < (pPathSync->m_aSyncPointList[k] + bound))
        {
            const core::COrientedVariant* pOvar = pBaseIncluded[baseIncludedItr];
            ssPoint.m_baseVariantsIncluded.push_back(pOvar);
            baseIncludedItr++;
        }
        
        while(calledIncludedItr < pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() < (pPathSync->m_aSyncPointList[k] + bound))
        {
            const core::COrientedVariant* pOvar = pCalledIncluded[calledIncludedItr];
            ssPoint.m_calledVariantsIncluded.push_back(pOvar);
            calledIncludedItr++;
        }
        
        while(baseExcludedItr < pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos < (pPathSync->m_aSyncPointList[k] + bound))
        {
            ssPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
            baseExcludedItr++;
        }
        
        while(calledExcludedItr < pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos < (pPathSync->m_aSyncPointList[k] + bound))
        {
            ssPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
            calledExcludedItr++;
        }
        
        a_rSyncPointList.push_back(ssPoint);
    }
    
    //Add Remaining variants to the last syncPoint
    core::CSyncPoint sPoint;
    sPoint.m_nStartPosition = pPathSync->m_aSyncPointList[pPathSync->m_aSyncPointList.size()-1];
    sPoint.m_nEndPosition = INT_MAX;
    sPoint.m_nIndex = static_cast<int>(pPathSync->m_aSyncPointList.size()-1);
    
    while(baseIncludedItr < pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsIncluded.push_back(pBaseIncluded[baseIncludedItr]);
        baseIncludedItr++;
    }
    while(calledIncludedItr < pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsIncluded.push_back(pCalledIncluded[calledIncludedItr]);
        calledIncludedItr++;
    }
    
    while(baseExcludedItr < pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
        baseExcludedItr++;
    }
    
    while(calledExcludedItr < pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
        calledExcludedItr++;
    }
    a_rSyncPointList.push_back(sPoint);
    
}

void CMendelianDecider::CheckFor0PathFor00(SChrIdTriplet& a_rTriplet,
                                            bool a_bIsFatherChild,
                                            std::vector<const CVariant*>& a_rVarList,
                                            std::vector<const CVariant*>& a_rViolationList,
                                            std::vector<const CVariant*>& a_rCompliantList)
{
    //Get sync point list
    std::vector<core::CSyncPoint> a_rSyncPointList;
    GetSyncPointList(a_rTriplet, a_bIsFatherChild, a_rSyncPointList, true);
    
    int varlistItr = 0;
    for(unsigned int k = 0; k < a_rSyncPointList.size() && varlistItr < (int)a_rVarList.size(); k++)
    {
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < (int)a_rVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_rVarList[varlistItr]->m_nStartPos)
        {
            a_rViolationList.push_back(a_rVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(unsigned int m = 0; m < a_rSyncPointList[k].m_calledVariantsIncluded.size() && varlistItr != (int)a_rVarList.size(); m++)
        {
            if(a_rSyncPointList[k].m_calledVariantsIncluded[m]->GetVariant().m_nId == a_rVarList[varlistItr]->m_nId)
            {
                tmpVarList.push_back(a_rVarList[varlistItr]);
                varlistItr++;
                bDoCheck = true;
            }
        }
        
        if(true == bDoCheck)
        {
            bool bIsCompliant = true;
            
            for(unsigned int m = 0; m < a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                {
                    if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                    {
                        bIsCompliant = false;
                        break;
                    }
                }
            }
            
            if(bIsCompliant == true)
            {
                for(unsigned int aa = 0; aa < tmpVarList.size(); aa++)
                    a_rCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(unsigned int aa = 0; aa < tmpVarList.size(); aa++)
                    a_rViolationList.push_back(tmpVarList[aa]);
            }
        }
    }
    
}


void CMendelianDecider::CheckFor0Path(SChrIdTriplet& a_rTriplet,
                                       bool a_bIsFatherChild,
                                       std::vector<const CVariant *> &a_pVarList,
                                       std::vector<const CVariant *> &a_pViolantList,
                                       std::vector<const CVariant *> &a_pCompliantList,
                                       std::vector<EMendelianDecision>& a_rParentDecisions,
                                       bool a_bIsUpdateDecisionList)
{
    
    //Get sync point list
    std::vector<core::CSyncPoint> a_rSyncPointList;
    GetSyncPointList(a_rTriplet, a_bIsFatherChild, a_rSyncPointList);
    
    
    int varlistItr = 0;
    for(unsigned int k = 0; k < a_rSyncPointList.size() && varlistItr < (int)a_pVarList.size(); k++)
    {
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < (int)a_pVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_pVarList[varlistItr]->m_nStartPos)
        {
            a_pViolantList.push_back(a_pVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(unsigned int m = 0; m < a_rSyncPointList[k].m_calledVariantsExcluded.size() && varlistItr != (int)a_pVarList.size(); m++)
        {
            if(a_rSyncPointList[k].m_calledVariantsExcluded[m]->m_nId == a_pVarList[varlistItr]->m_nId)
            {
                tmpVarList.push_back(a_pVarList[varlistItr]);
                varlistItr++;
                bDoCheck = true;
            }
        }
        
        if(true == bDoCheck)
        {
            bool bIsCompliant = true;
            
            for(unsigned int m = 0; m < a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                {
                    
                    if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                    {
                        if(a_bIsUpdateDecisionList)
                        {
                            //We are marking decision of mother/father variant as violation
                            if(true == a_bIsFatherChild)
                                a_rParentDecisions[pVar->m_nId] = eViolation;
                        }
                        
                        bIsCompliant = false;
                        break;
                    }
                    
                    else
                    {
                        if(a_bIsUpdateDecisionList)
                        {
                            //We are marking decision of mother/father variant as compliant
                            if(true == a_bIsFatherChild)
                                a_rParentDecisions[pVar->m_nId] = eCompliant;
                        }
                    }
                }
            }
            
            if(bIsCompliant == true)
            {
                for(unsigned int aa = 0; aa < tmpVarList.size(); aa++)
                    a_pCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(unsigned int aa = 0; aa < tmpVarList.size(); aa++)
                    a_pViolantList.push_back(tmpVarList[aa]);
            }
        }
    }
    
}

void CMendelianDecider::CheckFor00Child(SChrIdTriplet& a_rTriplet,
                                         std::vector<const CVariant*>& a_rOvarList,
                                         std::vector<const CVariant*>& a_rViolationList,
                                         std::vector<const CVariant*>& a_rCompliantList,
                                         bool a_bIsGTMatch,
                                        std::vector<EMendelianDecision>& a_rMotherDecisions,
                                        std::vector<EMendelianDecision>& a_rFatherDecisions)
{
    std::vector<const CVariant*> fatherCompliants;
    std::vector<const CVariant*> motherCompliants;
    std::vector<const CVariant*> fatherViolants;
    std::vector<const CVariant*> motherViolants;
    
    std::sort(a_rOvarList.begin(), a_rOvarList.end(), variantCompare);
    
    if(false == a_bIsGTMatch)
    {
        CheckFor0Path(a_rTriplet, true,  a_rOvarList, fatherViolants, fatherCompliants, a_rFatherDecisions, false);
        CheckFor0Path(a_rTriplet, false, a_rOvarList, motherViolants, motherCompliants, a_rMotherDecisions, false);
    }
    else
    {
        CheckFor0PathFor00(a_rTriplet, true,  a_rOvarList, fatherViolants, fatherCompliants);
        CheckFor0PathFor00(a_rTriplet, false, a_rOvarList, motherViolants, motherCompliants);
    }
    
    //Intersect father and mother compliant variants and write to a_rCompliantList
    for (std::vector<const CVariant*>::iterator i = fatherCompliants.begin(); i != fatherCompliants.end(); ++i)
    {
        if (std::find(motherCompliants.begin(), motherCompliants.end(), *i) != motherCompliants.end())
        {
            a_rCompliantList.push_back(*i);
        }
    }
    
    //Write all other variants to a_rViolationList
    for (std::vector<const CVariant*>::iterator i = a_rOvarList.begin(); i != a_rOvarList.end(); ++i)
    {
        if (std::find(a_rCompliantList.begin(), a_rCompliantList.end(), *i) == a_rCompliantList.end())
        {
            a_rViolationList.push_back(*i);
        }
    }
    
}

void CMendelianDecider::AssignDecisionToParentVars(EMendelianVcfName a_checkSide, SChrIdTriplet& a_rTriplet, std::vector<EMendelianDecision>& a_rParentDecisions, std::vector<EMendelianDecision>& a_rChildDecisions)
{
    std::vector<const CVariant*> varListToCheckParent = a_checkSide == eMOTHER ? m_provider.GetVariantList(eMOTHER, a_rTriplet.m_nMid) : m_provider.GetVariantList(eFATHER, a_rTriplet.m_nFid);
    std::vector<const CVariant*> varListToCheckChild = m_provider.GetVariantList(eCHILD, a_rTriplet.m_nCid);
    
    assert(varListToCheckParent.size() == a_rParentDecisions.size());
    
    //Generate the sync point list
    std::vector<int> syncPoints(a_checkSide == eFATHER ? m_aBestPathsFatherChildGT[a_rTriplet.m_nTripleIndex].m_aSyncPointList : m_aBestPathsMotherChildGT[a_rTriplet.m_nTripleIndex].m_aSyncPointList);
    std::vector<int> syncPointsAM(a_checkSide == eFATHER ? m_aBestPathsFatherChildAM[a_rTriplet.m_nTripleIndex].m_aSyncPointList : m_aBestPathsMotherChildAM[a_rTriplet.m_nTripleIndex].m_aSyncPointList);
    
    //Concat AM + GT sync points
    syncPoints.insert( syncPoints.end(), syncPointsAM.begin(), syncPointsAM.end());
    //Sort sync points
    std::sort(syncPoints.begin(), syncPoints.end());

    //Iterator to synchronization point list
    unsigned int itrSyncPList = 0;
    unsigned int itrChildVars = 0;
    
    for(unsigned int k = 0; k < varListToCheckParent.size(); k++)
    {
        //Skipped assigned variants
        if(a_rParentDecisions[varListToCheckParent[k]->m_nId] != eUnknown)
            continue;
        
        //If parent variant is genotype match, this side has no problem
        if(varListToCheckParent[k]->m_variantStatus == eGENOTYPE_MATCH || varListToCheckParent[k]->m_variantStatus == eALLELE_MATCH)
        {
            a_rParentDecisions[varListToCheckParent[k]->m_nId] = eCompliant;
            continue;
        }
        
        //If variant is no match and has no reference allele, mark it as violation
        if(varListToCheckParent[k]->m_variantStatus == eNO_MATCH && varListToCheckParent[k]->m_genotype[0] != 0 && varListToCheckParent[k]->m_genotype[1] != 0)
        {
            a_rParentDecisions[varListToCheckParent[k]->m_nId] = eViolation;
            continue;
        }
        
        //Skip irrelevant sync points
        while(itrSyncPList < syncPoints.size() && syncPoints[itrSyncPList] <= varListToCheckParent[k]->m_nStartPos)
            itrSyncPList++;
        
        //Terminate if the sync point list is ended
        if(itrSyncPList == syncPoints.size())
            break;
        
        //Skip irrelevant child variants
        while (itrChildVars < varListToCheckChild.size() && varListToCheckChild[itrChildVars]->m_nEndPos < varListToCheckParent[k]->m_nStartPos)
            itrChildVars++;

        //Terminate if the child variant list is ended
        if(itrChildVars == varListToCheckChild.size())
            break;
        
        bool bIsViolationFound = false;
        
        std::vector<const CVariant*> betweenSyncChildVars;
        unsigned int secondChildItr = itrChildVars;
        while (secondChildItr < varListToCheckChild.size() && varListToCheckChild[secondChildItr]->m_nEndPos < syncPoints[itrSyncPList])
            betweenSyncChildVars.push_back(varListToCheckChild[secondChildItr++]);
        
        //If there is no corresponding child variant between the sync points
        if(betweenSyncChildVars.size() == 0 && varListToCheckParent[k]->m_genotype[0] != 0 && varListToCheckParent[k]->m_genotype[1] != 0)
        {
            a_rParentDecisions[varListToCheckParent[k]->m_nId] = eViolation;
            bIsViolationFound = true;
            continue;
        }
        
        for(unsigned int m = 0; m < betweenSyncChildVars.size(); m++)
        {
            if(a_rChildDecisions[betweenSyncChildVars[m]->m_nId] == eViolation)
            {
                a_rParentDecisions[varListToCheckParent[k]->m_nId] = eViolation;
                bIsViolationFound = true;
                break;
            }
        }
        
        if(!bIsViolationFound)
            a_rParentDecisions[varListToCheckParent[k]->m_nId] = eCompliant;
    }
}

void CMendelianDecider::ReportChildChromosomeData(SChrIdTriplet& a_rTriplet, std::vector<const CVariant*>& a_rCompliants, std::vector<const CVariant*>& a_rViolations)
{
    int compliantSNPcount = 0;
    int compliantINDELcount = 0;
    int violationSNPcount = 0;
    int violationINDELcount = 0;
    
    for(const CVariant* pVar : a_rCompliants)
    {
        if(pVar->GetVariantType() == eSNP)
            compliantSNPcount++;
        else
            compliantINDELcount++;
    }
    
    for(const CVariant* pVar : a_rViolations)
    {
        if(pVar->GetVariantType() == eSNP)
            violationSNPcount++;
        else
            violationINDELcount++;
    }
    
    m_resultLog.LogShortReport(a_rTriplet.m_chrName, compliantSNPcount, violationSNPcount, compliantINDELcount, violationINDELcount);
}

void CMendelianDecider::MergeFunc(SChrIdTriplet& a_triplet,
                                  std::vector<EMendelianDecision>& a_rMotherDecisions,
                                  std::vector<EMendelianDecision>& a_rFatherDecisions,
                                  std::vector<EMendelianDecision>& a_rChildDecisions)
{
    std::vector<const CVariant*> compliants;
    std::vector<const CVariant*> violations;
    
    std::vector<const CVariant*> SameAlleleMatchViolationVars;
    std::vector<const CVariant*> MendelianCompliantVars;
    
    std::vector<const CVariant*> motherChildOnly;
    std::vector<const CVariant*> fatherChildOnly;
    
    std::vector<const CVariant*> check0atMotherSide;
    std::vector<const CVariant*> check0atFatherSide;
    std::vector<const CVariant*> check00Child;
    std::vector<const CVariant*> childUniqueList;
    std::vector<const CVariant*> check00ChildGTMatched;
    
    //Sort variants according to variant ids
    m_aBestPathsFatherChildGT[a_triplet.m_nTripleIndex].SortIncludedVariants();
    m_aBestPathsFatherChildAM[a_triplet.m_nTripleIndex].SortIncludedVariants();
    m_aBestPathsMotherChildGT[a_triplet.m_nTripleIndex].SortIncludedVariants();
    m_aBestPathsMotherChildAM[a_triplet.m_nTripleIndex].SortIncludedVariants();
    
    //Included lists of child
    CVariantIterator FatherChildVariants(m_aBestPathsFatherChildGT[a_triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants(),
                                         m_aBestPathsFatherChildAM[a_triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants());
    
    CVariantIterator MotherChildVariants(m_aBestPathsMotherChildGT[a_triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants(),
                                         m_aBestPathsMotherChildAM[a_triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants());
    
    //Check if the two list have common variants
    if(FatherChildVariants.hasNext() == false && MotherChildVariants.hasNext() == false)
        return;
    
    const core::COrientedVariant* varMC = (MotherChildVariants.hasNext() ? MotherChildVariants.Next() : NULL);
    const core::COrientedVariant* varFC = (FatherChildVariants.hasNext() ? FatherChildVariants.Next() : NULL);
    
    while(varMC != NULL && varFC != NULL)
    {
        //If Mother variant and Father variant is Common check for same allele match condition
        if(varMC->GetVariant().m_nId == varFC->GetVariant().m_nId)
        {
            if(varMC->GetVariant().m_genotype[0] == 0 && varMC->GetVariant().m_genotype[1] == 0)
                check00ChildGTMatched.push_back(&varMC->GetVariant());
            
            else if(varMC->GetVariant().m_bIsHeterozygous)
            {
                //ELIMINATE SAME ALLELE MATCHING(EXCEPTION 1)
                if(varMC->GetAlleleIndex() == varFC->GetAlleleIndex())
                {
                    if(varMC->GetVariant().m_variantStatus == eGENOTYPE_MATCH
                       ||
                       varFC->GetVariant().m_variantStatus == eGENOTYPE_MATCH)
                        MendelianCompliantVars.push_back(&varMC->GetVariant());
                    else
                        SameAlleleMatchViolationVars.push_back(&varMC->GetVariant());
                }
                else
                    MendelianCompliantVars.push_back(&varMC->GetVariant());
            }
            else
                MendelianCompliantVars.push_back(&varMC->GetVariant());
            
            if(MotherChildVariants.hasNext() && FatherChildVariants.hasNext())
            {
                varMC = MotherChildVariants.Next();
                varFC = FatherChildVariants.Next();
            }
            
            else
            {
                varMC = NULL;
                varFC = NULL;
                break;
            }
        }
        
        //If we have variant match with Father side only, we filter 0/x variants and rest of them are marked as violation
        else if(varMC->GetVariant().m_nId > varFC->GetVariant().m_nId)
        {
            if(varFC->GetVariant().m_genotype[0] == 0 || varFC->GetVariant().m_genotype[1] == 0)
                check0atMotherSide.push_back(&varFC->GetVariant());
            else
                fatherChildOnly.push_back(&varFC->GetVariant());
            
            if(FatherChildVariants.hasNext())
                varFC = FatherChildVariants.Next();
            else
            {
                varFC = NULL;
                break;
            }
        }
        
        //If we have variant match with Mother side only, we filter 0/x variants and rest of them are marked as violation
        else
        {
            if(varMC->GetVariant().m_genotype[0] == 0 || varMC->GetVariant().m_genotype[1] == 0)
                check0atFatherSide.push_back(&varMC->GetVariant());
            else
                motherChildOnly.push_back(&varMC->GetVariant());

            if(MotherChildVariants.hasNext())
                varMC = MotherChildVariants.Next();
            else
            {
                varMC = NULL;
                break;
            }
        }
        
    }
    
    //Process remaining vars in FatherChild explained as above
    while(true)
    {
        if(varFC == NULL)
        {
            if(FatherChildVariants.hasNext())
                varFC = FatherChildVariants.Next();
            else
                break;
        }
        
        if(varFC->GetVariant().m_genotype[0] == 0 || varFC->GetVariant().m_genotype[1] == 0)
            check0atMotherSide.push_back(&varFC->GetVariant());
        else
            fatherChildOnly.push_back(&varFC->GetVariant());

        if(!FatherChildVariants.hasNext())
            break;
        else
            varFC = FatherChildVariants.Next();
    }
    
    //Process remaining vars in MotherChild explined as above
    while(true)
    {
        if(varMC == NULL)
        {
            if(MotherChildVariants.hasNext())
                varMC = MotherChildVariants.Next();
            else
                break;
        }
        
        if(varMC->GetVariant().m_genotype[0] == 0 || varMC->GetVariant().m_genotype[1] == 0)
            check0atFatherSide.push_back(&varMC->GetVariant());
        else
            motherChildOnly.push_back(&varMC->GetVariant());

        if(!MotherChildVariants.hasNext())
            break;
        else
            varMC = MotherChildVariants.Next();
    }
    
    std::vector<const CVariant*> compliantVarsFrom0CheckMother;
    std::vector<const CVariant*> compliantVarsFrom0CheckFather;
    std::vector<const CVariant*> violationVarsFrom0CheckMother;
    std::vector<const CVariant*> violationVarsFrom0CheckFather;
    
    //Check for 0/x child variant set at father side
    CheckFor0Path(a_triplet, true, check0atFatherSide, violationVarsFrom0CheckFather, compliantVarsFrom0CheckFather, a_rFatherDecisions);
    //Check for 0/x child variant set at the mother side
    CheckFor0Path(a_triplet, false, check0atMotherSide, violationVarsFrom0CheckMother, compliantVarsFrom0CheckMother, a_rMotherDecisions);
    
    //Gather all compliant variants of child we found so far
    compliants.insert(std::end(compliants), std::begin(MendelianCompliantVars), std::end(MendelianCompliantVars));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckFather), std::end(compliantVarsFrom0CheckFather));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckMother), std::end(compliantVarsFrom0CheckMother));
    std::sort(compliants.begin(), compliants.end(), variantCompare);
    
    //Gather all violation variants of child we found so far
    violations.insert(std::end(violations), std::begin(SameAlleleMatchViolationVars), std::end(SameAlleleMatchViolationVars));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckFather), std::end(violationVarsFrom0CheckFather));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckMother), std::end(violationVarsFrom0CheckMother));
    violations.insert(std::end(violations), std::begin(fatherChildOnly), std::end(fatherChildOnly));
    violations.insert(std::end(violations), std::begin(motherChildOnly), std::end(motherChildOnly));
    
    //Find Child Unique variants
    std::vector<const CVariant*> childVariants = m_provider.GetSortedVariantListByID(eCHILD, a_triplet.m_nCid);
    std::vector<int>childProcessedArray(childVariants.size());
    for(unsigned int k = 0; k < childProcessedArray.size(); k++)
        childProcessedArray[k] = 0;
    
    //Mark mendelian compliant vars as processed
    for(unsigned int k = 0; k < compliants.size(); k++)
    {
        childProcessedArray[compliants[k]->m_nId]++;
    }
    
    //Mark mendelian violation vars as processed
    for(unsigned int k = 0; k < violations.size(); k++)
    {
        childProcessedArray[violations[k]->m_nId]++;
    }
    
    for(unsigned int childItr = 0; childItr < childProcessedArray.size(); childItr++)
    {
        if(childProcessedArray[childItr] == 0)
        {
            if(childVariants[childItr]->m_genotype[0] == 0 && childVariants[childItr]->m_genotype[1] == 0)
            {
                check00Child.push_back(childVariants[childItr]);
            }
            else
            {
                const CVariant* pVar = childVariants[childItr];
                childUniqueList.push_back(pVar);
            }
        }
    }
    
    //Add the new violations we found to violation list
    violations.insert(std::end(violations), std::begin(childUniqueList), std::end(childUniqueList));
    std::sort(violations.begin(), violations.end(), variantCompare);
    
    //We looked up all child variants. Now, we will look at parent variants where there is no corresponding child variant exist in the child.vcf (check for hidden 0/0 child variants)
    
    //Fill the child decision array
    std::vector<const CVariant*>::iterator compliantsIterator = compliants.begin();
    std::vector<const CVariant*>::iterator violationsIterator = violations.begin();
    
    int counterCompliant = 0;
    int counterViolation = 0;
    
    for(unsigned int k = 0; k < childVariants.size(); k++)
    {
        if(compliantsIterator != compliants.end() && childVariants[k]->m_nId == (*compliantsIterator)->m_nId)
        {
            a_rChildDecisions[childVariants[k]->m_nId] = EMendelianDecision::eCompliant;
            compliantsIterator++;
            counterCompliant++;
        }
        else if(violationsIterator != violations.end() && childVariants[k]->m_nId == (*violationsIterator)->m_nId)
        {
            a_rChildDecisions[childVariants[k]->m_nId] = EMendelianDecision::eViolation;
            violationsIterator++;
            counterViolation++;
        }
        else
            a_rChildDecisions[childVariants[k]->m_nId] = EMendelianDecision::eUnknown;
    }
    
    //Excluded Mother variant Check - If we can find 0/0 hidden child site correspond to mother variant
    std::vector<const CVariant*> uniqueMotherVars = m_provider.GetVariantList(m_provider.GetVariantList(eMOTHER, a_triplet.m_nMid,  m_aBestPathsMotherChildGT[a_triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded()),
                                                                              m_aBestPathsMotherChildAM[a_triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
    std::vector<bool> motherDecisions(uniqueMotherVars.size());
    CheckUniqueVars(eMOTHER, a_triplet, uniqueMotherVars, motherDecisions, a_rMotherDecisions, a_rChildDecisions);
    
    //Excluded Father variant Check - If we can find 0/0 hidden child site correspond to father variant
    std::vector<const CVariant*> uniqueFatherVars = m_provider.GetVariantList(m_provider.GetVariantList(eFATHER, a_triplet.m_nFid,  m_aBestPathsFatherChildGT[a_triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded()),
                                                                              m_aBestPathsFatherChildAM[a_triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
    std::vector<bool> fatherDecisions(uniqueFatherVars.size());
    CheckUniqueVars(eFATHER, a_triplet, uniqueFatherVars, fatherDecisions, a_rFatherDecisions, a_rChildDecisions);
    
    //Fill the mother decision array
    for(unsigned int k = 0; k < motherDecisions.size(); k++)
        a_rMotherDecisions[uniqueMotherVars[k]->m_nId] = (motherDecisions[k] ? eCompliant : eViolation);
    
    //Fill the father decision array
    for(unsigned int k = 0; k < fatherDecisions.size(); k++)
        a_rFatherDecisions[uniqueFatherVars[k]->m_nId] = (fatherDecisions[k] ? eCompliant : eViolation);
    
    //If NoCall Mode Is not enabled, mark all decisions of nocall childs as NoCallChild and all nocall parents as NoCallParent
    if(m_nocallMode != eNone)
    {
        std::vector<const CVariant*> motherVariants = m_provider.GetVariantList(eMOTHER, a_triplet.m_nMid);
        std::vector<const CVariant*> fatherVariants = m_provider.GetVariantList(eFATHER, a_triplet.m_nFid);
        
        for(unsigned int k = 0; k < motherVariants.size(); k ++)
        {
            if(motherVariants[k]->m_bIsNoCall)
                a_rMotherDecisions[motherVariants[k]->m_nId] = eNoCallParent;
        }
        
        for(unsigned int k = 0; k < fatherVariants.size(); k ++)
        {
            if(fatherVariants[k]->m_bIsNoCall)
                a_rFatherDecisions[fatherVariants[k]->m_nId] = eNoCallParent;
        }
        
        for(unsigned int k = 0; k < childVariants.size(); k ++)
        {
            if(childVariants[k]->m_bIsNoCall)
                a_rChildDecisions[childVariants[k]->m_nId] = eNoCallChild;
        }
    }
    
    //Assign the remaining unassigned parent variants
    AssignDecisionToParentVars(eMOTHER, a_triplet, a_rMotherDecisions, a_rChildDecisions);
    AssignDecisionToParentVars(eFATHER, a_triplet, a_rFatherDecisions, a_rChildDecisions);
    
    ReportChildChromosomeData(a_triplet, compliants, violations);
    
    std::cerr << "===================== STATISTICS " << a_triplet.m_chrName << " ===================" << std::endl;
    std::cerr << "Total Compliants:" << compliants.size() << std::endl;
    std::cerr << "Total Violations:" << violations.size() << std::endl;
    std::cerr << "Child Var Size:" << childVariants.size()<< std::endl;
    std::cerr << "=====================================================" << std::endl << std::endl;
    
}
