//
//  CPhaseEvaluator.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 11/16/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CPhaseEvaluator.h"


CPhaseEvaluator::CPhaseEvaluator(CVariantProvider *a_pVarProvider)
{
    m_pVariantProvider = a_pVarProvider;
    m_nMisPhasingCount = 0;
    m_UnphaseableCount = 0;
    m_nCorrectPhasingCount = 0;
    
    m_bBasePhase = false;
    m_bCallPhase = false;
    m_bBaseIsPhased = false;
    m_bCallIsPhased = false;
}


bool groupInPhase(std::vector<SVariantSummary> group)
{
    if (group.size() < 1)
        return true;
    if (!group[0].isPhased())
        return false;
    bool phase = group[0].phase();
    
    for (SVariantSummary v : group)
    {
        if (!v.isPhased() || v.phase() != phase)
            return false;
    }
    return true;
}


SPhasingResult CPhaseEvaluator::CountMisphasings(const CPath& a_rBestPath, int a_nChrId)
{
    std::vector<CVariant> excludedVarsBase = m_pVariantProvider->GetVariantList(eBASE, a_nChrId, a_rBestPath.m_baseSemiPath.GetExcluded());
    std::vector<CVariant> excludedVarsCall = m_pVariantProvider->GetVariantList(eCALLED, a_nChrId, a_rBestPath.m_calledSemiPath.GetExcluded());
    
    //std::vector<COrientedVariant> includedVarsBase = m_pVariantProvider->Get
    
    
    CCallIterator baseline(a_rBestPath.m_baseSemiPath.GetIncludedVariants(), excludedVarsBase);
    CCallIterator calls(a_rBestPath.m_calledSemiPath.GetIncludedVariants(), excludedVarsCall);
    
    std::vector<int> syncpoints = a_rBestPath.m_aSyncPointList;
    
    int misPhasings = 0;
    int unphaseable = 0;
    int correctPhasings = 0;
    bool baseIsPhased = false;
    bool basePhase = false;
    bool callIsPhased = false;
    bool callPhase = false;
    SVariantSummary call;
    SVariantSummary base;
    
    // Rather than assuming baseline is all phased, we'll be a bit more careful.
    for(int point : syncpoints)
    {
        // Collect all baseline variants in the sync region
        std::vector<SVariantSummary> baselineSection;
        do
        {
            if (!base.isNull() && base.startPos() < point)
            {
                baselineSection.push_back(base);
                base = baseline.hasNext() ? baseline.next() : SVariantSummary();
            }
            
            if (base.isNull())
            {
                base = baseline.hasNext() ? baseline.next() : SVariantSummary();
            }
        }
        while(!base.isNull() && base.startPos() < point);
        
        // Collect all called variants in the sync region
        std::vector<SVariantSummary> callSection;
        do
        {
            if (!call.isNull() && call.startPos() < point)
            {
                callSection.push_back(call);
                call = calls.hasNext() ? calls.next() : SVariantSummary();
            }
            
            if (call.isNull())
            {
                call = calls.hasNext() ? calls.next() : SVariantSummary();
            }
        }
        while(!call.isNull() && call.startPos() < point);
        
        // We need all of these to be in the same phasing otherwise we can't tell which calls are swapped
        if(!groupInPhase(baselineSection))
        {
            for (SVariantSummary summary : callSection)
            {
                if (summary.isPhased())
                    unphaseable++;
            }
            
            baseIsPhased = false;
            callIsPhased = false;
            continue;
        }
        
        // When the baseline calls have flipped orientation
        bool transition = false;
        if (!baseIsPhased)
            callIsPhased = false;
    
        for (SVariantSummary baseSummary : baselineSection)
        {
            if (baseSummary.isPhased())
            {
                if (basePhase != baseSummary.phase())
                    transition = true;
                baseIsPhased = true;
            }
            
            basePhase = baseSummary.phase();
        }
        
        for(SVariantSummary currentCall : callSection)
        {
            if (!currentCall.isPhased())
                callIsPhased = false;
            else if (!callIsPhased)
            {
                //Start phasing run
                callIsPhased = true;
                callPhase = currentCall.phase();
            }
            else
            {
                //Continue phasing
                const bool callTransition = (currentCall.phase() != callPhase);
                if (currentCall.included() && !(callTransition == transition))
                {
                    misPhasings++;
                    callPhase = currentCall.phase();
                }
                else if (currentCall.included())
                {
                    correctPhasings++;
                    callPhase = currentCall.phase();
                }
            }
            transition = false;
        }
    }
    
    return SPhasingResult(misPhasings, correctPhasings, unphaseable);

}


