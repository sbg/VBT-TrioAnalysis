/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
 *
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_PHASE_EVALUATOR_H_
#define _C_PHASE_EVALUATOR_H_

#include "CPath.h"
#include "CVariantProvider.h"

struct SPhasingResult
{
    SPhasingResult(int a_nMisPhasingCount, int a_nCorrectPhasingCount, int a_nUnphaseableCount)
    {
        m_nMisPhasingCount = a_nMisPhasingCount;
        m_nCorrectPhasingCount = a_nCorrectPhasingCount;
        m_UnphaseableCount = a_nUnphaseableCount;
    }
    
    int m_nMisPhasingCount;
    int m_nCorrectPhasingCount;
    int m_UnphaseableCount;
};

struct SVariantSummary
{
    const CVariant* m_pVariant;
    bool m_bIncluded;
    bool m_bPhase;
    
    bool m_bIsNull;
    
    SVariantSummary()
    {
        m_pVariant = 0;
        m_bIncluded = false;
        m_bPhase = false;
        m_bIsNull = true;
    }
    
    SVariantSummary(const CVariant& v, bool include, bool alternate)
    {
        m_pVariant = &v;
        m_bIncluded = include;
        m_bPhase = alternate;
        m_bIsNull = false;
    }
    
    SVariantSummary(const SVariantSummary& a_rObj)
    {
        m_pVariant = a_rObj.m_pVariant;
        m_bIncluded = a_rObj.m_bIncluded;
        m_bPhase = a_rObj.m_bPhase;
        m_bIsNull = a_rObj.m_bIsNull;
    }
    
    bool isPhased()
    {
        return m_pVariant->IsPhased();
    }
    int startPos()
    {
        return m_pVariant->GetStart();
    }
    bool included()
    {
        return m_bIncluded;
    }
    bool phase()
    {
        assert(isPhased());
        return m_bPhase;
    }
    
    bool isNull()
    {
        return m_bIsNull;
    }
    
    
};


//Combines included/excluded calls into one stream for purposes of bridging phasing across fp
class CCallIterator
{
    std::vector<COrientedVariant*>::iterator it_Included;
    std::vector<CVariant>::iterator it_Excluded;
    
    std::vector<COrientedVariant*> m_aIncluded;
    std::vector<CVariant> m_aExcluded;

    public:

    CCallIterator(std::vector<COrientedVariant*> included, std::vector<CVariant> excluded)
    {
        m_aIncluded = included;
        m_aExcluded = excluded;
        
        it_Included = m_aIncluded.begin();
        it_Excluded = m_aExcluded.begin();
    }
    
    bool hasNext()
    {
        return (it_Included != m_aIncluded.end() || it_Excluded != m_aExcluded.end());
    }
    
    SVariantSummary next()
    {
        SVariantSummary result;
        if(it_Included == m_aIncluded.end() || (it_Excluded != m_aExcluded.end() && it_Excluded->GetStart() < (*it_Included)->GetVariant().GetStart()))
        {
            result = SVariantSummary(*it_Excluded, false, false);
            it_Excluded++;
        }
        else
        {
            result = SVariantSummary((*it_Included)->GetVariant(), true, (*it_Included)->IsOrderOfGenotype());
            it_Included++;
        }
        
        return result;
    }
};


class CPhaseEvaluator
{
    public:
    
    CPhaseEvaluator(CVariantProvider *a_pVarProvider);

    SPhasingResult CountMisphasings(CPath& a_rPath, int a_nChrId);
    
    private:
    
    //Access to variant provider
    CVariantProvider* m_pVariantProvider;
    
    int m_nMisPhasingCount;
    int m_nCorrectPhasingCount;
    int m_UnphaseableCount;
    
    bool m_bBaseIsPhased;
    bool m_bBasePhase;
    bool m_bCallIsPhased;
    bool m_bCallPhase;
};



#endif // _C_PHASE_EVALUATOR_H_
