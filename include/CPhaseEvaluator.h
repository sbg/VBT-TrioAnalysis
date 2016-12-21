/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
 *
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_PHASE_EVALUATOR_H_
#define _C_PHASE_EVALUATOR_H_

class CVariantProvider;
class CPath;

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



#endif //_C_PHASE_EVALUATOR_H_
