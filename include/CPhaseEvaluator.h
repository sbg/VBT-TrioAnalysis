//
//  CPhaseEvaluator.h
//  VCFComparison
//
//  Created by Berke.Toptas on 11/16/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_PHASE_EVALUATOR_H_
#define _C_PHASE_EVALUATOR_H_

#include "CPath.h"

struct SPhasingResult
{
    int m_nMisPhasingCount;
    int m_nCorrectPhasingCount;
    int m_UnphaseableCount;
};

class CPhaseEvaluator
{
    public:
    
    CPhaseEvaluator();

    SPhasingResult CountMisphasings(const CPath& a_rPath);
    
    private:
    
    int m_nMisPhasingCount;
    int m_nCorrectPhasingCount;
    int m_UnphaseableCount;
    
    bool m_bBaseIsPhased;
    bool m_bBasePhase;
    bool m_bCallIsPhased;
    bool m_bCallPhase;
};



#endif // _C_PHASE_EVALUATOR_H_
