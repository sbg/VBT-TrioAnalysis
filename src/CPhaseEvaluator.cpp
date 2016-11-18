//
//  CPhaseEvaluator.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 11/16/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CPhaseEvaluator.h"



CPhaseEvaluator::CPhaseEvaluator()
{

    m_nMisPhasingCount = 0;
    m_UnphaseableCount = 0;
    m_nCorrectPhasingCount = 0;
    
    m_bBasePhase = false;
    m_bCallPhase = false;
    m_bBaseIsPhased = false;
    m_bCallIsPhased = false;
}

SPhasingResult CPhaseEvaluator::CountMisphasings(const CPath& a_rPath)
{
    SPhasingResult res;

    //TODO: Implement!!


    return res;
}


