//
//  CVariantIterator.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 11/12/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#include "CVariantIterator.h"
#include "COrientedVariant.h"

using namespace mendelian;

CVariantIterator::CVariantIterator(const std::vector<const core::COrientedVariant*>& includedGT, const std::vector<const core::COrientedVariant*>& includedAM)
: m_aIncludedGT(includedGT),
m_aIncludedAM(includedAM)
{
    it_IncludedGT = m_aIncludedGT.begin();
    it_IncludedAM = m_aIncludedAM.begin();
}

bool CVariantIterator::hasNext()
{
    bool hasNext = (it_IncludedGT != m_aIncludedGT.end() || it_IncludedAM != m_aIncludedAM.end());
    return hasNext;
}

const core::COrientedVariant*  CVariantIterator::Next()
{
    const core::COrientedVariant* toRet;
    
    if(it_IncludedAM == m_aIncludedAM.end())
    {
        toRet = *it_IncludedGT;
        it_IncludedGT++;
    }
    
    else if(it_IncludedGT != m_aIncludedGT.end() && (*it_IncludedGT)->GetVariant().GetOriginalPos() < (*it_IncludedAM)->GetVariant().GetOriginalPos())
    {
        toRet = *it_IncludedGT;
        it_IncludedGT++;
    }
    
    else if(it_IncludedGT != m_aIncludedGT.end() && (*it_IncludedGT)->GetVariant().GetOriginalPos() == (*it_IncludedAM)->GetVariant().GetOriginalPos())
    {
        toRet = (*it_IncludedGT)->GetVariant().m_nId < (*it_IncludedAM)->GetVariant().m_nId ? *(it_IncludedGT++) : *(it_IncludedAM++);
    }
    
    else
    {
        toRet = *it_IncludedAM;
        it_IncludedAM++;
    }
    
    return toRet;
}
