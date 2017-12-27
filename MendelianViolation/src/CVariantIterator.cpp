/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  CVariantIterator.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 1/13/17.
 *
 */

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
