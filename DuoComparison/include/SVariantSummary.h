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
 *  SVariantSummary.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 12/21/16.
 *
 */

#ifndef _S_VARIANT_SUMMARY_H_
#define _S_VARIANT_SUMMARY_H_

#include "CVariant.h"

namespace duocomparison
{

/**
 * @brief Stores Phasing information and Comparison decision information for a CVariant
 *
 */
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
    
    int originalStartPos()
    {
        return m_pVariant->GetOriginalPos();
    }
    
    bool included()
    {
        return m_bIncluded;
    }
    
    bool phase()
    {
        return m_bPhase;
    }
    
    bool isNull()
    {
        return m_bIsNull;
    }
    
};

}

#endif //_S_VARIANY_SUMMARY_H_
