//
//  SVariantSummary.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/21/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

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
