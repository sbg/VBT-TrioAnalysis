/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/OrientedVariant.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 * Copyright (c) 2016 Seven Bridges Genomics
 *               2017 SBGD Inc.
 *
 */

#include "COrientedVariant.h"
#include <iostream>

using namespace core;


COrientedVariant::COrientedVariant()
{
    m_nAlleleIndex = -1;
    m_nOtherAlleleIndex = -1;
    m_bIsNull = true;
}

COrientedVariant::COrientedVariant(const CVariant& a_rObj, bool a_bIsOrder)
{
    m_variant = &a_rObj;
    m_bIsNull = a_rObj.IsNull();
    
    if(a_rObj.IsHeterozygous())
    {
        if(a_bIsOrder == true)
        {
            m_nAlleleIndex = 0;
            m_nOtherAlleleIndex = 1;
        }
    
        else
        {
            m_nAlleleIndex = 1;
            m_nOtherAlleleIndex = 0;
        }
    }
    else
    {
        m_nAlleleIndex = 0;
        m_nOtherAlleleIndex = 0;
    }

    m_bIsOrderOfGenotype = a_bIsOrder;
}

COrientedVariant::COrientedVariant(const CVariant& a_rObj, int a_nAlleleIndex)
{
    m_variant = &a_rObj;
    m_bIsNull = a_rObj.IsNull();
    m_bIsOrderOfGenotype = true;
    m_nAlleleIndex = a_nAlleleIndex;
    m_nOtherAlleleIndex = a_nAlleleIndex;
}


COrientedVariant::COrientedVariant(const COrientedVariant& a_rObj)
{
    m_variant = a_rObj.m_variant;
    m_nAlleleIndex = a_rObj.m_nAlleleIndex;
    m_nOtherAlleleIndex = a_rObj.m_nOtherAlleleIndex;
    m_bIsOrderOfGenotype = a_rObj.m_bIsOrderOfGenotype;
    m_bIsNull = a_rObj.m_bIsNull;
}


int COrientedVariant::CompareTo(const COrientedVariant& a_rObj) const
{
    //decide by variant id
    int id = (m_variant->GetId() < a_rObj.m_variant->GetId()) ? -1 : ((m_variant->GetId() == a_rObj.m_variant->GetId()) ? 0 : 1);
    if(id != 0)
        return id;
    
    //decide by genotype order
    int genotype = (m_bIsOrderOfGenotype == a_rObj.m_bIsOrderOfGenotype) ? 0 : (m_bIsOrderOfGenotype ? 1 : -1);
    if(genotype != 0)
        return genotype;
    
    // Decide by allele index
    int alleleId = (m_nAlleleIndex < a_rObj.m_nAlleleIndex) ? -1 : ((m_nAlleleIndex == a_rObj.m_nAlleleIndex) ? 0 : 1);
    if (alleleId != 0)
        return alleleId;
    
    //Decide by other allele index
    return (m_nOtherAlleleIndex < a_rObj.m_nOtherAlleleIndex) ? -1 : ((m_nOtherAlleleIndex == a_rObj.m_nOtherAlleleIndex) ? 0 : 1);
}

const SAllele& COrientedVariant::GetAllele() const
{
    return m_variant->m_alleles[m_nAlleleIndex];
}

int COrientedVariant::GetStartPos() const
{
    return m_variant->GetStart();
}

int COrientedVariant::GetEndPos() const
{
    return m_variant->GetEnd();
}

int COrientedVariant::GetAlleleIndex() const
{
    return m_nAlleleIndex;
}

COrientedVariant COrientedVariant::Other() const
{
    COrientedVariant oVar2;
    
    oVar2.m_bIsOrderOfGenotype = !m_bIsOrderOfGenotype;
    oVar2.m_nAlleleIndex = m_nOtherAlleleIndex;
    oVar2.m_nOtherAlleleIndex = m_nAlleleIndex;
    oVar2.m_variant = m_variant;
    oVar2.m_bIsNull = m_bIsNull;
    return oVar2;
}

bool COrientedVariant::IsNull() const
{
    return m_bIsNull;
}

bool COrientedVariant::IsOrderOfGenotype() const
{
    return m_bIsOrderOfGenotype;
}

void COrientedVariant::SetToNull()
{
    m_bIsNull = true;
}

const CVariant& COrientedVariant::GetVariant() const
{
    return *m_variant;
}


void COrientedVariant::Print() const
{
    std::cout << m_nAlleleIndex << ":" << m_nOtherAlleleIndex << " " << (m_bIsOrderOfGenotype ? "true" : "false");
    std::cout << "    " << m_variant->ToString() << "  GT:" << m_variant->m_genotype[0] << "/" << m_variant->m_genotype[1]  <<std::endl;
}














