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

#ifndef _C_ORIENTED_VARIANT_H_
#define _C_ORIENTED_VARIANT_H_

#include "CVariant.h"

namespace core
{

/**
 * @brief A Container that stores the Phasing information of the variant
 *
 * COrientedVariant is an extended version of CVariant which stores an additional phasing data
 * to be used while replaying variant to the corresponding haplotype sequence
 */
class COrientedVariant
{
    public:
    COrientedVariant();
    
    ///Initialize Oriented with a variant and orientation selection
    COrientedVariant(const CVariant& a_rObj, bool a_bIsOrderOfGenotype);
    
    ///Create homozygous oriented variant from given allele index
    COrientedVariant(const CVariant& a_rObj, int a_nAlleleIndex);
    
    ///Copy constructor
    COrientedVariant(const COrientedVariant& a_rObj);    
    
    ///Get the allele alt string
    const SAllele& GetAllele() const;
    
    ///Compare variants according to start/end position
    int CompareTo(const COrientedVariant& a_rObj) const;
    
    ///Gets the start position of the variant
    int GetStartPos() const;
    
    ///Gets the index of the allele
    int GetAlleleIndex() const;
    
    ///Gets the end position of the allele
    int GetEndPos() const;
    
    ///Get the Allele in the reverse side of chosen orientation
    COrientedVariant Other() const;
    
    ///Return if the variant is null
    bool IsNull() const;
    
    ///Return true if hap A is used and return false if hap B is used
    bool IsOrderOfGenotype() const;
    
    ///Return the variant
    const CVariant& GetVariant() const;
    
    ///Set variant to null
    void SetToNull();
    
    ///[For Test purpose] print the oriented variant
    void Print() const;

    private:
    ///Index of the selected allele of this variant
    int m_nAlleleIndex;
    ///Index of the other allele of this variant
    int m_nOtherAlleleIndex;
    ///Pointer Access to variant
    const CVariant* m_variant;
    ///If the selected allele is first number or not (eg.  a/b   a-> true b-> false)
    bool m_bIsOrderOfGenotype;
    ///If oriented variant is null
    bool m_bIsNull;

};

}

#endif //_C_ORIENTED_VARIANT_H_
