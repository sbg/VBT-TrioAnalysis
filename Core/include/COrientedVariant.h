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
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/OrientedVariant.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 * Licensed under the Simplified BSD License: https://github.com/RealTimeGenomics/rtg-tools/blob/master/LICENSE.txt
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
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
