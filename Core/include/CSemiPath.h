/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/HalfPath.java
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

#ifndef _C_SEMI_PATH_H_
#define _C_SEMI_PATH_H_

#include "CHaplotypeSequence.h"
#include "EVcfName.h"
#include <vector>

namespace core
{

/**
 * @brief A Container that stores the current state of 1 vcf side during variant replay
 *
 * CSemiPath stores the variant replay information of single vcf. Each CSemipath contains two haplotype 
 * (since human is diploid). Variants included/excluded so far is stored at CSemipath level.
 *
 */
class CSemiPath
{
    public:

    CSemiPath();
    CSemiPath(const char* a_aRefSequence, int a_nRefSize, EVcfName a_uVcfName);
    CSemiPath(const CSemiPath& a_rObj);

    ///Include the given variant to this semipath
    void IncludeVariant(const COrientedVariant& a_rVariant, int a_nVariantIndex);
    
    ///Exclude the given variant to this semipath
    void ExcludeVariant(const CVariant& a_rVariant, int a_nVariantIndex);

    ///Compare and returns the template position diffence of hapA and hapB
    int CompareHaplotypePositions() const;

    ///Return vcf name of the semi path
    EVcfName GetVcfName() const;
    
    ///Return the list of excluded variants
    const std::vector<int>& GetExcluded() const;
    
    ///Gets the end position the semipath (max of haplotypeA and haplotypeB)
    int GetPosition() const;
    
    ///Return the end position of last variant added
    int GetVariantEndPosition() const;
    
    ///Return the end position of last included variant
    int GetIncludedVariantEndPosition() const;
    
    ///Return the index of last variant added
    int GetVariantIndex() const;
    
    ///Set the index of last variant added
    void SetVariantIndex(int a_nVariantIndex);
    
    ///Return pointer to included variants
    const std::vector<const COrientedVariant*>& GetIncludedVariants() const;
    
    ///Check whether this half path is fully on the template (i.e. no haplotypes are within a variant)
    bool IsOnTemplate() const;
    
    ///Detects overlapping variants
    bool IsNew(const COrientedVariant& a_rVar) const;
    
    ///Compare this half path with the given half path
    int CompareTo(const CSemiPath& a_rObj) const;
    
    ///Check whether this half path is equal to the given half path
    bool IsEqual(const CSemiPath& a_rObj) const;

    ///Checks if both haplotype A and B is finished
    bool HasFinished() const;

    ///Force semipath to move the given position
    void MoveForward(int a_nPosition);
    
    ///Chech whether this half path matches with the given half path
    bool Matches(const CSemiPath& a_rOther);

    /**
     *Test whether a deficit of variant bases are upstream in the queue in order to perform a step.
     *return false indicates that no variants need to be immediately enqueued
     */
    bool WantsFutureVariantBases() const;
  
    bool FinishedHaplotypeA() const;
    bool FinishedHaplotypeB() const;

    ///Return the next base on the B haplotype
    char NextHaplotypeBBase() const;
    
    ///Return the next base on the A haplotype
    char NextHaplotypeABase() const;

    void StepHaplotypeA();
    void StepHaplotypeB();

    ///Clear the included variants
    void ClearIncludedVariants();
    ///Set the included variants
    void AddIncludedVariants(std::vector<const COrientedVariant*>& a_rIncludedVarList);
    
    ///Clear the included variants
    void ClearExcludedVariants();
    ///Set the included variants
    void AddExcludedVariants(std::vector<int>& a_rExcludedVarList);

    ///Sorts included variants according to variant ids
    void SortIncludedVariants();
    
    ///[TEST Purpose]Print semipath
    void Print() const;
    
    private:

    ///name of the semipath
    EVcfName m_uVcfName;
    ///Index of last variant added
    int m_nVariantIndex;
    ///End of last variant added
    int m_nVariantEndPosition;
    ///Last variant included
    int m_nIncludedVariantEndPosition;

    std::vector<const COrientedVariant*> m_aIncludedVariants;
    std::vector<int> m_aExcludedVariants;

    CHaplotypeSequence m_haplotypeA;
    CHaplotypeSequence m_haplotypeB;

    bool m_bFinishedHapA;
    bool m_bFinishedHapB;

};

}

#endif // _C_SEMI_PATH_H_
