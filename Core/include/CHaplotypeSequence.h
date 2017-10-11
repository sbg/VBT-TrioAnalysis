/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/HaplotypePlayback.java
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

#ifndef _C_HAPLOTYPE_SEQUENCE_H_
#define _C_HAPLOTYPE_SEQUENCE_H_

#include <deque>
#include "CVariant.h"
#include "COrientedVariant.h"

namespace core
{

const int g_nINVALID = -1;

/**
 * @brief Container that stores information of each haplotype during the variant replay processing
 *
 * CHaplotypeSequence is the smallest container of Path-SemiPath structure. Each object stores one haplotype information. 
 * The current position during variant replay and the variant in play is stored in this object
 */
class CHaplotypeSequence
{    
  public:
    ///Empty Object
    CHaplotypeSequence();
    
    /**
     * @brief Default constructor
     *
     * @param a_aRefSequence Pointer to the reference sequence belong to the processed contig
     * @param a_nRefSize Length of the reference sequence given with a_aRefSequence
     */
    CHaplotypeSequence(const char* a_aRefSequence, int a_nRefSize);

    ///Copy constructor
    CHaplotypeSequence(const CHaplotypeSequence& a_rObj);

    ///Adds the variant allele specified phasing to the haplotype
    void AddVariant(const COrientedVariant& a_rVariant);

    ///Get the m_nTemplatePosition
    int GetTemplatePosition() const;
    
    ///Checks if the given haplotype sequence is equal with this
    bool IsEqual(const CHaplotypeSequence& a_rObj) const;

    ///Compare given haplotype sequence with this
    int CompareTo(const CHaplotypeSequence& a_rObj) const;
    
    ///Test if the haplotype is currently within a variant
    bool IsOnTemplate() const;

    ///Test if there are more nucleotides available
    bool HasNext() const;

    ///Step to the next nucleotide
    void Next();

    ///Returns the nucleotide on m_nTemplatePosition
    char NextBase() const;

    /**
     * Force the template position to the first template position at or beyond "a_nPosition" and the current template position which is not
     * in a variant. Force the state of any otherwise unmarked variants as UNKNOWN (a_nPosition is 0 based).
     */
    void MoveForward(int a_nPosition);

    /**
     *Test whether a deficit of variant bases are upstream in the queue in order to perform a step.
     *Return false indicates that no variants need to be immediately enqueued.
     */
    bool WantsFutureVariantBases() const;
 
    ///Detects variant overlaps
    bool IsNew(const COrientedVariant& a_rVar) const;
    
    ///[TEST Purpose] print the haplotype
    void Print() const;
    
  private:
  /// Sorted list of variants yet to be processed
  std::deque<COrientedVariant> m_aVariants;

  ///Reference nucleotid sequence
  const char* m_aRefSequence;
  int m_nRefSequenceLength;

  /// Position in template (start of current variant if one is active). 0 based
  int m_nTemplatePosition;

  /// Position in variant. INVALID if not currently in variant. 0 based
  int m_nPositionInVariant;

  int m_nLastVariantEnd;

  /// Variant that currently in or next one.
  COrientedVariant m_nextVariant;

};

}

#endif //_C_HAPLOTYPE_SEQUENCE_H_
