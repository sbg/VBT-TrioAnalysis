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
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/HaplotypePlayback.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 * Licensed under the Simplified BSD License: https://github.com/RealTimeGenomics/rtg-tools/blob/master/LICENSE.txt
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
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
