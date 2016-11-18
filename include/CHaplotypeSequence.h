/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_HAPLOTYPE_SEQUENCE_H_
#define _C_HAPLOTYPE_SEQUENCE_H_

#include <deque>
#include "CVariant.h"
#include "COrientedVariant.h"

const int g_nINVALID = -1;

class CHaplotypeSequence
{    
  public:
    //Empty Object
    CHaplotypeSequence();
    
    //Default constructor
    CHaplotypeSequence(const char* a_aRefSequence, int a_nRefSize);

    //Copy constructor
    CHaplotypeSequence(const CHaplotypeSequence& a_rObj);

    //Adds the variant allele specified with the a_nAlleleId to the haplotype
    void AddVariant(const COrientedVariant& a_rVariant);

    //Get the m_nTemplatePosition
    int GetTemplatePosition() const;
    
    //Checks if the given haplotype sequence is equal with this
    bool IsEqual(const CHaplotypeSequence& a_rObj) const;

    //Compare given haplotype sequence with this
    int CompareTo(const CHaplotypeSequence& a_rObj) const;
    
    //Test if the haplotype is currently within a variant
    bool IsOnTemplate() const;

    //Test if there are more nucleotides available
    bool HasNext() const;

    //Step to the next nucleotide
    void Next();

    //Returns the nucleotide on m_nTemplatePosition
    char NextBase() const;

    //Force the template position to the first template position at or beyond "a_nPosition" and the current template position which is not
    // in a variant. Force the state of any otherwise unmarked variants as UNKNOWN (a_nPosition is 0 based).
    void MoveForward(int a_nPosition);

    //Test whether a deficit of variant bases are upstream in the queue in order to perform a step.
    //return false indicates that no variants need to be immediately enqueued
    bool WantsFutureVariantBases() const;
 
    //[TEST Purpose] print the haplotype
    void Print() const;
    
  private:
  // Sorted list of variants yet to be processed
  std::deque<COrientedVariant> m_aVariants;

  //Reference nucleotid sequence
  const char* m_aRefSequence;
  int m_nRefSequenceLength;

  // Position in template (start of current variant if one is active). 0 based
  int m_nTemplatePosition;

  // Position in variant. INVALID if not currently in variant. 0 based
  int m_nPositionInVariant;

  int m_nLastVariantEnd;

  // Variant that currently in or next one.
  COrientedVariant m_nextVariant;

};


#endif //_C_HAPLOTYPE_SEQUENCE_H_
