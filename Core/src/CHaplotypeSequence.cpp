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

#include "CHaplotypeSequence.h"
#include <iostream>

using namespace core;

CHaplotypeSequence::CHaplotypeSequence()
{}

CHaplotypeSequence::CHaplotypeSequence(const char* a_aRefSequence, int a_nRefSize) 
: m_aRefSequence(a_aRefSequence),
  m_nRefSequenceLength(a_nRefSize)
{    
    m_nTemplatePosition = -1;
    m_nLastVariantEnd = -1;
    m_nPositionInVariant = -1;
}

CHaplotypeSequence::CHaplotypeSequence(const CHaplotypeSequence& a_rObj)
: m_aVariants(a_rObj.m_aVariants),
  m_aRefSequence(a_rObj.m_aRefSequence),
  m_nRefSequenceLength(a_rObj.m_nRefSequenceLength),
  m_nextVariant(a_rObj.m_nextVariant)
{
    m_nPositionInVariant = a_rObj.m_nPositionInVariant;
    m_nLastVariantEnd = a_rObj.m_nLastVariantEnd;
    m_nTemplatePosition = a_rObj.m_nTemplatePosition;
}

void CHaplotypeSequence::AddVariant(const COrientedVariant& a_rVariant)
{
    
    const SAllele a = a_rVariant.GetAllele();
    
    if (a.m_nStartPos == a.m_nEndPos && a.m_sequence.length() == 0)
    {
        // Adding the opposite side of a pure insert is redundant
        return;
    }
    
    //Allele is null
    if(a.m_nStartPos == -1 || a.m_bIsIgnored == true)
        return;
    
    m_nLastVariantEnd = a_rVariant.GetAllele().m_nEndPos;
    
    if(true == m_nextVariant.IsNull())
    {
        m_nextVariant = a_rVariant;
    }    
    else
    {
        m_aVariants.push_back(a_rVariant);
    }
}

bool CHaplotypeSequence::IsEqual(const CHaplotypeSequence& a_rObj) const
{
    return (CompareTo(a_rObj) == 0);
}

int CHaplotypeSequence::CompareTo(const CHaplotypeSequence& a_rObj) const
{
    
    //Check by Position
    int position = m_nTemplatePosition - a_rObj.m_nTemplatePosition;
    if (position != 0)
        return position;
    
    //Check if next variant exists
    else if (m_nextVariant.IsNull())
    {
        if (a_rObj.m_nextVariant.IsNull())
            return 0;
        else
            return -1;
    }
    else if (a_rObj.m_nextVariant.IsNull())
    {
        return 1;
    }
    
    // Check by next variant position
    int current = m_nextVariant.CompareTo(a_rObj.m_nextVariant);
    if (current != 0)
        return current;
    
    // Check by position in variant
    int varPos = m_nPositionInVariant - a_rObj.m_nPositionInVariant;
    if (varPos != 0)
        return varPos;

    
    std::deque<COrientedVariant>::const_iterator itThis = m_aVariants.begin();
    std::deque<COrientedVariant>::const_iterator itObj = a_rObj.m_aVariants.begin();
    
    while((itThis) != m_aVariants.end())
    {
        if((itObj) == a_rObj.m_aVariants.end())
            return 1;
        
        int future = itThis->CompareTo(*itObj);
        if(future != 0)
            return future;
        
        ++itThis;
        ++itObj;
    }
    
    if(itObj != a_rObj.m_aVariants.end())
        return -1;
        
    return 0;
}

int CHaplotypeSequence::GetTemplatePosition() const
{
    return m_nTemplatePosition;
}


bool CHaplotypeSequence::IsOnTemplate() const
{
    return m_nPositionInVariant == g_nINVALID;
}

bool CHaplotypeSequence::HasNext() const
{
   return m_nTemplatePosition < m_nRefSequenceLength - 1;
}

void CHaplotypeSequence::MoveForward(int a_nPosition)
{
    if (!IsOnTemplate()) 
    {
        std::cout << "Attempt to move forward while still in a variant" << std::endl;
        assert(IsOnTemplate());
    }

    m_nTemplatePosition = a_nPosition - 1;
    Next();
}

char CHaplotypeSequence::NextBase() const
{
    if(m_nPositionInVariant == g_nINVALID)
        return m_nRefSequenceLength > m_nTemplatePosition ? m_aRefSequence[m_nTemplatePosition] : 0;
    else
        return m_nextVariant.GetAllele().m_sequence[m_nPositionInVariant];
}

void CHaplotypeSequence::Next()
{
    if(IsOnTemplate())
    {
        m_nTemplatePosition++;
        if(!m_nextVariant.IsNull() && m_nextVariant.GetAllele().m_nStartPos == m_nTemplatePosition)
        {
            // Position to consume the variant
            m_nPositionInVariant = 0;
        }
    }

    else
    {
        assert(m_nPositionInVariant != -1);
        m_nPositionInVariant++;
    }
    
    if(!(!m_nextVariant.IsNull() || m_nPositionInVariant == -1))
    {
        std::cout << "Next Variant is NULL" << std::endl;
    }
    assert(!m_nextVariant.IsNull() || m_nPositionInVariant == -1);

    if(!m_nextVariant.IsNull())
    {
        while(true)
        {
            
            if(m_nPositionInVariant != static_cast<int>(m_nextVariant.GetAllele().m_sequence.length()))
            {
                // Haven't reached the end of the current variant.
                break;
            }
            else
            {
                // Finished variant, so position for next baseStart consuming next variant from the queue
                m_nTemplatePosition = m_nextVariant.GetAllele().m_nEndPos;
                m_nPositionInVariant = g_nINVALID;
               
                if(!m_aVariants.empty())
                {
                    m_nextVariant = m_aVariants.front();
                    m_aVariants.pop_front();
                }
                else
                {
                    //Set next variant to null
                    m_nextVariant.SetToNull();
                    break; 
                }

                if(m_nTemplatePosition < m_nextVariant.GetAllele().m_nStartPos)
                    break;

                m_nPositionInVariant = 0;

                if(m_nTemplatePosition != m_nextVariant.GetAllele().m_nStartPos)
                {
                    std::cout << "templatePosition=" << m_nTemplatePosition << " varStartPosition=" << m_nextVariant.GetVariant().GetStart() << std::endl;
                    std::cout << "Out of order variants during replay" << std::endl;
                    assert(m_nTemplatePosition == m_nextVariant.GetAllele().m_nStartPos);
                }
            }
        }
    }
}


bool CHaplotypeSequence::IsNew(const COrientedVariant& a_rVar) const
{
    return a_rVar.GetAllele().m_bIsIgnored || a_rVar.GetAllele().m_nStartPos >= m_nLastVariantEnd;
}

bool CHaplotypeSequence::WantsFutureVariantBases() const
{
    if (m_nextVariant.IsNull())
        return true;
    
    if (m_nPositionInVariant != g_nINVALID && m_nPositionInVariant <  static_cast<int>(m_nextVariant.GetAllele().m_sequence.length()) - 1)
        return false;
    
    for(int k= 0; k < static_cast<int>(m_aVariants.size()); k++)
    {
        if(m_aVariants[k].GetAllele().m_sequence.length() > 0)
            return false;
    }
    
    return true;
}

void CHaplotypeSequence::Print() const
{
    std::cout <<"Template Pos:" << m_nTemplatePosition << " Pos in Variant:" << m_nPositionInVariant << " Last Var end:" << m_nLastVariantEnd << std::endl;
    std::cout <<"Variant Cnt:" << m_aVariants.size() << std::endl;
    std::cout <<"Next Variant is " <<(m_nextVariant.IsNull() ? "null" : "not null") << std::endl;
    std::cout <<"Next Variant start:" << (m_nextVariant.IsNull() ? -1 : m_nextVariant.GetAllele().m_nStartPos) << std::endl;
}

