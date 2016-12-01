#include "CHaplotypeSequence.h"
#include <iostream>

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
: m_nextVariant(a_rObj.m_nextVariant),
  m_nRefSequenceLength(a_rObj.m_nRefSequenceLength),
  m_aRefSequence(a_rObj.m_aRefSequence),
  m_aVariants(a_rObj.m_aVariants)
{
    m_nPositionInVariant = a_rObj.m_nPositionInVariant;
    m_nLastVariantEnd = a_rObj.m_nLastVariantEnd;
    m_nTemplatePosition = a_rObj.m_nTemplatePosition;
}

void CHaplotypeSequence::AddVariant(const COrientedVariant& a_rVariant)
{
    assert(a_rVariant.GetStartPos() > m_nTemplatePosition);
    
    const SAllele a = a_rVariant.GetAllele();

    if (a.m_nStartPos == a.m_nEndPos && a.m_sequence.length() == 0)
    {
        // Adding the opposite side of a pure insert is redundant
        //std::cout << "Alt string is empty" << std::endl;
        return;
    }
    
    //Allele is null
    if(a.m_nStartPos == -1)
        return;
    
    m_nLastVariantEnd = a_rVariant.GetEndPos();
    
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
      return;
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
    
    assert(!m_nextVariant.IsNull() || m_nPositionInVariant == -1);

    if(!m_nextVariant.IsNull())
    {
        while(1)
        {
            if(m_nPositionInVariant != m_nextVariant.GetAllele().m_sequence.length())
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
                //std::cout << "templatePosition=" << m_nTemplatePosition << " varStartPosition=" << m_nextVariant.GetVariant().GetStart() << std::endl;

                if(m_nTemplatePosition != m_nextVariant.GetAllele().m_nStartPos)
                {
                    std::cout << "templatePosition=" << m_nTemplatePosition << " varStartPosition=" << m_nextVariant.GetVariant().GetStart() << std::endl;
                    std::cout << "Out of order variants during replay" << std::endl;
                }
            }
        }
    }
}


bool CHaplotypeSequence::IsNew(const COrientedVariant& a_rVar) const
{
    return a_rVar.IsNull() || a_rVar.GetAllele().m_nStartPos > m_nLastVariantEnd;
}

bool CHaplotypeSequence::WantsFutureVariantBases() const
{
    if (m_nextVariant.IsNull())
        return true;
    
    if (m_nPositionInVariant != g_nINVALID && m_nPositionInVariant < m_nextVariant.GetAllele().m_sequence.length() - 1)
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

