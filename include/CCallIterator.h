//
//  CCallIterator.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/21/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//
//Combines included/excluded calls into one stream

#ifndef _C_CALL_ITERATOR_H_
#define _C_CALL_ITERATOR_H_

#include <vector>
#include "COrientedVariant.h"
#include "SVariantSummary.h"


class CCallIterator
{
public:
    
    CCallIterator(std::vector<const COrientedVariant*>& included, std::vector<CVariant*>& excluded)
    : m_aExcluded(excluded),
      m_aIncluded(included)
    {
        it_Included = m_aIncluded.begin();
        it_Excluded = m_aExcluded.begin();
    }
    
    bool hasNext()
    {
        return (it_Included != m_aIncluded.end() || it_Excluded != m_aExcluded.end());
    }
    
    SVariantSummary next()
    {
        SVariantSummary result;
        if(it_Included == m_aIncluded.end() || (it_Excluded != m_aExcluded.end() && (*it_Excluded)->GetStart() < (*it_Included)->GetVariant().GetStart()))
        {
            result = SVariantSummary(*(*it_Excluded), false, false);
            it_Excluded++;
        }
        else
        {
            result = SVariantSummary((*it_Included)->GetVariant(), true, (*it_Included)->IsOrderOfGenotype());
            it_Included++;
        }
        
        return result;
    }
    
private:
    
    std::vector<const COrientedVariant*>::iterator it_Included;
    std::vector<CVariant*>::iterator it_Excluded;
    
    std::vector<const COrientedVariant*>& m_aIncluded;
    std::vector<CVariant*>& m_aExcluded;
    
};


#endif // _C_CALL_ITERATOR_H_
