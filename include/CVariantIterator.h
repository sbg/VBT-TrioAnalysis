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


class CVariantIterator
{
public:
    
    CVariantIterator(std::vector<const COrientedVariant*>& included, std::vector<CVariant*>& excluded, std::vector<CVariant>& notAssessed)
    : m_aExcluded(excluded),
      m_aIncluded(included),
      m_aNotAssessed(notAssessed)
    {
        it_Included = m_aIncluded.begin();
        it_Excluded = m_aExcluded.begin();
        it_NotAssessed = m_aNotAssessed.begin();
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
    
    int CompareTrio(int included, int excluded, int notAssessed)
    {
        if(included >= excluded && included >= notAssessed)
            return 0;
        else if(excluded >= included && excluded >= notAssessed)
            return 1;
        else
            return 2;
    }
    
    std::vector<const COrientedVariant*>::iterator it_Included;
    std::vector<CVariant*>::iterator it_Excluded;
    std::vector<CVariant>::iterator it_NotAssessed;
    
    std::vector<const COrientedVariant*>& m_aIncluded;
    std::vector<CVariant*>& m_aExcluded;
    
    std::vector<CVariant>& m_aNotAssessed;
    
};


#endif // _C_CALL_ITERATOR_H_
