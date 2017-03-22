//
//  CCallIteratorV2.h
//  VCFComparison
//
//  Created by Berke.Toptas on 1/13/17.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//
//Combines included/excluded calls into one stream

#ifndef _C_CALL_ITERATOR_V2_H_
#define _C_CALL_ITERATOR_V2_H_

#include <vector>
#include "COrientedVariant.h"
#include "SVariantSummary.h"


class CVariantIteratorV2
{
public:
    
    CVariantIteratorV2(std::vector<const COrientedVariant*>& included, std::vector<const CVariant*>& excluded, std::vector<CVariant>& notAssessed)
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
        bool hasNext = (it_Included != m_aIncluded.end() || it_Excluded != m_aExcluded.end() || it_NotAssessed != m_aNotAssessed.end());
        return hasNext;
    }
    
    bool FillNext(std::vector<SVariantSummary>& a_rResultsVec)
    {
        SVariantSummary result;
        
        int nLastPosition = -1;
        
        while(true)
        {
            int includeStart = it_Included != m_aIncluded.end() ? (*it_Included)->GetVariant().GetOriginalPos() : INT_MAX;
            int excludeStart = it_Excluded != m_aExcluded.end() ? (*it_Excluded)->GetOriginalPos() : INT_MAX;
            int notAssessedStart = it_NotAssessed != m_aNotAssessed.end() ? it_NotAssessed->GetOriginalPos() : INT_MAX;

            //All lists are finished return false
            if(includeStart == INT_MAX && excludeStart == INT_MAX && notAssessedStart == INT_MAX)
                return false;
            
            //All variants are filled with the same position
            if(nLastPosition != -1 && nLastPosition != includeStart && nLastPosition != excludeStart && nLastPosition != notAssessedStart)
                return true;
            
            int index = CompareTrio(includeStart, excludeStart, notAssessedStart);
            
            switch(index)
            {
                case 1:
                    result = SVariantSummary((*it_Included)->GetVariant(), true, (*it_Included)->IsOrderOfGenotype());
                    a_rResultsVec.push_back(result);
                    nLastPosition = includeStart;
                    it_Included++;
                    break;
                case 2:
                    result = SVariantSummary(*(*it_Excluded), false, false);
                    a_rResultsVec.push_back(result);
                    nLastPosition = excludeStart;
                    it_Excluded++;
                    break;
                case 3:
                    result = SVariantSummary(*it_NotAssessed, false, false);
                    it_NotAssessed++;
                    a_rResultsVec.push_back(result);
                    nLastPosition = notAssessedStart;
                    break;
            }

        }
    
        return true;
    }
    
private:
    
    int CompareTrio(int included, int excluded, int notAssessed)
    {
        if(included <= excluded && included <= notAssessed)
            return 1;
        else if(excluded <= included && excluded <= notAssessed)
            return 2;
        else
            return 3;
    }
    
    std::vector<const COrientedVariant*>::iterator it_Included;
    std::vector<const CVariant*>::iterator it_Excluded;
    std::vector<CVariant>::iterator it_NotAssessed;
    
    std::vector<const COrientedVariant*>& m_aIncluded;
    std::vector<const CVariant*>& m_aExcluded;
    
    std::vector<CVariant>& m_aNotAssessed;
    
};


#endif // _C_CALL_ITERATOR_V2_H_
