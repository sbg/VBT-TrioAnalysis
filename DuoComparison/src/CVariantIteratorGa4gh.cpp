//
//  CVariantIteratorGa4gh.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/6/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CVariantIteratorGa4gh.h"

using namespace duocomparison;

///Constructor
CVariantIteratorGa4gh::CVariantIteratorGa4gh(std::vector<const core::COrientedVariant*>& included, std::vector<const CVariant*>& excluded, std::vector<const CVariant*>& notAssessed)
: m_aIncluded(included),
m_aExcluded(excluded),
m_aNotAssessed(notAssessed)
{
    it_Included = m_aIncluded.begin();
    it_Excluded = m_aExcluded.begin();
    it_NotAssessed = m_aNotAssessed.begin();
}

///Return true if the iterator is not at the end of variant list
bool CVariantIteratorGa4gh::hasNext()
{
    bool hasNext = (it_Included != m_aIncluded.end() || it_Excluded != m_aExcluded.end() || it_NotAssessed != m_aNotAssessed.end());
    return hasNext;
}

///Return next variant(s) that the iterator point. Multiple variant is return in case they are at the same position
bool CVariantIteratorGa4gh::FillNext(std::vector<SVariantSummary>& a_rResultsVec)
{
    SVariantSummary result;
    
    int nLastPosition = -1;
    
    while(true)
    {
        int includeStart = it_Included != m_aIncluded.end() ? (*it_Included)->GetVariant().GetOriginalPos() : INT_MAX;
        int excludeStart = it_Excluded != m_aExcluded.end() ? (*it_Excluded)->GetOriginalPos() : INT_MAX;
        int notAssessedStart = it_NotAssessed != m_aNotAssessed.end() ? (*it_NotAssessed)->GetOriginalPos() : INT_MAX;
        
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
                result = SVariantSummary(*(*it_NotAssessed), false, false);
                it_NotAssessed++;
                a_rResultsVec.push_back(result);
                nLastPosition = notAssessedStart;
                break;
        }
        
    }
    
    return true;
}


int CVariantIteratorGa4gh::CompareTrio(int included, int excluded, int notAssessed)
{
    if(included <= excluded && included <= notAssessed)
        return 1;
    else if(excluded <= included && excluded <= notAssessed)
        return 2;
    else
        return 3;
}
