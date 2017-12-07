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

namespace duocomparison
{

/**
 * @brief Groups indexes of common chromosomes for truth and query vcfs
 *
 */
class CVariantIteratorGa4gh
{
public:
    
    ///Constructor
    CVariantIteratorGa4gh(std::vector<const core::COrientedVariant*>& included, std::vector<const CVariant*>& excluded, std::vector<const CVariant*>& notAssessed);
    
    ///Return true if the iterator is not at the end of variant list
    bool hasNext();
    
    ///Return next variant(s) that the iterator point. Multiple variant is return in case they are at the same position
    bool FillNext(std::vector<SVariantSummary>& a_rResultsVec);
    
private:
    
    int CompareTrio(int included, int excluded, int notAssessed);
    
    std::vector<const core::COrientedVariant*>::iterator it_Included;
    std::vector<const CVariant*>::iterator it_Excluded;
    std::vector<const CVariant*>::iterator it_NotAssessed;
    
    std::vector<const core::COrientedVariant*>& m_aIncluded;
    std::vector<const CVariant*>& m_aExcluded;
    std::vector<const CVariant*>& m_aNotAssessed;
    
};

}

#endif // _C_CALL_ITERATOR_V2_H_
