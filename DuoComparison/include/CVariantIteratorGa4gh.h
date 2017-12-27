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
 *  CCallIteratorGa4gh.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 1/13/17.
 *
 */

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
