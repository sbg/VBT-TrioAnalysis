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
 *  CTrioVcfMerge.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 11/8/17.
 *
 */

#ifndef _C_TRIO_VCF_MERGE_H_
#define _C_TRIO_VCF_MERGE_H_

#include <vector>

class CVariant;

namespace mendelian
{

/**
 * @brief Merges given three variant set belongs to the same chromosome
 *
 */
class CTrioVcfMerge
{
    
public:
    
    ///Constructor. Set the variant sets
    CTrioVcfMerge(std::vector<const CVariant*>& a_rMotherVariants, std::vector<const CVariant*>& a_rFatherVariants, std::vector<const CVariant*>& a_rChildVariants);
    
    ///Return the next record variants. If a sample is 0/0 null pointer will be returned for that sample. false will be return when iterator comes to the end
    bool GetNext(const CVariant** a_rMotherVar, const CVariant** a_rFatherVariant,  const CVariant** a_rChildVariant);
    
    
private:
    
    ///If buffers are empty, refill buffers
    bool RefillBuffers();
    
    std::vector<const CVariant*> motherVariantBuffer;
    std::vector<const CVariant*> fatherVariantBuffer;
    std::vector<const CVariant*> childVariantBuffer;

    //Keep track of the smallest original variant position at each iteration
    int smallestVariantPosition;
    
    //Reference to the variant set to be merged
    std::vector<const CVariant*>& m_motherVariants;
    std::vector<const CVariant*>& m_fatherVariants;
    std::vector<const CVariant*>& m_childVariants;
    
    //Sample iterators
    int motherItr, fatherItr, childItr;
};

}

#endif /* _C_TRIO_VCF_MERGE_H_ */
