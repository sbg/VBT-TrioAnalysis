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
 *  CVariantIterator.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 1/13/17.
 *
 */

#ifndef _C_CALL_ITERATOR_H_
#define _C_CALL_ITERATOR_H_

#include <vector>

namespace core
{
    class COrientedVariant;
}

namespace mendelian
{

/**
 * @brief A container that stores allele match and genotype match variants and implements Next function to iterate between these two set according to variant order
 *
 */
class CVariantIterator
{
public:
    
    CVariantIterator(const std::vector<const core::COrientedVariant*>& includedGT, const std::vector<const core::COrientedVariant*>& includedAM);
    
    bool hasNext();
    
    const core::COrientedVariant* Next();

    
private:
    
    std::vector<const core::COrientedVariant*>::const_iterator it_IncludedGT;
    std::vector<const core::COrientedVariant*>::const_iterator it_IncludedAM;
    
    const std::vector<const core::COrientedVariant*>& m_aIncludedGT;
    const std::vector<const core::COrientedVariant*>& m_aIncludedAM;
    
};

}

#endif // _C_CALL_ITERATOR_H_
