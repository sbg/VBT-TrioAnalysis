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
 *  CUtils.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 12/4/17.
 *
 */

#ifndef _C_UTILS_H_
#define _C_UTILS_H_

#include <string>

class CVariant;
namespace core
{
    class COrientedVariant;
}


class CUtils
{
    
public:

    ///Checks if given two range is overlapping. End position of the range is exclusive
    static bool IsOverlap(int left1, int right1, int left2, int right2);

    ///Checks whether given variant is a structural variant type (A complex type)
    static bool IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength);

    ///Compare 2 variants according to variant start end positions and variant id
    static bool CompareVariants(const CVariant& var1, const CVariant& var2);

    ///Checks if the variant genotype is 0/0 (homref)
    static bool IsHomRef(const CVariant& a_rVariant);
    
    ///Compare variants according to id for sort operation
    static bool CompareVariantsById(const CVariant* v1, const CVariant* v2);
    
    ///Compare  oriented variants according to id for sort operation
    static bool CompareOrientedVariantsById(const core::COrientedVariant* v1, const core::COrientedVariant* v2);
    
    ///Return true if the given file path exists
    static bool IsFileExists (const std::string& name);
    
    ///Return true if the given directory path exists
    static bool IsDirectoryExists(const std::string& path);
    
};


#endif /* _C_UTILS_H_ */
