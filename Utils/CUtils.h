//
//  CUtils.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/4/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_UTILS_H_
#define _C_UTILS_H_

#include <string>

class CVariant;


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
    
    ///Return true if the given file path exists
    static bool IsFileExists (const std::string& name);

    
};


#endif /* _C_UTILS_H_ */
