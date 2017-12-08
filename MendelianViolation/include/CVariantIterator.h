//
//  CVariantIterator.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 1/13/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

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
