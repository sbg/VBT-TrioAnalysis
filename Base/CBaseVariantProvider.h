//
//  CBaseVariantProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/7/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_BASE_VARIANT_PROVIDER_H_
#define _C_BASE_VARIANT_PROVIDER_H_

#include "EVariantMatch.h"
#include "CVcfReader.h"
#include "CFastaParser.h"
#include <vector>

class CVariant;

namespace core
{
    class COrientedVariant;
}

class CBaseVariantProvider
{
    
public:
    
    ///Set the status of each variant in the given list
    void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
    
    ///Set the status of each variant in the given list
    void SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
            
    //Read contig given by the chromosome id
    bool ReadContig(std::string a_chrId, SContig& a_rContig);

    
protected:

    //REFERENCE FASTA
    CFastaParser m_referenceFasta;


};



#endif /* _C_BASE_VARIANT_PROVIDER_H_ */