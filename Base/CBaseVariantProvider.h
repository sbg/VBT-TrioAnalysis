//
//  CBaseVariantProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/7/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
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
            
    ///Read contig given by the chromosome id
    bool ReadContig(std::string a_chrId, SContig& a_rContig);
    
protected:

    ///Find the optimal Trimming for variant list that have more than 1 trimming options
    void FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, std::vector<std::vector<CVariant>>* a_pAllVarList, const SConfig& a_rConfig);
    
    //REFERENCE FASTA
    CFastaParser m_referenceFasta;


};



#endif /* _C_BASE_VARIANT_PROVIDER_H_ */
