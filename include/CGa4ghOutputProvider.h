//
//  CGa4ghOutputProvider.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/20/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_GA4GH_OUTPUT_PROVIDER_H_
#define _C_GA4GH_OUTPUT_PROVIDER_H_

#include "CVcfWriter.h"

class CPath;
class CVariantProvider;

class CGa4ghOutputProvider
{

public:
 
    //Fill the header of output vcf according to ga4gh standards
    void FillHeader();
    
    //Add records of best path to the vcf file (Contain single chromosome)
    void AddRecords(CPath* a_pBestPath);
    
    
private:
    
    
    
    //Return true if base and called variants can be merged
    bool CanMerge(const CVariant& a_rVariantBase, const CVariant& a_rVariantCalled);
    
    //Merge the two variant and fills the outputRec
    void MergeVariants(const CVariant& a_rVariantBase,
                       const CVariant& a_rVariantCalled,
                       const std::string& a_rMatchTypeBase,
                       const std::string& a_rMatchTypeCalled,
                       const::std::string& a_rDecision,
                       SVcfRecord& a_rOutputRec);
    
    //Write the content of variant into the output record
    void VariantToVcfRecord(const CVariant& a_rVariant, SVcfRecord& a_rOutputRec, bool a_bIsBase, const std::string& a_rMatchType, const::std::string& a_rDecision);
    
    //Vcf writer instance for output
    CVcfWriter m_vcfWriter;

    
    
    //Pointer to the variant provider
    CVariantProvider* m_pVariantProvider;
    
    //Pointer to the best path list
    CPath* m_pBestPaths;
    
};




#endif //_C_GA4GH_OUTPUT_PROVIDER_H_
