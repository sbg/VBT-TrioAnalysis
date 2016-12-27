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
//#include "CVcfReader.h"

class CPath;
class CVariantProvider;

class CGa4ghOutputProvider
{

public:
 
    //Set access to variant provider
    void SetVariantProvider(CVariantProvider* a_pProvider);
    
    //Set access to best path list
    void SetBestPaths(CPath* a_pBestPathList);
    
    //Set the output vcf path
    void SetVcfPath(const std::string& a_rVcfPath);
    
    //Generates the vcf file by merging all chromosomes
    void GenerateGa4ghVcf();
    
private:
    
    //Fill the header of output vcf according to ga4gh standards
    void FillHeader();
    
    //Add records of best path to the vcf file (Contain single chromosome)
    void AddRecords(const CPath& a_rBestPath, int a_nChrId);
    
    //Return true if base and called variants can be merged
    bool CanMerge(const CVariant* a_pVariantBase, const CVariant* a_pVariantCalled) const;
    
    //Merge the two variant and fills the outputRec
    void MergeVariants(const CVariant* a_rVariantBase,
                       const CVariant* a_rVariantCalled,
                       const std::string& a_rMatchType,
                       const std::string& a_rDecisionBase,
                       const std::string& a_rDecisionCalled,
                       SVcfRecord& a_rOutputRec);
    
    //Write the content of variant into the output record
    void VariantToVcfRecord(const CVariant* a_rVariant, SVcfRecord& a_rOutputRec, bool a_bIsBase, const std::string& a_rMatchType, const::std::string& a_rDecision);
    
    //Vcf writer instance for output
    CVcfWriter m_vcfWriter;
    
    //Vcf reader of baseline
//    CVcfReader m_baseReader;
    //Vcf reader of call
//    CVcfReader m_calledReader;
    

    //Pointer to the variant provider
    CVariantProvider* m_pVariantProvider;
    
    //Pointer to the best path list
    CPath* m_pBestPaths;
    
    //Path of the vcf file to be generated
    std::string m_vcfPath;
    
};


#endif //_C_GA4GH_OUTPUT_PROVIDER_H_
