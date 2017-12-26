//
//  CGa4ghOutputProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/20/16.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#ifndef _C_GA4GH_OUTPUT_PROVIDER_H_
#define _C_GA4GH_OUTPUT_PROVIDER_H_

#include "CVcfWriter.h"
#include "CVcfReader.h"
#include "SChrIdTuple.h"
#include "SVariantSummary.h"

namespace core
{
class CPath;
}

namespace duocomparison
{

class CVariantProvider;

/**
 * @brief Output merged baseline and called vcf annotated with comparison decisions in ga4gh benchmarking output standards
 */
class CGa4ghOutputProvider
{

public:
 
    ///Set access to variant provider
    void SetVariantProvider(CVariantProvider* a_pProvider);
    
    ///Set access to best path list
    void SetBestPaths(std::vector<core::CPath>& a_rBestPathList, std::vector<core::CPath>& a_rBestAlleleMatchPathList);
    
    ///Set the output vcf path
    void SetVcfPath(const std::string& a_rVcfPath);
    
    ///Set contigs [id, name and length] to write output header
    void SetContigList(const std::vector<SVcfContig>& a_rContigs);
    
    ///Generates the vcf file by merging all chromosomes
    void GenerateGa4ghVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes);
    
private:
    
    //Return the match string
    std::string GetMatchStr(EVariantMatch a_match);
    
    //Fill the header of output vcf according to ga4gh standards
    void FillHeader();
    
    //Add records of best path to the vcf file (Contain single chromosome)
    void AddRecords(const core::CPath& a_rBestPath, SChrIdTuple a_rTuple);

    //Generate and two sample record with one of them is empty
    void AddSingleSampleRecord(const SVariantSummary& a_rVariant, bool a_bIsBase);
    
    //Return true if base and called variants can be merged
    bool CanMerge(const CVariant* a_pVariantBase, const CVariant* a_pVariantCalled) const;
    
    //Merge the two variant and fills the outputRec
    void MergeVariants(const CVariant* a_rVariantBase,
                       const CVariant* a_rVariantCalled,
                       const std::string& a_rMatchTypeBase,
                       const std::string& a_rMatchTypeCalled,
                       const std::string& a_rDecisionBase,
                       const std::string& a_rDecisionCalled,
                       SVcfRecord& a_rOutputRec);
    
    //Write the content of variant into the output record
    void VariantToVcfRecord(const CVariant* a_rVariant,
                            SVcfRecord& a_rOutputRec,
                            bool a_bIsBase,
                            const std::string& a_rMatchType,
                            const::std::string& a_rDecision);
    
    //Vcf writer instance for output
    CVcfWriter m_vcfWriter;
    
    //Vcf reader of baseline
    CVcfReader m_baseReader;
    //Vcf reader of call
    CVcfReader m_calledReader;
    
    //Pointer to the variant provider
    CVariantProvider* m_pVariantProvider;
    
    //Pointer to the best path list
    std::vector<core::CPath> m_aBestPaths;
    std::vector<core::CPath> m_aBestAlleleMatchPaths;
    
    //Path of the vcf file to be generated
    std::string m_vcfPath;
    
    //Contig id list to write output header
    std::vector<SVcfContig> m_contigs;
    
};

}

#endif //_C_GA4GH_OUTPUT_PROVIDER_H_
