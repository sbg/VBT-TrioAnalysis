//
//  CSplitOutputProvider.h
//  VCFComparison
//
//  Created by Berke.Toptas on 4/7/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_SPLIT_OUTPUT_PROVIDER_H_
#define _C_SPLIT_OUTPUT_PROVIDER_H_

#include "CVcfWriter.h"
#include "CVcfReader.h"
#include "SChrIdTuple.h"

class COrientedVariant;
class CVariantProvider;
class CPath;

class CSplitOutputProvider
{

public:
    
    //Default constructor
    CSplitOutputProvider();
    
    //Set access to variant provider
    void SetVariantProvider(CVariantProvider* a_pProvider);
    
    //Set access to best path list
    void SetBestPaths(std::vector<CPath>& a_rBestPathList);
    
    //Set the output vcfs path FOLDER
    void SetVcfPath(const std::string& a_rVcfPath);
    
    //Set contigs [id, name and length] to write output header
    void SetContigList(const std::vector<SVcfContig>& a_rContigs);
    
    //Generates 4 vcf files splitting each variant decisions
    void GenerateSplitVcfs(const std::vector<SChrIdTuple>& a_rCommonChromosomes);

    
private:

    //Add the given variant list to the given writer vcf
    void AddRecords(CVcfWriter* a_pWriter, const std::vector<const COrientedVariant*>& a_pOvarList);
    void AddRecords(CVcfWriter* a_pWriter, const std::vector<const CVariant*>& a_pVarList);
    
    //Convert CVariant/COrientedVariant to vcf record
    void VariantToVcfRecord(const COrientedVariant* a_pOvar, SVcfRecord& a_rVcfRecord);
    void VariantToVcfRecord(const CVariant* a_pVar, SVcfRecord& a_rVcfRecord);
    
    //Fill the header of given writers
    //@a_bIsBaseSide : if we take filter names from base vcf or called vcf
    void FillHeader(CVcfWriter *a_pWriter, bool a_bIsBaseSide);
    
    //Generates the vcf files of each category
    void GenerateTpBaseVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes);
    void GenerateTpCalledVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes);
    void GenerateFnVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes);
    void GenerateFpVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes);
    
    //Path of output folder where we place vcf files
    std::string m_vcfsFolder;

    //Vcf Writer objects of each decision category
    CVcfWriter m_TPBaseWriter;
    CVcfWriter m_TPCalledWriter;
    CVcfWriter m_FPWriter;
    CVcfWriter m_FNWriter;
    
    //Access to variant provider
    CVariantProvider* m_pProvider;
    
    //Access to best paths
    std::vector<CPath> m_aBestPaths;
    
    //Contig list which will be used to fill header part of output vcf
    std::vector<SVcfContig> m_contigs;
    
};

#endif // _C_SPLIT_OUTPUT_PROVIDER_H_
