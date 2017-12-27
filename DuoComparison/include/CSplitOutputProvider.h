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
 *  CSplitOutputProvider.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 4/7/17.
 *
 */

#ifndef _C_SPLIT_OUTPUT_PROVIDER_H_
#define _C_SPLIT_OUTPUT_PROVIDER_H_

#include "CVcfWriter.h"
#include "CVcfReader.h"
#include "SChrIdTuple.h"

namespace core
{
class COrientedVariant;
class CPath;
}

namespace duocomparison
{
class CVariantProvider;

/**
 * @brief Outputs TPBase, TPCalled, FP and FN vcf files as output after comparison operation
 *
 */
class CSplitOutputProvider
{

public:
    
    ///Default constructor
    CSplitOutputProvider();
    
    ///Set access to variant provider
    void SetVariantProvider(CVariantProvider* a_pProvider);
    
    ///Set access to best path list
    void SetBestPaths(std::vector<core::CPath>& a_rBestPathList);
    
    ///Set the output vcfs path FOLDER
    void SetVcfPath(const std::string& a_rVcfPath);
    
    ///Set contigs [id, name and length] to write output header
    void SetContigList(const std::vector<SVcfContig>& a_rContigs);
    
    ///Generates 4 vcf files splitting each variant decisions for given common chromosome list
    void GenerateSplitVcfs(const std::vector<SChrIdTuple>& a_rCommonChromosomes);

    
private:

    //Add the given variant list to the given writer vcf
    void AddRecords(CVcfWriter* a_pWriter, const std::vector<const core::COrientedVariant*>& a_pOvarList);
    void AddRecords(CVcfWriter* a_pWriter, const std::vector<const CVariant*>& a_pVarList);
    
    //Convert CVariant/COrientedVariant to vcf record
    void VariantToVcfRecord(const core::COrientedVariant* a_pOvar, SVcfRecord& a_rVcfRecord);
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
    std::vector<core::CPath> m_aBestPaths;
    
    //Contig list which will be used to fill header part of output vcf
    std::vector<SVcfContig> m_contigs;
    
};
    
}
#endif // _C_SPLIT_OUTPUT_PROVIDER_H_
