//
//  CSplitOutputProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 4/7/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CSplitOutputProvider.h"
#include "CPath.h"
#include "CVariantProvider.h"


CSplitOutputProvider::CSplitOutputProvider()
{

}

void CSplitOutputProvider::SetVariantProvider(CVariantProvider* a_pProvider)
{
    m_pProvider = a_pProvider;
}

void CSplitOutputProvider::SetBestPaths(CPath* a_pBestPathList)
{
    m_pBestPaths = a_pBestPathList;
}


void CSplitOutputProvider::SetVcfPath(const std::string& a_rVcfPath)
{
    m_vcfsFolder = a_rVcfPath;
}

void CSplitOutputProvider::GenerateSplitVcfs()
{
    GenerateTpCalledVcf();
    GenerateTpBaseVcf();
    GenerateFnVcf();
    GenerateFpVcf();
}



void CSplitOutputProvider::GenerateTpBaseVcf()
{
    std::string filePath = m_vcfsFolder + "/TPBase.vcf";
    m_TPBaseWriter.CreateVcf(filePath.c_str());
    
    //Fill the header section
    FillHeader(&m_TPBaseWriter, true);
    
    //Process each chromosome
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        const std::vector<const COrientedVariant*> ovarList = m_pBestPaths[k].m_baseSemiPath.GetIncludedVariants();
        AddRecords(&m_TPBaseWriter, ovarList);
    }
    
    m_TPBaseWriter.CloseVcf();
}


void CSplitOutputProvider::GenerateTpCalledVcf()
{
    std::string filePath = m_vcfsFolder + "/TPCalled.vcf";
    m_TPCalledWriter.CreateVcf(filePath.c_str());
    
    //Fill the header section
    FillHeader(&m_TPCalledWriter, false);
    
    //Process each chromosome
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        const std::vector<const COrientedVariant*> ovarList = m_pBestPaths[k].m_calledSemiPath.GetIncludedVariants();
        AddRecords(&m_TPCalledWriter, ovarList);
    }
    
    m_TPCalledWriter.CloseVcf();
}

void CSplitOutputProvider::GenerateFnVcf()
{
    std::string filePath = m_vcfsFolder + "/FN.vcf";
    m_FNWriter.CreateVcf(filePath.c_str());
    
    //Fill the header section
    FillHeader(&m_FNWriter, true);
    
    //Process each chromosome
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        const std::vector<const CVariant*> varList = m_pProvider->GetVariantList(eBASE, k, m_pBestPaths[k].m_baseSemiPath.GetExcluded());
        AddRecords(&m_FNWriter, varList);
    }
    
    m_FNWriter.CloseVcf();
}

void CSplitOutputProvider::GenerateFpVcf()
{
    std::string filePath = m_vcfsFolder + "/FP.vcf";
    m_FPWriter.CreateVcf(filePath.c_str());
    
    //Fill the header section
    FillHeader(&m_FPWriter, false);
    
    //Process each chromosome
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        const std::vector<const CVariant*> varList = m_pProvider->GetVariantList(eCALLED, k, m_pBestPaths[k].m_calledSemiPath.GetExcluded());
        AddRecords(&m_FPWriter, varList);
    }
    
    m_FPWriter.CloseVcf();
}


void CSplitOutputProvider::VariantToVcfRecord(const CVariant* a_pVariant, SVcfRecord& a_rOutputRec)
{
    //Fill basic variant data
    a_rOutputRec.m_nChrId = a_pVariant->m_nChrId;
    a_rOutputRec.m_nPosition = a_pVariant->m_nOriginalPos;
    a_rOutputRec.m_alleles = a_pVariant->m_allelesStr;
    a_rOutputRec.m_aFilterString = a_pVariant->m_filterString;
    
    //Fill genotype of sample data
    SPerSampleData data;
    data.m_bIsPhased = a_pVariant->m_bIsPhased;
    data.m_nHaplotypeCount = a_pVariant->m_nZygotCount;
    for(int k = 0; k < data.m_nHaplotypeCount; k++)
        data.m_aGenotype[k] = a_pVariant->m_genotype[k];
    
    a_rOutputRec.m_aSampleData.push_back(data);
}


void CSplitOutputProvider::AddRecords(CVcfWriter* a_pWriter, const std::vector<const COrientedVariant*>& a_pOvarList)
{
    for(const COrientedVariant* pOvar : a_pOvarList)
    {
        SVcfRecord record;
        VariantToVcfRecord(&pOvar->GetVariant(), record);
        a_pWriter->AddRecord(record);
    }
}


void CSplitOutputProvider::AddRecords(CVcfWriter* a_pWriter, const std::vector<const CVariant*>& a_pVarList)
{
    for(const CVariant* pVar : a_pVarList)
    {
        SVcfRecord record;
        VariantToVcfRecord(pVar, record);
        a_pWriter->AddRecord(record);
    }
}



void CSplitOutputProvider::FillHeader(CVcfWriter *a_pWriter, bool a_bIsBaseSide)
{
    //INIT VCF HEADER
    a_pWriter->InitHeader();
    a_pWriter->AddHeaderLine("##source= SBG Vcf Comparison Tool Ver. 1.0 (Beta), 2016");
    
    //ADD REQUIRED FORMATS BY GA4GH
    a_pWriter->AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    a_pWriter->AddHeaderLine("##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
    a_pWriter->AddHeaderLine("##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
    
    //ADD FILTER COLUMNS FROM CALL FILE
    std::vector<std::string> filterNames;
    std::vector<std::string> filterDescriptions;
    m_pProvider->GetFilterInfo((a_bIsBaseSide ? eBASE : eCALLED), filterNames, filterDescriptions);
    for(int k = 1; k < (int)filterNames.size(); k++)
        a_pWriter->AddHeaderLine("##FILTER=<ID=" + filterNames[k] + ",Description=" + filterDescriptions[k] + ">");
    
    //ADD CONTIG IDs
    a_pWriter->AddHeaderLine("##contig=<ID=chr1,length=249250621>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr2,length=243199373>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr3,length=198022430>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr4,length=191154276>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr5,length=180915260>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr6,length=171115067>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr7,length=159138663>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr8,length=146364022>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr9,length=141213431>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr10,length=135534747>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr11,length=135006516>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr12,length=133851895>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr13,length=115169878>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr14,length=107349540>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr15,length=102531392>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr16,length=90354753>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr17,length=81195210>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr18,length=78077248>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr19,length=59128983>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr20,length=63025520>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr21,length=48129895>");
    a_pWriter->AddHeaderLine("##contig=<ID=chr22,length=51304566>");
    a_pWriter->AddHeaderLine("##contig=<ID=chrM,length=16571>");
    a_pWriter->AddHeaderLine("##contig=<ID=chrX,length=155270560>");
    a_pWriter->AddHeaderLine("##contig=<ID=chrY,length=59373566>");
    
    //ADD REQUIRED SAMPLES
    a_pWriter->AddSampleName("S");
    
    //CLOSE HEADER
    a_pWriter->WriteHeaderToVcf();

}




