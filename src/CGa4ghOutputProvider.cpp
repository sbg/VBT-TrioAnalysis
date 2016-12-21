//
//  CGa4ghOutputProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/20/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CGa4ghOutputProvider.h"
#include <sstream>
#include "CPath.h"

void CGa4ghOutputProvider::FillHeader()
{
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= SBG Vcf Comparison Tool Ver. 1.0 (Beta), 2016");
    
    //ADD REQUIRED FORMATS BY GA4GH
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
    
    //ADD CONTIG IDs
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr1,length=249250621>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr2,length=243199373>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr3,length=198022430>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr4,length=191154276>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr5,length=180915260>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr6,length=171115067>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr7,length=159138663>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr8,length=146364022>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr9,length=141213431>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr10,length=135534747>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr11,length=135006516>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr12,length=133851895>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr13,length=115169878>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr14,length=107349540>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr15,length=102531392>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr16,length=90354753>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr17,length=81195210>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr18,length=78077248>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr19,length=59128983>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr20,length=63025520>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr21,length=48129895>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr22,length=51304566>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrM,length=16571>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrX,length=155270560>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrY,length=59373566>");

    //ADD REQUIRED SAMPLES
    m_vcfWriter.AddSampleName("TRUTH");
    m_vcfWriter.AddSampleName("QUERY");
    
    //CLOSE HEADER
    m_vcfWriter.WriteHeaderToVcf();
}


void CGa4ghOutputProvider::AddRecords(CPath* a_pBestPath)
{
    std::vector<COrientedVariant*> basePaths;
    
    

    
}


void CGa4ghOutputProvider::VariantToVcfRecord(const CVariant& a_rVariant, SVcfRecord& a_rOutputRec, bool a_bIsBase, const std::string& a_rMatchType, const::std::string& a_rDecision)
{
    //Fill basic variant data
    a_rOutputRec.m_chrName = a_rVariant.m_chrName;
    a_rOutputRec.m_nPosition = a_rVariant.m_nStartPos;
    a_rOutputRec.m_alleles = a_rVariant.m_allelesStr;
    a_rOutputRec.m_filterString = a_rVariant.m_filterString;
    
    //If this is called variant create an empty sample data for TRUTH column
    if(a_bIsBase == false)
        a_rOutputRec.m_aSampleData.push_back(SPerSampleData());
    
    //Fill genotype of sample data
    SPerSampleData data;
    data.m_bIsPhased = a_rVariant.m_bIsPhased;
    data.m_nHaplotypeCount = a_rVariant.m_nZygotCount;
    for(int k = 0; k < data.m_nHaplotypeCount; k++)
        data.m_aGenotype[k] = a_rVariant.m_genotype[k];
    a_rOutputRec.m_aSampleData.push_back(data);
}


void CGa4ghOutputProvider::MergeVariants(const CVariant& a_rVariantBase,
                                         const CVariant& a_rVariantCalled,
                                         const std::string& a_rMatchType,
                                         const std::string& a_rDecisionBase,
                                         const std::string& a_rDecisionCalled,
                                         SVcfRecord& a_rOutputRec)
{
    //Fill basic variant data
    a_rOutputRec.m_chrName = a_rVariantBase.m_chrName;
    a_rOutputRec.m_nPosition = a_rVariantBase.m_nStartPos;
    a_rOutputRec.m_alleles = a_rVariantBase.m_allelesStr;
    a_rOutputRec.m_filterString = a_rVariantBase.m_filterString;
    
    //Fill base sample (TRUTH)
    SPerSampleData data;
    data.m_bIsPhased = a_rVariantBase.m_bIsPhased;
    data.m_nHaplotypeCount = a_rVariantBase.m_nZygotCount;
    for(int k = 0; k < data.m_nHaplotypeCount; k++)
        data.m_aGenotype[k] = a_rVariantBase.m_genotype[k];
    data.m_decisionBD = a_rDecisionBase;
    data.m_matchTypeBK = a_rMatchType;
    
    a_rOutputRec.m_aSampleData.push_back(data);
    
    //Fill called sample (QUERY)
    SPerSampleData data2;
    data2.m_bIsPhased = a_rVariantCalled.m_bIsPhased;
    data2.m_nHaplotypeCount = a_rVariantCalled.m_nZygotCount;
    
    std::stringstream ss(a_rVariantBase.m_allelesStr);
    std::vector<std::string> baseVariants;
        
    while(ss.good())
    {
        std::string substr;
        getline(ss, substr, ',');
        baseVariants.push_back(substr);
    }
    
    for(int k=0; k < a_rVariantCalled.m_nZygotCount; k++)
    {
        for(int p = 0; p < baseVariants.size(); p++)
        {
            if(0 == a_rVariantCalled.m_alleles[k].m_sequence.compare(baseVariants[p]))
            {
                data2.m_aGenotype[k] = p;
            }
        }
    }

    
    data2.m_decisionBD = a_rDecisionCalled;
    data2.m_matchTypeBK = a_rMatchType;
    a_rOutputRec.m_aSampleData.push_back(data2);
    
}


bool CGa4ghOutputProvider::CanMerge(const CVariant &a_rVariantBase, const CVariant &a_rVariantCalled)
{
    bool bIsPosEqual = a_rVariantBase.m_nStartPos == a_rVariantCalled.m_nStartPos;
    bool bIsRefEqual = a_rVariantBase.m_refSequence == a_rVariantCalled.m_refSequence;
    
    if(bIsPosEqual && bIsRefEqual)
    {
        
        std::stringstream ss(a_rVariantBase.m_allelesStr);
        std::vector<std::string> baseVariants;
    
        while(ss.good())
        {
            std::string substr;
            getline(ss, substr, ',');
            baseVariants.push_back(substr);
        }
    
        bool bIsAlleleExists;
    
        for(int k=0; k < a_rVariantCalled.m_nAlleleCount; k++)
        {
        
            bIsAlleleExists = false;
            for(int p = 0; p < baseVariants.size(); p++)
            {
                if(0 == a_rVariantCalled.m_alleles[k].m_sequence.compare(baseVariants[p]))
                {
                    bIsAlleleExists = true;
                    break;
                }
            }
        
            if(false == bIsAlleleExists)
                 return false;
        }
        
        return true;
    }
    else
        return false;
}

































