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
#include "CVariantIterator.h"
#include "CVariantProvider.h"
#include <iostream>

void CGa4ghOutputProvider::SetVariantProvider(CVariantProvider* a_pProvider)
{
    m_pVariantProvider = a_pProvider;
}

void CGa4ghOutputProvider::SetBestPaths(CPath* a_pBestPathList, CPath* a_pBestAlleleMatchPathList)
{
    m_pBestPaths = a_pBestPathList;
    m_pBestAlleleMatchPaths = a_pBestAlleleMatchPathList;
}

void CGa4ghOutputProvider::SetVcfPath(const std::string& a_rVcfPath)
{
    m_vcfPath = a_rVcfPath;
}

void CGa4ghOutputProvider::GenerateGa4ghVcf()
{
    m_vcfWriter.CreateVcf(m_vcfPath.c_str());
    FillHeader();
    
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
        AddRecords(m_pBestPaths[k], m_pBestAlleleMatchPaths[k], k);
    
    m_vcfWriter.CloseVcf();
}

void CGa4ghOutputProvider::FillHeader()
{
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= SBG Vcf Comparison Tool Ver. 1.0 (Beta), 2016");
    
    //ADD REQUIRED FORMATS BY GA4GH
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
    
    //ADD FILTER COLUMNS FROM CALL FILE
    std::vector<std::string> filterNames;
    std::vector<std::string> filterDescriptions;
    m_pVariantProvider->GetFilterInfo(eCALLED, filterNames, filterDescriptions);
    for(int k = 1; k < filterNames.size(); k++)
        m_vcfWriter.AddHeaderLine("##FILTER=<ID=" + filterNames[k] + ",Description=" + filterDescriptions[k] + ">");
    
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


void CGa4ghOutputProvider::AddRecords(const CPath& a_rBestPath, const CPath& a_rBestAlleleMatchPath, int a_nChrId)
{
    
    //Best Path included/excluded variants
    std::vector<const CVariant*> excludedVarsBase = m_pVariantProvider->GetVariantList(eBASE, a_nChrId, a_rBestPath.m_baseSemiPath.GetExcluded());
    std::vector<const CVariant*> excludedVarsCall = m_pVariantProvider->GetVariantList(eCALLED, a_nChrId,a_rBestPath.m_calledSemiPath.GetExcluded());
    
    std::vector<const COrientedVariant*> includedVarsBase = a_rBestPath.m_baseSemiPath.GetIncludedVariants();
    std::vector<const COrientedVariant*> includedVarsCall = a_rBestPath.m_calledSemiPath.GetIncludedVariants();
    
    int cntGenMatch = 0;
    int cntNoMatch = 0;
    int cntAlleleMatch = 0;
    
    for(const COrientedVariant* pOvar : includedVarsBase)
    {
        if(pOvar->GetVariant().m_variantStatus == eGENOTYPE_MATCH)
            cntGenMatch++;
        else if(pOvar->GetVariant().m_variantStatus == eNO_MATCH)
            cntNoMatch++;
        else if(pOvar->GetVariant().m_variantStatus == eALLELE_MATCH)
            cntAlleleMatch++;
    }
    for(const CVariant* pVar : excludedVarsBase)
    {
        if(pVar->m_variantStatus == eALLELE_MATCH)
            cntAlleleMatch++;
    }
 
    int cntGenMatch2 = 0;
    int cntNoMatch2 = 0;
    int cntAlleleMatch2 = 0;
    
    for(const COrientedVariant* pOvar : includedVarsCall)
    {
        if(pOvar->GetVariant().m_variantStatus == eGENOTYPE_MATCH)
            cntGenMatch2++;
        else if(pOvar->GetVariant().m_variantStatus == eNO_MATCH)
            cntNoMatch2++;
        else if(pOvar->GetVariant().m_variantStatus == eALLELE_MATCH)
            cntAlleleMatch2++;
    }
    for(const CVariant* pVar : excludedVarsCall)
    {
        if(pVar->m_variantStatus == eALLELE_MATCH)
            cntAlleleMatch2++;
    }
    
    //Not Asessed variants
    std::vector<CVariant> notAssessedBase = m_pVariantProvider->GetNotAssessedVariantList(eBASE, a_nChrId);
    std::vector<CVariant> notAssessedCalled = m_pVariantProvider->GetNotAssessedVariantList(eCALLED, a_nChrId);
    
    CVariantIterator baseVariants(includedVarsBase, excludedVarsBase, notAssessedBase);
    CVariantIterator calledVariants(includedVarsCall, excludedVarsCall, notAssessedCalled);

    SVariantSummary varBase;
    SVariantSummary varCalled;
    
    varBase = baseVariants.hasNext() ? baseVariants.next() : SVariantSummary();
    varCalled = calledVariants.hasNext() ? calledVariants.next() : SVariantSummary();
    
    while(!varBase.isNull() || !varCalled.isNull())
    {
        int baseStart = varBase.isNull() ? 0 : varBase.m_pVariant->m_nStartPos - (varBase.m_pVariant->m_bIsFirstNucleotideTrimmed ? 1 : 0);
        int callStart = varCalled.isNull() ? 0 : varCalled.m_pVariant->m_nStartPos - (varCalled.m_pVariant->m_bIsFirstNucleotideTrimmed ? 1 : 0);
        
        if(!varBase.isNull() && !varCalled.isNull() && CanMerge(varBase.m_pVariant, varCalled.m_pVariant))
        {
            SVcfRecord record;
            std::string decisionBase = varBase.m_bIncluded ? "TP" : (varBase.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FN");
            std::string decisionCalled = varCalled.m_bIncluded ? "TP" : (varCalled.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FP");
            std::string matchBase = GetMatchStr(varBase.m_pVariant->m_variantStatus);
            std::string matchCalled = GetMatchStr(varCalled.m_pVariant->m_variantStatus);
            MergeVariants(varBase.m_pVariant, varCalled.m_pVariant, matchBase, matchCalled, decisionBase, decisionCalled, record);
            m_vcfWriter.AddRecord(record);
            
            //We used both variant. Get the next variant from both list
            varBase = baseVariants.hasNext() ? baseVariants.next() : SVariantSummary();
            varCalled = calledVariants.hasNext() ? calledVariants.next() : SVariantSummary();
            continue;
        }
        
        else if(varCalled.isNull() || (!varBase.isNull() && baseStart < callStart))
        {
            SVcfRecord record;
            std::string decision = varBase.m_bIncluded ? "TP" : (varBase.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FN");
            std::string match = GetMatchStr(varBase.m_pVariant->m_variantStatus);
            VariantToVcfRecord(varBase.m_pVariant, record, true, match, decision);
            m_vcfWriter.AddRecord(record);
            varBase = baseVariants.hasNext() ? baseVariants.next() : SVariantSummary();
            continue;
        }
        
        else
        {
            SVcfRecord record;
            record.m_aSampleData.push_back(SPerSampleData());
            std::string decision = varCalled.m_bIncluded ? "TP" : (varCalled.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FP");
            std::string match = GetMatchStr(varCalled.m_pVariant->m_variantStatus);
            VariantToVcfRecord(varCalled.m_pVariant, record, false, match, decision);
            m_vcfWriter.AddRecord(record);
            varCalled = calledVariants.hasNext() ? calledVariants.next() : SVariantSummary();
            continue;
        }
    }
}

void CGa4ghOutputProvider::VariantToVcfRecord(const CVariant* a_pVariant, SVcfRecord& a_rOutputRec, bool a_bIsBase, const std::string& a_rMatchType, const::std::string& a_rDecision)
{
    //Fill basic variant data
    a_rOutputRec.m_chrName = a_pVariant->m_chrName;
    a_rOutputRec.m_nPosition = a_pVariant->m_nStartPos - (a_pVariant->m_bIsFirstNucleotideTrimmed ? 1 : 0);
    a_rOutputRec.m_alleles = a_pVariant->m_allelesStr;
    if(!a_bIsBase)
        a_rOutputRec.m_aFilterString = a_pVariant->m_filterString;
    
    //Fill genotype of sample data
    SPerSampleData data;
    data.m_bIsPhased = a_pVariant->m_bIsPhased;
    data.m_nHaplotypeCount = a_pVariant->m_nZygotCount;
    for(int k = 0; k < data.m_nHaplotypeCount; k++)
        data.m_aGenotype[k] = a_pVariant->m_genotype[k];
    data.m_decisionBD = a_rDecision;
    if(a_rMatchType != "nm")
        data.m_matchTypeBK = a_rMatchType;
    a_rOutputRec.m_aSampleData.push_back(data);
}

void CGa4ghOutputProvider::MergeVariants(const CVariant* a_pVariantBase,
                                         const CVariant* a_pVariantCalled,
                                         const std::string& a_rMatchTypeBase,
                                         const std::string& a_rMatchTypeCalled,
                                         const std::string& a_rDecisionBase,
                                         const std::string& a_rDecisionCalled,
                                         SVcfRecord& a_rOutputRec)
{

    //Fill basic variant data
    a_rOutputRec.m_chrName = a_pVariantBase->m_chrName;
    a_rOutputRec.m_nPosition = a_pVariantBase->m_nStartPos - (a_pVariantBase->m_bIsFirstNucleotideTrimmed ? 1 : 0);
    a_rOutputRec.m_alleles = a_pVariantBase->m_allelesStr;
    a_rOutputRec.m_aFilterString = a_pVariantCalled->m_filterString;
    
    //Fill base sample (TRUTH)
    SPerSampleData data;
    data.m_bIsPhased = a_pVariantBase->m_bIsPhased;
    data.m_nHaplotypeCount = a_pVariantBase->m_nZygotCount;
    for(int k = 0; k < data.m_nHaplotypeCount; k++)
        data.m_aGenotype[k] = a_pVariantBase->m_genotype[k];
    data.m_decisionBD = a_rDecisionBase;
    if(a_rMatchTypeBase != "nm")
        data.m_matchTypeBK = a_rMatchTypeBase;
    
    a_rOutputRec.m_aSampleData.push_back(data);
    
    //Fill called sample (QUERY)
    SPerSampleData data2;
    data2.m_bIsPhased = a_pVariantCalled->m_bIsPhased;
    data2.m_nHaplotypeCount = a_pVariantCalled->m_nZygotCount;
    
    std::stringstream ss(a_pVariantBase->m_allelesStr);
    std::vector<std::string> baseVariants;
    std::string substr;
    
    while(getline(ss, substr, ','))
    {
        baseVariants.push_back(substr);
    }
    for(int k=0; k < a_pVariantCalled->m_nZygotCount; k++)
    {
        for(int p = 0; p < baseVariants.size(); p++)
        {
            std::string allele = a_pVariantCalled->m_bIsFirstNucleotideTrimmed ? (a_pVariantCalled->GetRefSeq()[0] + a_pVariantCalled->m_alleles[k].m_sequence) : a_pVariantCalled->m_alleles[k].m_sequence;
            if(0 == baseVariants[p].compare(allele))
            {
                data2.m_aGenotype[k] = p;
                break;
            }
        }
    }

    
    data2.m_decisionBD = a_rDecisionCalled;
    if(a_rMatchTypeCalled != "nm")
        data2.m_matchTypeBK = a_rMatchTypeCalled;
    a_rOutputRec.m_aSampleData.push_back(data2);
    
}


bool CGa4ghOutputProvider::CanMerge(const CVariant* a_pVariantBase, const CVariant* a_pVariantCalled) const
{
    int baseStart = a_pVariantBase->m_nStartPos - (a_pVariantBase->m_bIsFirstNucleotideTrimmed ? 1 : 0);
    int calledStart = a_pVariantCalled->m_nStartPos - (a_pVariantCalled->m_bIsFirstNucleotideTrimmed ? 1 : 0);

    bool bIsPosEqual = baseStart == calledStart;
    bool bIsRefEqual = a_pVariantBase->m_refSequence == a_pVariantCalled->m_refSequence;
    
    if(bIsPosEqual && bIsRefEqual)
    {
        
        std::stringstream ss(a_pVariantBase->m_allelesStr);
        std::vector<std::string> baseVariants;
    
        std::string substr;
        while(getline(ss, substr, ','))
        {
            baseVariants.push_back(substr);
        }
    
        bool bIsAlleleExists;
    
        for(int k=0; k < a_pVariantCalled->m_nAlleleCount; k++)
        {
            bIsAlleleExists = false;
            for(int p = 0; p < baseVariants.size(); p++)
            {
                std::string allele = a_pVariantCalled->m_bIsFirstNucleotideTrimmed ? (a_pVariantCalled->GetRefSeq()[0] + a_pVariantCalled->m_alleles[k].m_sequence) : a_pVariantCalled->m_alleles[k].m_sequence;
                if(0 == baseVariants[p].compare(allele))
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

std::string CGa4ghOutputProvider::GetMatchStr(EVariantMatch a_match)
{
    const std::string GENOTYPE_MATCH = "gm";
    const std::string ALLELE_MATCH = "am";
    const std::string NO_MATCH = "nm";
    
    switch (a_match) {
        case eALLELE_MATCH:
            return ALLELE_MATCH;
        case eGENOTYPE_MATCH:
            return GENOTYPE_MATCH;
        default:
            return NO_MATCH;
    }
}
































