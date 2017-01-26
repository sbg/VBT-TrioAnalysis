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
#include "CVariantIteratorV2.h"
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
    m_vcfPath = a_rVcfPath + "/Ga4ghOutput.vcf";
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
    
    //Not Asessed variants
    std::vector<CVariant>& notAssessedBase = m_pVariantProvider->GetNotAssessedVariantList(eBASE, a_nChrId);
    std::vector<CVariant>& notAssessedCalled = m_pVariantProvider->GetNotAssessedVariantList(eCALLED, a_nChrId);

    //Variant Iterators
    CVariantIteratorV2 baseVariants(includedVarsBase, excludedVarsBase, notAssessedBase);
    CVariantIteratorV2 calledVariants(includedVarsCall, excludedVarsCall, notAssessedCalled);
    
    std::vector<SVariantSummary> nextVarBaseList;
    std::vector<SVariantSummary> nextVarCalledList;

    baseVariants.FillNext(nextVarBaseList);
    calledVariants.FillNext(nextVarCalledList);
    
    while(true)
    {
        if(nextVarBaseList.size() == 0 && nextVarCalledList.size() == 0)
            break;
        
        int basePosition = nextVarBaseList.size() > 0 ? nextVarBaseList[0].originalStartPos() : INT_MAX;
        int calledPosition = nextVarCalledList.size() > 0 ? nextVarCalledList[0].originalStartPos() : INT_MAX;

        if(basePosition == calledPosition)
        {
            //Look for a match between base and called list
            for(int i = 0; i < nextVarBaseList.size(); i++)
            {
                for(int j = 0; j < nextVarCalledList.size(); j++)
                {
                    if(CanMerge(nextVarBaseList[i].m_pVariant, nextVarCalledList[j].m_pVariant))
                    {
                        SVcfRecord record;
                        std::string decisionBase = nextVarBaseList[i].m_bIncluded ? "TP" : (nextVarBaseList[i].m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FN");
                        std::string decisionCalled = nextVarCalledList[j].m_bIncluded ? "TP" : (nextVarCalledList[j].m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FP");
                        std::string matchBase = GetMatchStr(nextVarBaseList[i].m_pVariant->m_variantStatus);
                        std::string matchCalled = GetMatchStr(nextVarCalledList[j].m_pVariant->m_variantStatus);
                        MergeVariants(nextVarBaseList[i].m_pVariant, nextVarCalledList[j].m_pVariant, matchBase, matchCalled, decisionBase, decisionCalled, record);
                        m_vcfWriter.AddRecord(record);
                        
                        nextVarBaseList.erase(nextVarBaseList.begin() + i);
                        nextVarCalledList.erase(nextVarCalledList.begin() + j);
                        j--;
                        i--;
                        break;
                    }
                }
            }

            //If there are base variants exists which doesnt merge already
            for (SVariantSummary var : nextVarCalledList)
            {
                SVcfRecord record;
                record.m_aSampleData.push_back(SPerSampleData());
                std::string decision = var.m_bIncluded ? "TP" : (var.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FP");
                std::string match = GetMatchStr(var.m_pVariant->m_variantStatus);
                VariantToVcfRecord(var.m_pVariant, record, false, match, decision);
                m_vcfWriter.AddRecord(record);
            }
            
            nextVarCalledList.clear();
            calledVariants.FillNext(nextVarCalledList);


            //If there are called variants exists which doesnt merge already
            for(SVariantSummary var : nextVarBaseList)
            {
                SVcfRecord record;
                std::string decision = var.m_bIncluded ? "TP" : (var.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FN");
                std::string match = GetMatchStr(var.m_pVariant->m_variantStatus);
                VariantToVcfRecord(var.m_pVariant, record, true, match, decision);
                m_vcfWriter.AddRecord(record);
            }
            
            nextVarBaseList.clear();
            baseVariants.FillNext(nextVarBaseList);
        }
        
        else if (basePosition > calledPosition)
        {
            for (SVariantSummary var : nextVarCalledList)
            {
                SVcfRecord record;
                record.m_aSampleData.push_back(SPerSampleData());
                std::string decision = var.m_bIncluded ? "TP" : (var.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FP");
                std::string match = GetMatchStr(var.m_pVariant->m_variantStatus);
                VariantToVcfRecord(var.m_pVariant, record, false, match, decision);
                m_vcfWriter.AddRecord(record);
            }
            
            nextVarCalledList.clear();
            calledVariants.FillNext(nextVarCalledList);

        }
        
        else
        {
            for(SVariantSummary var : nextVarBaseList)
            {
                SVcfRecord record;
                std::string decision = var.m_bIncluded ? "TP" : (var.m_pVariant->m_variantStatus == eNOT_ASSESSED ? "N" : "FN");
                std::string match = GetMatchStr(var.m_pVariant->m_variantStatus);
                VariantToVcfRecord(var.m_pVariant, record, true, match, decision);
                m_vcfWriter.AddRecord(record);
            }
            
            nextVarBaseList.clear();
            baseVariants.FillNext(nextVarBaseList);

        }
    }
    
}



void CGa4ghOutputProvider::VariantToVcfRecord(const CVariant* a_pVariant, SVcfRecord& a_rOutputRec, bool a_bIsBase, const std::string& a_rMatchType, const::std::string& a_rDecision)
{
    //Fill basic variant data
    a_rOutputRec.m_chrName = a_pVariant->m_chrName;
    a_rOutputRec.m_nPosition = a_pVariant->m_nOriginalPos;
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
    a_rOutputRec.m_chrName = a_pVariantCalled->m_chrName;
    a_rOutputRec.m_nPosition = a_pVariantCalled->m_nOriginalPos;
    a_rOutputRec.m_alleles = a_pVariantCalled->m_allelesStr;
    a_rOutputRec.m_aFilterString = a_pVariantCalled->m_filterString;
    
    //Fill base sample (TRUTH)
    SPerSampleData data;
    data.m_bIsPhased = a_pVariantBase->m_bIsPhased;
    data.m_nHaplotypeCount = a_pVariantBase->m_nZygotCount;
    
    std::stringstream ss(a_pVariantCalled->m_allelesStr);
    std::vector<std::string> calledVariants;
    std::string substr;
    
    while(getline(ss, substr, ','))
    {
        calledVariants.push_back(substr);
    }
    
    for(int k=0; k < a_pVariantBase->m_nZygotCount; k++)
    {
        int bHasFound = false;
        for(int p = 0; p < calledVariants.size(); p++)
        {
            std::string allele = a_pVariantBase->m_bIsFirstNucleotideTrimmed ? (a_pVariantBase->GetRefSeq()[0] + a_pVariantBase->m_alleles[k].m_sequence) : a_pVariantBase->m_alleles[k].m_sequence;
            if(0 == calledVariants[p].compare(allele))
            {
                data.m_aGenotype[k] = p;
                bHasFound = true;
                break;
            }
        }
        
        if(!bHasFound)
        {
            std::string allele = a_pVariantBase->m_bIsFirstNucleotideTrimmed ? (a_pVariantBase->GetRefSeq()[0] + a_pVariantBase->m_alleles[k].m_sequence) : a_pVariantBase->m_alleles[k].m_sequence;
            calledVariants.push_back(allele);
            a_rOutputRec.m_alleles += ("," + allele);
            data.m_aGenotype[k] = static_cast<int>(calledVariants.size()) - 1;
        }
        
        
    }
    
    data.m_decisionBD = a_rDecisionBase;
    if(a_rMatchTypeBase != "nm")
        data.m_matchTypeBK = a_rMatchTypeBase;
    
    a_rOutputRec.m_aSampleData.push_back(data);
    
    //Fill called sample (QUERY)
    SPerSampleData data2;
    data2.m_bIsPhased = a_pVariantCalled->m_bIsPhased;
    data2.m_nHaplotypeCount = a_pVariantCalled->m_nZygotCount;
    
    for(int k = 0; k < data2.m_nHaplotypeCount; k++)
        data2.m_aGenotype[k] = a_pVariantCalled->m_genotype[k];
    
    data2.m_decisionBD = a_rDecisionCalled;
    if(a_rMatchTypeCalled != "nm")
        data2.m_matchTypeBK = a_rMatchTypeCalled;
    a_rOutputRec.m_aSampleData.push_back(data2);
    
}


bool CGa4ghOutputProvider::CanMerge(const CVariant* a_pVariantBase, const CVariant* a_pVariantCalled) const
{
    bool bIsPosEqual = a_pVariantBase->m_nOriginalPos == a_pVariantCalled->m_nOriginalPos;
    bool bIsRefEqual = a_pVariantBase->m_refSequence == a_pVariantCalled->m_refSequence;
    
    if(bIsPosEqual && bIsRefEqual)
        return true;
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











