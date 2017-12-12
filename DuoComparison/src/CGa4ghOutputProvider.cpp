//
//  CGa4ghOutputProvider.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/20/16.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#include "CGa4ghOutputProvider.h"
#include <sstream>
#include "CPath.h"
#include "CVariantIteratorGa4gh.h"
#include "CVariantProvider.h"
#include "Constants.h"
#include <iostream>
#include <algorithm>
#include <Utils/CUtils.h>

using namespace duocomparison;

void CGa4ghOutputProvider::SetVariantProvider(CVariantProvider* a_pProvider)
{
    m_pVariantProvider = a_pProvider;
}

void CGa4ghOutputProvider::SetBestPaths(std::vector<core::CPath>& a_rBestPathList, std::vector<core::CPath>& a_rBestAlleleMatchPathList)
{
    m_aBestPaths = a_rBestPathList;
    m_aBestAlleleMatchPaths = a_rBestAlleleMatchPathList;
}

void CGa4ghOutputProvider::SetVcfPath(const std::string& a_rVcfPath)
{
    m_vcfPath = a_rVcfPath + "/Ga4ghOutput.vcf";
}

void CGa4ghOutputProvider::SetContigList(const std::vector<SVcfContig>& a_rContigs)
{
    m_contigs = a_rContigs;
}

void CGa4ghOutputProvider::GenerateGa4ghVcf(const std::vector<SChrIdTuple>& a_rCommonChromosomes)
{
    
    m_vcfWriter.CreateVcf(m_vcfPath.c_str());
    FillHeader();
    
    for(SChrIdTuple tuple : a_rCommonChromosomes)
    {
        std::cout << "Processing Chromosome " << tuple.m_chrName << std::endl;
        AddRecords(m_aBestPaths[tuple.m_nTupleIndex], tuple);
    }
    
    m_vcfWriter.CloseVcf();
}

void CGa4ghOutputProvider::FillHeader()
{
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= VBT Variant Comparison Tool " + VBT_VERSION);
    
    //ADD REQUIRED FORMATS BY GA4GH
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">");
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">");
    
    //ADD FILTER COLUMNS FROM CALL FILE
    std::vector<std::string> filterNames;
    std::vector<std::string> filterDescriptions;
    m_pVariantProvider->GetFilterInfo(eCALLED, filterNames, filterDescriptions);
    for(unsigned int k = 1; k < filterNames.size(); k++)
        m_vcfWriter.AddHeaderLine("##FILTER=<ID=" + filterNames[k] + ",Description=" + filterDescriptions[k] + ">");
    
    //ADD CONTIG IDs
    for(unsigned int k = 0; k < m_contigs.size(); k++)
        m_vcfWriter.AddHeaderLine("##contig=<ID=" + m_contigs[k].name + ",length=" + std::to_string(m_contigs[k].length) + ">");
    
    //ADD REQUIRED SAMPLES
    m_vcfWriter.AddSampleName("TRUTH");
    m_vcfWriter.AddSampleName("QUERY");
    
    //CLOSE HEADER
    m_vcfWriter.WriteHeaderToVcf();
}

void CGa4ghOutputProvider::AddRecords(const core::CPath& a_rBestPath, SChrIdTuple a_rTuple)
{
    //Best Path excluded variants
    std::vector<const CVariant*> excludedVarsBase = m_pVariantProvider->GetVariantList(eBASE, a_rTuple.m_nBaseId, a_rBestPath.m_baseSemiPath.GetExcluded());
    std::vector<const CVariant*> excludedVarsCall = m_pVariantProvider->GetVariantList(eCALLED, a_rTuple.m_nCalledId, a_rBestPath.m_calledSemiPath.GetExcluded());
    
    //Sort excluded variants by original positions
    std::sort(excludedVarsBase.begin(), excludedVarsBase.end(), CUtils::CompareVariantsById);
    std::sort(excludedVarsCall.begin(), excludedVarsCall.end(), CUtils::CompareVariantsById);
    
    //Best Path included variants
    std::vector<const core::COrientedVariant*> includedVarsBase = a_rBestPath.m_baseSemiPath.GetIncludedVariants();
    std::vector<const core::COrientedVariant*> includedVarsCall = a_rBestPath.m_calledSemiPath.GetIncludedVariants();
    
    //Sort included variants by original positions
    std::sort(includedVarsBase.begin(), includedVarsBase.end(), CUtils::CompareOrientedVariantsById);
    std::sort(includedVarsCall.begin(), includedVarsCall.end(), CUtils::CompareOrientedVariantsById);
    
    //Not Assessed variants
    std::vector<const CVariant*> notAssessedBase = m_pVariantProvider->GetNotAssessedVariantList(eBASE, a_rTuple.m_nBaseId);
    std::vector<const CVariant*> notAssessedCalled = m_pVariantProvider->GetNotAssessedVariantList(eCALLED, a_rTuple.m_nCalledId);

    //Get Skipped variants in Complex Regions (Will be marked as Non Assessed Variant)
    std::vector<const CVariant*> skippedComplexBase = m_pVariantProvider->GetSkippedComplexVariantList(eBASE, a_rTuple.m_nBaseId);
    std::vector<const CVariant*> skippedComplexCalled = m_pVariantProvider->GetSkippedComplexVariantList(eCALLED, a_rTuple.m_nCalledId);
    
    //Merge skipped complex and non assessed variants
    notAssessedBase.insert(std::end(notAssessedBase), std::begin(skippedComplexBase), std::end(skippedComplexBase));
    notAssessedCalled.insert(std::end(notAssessedCalled), std::begin(skippedComplexCalled), std::end(skippedComplexCalled));
    
    //Sort non assessed variants by original positions
    std::sort(notAssessedBase.begin(), notAssessedBase.end(), CUtils::CompareVariantsById);
    std::sort(notAssessedCalled.begin(), notAssessedCalled.end(), CUtils::CompareVariantsById);
    
    //Variant Iterators
    CVariantIteratorGa4gh baseVariants(includedVarsBase, excludedVarsBase, notAssessedBase);
    CVariantIteratorGa4gh calledVariants(includedVarsCall, excludedVarsCall, notAssessedCalled);
    
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
            for(unsigned int i = 0; i < nextVarBaseList.size(); i++)
            {
                for(unsigned j = 0; j < nextVarCalledList.size(); j++)
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
                record.m_aSampleData.push_back(SPerSampleData());
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
                record.m_aSampleData.push_back(SPerSampleData());
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
    data.m_bIsNoCallVariant = a_pVariant->m_bIsNoCall;
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
    data.m_bIsNoCallVariant = a_pVariantBase->m_bIsNoCall;
    
    std::stringstream ss(a_pVariantCalled->m_allelesStr);
    std::stringstream ss2(a_pVariantBase->m_allelesStr);
    std::vector<std::string> calledVariants;
    std::vector<std::string> baseVariants;
    std::string substr;
    
    while(getline(ss, substr, ','))
    {
        calledVariants.push_back(substr);
    }

    while(getline(ss2, substr, ','))
    {
        baseVariants.push_back(substr);
    }

    
    for(int k=0; k < (int)a_pVariantBase->m_nZygotCount; k++)
    {
        std::string allele = baseVariants[(unsigned int)a_pVariantBase->m_genotype[k]];
        
        int bHasFound = false;
        for(unsigned int p = 0; p < calledVariants.size(); p++)
        {
            if(0 == calledVariants[p].compare(allele))
            {
                data.m_aGenotype[k] = (int)p;
                bHasFound = true;
                break;
            }
        }
        
        if(!bHasFound)
        {
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
    data2.m_bIsNoCallVariant = a_pVariantCalled->m_bIsNoCall;
    
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











