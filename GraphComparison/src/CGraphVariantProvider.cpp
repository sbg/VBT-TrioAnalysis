//
//  CGraphVariantProvider.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/11/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CGraphVariantProvider.h"
#include "CSimpleBEDParser.h"
#include "CVcfWriter.h"
#include <iostream>
#include <algorithm>

using namespace graphcomparison;

void CGraphVariantProvider::SetParameters(const std::string& a_rBaseVcf,
                                          const std::string& a_rCalledVcf,
                                          const std::string& a_rReference,
                                          const std::string& a_rBedFilePath,
                                          bool a_bIsPassFilterEnabled,
                                          int a_nMaxBasePairLength)
{
    m_baseVcfPath = a_rBaseVcf;
    m_calledVcfPath = a_rCalledVcf;
    m_fastaPath = a_rReference;
    m_bPassFilterEnabled = a_bIsPassFilterEnabled;
    m_nMaxBasePairLength = a_nMaxBasePairLength;
    m_bedFilePath = a_rBedFilePath;
    m_bIsBEDFileEnabled = (a_rBedFilePath != "");
    
}

void CGraphVariantProvider::SetChromosomeIdTuples()
{
    int tupleIndex = 0;
    
    for(auto baseItr = m_baseVCF.m_chrIndexMap.begin(); baseItr != m_baseVCF.m_chrIndexMap.end(); baseItr++)
    {
        for(auto calledItr = m_calledVCF.m_chrIndexMap.begin(); calledItr != m_calledVCF.m_chrIndexMap.end(); calledItr++)
        {
            if(baseItr->first == calledItr->first)
            {
                if(m_aBaseVariantList[baseItr->second].size() > LEAST_VARIANT_THRESHOLD
                   &&
                   m_aCalledVariantList[calledItr->second].size() > LEAST_VARIANT_THRESHOLD)
                {
                    m_aCommonChrTupleList.push_back(duocomparison::SChrIdTuple(baseItr->second, calledItr->second, baseItr->first, tupleIndex++));
                }
            }
        }
    }
    
    std::sort(m_aCommonChrTupleList.begin(), m_aCommonChrTupleList.end(), [](const duocomparison::SChrIdTuple& t1, const duocomparison::SChrIdTuple& t2){return t1.m_nBaseId < t2.m_nBaseId;});
    
}

std::vector<duocomparison::SChrIdTuple>& CGraphVariantProvider::GetChromosomeIdTuples()
{
    return m_aCommonChrTupleList;
}

bool CGraphVariantProvider::InitializeReaders()
{
    bool bIsSuccess;
    
    //OPEN VCF FILES
    bIsSuccess = m_baseVCF.Open(m_baseVcfPath.c_str());
    if(!bIsSuccess)
    {
        std::cerr << "Baseline VCF file is unable to open!: " << m_baseVcfPath << std::endl;
        return false;
    }
    
    bIsSuccess = m_calledVCF.Open(m_calledVcfPath.c_str());
    if(!bIsSuccess)
    {
        std::cerr << "Called VCF file is unable to open!: " << m_calledVcfPath << std::endl;
        return false;
    }
    
    // OPEN FASTA FILE
    bIsSuccess = m_fastaParser.OpenFastaFile(m_fastaPath.c_str());
    if(!bIsSuccess)
    {
        std::cerr << "FASTA file is unable to open!: " << m_fastaPath << std::endl;
        return false;
    }
    
    if(bIsSuccess)
    {
        FillVariantLists();
        //FillOrientedVariantLists();
        
        //Set the common chromosome list for parsing
        SetChromosomeIdTuples();
        
        for(duocomparison::SChrIdTuple tuple : m_aCommonChrTupleList)
        {
            if(m_aBaseVariantList[tuple.m_nBaseId].size() >= LEAST_VARIANT_THRESHOLD && m_aCalledVariantList[tuple.m_nCalledId].size() < LEAST_VARIANT_THRESHOLD)
            {
                std::cout << "Called VCF does not contain Chromosome " << tuple.m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aCalledVariantList[tuple.m_nCalledId].clear();
            }
            else if(m_aBaseVariantList[tuple.m_nBaseId].size() < LEAST_VARIANT_THRESHOLD && m_aCalledVariantList[tuple.m_nCalledId].size() >= LEAST_VARIANT_THRESHOLD)
            {
                std::cout << "Baseline VCF does not contain Chromosome " << tuple.m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aBaseVariantList[tuple.m_nBaseId].clear();
            }
        }
    }
    return bIsSuccess;
    
}

bool CompareVariants(const CVariant& var1, const CVariant& var2)
{
    if(var1.m_nStartPos != var2.m_nStartPos)
        return var1.m_nStartPos < var2.m_nStartPos;
    else if(var1.m_nOriginalPos != var2.m_nOriginalPos)
        return var1.m_nOriginalPos < var2.m_nOriginalPos;
    else
        return var1.m_nId < var2.m_nId;
}

void CGraphVariantProvider::FillVariantLists()
{
    //BED region parameters if there is a defined bed file
    CSimpleBEDParser bedParser;
    SBedRegion bedRegion;
    bool hasNextRegion = true;
    
    if(m_bIsBEDFileEnabled)
    {
        bedParser.InitBEDFile(m_bedFilePath);
        hasNextRegion = bedParser.GetNextRegion(bedRegion);
    }
    
    if(m_bIsBEDFileEnabled && hasNextRegion == false)
        return;
    
    CVariant variant;
    int id = 0;
    std::string preChrId = "";
    
    //Initialize variantLists
    m_aBaseVariantList = std::deque<std::deque<CVariant>>(m_baseVCF.GetContigSize());
    m_aCalledVariantList = std::deque<std::deque<CVariant>>(m_calledVCF.GetContigSize());
    
    while(m_baseVCF.GetNextRecord(&variant, id))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cerr << "Reading chromosome " << preChrId << " of base vcf..." << std::endl;
            variant.m_nId = 0;
            id = 0;
        }
        
        //Pass to the next region
        if(m_bIsBEDFileEnabled &&
           ((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
           ||
           m_baseVCF.m_chrIndexMap[variant.m_chrName] > m_baseVCF.m_chrIndexMap[bedRegion.m_chrName]))
        {
            hasNextRegion = bedParser.GetNextRegion(bedRegion);
            if(false == hasNextRegion)
                break;
        }

        //Variant Could not pass from BED region
        if(m_bIsBEDFileEnabled &&
           (bedRegion.m_chrName != variant.m_chrName || std::min(bedRegion.m_nEndPos, variant.m_nOriginalPos + static_cast<int>(variant.m_refSequence.length())) - std::max(bedRegion.m_nStartPos, variant.m_nOriginalPos) < 0))
            continue;
        //Variant Count not pass from PASS filtering
        else if(m_bPassFilterEnabled && variant.m_bIsFilterPASS == false)
            continue;
        //Variant Length is higher than the max base-pair limit
        else if(variant.m_nEndPos - variant.m_nStartPos > m_nMaxBasePairLength)
            continue;
        else if(variant.m_alleles[0].m_sequence == "*")
            continue;
        //Add variant to the variant list
        else
        {
            m_aBaseVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    for(int k = 0; k < (int)m_baseVCF.GetContigSize(); k++)
        std::sort(m_aBaseVariantList[k].begin(), m_aBaseVariantList[k].end(), CompareVariants);
    
    if(m_bIsBEDFileEnabled)
    {
        bedParser.ResetIterator();
        bedParser.GetNextRegion(bedRegion);
    }
    
    preChrId = "";
    id = 0;
    
    while(m_calledVCF.GetNextRecord(&variant, id))
    {
        
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cerr << "Reading chromosome " << preChrId << " of called vcf..." << std::endl;
            variant.m_nId = 0;
            id = 0;
        }
        
        //Pass to the next region
        if(m_bIsBEDFileEnabled &&
           ((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
           ||
           m_calledVCF.m_chrIndexMap[variant.m_chrName] > m_calledVCF.m_chrIndexMap[bedRegion.m_chrName]))
        {
            hasNextRegion = bedParser.GetNextRegion(bedRegion);
            if(false == hasNextRegion)
                break;
        }
        
        //Variant Could not pass from BED region
        if(m_bIsBEDFileEnabled &&
           (bedRegion.m_chrName != variant.m_chrName || std::min(bedRegion.m_nEndPos, variant.m_nOriginalPos + static_cast<int>(variant.m_refSequence.length())) - std::max(bedRegion.m_nStartPos, variant.m_nOriginalPos) < 0))
            continue;
        //Variant Count not pass from PASS filtering
        else if(m_bPassFilterEnabled && variant.m_bIsFilterPASS == false)
            continue;
        //Variant Length is higher than the max base-pair limit
        else if(variant.m_nEndPos - variant.m_nStartPos > m_nMaxBasePairLength)
            continue;
        else if(variant.m_alleles[0].m_sequence == "*")
            continue;
        else
        {
            m_aCalledVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    for(int k = 0; k < (int)m_calledVCF.GetContigSize(); k++)
        std::sort(m_aCalledVariantList[k].begin(), m_aCalledVariantList[k].end(), CompareVariants);
}

void CGraphVariantProvider::FillOrientedVariantList(const duocomparison::SChrIdTuple& a_rTuple, std::deque<core::COrientedVariant>& a_rBaseOrientedVars, std::deque<core::COrientedVariant>& a_rCalledOrientedVars)
{
    
    for(int j=0; j < (int)m_aBaseVariantList[a_rTuple.m_nBaseId].size(); j++)
    {
        a_rBaseOrientedVars.push_back(core::COrientedVariant(m_aBaseVariantList[a_rTuple.m_nBaseId][j], 0));
        a_rBaseOrientedVars.push_back(core::COrientedVariant(m_aBaseVariantList[a_rTuple.m_nBaseId][j], 1));
    }

    for(int j=0; j < (int)m_aCalledVariantList[a_rTuple.m_nCalledId].size(); j++)
    {
        a_rCalledOrientedVars.push_back(core::COrientedVariant(m_aCalledVariantList[a_rTuple.m_nCalledId][j], 0));
        a_rCalledOrientedVars.push_back(core::COrientedVariant(m_aCalledVariantList[a_rTuple.m_nCalledId][j], 1));
    }
    
}

std::vector<const CVariant*> CGraphVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes)
{
    std::vector<const CVariant*> result;
    
    switch (a_uFrom)
    {
        case eBASE:
            for(int k=0; k< (int)a_VariantIndexes.size(); k++)
                result.push_back(&(m_aBaseVariantList[a_nChrNo][a_VariantIndexes[k]]));
            break;
            
        case eCALLED:
            for(int k=0; k< (int)a_VariantIndexes.size(); k++)
                result.push_back(&(m_aCalledVariantList[a_nChrNo][a_VariantIndexes[k]]));
            
        default:
            break;
    }
    
    return result;
}

std::vector<const CVariant*> CGraphVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo)
{
    unsigned long size = a_uFrom == eBASE ? m_aBaseVariantList[a_nChrNo].size() : m_aCalledVariantList[a_nChrNo].size();
    std::vector<const CVariant*> result(size);
    
    switch (a_uFrom)
    {
        case eBASE:
            for(int k=0; k < (int)size; k++)
                result[k] = &(m_aBaseVariantList[a_nChrNo][k]);
            break;
            
        case eCALLED:
            for(int k=0; k< (int)size; k++)
                result[k] = &(m_aCalledVariantList[a_nChrNo][k]);
            
        default:
            break;
    }
    
    return result;
}

std::vector<const core::COrientedVariant*> CGraphVariantProvider::GetOrientedVariantList(std::deque<core::COrientedVariant>& a_rOvarList)
{
    unsigned long size = a_rOvarList.size();
    std::vector<const core::COrientedVariant*> result(size);
    
    for(int k=0; k< (int)size; k++)
        result[k] = &(a_rOvarList[k]);
            
    return result;
}

std::vector<const core::COrientedVariant*> CGraphVariantProvider::GetOrientedVariantList(std::deque<core::COrientedVariant>& a_rOvarList, const std::vector<const CVariant*>& a_rVariants)
{
    std::vector<const core::COrientedVariant*> ovarList(a_rVariants.size() * 2);
    int itr = 0;
    
    for(unsigned int k = 0; k < a_rVariants.size(); k++)
    {
        ovarList[itr++] = &(a_rOvarList[2*(a_rVariants[k]->m_nId)]);
        ovarList[itr++] = &(a_rOvarList[2*(a_rVariants[k]->m_nId) +1]);
    }
    
    return ovarList;
}

std::vector<int> CGraphVariantProvider::GetExcludedIndexes(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_rExcludedIndexList)
{
    std::vector<int> excludedVariantIndexes;

    for(int ind : a_rExcludedIndexList)
        excludedVariantIndexes.push_back(a_rVariantList[ind]->m_nId);

    return excludedVariantIndexes;
}

void CGraphVariantProvider::GetContig(std::string a_chrName, SContig& a_rCtg)
{
    a_rCtg.Clean();
    m_fastaParser.FetchNewChromosome(a_chrName, a_rCtg);
}

void CGraphVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
        pVar->m_variantStatus = a_status;
}

void CGraphVariantProvider::SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const core::COrientedVariant* pOVar : a_rVariantList)
        pOVar->GetVariant().m_variantStatus = a_status;
}

std::vector<int> CGraphVariantProvider::GetVariantIndexesByStatus(EVcfName a_uFrom, int a_nChrNo, EVariantMatch a_nStatus)
{
    std::vector<int> returnIndexes;
    
    switch (a_uFrom)
    {
        case eBASE:
            for(unsigned int k = 0; k < m_aBaseVariantList[a_nChrNo].size(); k++)
            {
                if(m_aBaseVariantList[a_nChrNo][k].m_variantStatus == a_nStatus)
                    returnIndexes.push_back(k);
            }
            break;
        case eCALLED:
            for(unsigned int k = 0; k < m_aCalledVariantList[a_nChrNo].size(); k++)
            {
                if(m_aCalledVariantList[a_nChrNo][k].m_variantStatus == a_nStatus)
                    returnIndexes.push_back(k);
            }
            break;
        default:
            break;
    }
    
    return returnIndexes;
}



