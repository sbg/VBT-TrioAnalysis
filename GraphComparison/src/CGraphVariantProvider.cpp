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
        FillOrientedVariantLists();
        
        //Set the common chromosome list for parsing
        SetChromosomeIdTuples();
        
        for(duocomparison::SChrIdTuple tuple : m_aCommonChrTupleList)
        {
            if(m_aBaseVariantList[tuple.m_nBaseId].size() >= LEAST_VARIANT_THRESHOLD && m_aCalledVariantList[tuple.m_nCalledId].size() < LEAST_VARIANT_THRESHOLD)
            {
                std::cout << "Called VCF does not contain Chromosome " << tuple.m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aCalledVariantList[tuple.m_nCalledId].clear();
                m_aCalledOrientedVariantList[tuple.m_nCalledId].clear();
            }
            else if(m_aBaseVariantList[tuple.m_nBaseId].size() < LEAST_VARIANT_THRESHOLD && m_aCalledVariantList[tuple.m_nCalledId].size() >= LEAST_VARIANT_THRESHOLD)
            {
                std::cout << "Baseline VCF does not contain Chromosome " << tuple.m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aBaseVariantList[tuple.m_nBaseId].clear();
                m_aBaseOrientedVariantList[tuple.m_nBaseId].clear();
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
    m_aBaseVariantList = std::vector<std::vector<CVariant>>(m_baseVCF.GetContigSize());
    m_aCalledVariantList = std::vector<std::vector<CVariant>>(m_calledVCF.GetContigSize());
    
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
        else
        {
            m_aCalledVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    for(int k = 0; k < (int)m_calledVCF.GetContigSize(); k++)
        std::sort(m_aCalledVariantList[k].begin(), m_aCalledVariantList[k].end(), CompareVariants);
}

void CGraphVariantProvider::FillOrientedVariantLists()
{
    //Initialize OrientedVariantLists
    m_aBaseOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_baseVCF.GetContigSize());
    m_aCalledOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_calledVCF.GetContigSize());
    
    
    for(int i=0; i < (int)m_aBaseOrientedVariantList.size(); i++)
    {
        for(int j=0; j < (int)m_aBaseVariantList[i].size(); j++)
        {
            m_aBaseOrientedVariantList[i].push_back(core::COrientedVariant(m_aBaseVariantList[i][j], 0));
            m_aBaseOrientedVariantList[i].push_back(core::COrientedVariant(m_aBaseVariantList[i][j], 1));
        }
    }
    
    for(int i=0; i < (int)m_aCalledOrientedVariantList.size(); i++)
    {
        for(int j=0; j < (int)m_aCalledVariantList[i].size(); j++)
        {
            m_aCalledOrientedVariantList[i].push_back(core::COrientedVariant(m_aCalledVariantList[i][j], 0));
            m_aCalledOrientedVariantList[i].push_back(core::COrientedVariant(m_aCalledVariantList[i][j], 1));
        }
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

std::vector<const core::COrientedVariant*> CGraphVariantProvider::GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo)
{
    unsigned long size = a_uFrom == eBASE ? m_aBaseOrientedVariantList[a_nChrNo].size() : m_aCalledOrientedVariantList[a_nChrNo].size();
    std::vector<const core::COrientedVariant*> result(size);
    
    switch (a_uFrom)
    {
        case eBASE:
            for(int k=0; k < (int)size; k++)
                result[k] = &(m_aBaseOrientedVariantList[a_nChrNo][k]);
            break;
            
        case eCALLED:
            for(int k=0; k< (int)size; k++)
                result[k] = &(m_aCalledOrientedVariantList[a_nChrNo][k]);
            
        default:
            break;
    }
    
    return result;
}

std::vector<const core::COrientedVariant*> CGraphVariantProvider::GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int> a_rVariantIndexes)
{
    std::vector<const core::COrientedVariant*> ovarList(a_rVariantIndexes.size() * 2);
    int itr = 0;
    
    switch (a_uFrom)
    {
        case eBASE:
            for(int ind : a_rVariantIndexes)
            {
                ovarList[itr++] = &(m_aBaseOrientedVariantList[a_nChrNo][2*ind]);
                ovarList[itr++] = &(m_aBaseOrientedVariantList[a_nChrNo][2*ind +1]);
            }
            break;
            
        case eCALLED:
            for(int ind : a_rVariantIndexes)
            {
                ovarList[itr++] = &(m_aCalledOrientedVariantList[a_nChrNo][2*ind]);
                ovarList[itr++] = &(m_aCalledOrientedVariantList[a_nChrNo][2*ind +1]);
            }
            break;
            
        default:
            break;
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




