//
//  CVariantProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVariantProvider.h"
#include "COrientedVariant.h"
#include "CSimpleBEDParser.h"
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace duocomparison;

//Checks if the given two range is overlapping
bool isOverlap(int left1, int right1, int left2, int right2)
{
    //If the interval length is 0 (eg. 974791-974791) we need to check if the boundaries matches
    if(left1 == left2)
        return true;
    if(right1 == left1)
        return (right2 > left1 && left2 <= left1);
    else if(right2 == left2)
        return (right1 > left2 && left1 <= left2);
    else
        return std::min(right1, right2) - std::max(left1, left2) > 0;
}

CVariantProvider::CVariantProvider()
{
    m_bIsHomozygousOvarListInitialized = false;
}

CVariantProvider::~CVariantProvider()
{
    for(unsigned int k= 0; k < m_aContigList.size(); k++)
    {
        if(m_aContigList[k].m_nRefLength > 0 && m_aContigList[k].m_pRefSeq != 0)
            delete[] m_aContigList[k].m_pRefSeq;
    }
}

void CVariantProvider::SetChromosomeIdTuples()
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
                    m_aCommonChrTupleList.push_back(SChrIdTuple(baseItr->second, calledItr->second, baseItr->first, tupleIndex++));
                }
            }
        }
    }
    
    std::sort(m_aCommonChrTupleList.begin(), m_aCommonChrTupleList.end(), [](const SChrIdTuple& t1, const SChrIdTuple& t2){return t1.m_nBaseId < t2.m_nBaseId;});
}

bool CVariantProvider::InitializeReaders(const SConfig& a_rConfig)
{
    bool bIsSuccess;

    m_config = a_rConfig;
    
    //OPEN VCF FILES
    bIsSuccess = m_baseVCF.Open(a_rConfig.m_pBaseVcfFileName);
    if(!bIsSuccess)
    {
        std::cerr << "Baseline VCF file is unable to open!: " << a_rConfig.m_pBaseVcfFileName << std::endl;
        return false;
    }
    m_baseVCF.setID(0);

    bIsSuccess = m_calledVCF.Open(a_rConfig.m_pCalledVcfFileName);
    if(!bIsSuccess)
    {
        std::cerr << "Called VCF file is unable to open!: " << a_rConfig.m_pCalledVcfFileName << std::endl;
        return false;
    }
    m_calledVCF.setID(1);
    
    //SET SAMPLE NAME TO READ ONLY ONE SAMPLE FROM THE VCF
    if (true == m_config.m_bBaseSampleEnabled)
       bIsSuccess = m_baseVCF.SelectSample(m_config.m_pBaseSample);
    else
    {
        std::vector<std::string> sampleNames;
        m_baseVCF.GetSampleNames(sampleNames);
       bIsSuccess = m_baseVCF.SelectSample(sampleNames[0]);
    }
    
    if(!bIsSuccess)
    {
        std::cerr << "Baseline Sample name is incorrect!" << std::endl;
        return false;
    }
    
    if (true == m_config.m_bCalledSampleEnabled)
        bIsSuccess = m_calledVCF.SelectSample(m_config.m_pCalledSample);
    else
    {
        std::vector<std::string> sampleNames;
        m_calledVCF.GetSampleNames(sampleNames);
        bIsSuccess = m_calledVCF.SelectSample(sampleNames[0]);
    }
    
    if(!bIsSuccess)
    {
        std::cerr << "Called Sample name is incorrect!" << std::endl;
        return false;
    }

    
    // OPEN FASTA FILE
    bIsSuccess = m_fastaParser.OpenFastaFile(a_rConfig.m_pFastaFileName);
    if(!bIsSuccess)
    {
        std::cerr << "FASTA file is unable to open!: " << a_rConfig.m_pFastaFileName << std::endl;
        return false;
    }
    
    if(bIsSuccess)
    {
        if(m_config.m_bInitializeFromBed)
            FillVariantsFromBED();
        else
            FillVariantLists();
        FillOrientedVariantLists();
        
        //Set the common chromosome list for parsing
        SetChromosomeIdTuples();
        
        //Initialize contig list
        m_aContigList = std::vector<SContig>(m_aCommonChrTupleList.size());
    
        for(SChrIdTuple tuple : m_aCommonChrTupleList)
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
            else if(m_aBaseVariantList[tuple.m_nBaseId].size() > LEAST_VARIANT_THRESHOLD && m_aCalledVariantList[tuple.m_nCalledId].size() > LEAST_VARIANT_THRESHOLD)
            {
                //Read contig from FASTA file for the given chromosome
                std::cout << "Reading reference of chromosome " << tuple.m_chrName << " from the FASTA file" << std::endl;
                bool bIsSuccess2 = m_fastaParser.FetchNewChromosome(tuple.m_chrName, m_aContigList[tuple.m_nTupleIndex]);
                if(!bIsSuccess2)
                {
                    std::cerr << "Chromosome " << m_aCommonChrTupleList[tuple.m_nTupleIndex].m_chrName << "will be filtered out from the comparison since reference FASTA could not read or it does not contain given chromosome" << std::endl;
                    m_aCalledVariantList[tuple.m_nCalledId].clear();
                    m_aBaseVariantList[tuple.m_nBaseId].clear();
                    m_aBaseOrientedVariantList[tuple.m_nBaseId].clear();
                    m_aCalledOrientedVariantList[tuple.m_nCalledId].clear();
                }
                else
                {
                    //Trim the variants which are out of bound according to FASTA file
                    while(m_aBaseVariantList[tuple.m_nBaseId].back().GetEnd() > m_aContigList[tuple.m_nTupleIndex].m_nRefLength)
                        m_aBaseVariantList[tuple.m_nBaseId].pop_back();
                    while(m_aCalledVariantList[tuple.m_nCalledId].back().GetEnd() > m_aContigList[tuple.m_nTupleIndex].m_nRefLength)
                        m_aCalledVariantList[tuple.m_nCalledId].pop_back();
                }
            }
        }
        
    }
    
    return bIsSuccess;
}

bool CVariantProvider::CompareVariants(const CVariant& var1, const CVariant& var2)
{
    if(var1.m_nStartPos != var2.m_nStartPos)
        return var1.m_nStartPos < var2.m_nStartPos;
    else if(var1.m_nEndPos != var2.m_nEndPos)
        return var1.m_nEndPos < var2.m_nEndPos;
    else
        return var1.m_nId < var2.m_nId;
}

void CVariantProvider::FillVariantsFromBED()
{
    CVariant variant;
    int id = 0;
    std::string preChrId = "";

    //Initialize variantLists
    m_aBaseVariantList = std::vector<std::vector<CVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledVariantList = std::vector<std::vector<CVariant>>(m_calledVCF.GetContigs().size());
    m_aBaseNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_calledVCF.GetContigs().size());

    std::vector<CVariant> multiTrimmableVarListBase;
    std::vector<CVariant> multiTrimmableVarListCalled;
    
    CSimpleBEDParser bedParser;
    bedParser.InitBEDFile(m_config.m_pBedFileName);
    
    unsigned int regionIterator = 0;
    
    //std::string lastChrNameBed;
    
    while(m_baseVCF.GetNextRecord(&variant, id++, m_config))
    {
        
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            regionIterator = 0;
            std::cout << "Processing chromosome " << preChrId << " of base vcf" << std::endl;
        }
        
        //No Region exist for this chromosome
        if(bedParser.m_regionMap[variant.m_chrName].size() == 0)
            continue;
        
        //Skip to next region
        while(regionIterator < bedParser.m_regionMap[variant.m_chrName].size()
              &&
              variant.m_nOriginalPos >= bedParser.m_regionMap[variant.m_chrName][regionIterator].m_nEndPos)
        {
            regionIterator++;
        }
        
        //Skip if regions are finished for given chromosome
        if(regionIterator == bedParser.m_regionMap[variant.m_chrName].size())
            continue;
        
        //Variant Could not pass from BED region
        if(bedParser.m_regionMap[variant.m_chrName][regionIterator].m_nStartPos >= variant.m_nEndPos)
            continue;
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);

        else if(true == variant.m_bHaveMultipleTrimOption)
            multiTrimmableVarListBase.push_back(variant);
        
        else
            m_aBaseVariantList[variant.m_nChrId].push_back(variant);
        
    }
    
    FindOptimalTrimmings(multiTrimmableVarListBase, eBASE);
    AppendTrimmedVariants(multiTrimmableVarListBase, eBASE);

    
    for(unsigned int k = 0; k < m_baseVCF.GetContigs().size(); k++)
    {
        std::sort(m_aBaseNotAssessedVariantList[k].begin(), m_aBaseNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aBaseVariantList[k].begin(), m_aBaseVariantList[k].end(), CompareVariants);
    }
    
    preChrId = "";
    
    while(m_calledVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            regionIterator = 0;
            std::cout << "Processing chromosome " << preChrId << " of called vcf" << std::endl;
        }

        //No Region exist for this chromosome
        if(bedParser.m_regionMap[variant.m_chrName].size() == 0)
            continue;
        
        //Skip to next region
        while(regionIterator < bedParser.m_regionMap[variant.m_chrName].size()
              &&
              variant.m_nOriginalPos >= bedParser.m_regionMap[variant.m_chrName][regionIterator].m_nEndPos)
        {
            regionIterator++;
        }
        
        //Skip if regions are finished for given chromosome
        if(regionIterator == bedParser.m_regionMap[variant.m_chrName].size())
            continue;
        
        //Variant Could not pass from BED region
        if(bedParser.m_regionMap[variant.m_chrName][regionIterator].m_nStartPos >= variant.m_nEndPos)
            continue;
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        else if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(variant.m_bHaveMultipleTrimOption)
            multiTrimmableVarListCalled.push_back(variant);
        
        else
            m_aCalledVariantList[variant.m_nChrId].push_back(variant);

    }
    
    FindOptimalTrimmings(multiTrimmableVarListCalled, eCALLED);
    AppendTrimmedVariants(multiTrimmableVarListCalled, eCALLED);
    
    for(unsigned int k = 0; k < m_calledVCF.GetContigs().size(); k++)
    {
        std::sort(m_aCalledNotAssessedVariantList[k].begin(), m_aCalledNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aCalledVariantList[k].begin(), m_aCalledVariantList[k].end(), CompareVariants);
    }
        
}


void CVariantProvider::FillVariantLists()
{
    CVariant variant;
    int id = 0;
    std::string preChrId = "";
    
    //Initialize variantLists
    m_aBaseVariantList = std::vector<std::vector<CVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledVariantList = std::vector<std::vector<CVariant>>(m_calledVCF.GetContigs().size());
    m_aBaseNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_calledVCF.GetContigs().size());
    
    std::vector<CVariant> multiTrimmableVarListBase;
    std::vector<CVariant> multiTrimmableVarListCalled;

    while(m_baseVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of base vcf" << std::endl;
        }
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(true == variant.m_bHaveMultipleTrimOption)
            multiTrimmableVarListBase.push_back(variant);
        
        else
            m_aBaseVariantList[variant.m_nChrId].push_back(variant);
    }
    
    FindOptimalTrimmings(multiTrimmableVarListBase, eBASE);
    AppendTrimmedVariants(multiTrimmableVarListBase, eBASE);
    
    for(unsigned int k = 0; k < m_baseVCF.GetContigs().size(); k++)
    {
        std::sort(m_aBaseNotAssessedVariantList[k].begin(), m_aBaseNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aBaseVariantList[k].begin(), m_aBaseVariantList[k].end(), CompareVariants);
    }

    preChrId = "";
    
    while(m_calledVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of called vcf" << std::endl;
        }
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        else if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(true == variant.m_bHaveMultipleTrimOption)
            multiTrimmableVarListCalled.push_back(variant);
        
        else
            m_aCalledVariantList[variant.m_nChrId].push_back(variant);
    }
    
    FindOptimalTrimmings(multiTrimmableVarListCalled, eCALLED);
    AppendTrimmedVariants(multiTrimmableVarListCalled, eCALLED);
    
    for(unsigned int k = 0; k < m_calledVCF.GetContigs().size(); k++)
    {
        std::sort(m_aCalledNotAssessedVariantList[k].begin(), m_aCalledNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aCalledVariantList[k].begin(), m_aCalledVariantList[k].end(), CompareVariants);
    }
}

void CVariantProvider::FillOrientedVariantLists()
{
    //Initialize OrientedVariantLists
    m_aBaseOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_calledVCF.GetContigs().size());
    
    
    for(unsigned int i=0; i < m_aBaseOrientedVariantList.size(); i++)
    {
        for(unsigned int j=0; j < m_aBaseVariantList[i].size(); j++)
        {
            m_aBaseOrientedVariantList[i].push_back(core::COrientedVariant(m_aBaseVariantList[i][j], true));
            m_aBaseOrientedVariantList[i].push_back(core::COrientedVariant(m_aBaseVariantList[i][j], false));
        }
    }
    
    for(unsigned int i=0; i < m_aCalledOrientedVariantList.size(); i++)
    {
        for(unsigned int j=0; j < m_aCalledVariantList[i].size(); j++)
        {
            m_aCalledOrientedVariantList[i].push_back(core::COrientedVariant(m_aCalledVariantList[i][j], true));
            m_aCalledOrientedVariantList[i].push_back(core::COrientedVariant(m_aCalledVariantList[i][j], false));
        }
    }
}

void CVariantProvider::FillAlleleMatchVariantList(SChrIdTuple& a_tuple, std::vector<const CVariant*>& a_rBaseVariants, std::vector<const CVariant*>& a_rCalledVariants)
{
    //Initialize HomozygousOrientedVariantLists
    if(false == m_bIsHomozygousOvarListInitialized)
    {
        m_aBaseHomozygousOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_baseVCF.GetContigs().size());
        m_aCalledHomozygousOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_calledVCF.GetContigs().size());
        m_bIsHomozygousOvarListInitialized = true;
    }
    
    for(unsigned int j=0; j < a_rBaseVariants.size(); j++)
    {
        m_aBaseHomozygousOrientedVariantList[a_tuple.m_nBaseId].push_back(core::COrientedVariant(*a_rBaseVariants[j], 0));
        m_aBaseHomozygousOrientedVariantList[a_tuple.m_nBaseId].push_back(core::COrientedVariant(*a_rBaseVariants[j], 1));
    }
        
    for(unsigned int j=0; j < a_rCalledVariants.size(); j++)
    {
        m_aCalledHomozygousOrientedVariantList[a_tuple.m_nCalledId].push_back(core::COrientedVariant(*a_rCalledVariants[j], 0));
        m_aCalledHomozygousOrientedVariantList[a_tuple.m_nCalledId].push_back(core::COrientedVariant(*a_rCalledVariants[j], 1));
    }
}

std::vector<const CVariant*> CVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes)
{
    std::vector<const CVariant*> result;

    switch (a_uFrom)
    {
        case eBASE:
            for(unsigned int k=0; k < a_VariantIndexes.size(); k++)
                result.push_back(&(m_aBaseVariantList[a_nChrNo][a_VariantIndexes[k]]));
            break;
            
        case eCALLED:
            for(unsigned int k=0; k < a_VariantIndexes.size(); k++)
                result.push_back(&(m_aCalledVariantList[a_nChrNo][a_VariantIndexes[k]]));
            
        default:
            break;
    }

    return result;
}

std::vector<const CVariant*> CVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo)
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


std::vector<const CVariant*> CVariantProvider::GetVariantList(std::vector<const CVariant*>& a_varList, const std::vector<int>& a_VariantIndexes)
{
    std::vector<const CVariant*> result(a_VariantIndexes.size());
    
    for(unsigned int k = 0; k < a_VariantIndexes.size(); k++)
    {
        result[k] = a_varList[a_VariantIndexes[k]];
    }
    return result;
}


std::vector<const core::COrientedVariant*> CVariantProvider::GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, bool a_bIsGenotypeMatch)
{
    if(a_bIsGenotypeMatch)
    {
        unsigned long size = a_uFrom == eBASE ? m_aBaseOrientedVariantList[a_nChrNo].size() : m_aCalledOrientedVariantList[a_nChrNo].size();
        std::vector<const core::COrientedVariant*> result(size);
        
        switch (a_uFrom)
        {
            case eBASE:
                for(unsigned int k=0; k < size; k++)
                    result[k] = &(m_aBaseOrientedVariantList[a_nChrNo][k]);
                break;
                
            case eCALLED:
                for(unsigned int k=0; k < size; k++)
                    result[k] = &(m_aCalledOrientedVariantList[a_nChrNo][k]);
                
            default:
                break;
        }

        return result;
    }
    
    else
    {
        unsigned long size = a_uFrom == eBASE ? m_aBaseHomozygousOrientedVariantList[a_nChrNo].size() : m_aCalledHomozygousOrientedVariantList[a_nChrNo].size();
        std::vector<const core::COrientedVariant*> result(size);
        
        switch (a_uFrom)
        {
            case eBASE:
                for(unsigned int k=0; k < size; k++)
                    result[k] = &(m_aBaseHomozygousOrientedVariantList[a_nChrNo][k]);
                break;
                
            case eCALLED:
                for(unsigned int k=0; k < size; k++)
                    result[k] = &(m_aCalledHomozygousOrientedVariantList[a_nChrNo][k]);
                
            default:
                break;
        }

        return result;
    }
    
}


bool CVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
{
    std::size_t found;
    
    for(int k = 0; k < a_rVariant.m_nAlleleCount; k++)
    {
        std::string allele = a_rVariant.m_alleles[k].m_sequence;
        
        if((int)allele.size() > a_nMaxLength)
            return true;
        found = allele.find('[');
        if(found != std::string::npos)
            return true;
        found = allele.find('<');
        if(found != std::string::npos)
            return true;
        found = allele.find('*');
        if(found != std::string::npos)
            return true;
        found = allele.find('.');
        if(found != std::string::npos)
            return true;
    }
    return false;
}

std::vector<SChrIdTuple>& CVariantProvider::GetChromosomeIdTuples()
{
    return m_aCommonChrTupleList;
}

void CVariantProvider::GetFilterInfo(EVcfName a_vcfType, std::vector<std::string>& a_rFilterNames, std::vector<std::string>& a_rFilterDescriptions)
{
    switch (a_vcfType)
    {
        case eBASE:
            m_baseVCF.GetFilterInfo(a_rFilterNames, a_rFilterDescriptions);
            break;
        case eCALLED:
            m_calledVCF.GetFilterInfo(a_rFilterNames, a_rFilterDescriptions);
            break;
        default:
            break;
    }
}

std::vector<CVariant>& CVariantProvider::GetNotAssessedVariantList(EVcfName a_uFrom, int a_nChrNo)
{
    if(eBASE == a_uFrom)
        return m_aBaseNotAssessedVariantList[a_nChrNo];
    else
        return m_aCalledNotAssessedVariantList[a_nChrNo];
}



void CVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
        pVar->m_variantStatus = a_status;
}


void CVariantProvider::SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const core::COrientedVariant* pOvar : a_rVariantList)
        pOvar->GetVariant().m_variantStatus = a_status;
}


void CVariantProvider::GetContig(std::string a_chrName, SContig& a_rCtg)
{
    for(unsigned int k = 0; k < m_aContigList.size(); k++)
    {
        if(m_aContigList[k].m_chromosomeName == a_chrName)
        {
            a_rCtg.m_chromosomeName = a_chrName;
            a_rCtg.m_pRefSeq = m_aContigList[k].m_pRefSeq;
            a_rCtg.m_nRefLength = m_aContigList[k].m_nRefLength;
            break;
        }
    }
}

const std::vector<SVcfContig>& CVariantProvider::GetContigs() const
{
    return m_calledVCF.GetContigs();
}

bool CVariantProvider::IsHomRef(const CVariant& a_rVariant) const
{
    bool res = true;
    
    for(int k = 0; k < a_rVariant.m_nZygotCount; k++)
    {
        if(a_rVariant.m_genotype[k] != 0)
        {
            res = false;
            break;
        }
    }
    return res;
}

void CVariantProvider::FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom)
{
    if(a_rVariantList.size() == 0)
        return;
    
    std::vector<std::vector<CVariant>>* allVarList = a_uFrom == eBASE ? &m_aBaseVariantList : &m_aCalledVariantList;
    
    unsigned int varItr[2];
    varItr[0] = 0;
    varItr[1] = 0;
    
    int currentChrId = a_rVariantList[0].m_nChrId;
    
    std::vector<CVariant> affectedVariants;
    std::vector<CVariant> genotypeAffectedVariants;
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
    {
        if(a_rVariantList[k].m_chrName == "1" && a_rVariantList[k].m_nOriginalPos == 44776634)
        {
            int asd = 0;
            asd ++;
        }
        

        for(int i = 0; i < 2; i++)
        {
            if(a_rVariantList[k].m_alleles[i].m_bIsIgnored || a_rVariantList[k].m_alleles[i].m_bIsTrimmed)
                continue;
            
            //The maximum possible trimming nucleotide count from start and end of each allele
            unsigned int canTrimStart, canTrimEnd;
            a_rVariantList[k].GetMaxTrimStartEnd(i, canTrimStart, canTrimEnd);
            
            std::vector<CVariant> tmpoverlapVariants;
            
            if(currentChrId != a_rVariantList[k].m_nChrId)
            {
                varItr[0] = 0;
                varItr[1] = 0;
                currentChrId = a_rVariantList[k].m_nChrId;
            }
            
            while(varItr[i] < (*allVarList)[currentChrId].size() && a_rVariantList[k].m_alleles[i].m_nStartPos > (*allVarList)[currentChrId][varItr[i]].m_nEndPos)
                varItr[i]++;
            
            if(varItr[i] == (*allVarList)[currentChrId].size())
                continue;
            
            unsigned int secondItr = varItr[i];

            while(secondItr < (*allVarList)[currentChrId].size() && a_rVariantList[k].m_alleles[i].m_nEndPos >= (*allVarList)[currentChrId][secondItr].m_nStartPos)
            {
                if(isOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
                             a_rVariantList[k].m_alleles[i].m_nEndPos,
                             (*allVarList)[currentChrId][secondItr].m_nStartPos,
                             (*allVarList)[currentChrId][secondItr].m_nEndPos))
                    tmpoverlapVariants.push_back((*allVarList)[currentChrId][secondItr]);
                secondItr++;
            }
            
            //Trim variants as standard if there is no overlap
            if(tmpoverlapVariants.size() == 0)
                a_rVariantList[k].TrimVariant(i);
            
            else
            {
                affectedVariants.push_back(a_rVariantList[k]);
                for(unsigned int ovarItr = 0; ovarItr < tmpoverlapVariants.size(); ovarItr++)
                {
                    //Check each allele of overlapping variant
                    for(int tmpItr = 0; tmpItr < 2; tmpItr++)
                    {
                        //If the allele does not overlap, continue
                        if(!isOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
                                     a_rVariantList[k].m_alleles[i].m_nEndPos,
                                     tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos,
                                     tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos))
                            continue;
                        
                        unsigned int overlapStart = std::max(a_rVariantList[k].m_alleles[i].m_nStartPos, tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos);
                        unsigned int overlapEnd = std::min(a_rVariantList[k].m_alleles[i].m_nEndPos, tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos);
                        
                        //Trim from beginning
                        if(a_rVariantList[k].m_alleles[i].m_nStartPos == (int)overlapStart)
                        {
                            if(a_rVariantList[k].m_alleles[i].m_nStartPos + (int)canTrimStart >= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos)
                            {
                                int toClip = overlapEnd - a_rVariantList[k].m_alleles[i].m_nStartPos;
                                if(toClip > 0 && toClip <= (int)canTrimStart)
                                {
                                    a_rVariantList[k].TrimVariant(i, toClip, 0);
                                    canTrimStart -= toClip;
                                    continue;
                                }
                            }
                        }
                        
                        //Trim from beginning or end --
                        else if(a_rVariantList[k].m_alleles[i].m_nEndPos > (int) overlapEnd)
                        {
                            //Try to trim from end
                            if(a_rVariantList[k].m_alleles[i].m_nEndPos - (int)canTrimEnd <= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos)
                            {
                                int toClip = a_rVariantList[k].m_alleles[i].m_nEndPos - overlapStart;
                                if(toClip > 0 & toClip <= (int)canTrimEnd)
                                {
                                    a_rVariantList[k].TrimVariant(i, 0, toClip);
                                    canTrimEnd -= toClip;
                                    continue;
                                }
                            }
                            
                            //Try to trim from front
                            if(a_rVariantList[k].m_alleles[i].m_nStartPos + (int)canTrimStart >= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos)
                            {
                                int toClip = overlapEnd - a_rVariantList[k].m_alleles[i].m_nStartPos;
                                if(toClip > 0 && toClip <= (int)canTrimStart)
                                {
                                    a_rVariantList[k].TrimVariant(i, toClip, 0);
                                    canTrimStart -= toClip;
                                    continue;
                                }
                            }

                        }
                        
                        //Trim from end
                        else
                        {
                            if(a_rVariantList[k].m_alleles[i].m_nEndPos - (int)canTrimEnd < tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos)
                            {
                                int toClip = a_rVariantList[k].m_alleles[i].m_nEndPos - overlapStart;
                                if(toClip > 0 & toClip <= (int)canTrimEnd)
                                {
                                    a_rVariantList[k].TrimVariant(i, 0, toClip);
                                    canTrimEnd -= toClip;
                                    continue;
                                }
                            }
                        }
                    }
                }
                
                if(!a_rVariantList[k].m_alleles[i].m_bIsTrimmed)
                    a_rVariantList[k].TrimVariant(i);
                else
                    a_rVariantList[k].TrimVariant(i, canTrimStart, canTrimEnd);
            }
        }
        
        
        //Update Variant Range
        int minStart = INT_MAX;
        int maxEnd = -1;
        
        for(int i =0; i < a_rVariantList[k].m_nAlleleCount; i++)
        {
            if(!a_rVariantList[k].m_alleles[i].m_bIsIgnored)
            {
                maxEnd = std::max(maxEnd, static_cast<int>(a_rVariantList[k].m_alleles[i].m_nEndPos));
                minStart = std::min(minStart, static_cast<int>(a_rVariantList[k].m_alleles[i].m_nStartPos));
            }
        }
        
        a_rVariantList[k].m_nEndPos = maxEnd == -1 ? a_rVariantList[k].m_nOriginalPos : maxEnd;
        a_rVariantList[k].m_nStartPos = minStart == INT_MAX ? a_rVariantList[k].m_nOriginalPos : minStart;
        
    }
    
}

void CVariantProvider::AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom)
{
    std::vector<std::vector<CVariant>>& variantList = (a_uFrom == eCALLED) ? m_aCalledVariantList : m_aBaseVariantList;
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
        variantList[a_rVariantList[k].m_nChrId].push_back(a_rVariantList[k]);
}












