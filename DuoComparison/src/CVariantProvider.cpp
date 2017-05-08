//
//  CVariantProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVariantProvider.h"
#include <iostream>
#include "COrientedVariant.h"

CVariantProvider::~CVariantProvider()
{
    for(int k= 0; k < (int)m_aContigList.size(); k++)
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
        FillVariantLists();
        FillOrientedVariantLists();
        
        //Set the common chromosome list for parsing
        SetChromosomeIdTuples();
        
        //Initialize contig list
        m_aContigList = std::vector<SContig>(m_aCommonChrTupleList.size());
    
        for(SChrIdTuple tuple : m_aCommonChrTupleList)
        {
            if(m_aBaseVariantList[tuple.m_nBaseId].size() != 0 && m_aCalledVariantList[tuple.m_nCalledId].size() == 0)
            {
                std::cout << "Called VCF does not contain Chromosome " << m_aCommonChrTupleList[tuple.m_nTupleIndex].m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aCalledVariantList[tuple.m_nCalledId].clear();
                m_aCalledOrientedVariantList[tuple.m_nCalledId].clear();
            }
            else if(m_aBaseVariantList[tuple.m_nBaseId].size() == 0 && m_aCalledVariantList[tuple.m_nCalledId].size() != 0)
            {
                std::cout << "Baseline VCF does not contain Chromosome " << m_aCommonChrTupleList[tuple.m_nTupleIndex].m_chrName << ".The chromosome will be filtered out from the comparison!" << std::endl;
                m_aBaseVariantList[tuple.m_nBaseId].clear();
                m_aBaseOrientedVariantList[tuple.m_nBaseId].clear();
            }
            else if(m_aBaseVariantList[tuple.m_nBaseId].size() != 0 && m_aCalledVariantList[tuple.m_nCalledId].size() != 0)
            {
                //Read contig from FASTA file for the given chromosome
                std::cout << "Reading reference of chromosome " << m_aCommonChrTupleList[tuple.m_nTupleIndex].m_chrName << " from the FASTA file" << std::endl;
                bool bIsSuccess2 = m_fastaParser.FetchNewChromosome(m_aCommonChrTupleList[tuple.m_nTupleIndex].m_chrName, m_aContigList[tuple.m_nTupleIndex]);
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

    while(m_baseVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(m_calledVCF.GetContigId(variant.m_chrName) == -1)
            continue;
        
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of base vcf" << std::endl;
        }
        
        if(variant.m_nChrId < 8)
            continue;
        if(variant.m_nChrId > 12)
            break;

        if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aBaseNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        else
            PushVariant(variant, m_aBaseVariantList[variant.m_nChrId]);
    }
    
    preChrId = "";
    
    while(m_calledVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(m_baseVCF.GetContigId(variant.m_chrName) == -1)
            continue;
        
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of called vcf" << std::endl;
        }
        
        if(variant.m_nChrId < 8)
            continue;
        if(variant.m_nChrId > 12)
            break;
        
        else if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bSNPOnly && variant.GetVariantType() != eSNP)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(m_config.m_bINDELOnly && variant.GetVariantType() != eINDEL)
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            m_aCalledNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        else
            PushVariant(variant, m_aCalledVariantList[variant.m_nChrId]);
    }
}

void CVariantProvider::FillOrientedVariantLists()
{
    //Initialize OrientedVariantLists
    m_aBaseOrientedVariantList = std::vector<std::vector<COrientedVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledOrientedVariantList = std::vector<std::vector<COrientedVariant>>(m_calledVCF.GetContigs().size());
    
    
    for(int i=0; i < m_aBaseOrientedVariantList.size(); i++)
    {
        for(int j=0; j < (int)m_aBaseVariantList[i].size(); j++)
        {
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], true));
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], false));
        }
    }
    
    for(int i=0; i < m_aCalledOrientedVariantList.size(); i++)
    {
        for(int j=0; j < (int)m_aCalledVariantList[i].size(); j++)
        {
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], true));
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], false));
        }
    }
}

void CVariantProvider::FillAlleleMatchVariantList(SChrIdTuple& a_tuple, std::vector<const CVariant*>& a_rBaseVariants, std::vector<const CVariant*>& a_rCalledVariants)
{
    //Initialize HomozygousOrientedVariantLists
    m_aBaseHomozygousOrientedVariantList = std::vector<std::vector<COrientedVariant>>(m_baseVCF.GetContigs().size());
    m_aCalledHomozygousOrientedVariantList = std::vector<std::vector<COrientedVariant>>(m_calledVCF.GetContigs().size());
    
    for(int j=0; j < (int)a_rBaseVariants.size(); j++)
    {
        m_aBaseHomozygousOrientedVariantList[a_tuple.m_nBaseId].push_back(COrientedVariant(*a_rBaseVariants[j], 0));
        m_aBaseHomozygousOrientedVariantList[a_tuple.m_nBaseId].push_back(COrientedVariant(*a_rBaseVariants[j], 1));
    }
        
    for(int j=0; j < (int)a_rCalledVariants.size(); j++)
    {
        m_aCalledHomozygousOrientedVariantList[a_tuple.m_nCalledId].push_back(COrientedVariant(*a_rCalledVariants[j], 0));
        m_aCalledHomozygousOrientedVariantList[a_tuple.m_nCalledId].push_back(COrientedVariant(*a_rCalledVariants[j], 1));
    }
}

std::vector<const CVariant*> CVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes)
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
    
    for(int k = 0; k < (int)a_VariantIndexes.size(); k++)
    {
        result[k] = a_varList[a_VariantIndexes[k]];
    }
    return result;
}


std::vector<const COrientedVariant*> CVariantProvider::GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, bool a_bIsGenotypeMatch)
{
    if(a_bIsGenotypeMatch)
    {
        unsigned long size = a_uFrom == eBASE ? m_aBaseOrientedVariantList[a_nChrNo].size() : m_aCalledOrientedVariantList[a_nChrNo].size();
        std::vector<const COrientedVariant*> result(size);
        
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
    
    else
    {
        unsigned long size = a_uFrom == eBASE ? m_aBaseHomozygousOrientedVariantList[a_nChrNo].size() : m_aCalledHomozygousOrientedVariantList[a_nChrNo].size();
        std::vector<const COrientedVariant*> result(size);
        
        switch (a_uFrom)
        {
            case eBASE:
                for(int k=0; k < (int)size; k++)
                    result[k] = &(m_aBaseHomozygousOrientedVariantList[a_nChrNo][k]);
                break;
                
            case eCALLED:
                for(int k=0; k< (int)size; k++)
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


void CVariantProvider::PushVariant(CVariant& a_rVariant, std::vector<CVariant>& a_rVecToPush)
{
    int k = 0;
    int size = static_cast<int>(a_rVecToPush.size())-1;
    std::vector<CVariant>::iterator it;
    
    if(size == -1)
    {
        a_rVecToPush.push_back(a_rVariant);
        return;
    }
    
    while(k <= size && a_rVariant.GetStart() < a_rVecToPush[size-k].GetStart())
        k++;
    
    if(a_rVariant.GetStart() == a_rVecToPush[size-k].GetStart())
    {
        if(a_rVariant.m_bIsFirstNucleotideTrimmed && !a_rVecToPush[size-k].m_bIsFirstNucleotideTrimmed)
            k++;
        
        else if(a_rVariant.m_bIsFirstNucleotideTrimmed == a_rVecToPush[size-k].m_bIsFirstNucleotideTrimmed)
        {
            if(a_rVariant.GetEnd() < a_rVecToPush[size-k].GetEnd())
                k++;
        }
        
//        if(a_rVariant.GetEnd() == a_rVecToPush[size-k].GetEnd())
//        {
//            if(a_rVariant.m_alleles[0].m_sequence.length() < a_rVecToPush[size-k].m_alleles[0].m_sequence.length())
//                k++;
//        }
    }
    
    it = a_rVecToPush.begin() + (size - k) + 1;
    a_rVecToPush.insert(it, a_rVariant);
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


void CVariantProvider::SetVariantStatus(const std::vector<const COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const COrientedVariant* pOvar : a_rVariantList)
        pOvar->GetVariant().m_variantStatus = a_status;
}


void CVariantProvider::GetContig(std::string a_chrName, SContig& a_rCtg)
{
    for(int k = 0; k < (int)m_aContigList.size(); k++)
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











