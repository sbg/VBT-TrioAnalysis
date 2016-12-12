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
    for(int k= 0; k < CHROMOSOME_COUNT; k++)
    {
        if(m_aContigList[k].m_nRefLength > 0 && m_aContigList[k].m_pRefSeq != 0)
            delete m_aContigList[k].m_pRefSeq;
    }
}


bool CVariantProvider::InitializeReaders(const SConfig& a_rConfig)
{
    bool bIsSuccess;

    m_config = a_rConfig;
    
    bIsSuccess = m_baseVCF.Open(a_rConfig.m_pBaseVcfFileName);
    if(!bIsSuccess)
        std::cout << "Baseline VCF file is unable to open!: " << a_rConfig.m_pBaseVcfFileName << std::endl;
    m_baseVCF.setID(0);
    bIsSuccess = m_calledVCF.Open(a_rConfig.m_pCalledVcfFileName);
    if(!bIsSuccess)
        std::cout << "Called VCF file is unable to open!: " << a_rConfig.m_pCalledVcfFileName << std::endl;
    m_calledVCF.setID(1);
    
    bIsSuccess = m_fastaParser.OpenFastaFile(a_rConfig.m_pFastaFileName);
    if(!bIsSuccess)
        std::cout << "FASTA file is unable to open!: " << a_rConfig.m_pCalledVcfFileName << std::endl;
    
    if(bIsSuccess)
    {
        FillVariantLists();
        FillOrientedVariantLists();
    }
    
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        if(m_aBaseVariantList[k].size() != 0 && m_aCalledVariantList[k].size() == 0)
        {
            std::cout << "Called VCF does not contain Chromosome " << k+1 << ".The chromosome will be filtered out from the comparison!" << std::endl;
            m_aCalledVariantList[k].clear();
        }
        else if(m_aBaseVariantList[k].size() == 0 && m_aCalledVariantList[k].size() != 0)
        {
            std::cout << "Baseline VCF does not contain Chromosome " << k+1 << ".The chromosome will be filtered out from the comparison!" << std::endl;
            m_aBaseVariantList[k].clear();
        }
        else if(m_aBaseVariantList[k].size() != 0 && m_aCalledVariantList[k].size() != 0)
        {
            //Read contig from FASTA file for the given chromosome
            bool bIsSuccess2 = m_fastaParser.FetchNewChromosome(m_aBaseVariantList[k][0].m_chrName, m_aContigList[k]);
            if(!bIsSuccess2)
            {
                std::cout << "Chromosome " << k+1 << "will be filtered out from the comparison since reference FASTA could not read or it does not contain given chromosome" << std::endl;
                m_aCalledVariantList[k].clear();
                m_aBaseVariantList[k].clear();
            }
            else
            {
                //Trim the variants which are out of bound according to FASTA file
                while(m_aBaseVariantList[k].back().GetEnd() > m_aContigList[k].m_nRefLength)
                    m_aBaseVariantList[k].pop_back();
                while(m_aCalledVariantList[k].back().GetEnd() > m_aContigList[k].m_nRefLength)
                    m_aCalledVariantList[k].pop_back();
            }
        }
    }
    
    return bIsSuccess;
    
    //for(int k = 0; k < m_aCalledVariantList[21].size(); k++)
    //   std::cout << k << ": " << m_aCalledVariantList[21][k].ToString() << std::endl;
}


void CVariantProvider::FillVariantLists()
{
    
    CVariant variant;
    int id = 0;
    while(m_baseVCF.GetNextRecord(&variant, id++, m_config))
    {
        if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            continue;
        
        if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            continue;
        
        PushVariant(variant, m_aBaseVariantList[variant.m_nChrId-1]);
    }
    
    while(m_calledVCF.GetNextRecord(&variant, id++, m_config))
    {
        
        if(m_config.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            continue;
        
        if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            continue;
        
        PushVariant(variant, m_aCalledVariantList[variant.m_nChrId-1]);
    }
}

void CVariantProvider::FillOrientedVariantLists()
{
    for(int i=0; i < CHROMOSOME_COUNT; i++)
    {
        for(int j=0; j < m_aBaseVariantList[i].size(); j++)
        {
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], true));
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], false));
        }
        
        for(int j=0; j < m_aCalledVariantList[i].size(); j++)
        {
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], true));
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], false));
        }
    }
}

void CVariantProvider::GetContig(int a_nChrId, SContig& a_rContig) const
{
    a_rContig.m_chromosome = m_aContigList[a_nChrId].m_chromosome;
    a_rContig.m_nRefLength = m_aContigList[a_nChrId].m_nRefLength;
    a_rContig.m_nChrId = a_nChrId;
    a_rContig.m_pRefSeq = m_aContigList[a_nChrId].m_pRefSeq;
}

const CVariant* CVariantProvider::GetVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId) const
{
    switch(a_uFrom)
    {
        case eBASE:
            return &m_aBaseVariantList[a_nChrNo][a_nVariantId];
        case eCALLED:
            return &m_aCalledVariantList[a_nChrNo][a_nVariantId];
    }
}


const COrientedVariant* CVariantProvider::GetOrientedVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId, bool a_bOrientation) const
{
    int nExtra = a_bOrientation ? 0 : 1;
    
   if(a_uFrom == eBASE)
   {
       const COrientedVariant* o1 = &m_aBaseOrientedVariantList[a_nChrNo][a_nVariantId*2 + nExtra];
       return o1;
   }
   
   else if(a_uFrom == eCALLED)
   {
       const COrientedVariant* o2 = &m_aCalledOrientedVariantList[a_nChrNo][a_nVariantId*2 + nExtra];
       return o2;
   }
    
    return 0;
}

int CVariantProvider::GetVariantListSize(EVcfName a_uFrom, int a_nChrNo) const
{
    switch(a_uFrom)
    {
        case eBASE:
            return (int)m_aBaseVariantList[a_nChrNo].size();
        case eCALLED:
            return (int)m_aCalledVariantList[a_nChrNo].size();
    }
}



std::vector<CVariant> CVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo, std::vector<int> a_VariantIndexes)
{
    std::vector<CVariant> result;

    switch (a_uFrom)
    {
        case eBASE:
            for(int k=0; k< a_VariantIndexes.size(); k++)
                result.push_back(m_aBaseVariantList[a_nChrNo][a_VariantIndexes[k]]);
            break;
            
        case eCALLED:
            for(int k=0; k< a_VariantIndexes.size(); k++)
                result.push_back(m_aCalledVariantList[a_nChrNo][a_VariantIndexes[k]]);
            
        default:
            break;
    }

    return result;
}

bool CVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
{
    std::size_t found;
    
    for(int k = 0; k < a_rVariant.m_nAlleleCount; k++)
    {
        std::string allele = a_rVariant.m_alleles[k].m_sequence;
        
        if(allele.size() > a_nMaxLength)
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
        if(a_rVariant.GetEnd() < a_rVecToPush[size-k].GetEnd())
            k++;
        
//        if(a_rVariant.GetEnd() == a_rVecToPush[size-k].GetEnd())
//        {
//            if(a_rVariant.m_alleles[0].m_sequence.length() < a_rVecToPush[size-k].m_alleles[0].m_sequence.length())
//                k++;
//        }
    }
    
    it = a_rVecToPush.begin() + (size - k) + 1;
    a_rVecToPush.insert(it, a_rVariant);
}




