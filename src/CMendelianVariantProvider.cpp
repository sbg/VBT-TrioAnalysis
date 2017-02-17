//
//  CMendelianVariantProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianVariantProvider.h"
#include <iostream>


bool CMendelianVariantProvider::InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig)
{
    bool bIsSuccessFather = true;
    bool bIsSuccessMother = true;
    bool bIsSuccessChild = true;
    bool bIsSuccessFasta = true;
    
    
    m_motherChildConfig = a_rMotherChildConfig;
    m_fatherChildConfig = a_rFatherChildConfig;
    
    //Open FATHER vcf file
    bIsSuccessFather = m_FatherVcf.Open(a_rFatherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessFather)
        std::cout << "Father VCF file is unable to open!: " << a_rFatherChildConfig.m_pBaseVcfFileName << std::endl;
    
    //Set sample name of FATHER
    else if (true == m_fatherChildConfig.m_bBaseSampleEnabled)
        m_FatherVcf.SelectSample(m_fatherChildConfig.m_pBaseSample);
    else
    {
        std::vector<std::string> sampleNames;
        m_FatherVcf.GetSampleNames(sampleNames);
        m_FatherVcf.SelectSample(sampleNames[0]);
    }
    
    //Open MOTHER vcf file
    bIsSuccessMother = m_MotherVcf.Open(a_rMotherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessMother)
        std::cout << "Mother VCF file is unable to open!: " << a_rMotherChildConfig.m_pBaseVcfFileName << std::endl;
    
    //Set sample name of MOTHER
    else if (true == m_motherChildConfig.m_bBaseSampleEnabled)
        m_MotherVcf.SelectSample(m_motherChildConfig.m_pBaseSample);
    else
    {
        std::vector<std::string> sampleNames;
        m_MotherVcf.GetSampleNames(sampleNames);
        m_MotherVcf.SelectSample(sampleNames[0]);
    }
    
    //Open CHILD vcf file
    bIsSuccessChild = m_ChildVcf.Open(a_rFatherChildConfig.m_pCalledVcfFileName);
    if(!bIsSuccessChild)
        std::cout << "Child VCF file is unable to open!: " << a_rFatherChildConfig.m_pCalledVcfFileName << std::endl;
    
    //Set sample name of CHILD
    else if (true == m_fatherChildConfig.m_bCalledSampleEnabled)
        m_ChildVcf.SelectSample(m_fatherChildConfig.m_pCalledSample);
    else
    {
        std::vector<std::string> sampleNames;
        m_ChildVcf.GetSampleNames(sampleNames);
        m_ChildVcf.SelectSample(sampleNames[0]);
    }
    
    // OPEN FASTA FILE
    bIsSuccessFasta = m_fastaParser.OpenFastaFile(a_rFatherChildConfig.m_pFastaFileName);
    if(!bIsSuccessFasta)
        std::cout << "FASTA file is unable to open!: " << a_rFatherChildConfig.m_pFastaFileName << std::endl;
    
    if(bIsSuccessChild && bIsSuccessMother && bIsSuccessFather && bIsSuccessFasta)
    {
        //Fill the variants of 3 vcf file
        FillVariants();
        
        //Get the common chromosome ids and clear the uncommon variants from provider
        std::vector<int> commonChromosomeIds = GetCommonChromosomes(true);
        
        //Fill the oriented variants of 3 vcf for genotype matching
        FillGenotypeMatchOrientedVariants(commonChromosomeIds);
        
        //Fill the oriented variants of 3 vcf for allele matching
        FillAlleleMatchOrientedVariants(commonChromosomeIds);
        
        for(int k = 0; k < commonChromosomeIds.size(); k++)
        {
            //Read contig from FASTA file for the given chromosome
            std::cout << "Reading reference of chromosome " << m_aChildVariantList[commonChromosomeIds[k]][0].m_chrName << " from the FASTA file" << std::endl;
            bool bIsSuccess2 = m_fastaParser.FetchNewChromosome(m_aChildVariantList[commonChromosomeIds[k]][0].m_chrName, m_aContigList[commonChromosomeIds[k]]);
            if(!bIsSuccess2)
            {
                std::cout << "Chromosome " << k+1 << "will be filtered out from the comparison since reference FASTA could not read or it does not contain given chromosome" << std::endl;
                m_aFatherVariantList[commonChromosomeIds[k]].clear();
                m_aMotherVariantList[commonChromosomeIds[k]].clear();
                m_aChildVariantList[commonChromosomeIds[k]].clear();
                
                m_aFatherOrientedVariantList[commonChromosomeIds[k]].clear();
                m_aMotherOrientedVariantList[commonChromosomeIds[k]].clear();
                m_aChildOrientedVariantList[commonChromosomeIds[k]].clear();
            }
            
            else
            {
                //Trim the variants which are out of bound according to FASTA file
                while(m_aChildVariantList[commonChromosomeIds[k]].back().GetEnd() > m_aContigList[commonChromosomeIds[k]].m_nRefLength)
                    m_aChildVariantList[commonChromosomeIds[k]].pop_back();
                while(m_aFatherVariantList[commonChromosomeIds[k]].back().GetEnd() > m_aContigList[commonChromosomeIds[k]].m_nRefLength)
                    m_aFatherVariantList[commonChromosomeIds[k]].pop_back();
                while(m_aMotherVariantList[commonChromosomeIds[k]].back().GetEnd() > m_aContigList[commonChromosomeIds[k]].m_nRefLength)
                    m_aMotherVariantList[commonChromosomeIds[k]].pop_back();
            }
        }
        
    }

    return bIsSuccessChild && bIsSuccessMother && bIsSuccessFather && bIsSuccessFasta;
}


void CMendelianVariantProvider::FillVariants()
{
    CVariant variant;
    int id = 0;
    std::string preChrId = "";

    //READ VARIANTS OF FATHER
    while(m_FatherVcf.GetNextRecord(&variant, id++, m_fatherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of Parent[FATHER] vcf" << std::endl;
            //if(preChrId == "2")
            //    break;
        }
        
        if(variant.m_nChrId == -1)
            continue;
        else if(variant.m_nChrId > 23) //DISCARD Y and MT chromosome
            PushVariant(variant, m_aFatherNotAssessedVariantList[variant.m_nChrId-1]);

        else if(m_fatherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            PushVariant(variant, m_aFatherNotAssessedVariantList[variant.m_nChrId-1]);
        
        else if(IsStructuralVariant(variant, m_fatherChildConfig.m_nMaxVariantSize))
            PushVariant(variant, m_aFatherNotAssessedVariantList[variant.m_nChrId-1]);
        
        else
            PushVariant(variant, m_aFatherVariantList[variant.m_nChrId-1]);
    }

    preChrId = "";

    //READ VARIANTS OF MOTHER
    while(m_MotherVcf.GetNextRecord(&variant, id++, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of Parent[MOTHER] vcf" << std::endl;
            //if(preChrId == "2")
            //    break;

        }

        if(variant.m_nChrId == -1)
            continue;
        else if(variant.m_nChrId > 23) //DISCARD Y and MT chromosome
            PushVariant(variant, m_aMotherNotAssessedVariantList[variant.m_nChrId-1]);
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            PushVariant(variant, m_aMotherNotAssessedVariantList[variant.m_nChrId-1]);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            PushVariant(variant, m_aMotherNotAssessedVariantList[variant.m_nChrId-1]);
        
        else
            PushVariant(variant, m_aMotherVariantList[variant.m_nChrId-1]);
    }

    preChrId = "";
    
    //READ VARIANTS OF CHILD
    while(m_ChildVcf.GetNextRecord(&variant, id++, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of Child vcf" << std::endl;
            //if(preChrId == "2")
            //    break;

        }
        
        if(variant.m_nChrId == -1)
            continue;
        else if(variant.m_nChrId > 23) //DISCARD Y and MT chromosome
            PushVariant(variant, m_aChildNotAssessedVariantList[variant.m_nChrId-1]);
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            PushVariant(variant, m_aChildNotAssessedVariantList[variant.m_nChrId-1]);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            PushVariant(variant, m_aChildNotAssessedVariantList[variant.m_nChrId-1]);
        
        else
            PushVariant(variant, m_aChildVariantList[variant.m_nChrId-1]);
    }
}



void CMendelianVariantProvider::FillGenotypeMatchOrientedVariants(std::vector<int>& a_aCommonChromosomes)
{
    for(int i=0; i < a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aFatherVariantList[a_aCommonChromosomes[i]].size(); j++, k+=2)
        {
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i]][j], true);
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i]][j], false);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aChildVariantList[a_aCommonChromosomes[i]].size(); j++,k+=2)
        {
            m_aChildOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i]][j], true);
            m_aChildOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i]][j], false);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aMotherVariantList[a_aCommonChromosomes[i]].size(); j++, k+=2)
        {
            m_aMotherOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i]][j], true);
            m_aMotherOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i]][j], false);
        }
    }
}

void CMendelianVariantProvider::FillAlleleMatchOrientedVariants(std::vector<int>& a_aCommonChromosomes)
{
    for(int i=0; i < a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aFatherVariantList[a_aCommonChromosomes[i]].size(); j++, k+=2)
        {
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i]][j], 0);
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i]][j], 1);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aChildVariantList[a_aCommonChromosomes[i]].size(); j++,k+=2)
        {
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i]][j], 0);
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i]][j], 1);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]] = std::vector<COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i]].size() * 2);
        for(int j=0, k=0; j < m_aMotherVariantList[a_aCommonChromosomes[i]].size(); j++, k+=2)
        {
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k] = COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i]][j], 0);
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i]][k+1] = COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i]][j], 1);
        }
    }
}

bool CMendelianVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
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

void CMendelianVariantProvider::PushVariant(CVariant& a_rVariant, std::vector<CVariant>& a_rVecToPush)
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


std::vector<int> CMendelianVariantProvider::GetCommonChromosomes(bool a_bIsCalledInProviderInitialization)
{
    std::vector<int> commonChrIds;

    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        if(m_aChildVariantList[k].size() > LEAST_VARIANT_THRESHOLD && m_aFatherVariantList[k].size() > LEAST_VARIANT_THRESHOLD && m_aMotherVariantList[k].size() > LEAST_VARIANT_THRESHOLD)
            commonChrIds.push_back(k);
        else if(true == a_bIsCalledInProviderInitialization)
        {
            std::cout << "Warning! Chromosome " << k+1 << " is not contained by all three vcf files. Variants will be filtered out from comparison" << std::endl;
            
            //Clear redundant variants
            m_aChildVariantList[k].clear();
            m_aMotherVariantList[k].clear();
            m_aFatherVariantList[k].clear();
        }
    }
    
    return commonChrIds;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList;
    
    switch (a_uFrom) {
        case eCHILD:
            for(int k = 0; k < m_aChildVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][k]);
            break;
        case eFATHER:
            for(int k = 0; k < m_aFatherVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][k]);
            break;
        case eMOTHER:
            for(int k = 0; k < m_aMotherVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aMotherVariantList[a_nChrNo][k]);
            break;

        default:
            break;
    }
    
    return varList;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_nIndexList) const
{
    std::vector<const CVariant*> varList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            for(int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eFATHER:
            for(int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eMOTHER:
            for(int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aMotherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
            
        default:
            break;
    }
    
    return varList;
}

std::vector<const COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch) const
{
    std::vector<const COrientedVariant*> ovarList;
    
    if(false == a_bIsAlleleMatch)
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < m_aChildOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][k]);
                break;
            case eFATHER:
                for(int k = 0; k < m_aFatherOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][k]);
                break;
            case eMOTHER:
                for(int k = 0; k < m_aMotherOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aMotherOrientedVariantList[a_nChrNo][k]);
                break;
                
            default:
                break;
        }
    }
    
    else
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < m_aChildAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
            case eFATHER:
                for(int k = 0; k < m_aFatherAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
            case eMOTHER:
                for(int k = 0; k < m_aMotherAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aMotherAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
                
            default:
                break;
        }
    }
    
    return ovarList;
}

std::vector<const COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch, const std::vector<int>& a_nIndexList) const
{
    std::vector<const COrientedVariant*> ovarList;

    if(false == a_bIsAlleleMatch)
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][a_nIndexList[k]*2+1]);
                }

                break;
            case eFATHER:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2+1]);
                }
                break;
            case eMOTHER:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aMotherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aMotherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2+1]);
                }
                break;
                
            default:
                break;
        }
    }
    
    else
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2 +1]);
                }
                break;
            case eFATHER:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2 +1]);
                }
                break;
            case eMOTHER:
                for(int k = 0; k < a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aMotherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aMotherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2 +1]);
                }
                break;
                
            default:
                break;
        }
    
    
    }
    
    return ovarList;

}

void CMendelianVariantProvider::GetContig(int a_nChrId, SContig& a_rContig) const
{
    a_rContig.m_chromosome = m_aContigList[a_nChrId].m_chromosome;
    a_rContig.m_nRefLength = m_aContigList[a_nChrId].m_nRefLength;
    a_rContig.m_nChrId = a_nChrId;
    a_rContig.m_pRefSeq = m_aContigList[a_nChrId].m_pRefSeq;
}


void CMendelianVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
    {
        if(a_status == eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pVar->m_variantStatus != eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pVar-> m_variantStatus == eNOT_ASSESSED)
            pVar->m_variantStatus = a_status;
        else
            break;
    }
}


void CMendelianVariantProvider::SetVariantStatus(const std::vector<const COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const COrientedVariant* pOVar : a_rVariantList)
    {
        if(a_status == eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pOVar->GetVariant().m_variantStatus != eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pOVar->GetVariant().m_variantStatus == eNOT_ASSESSED)
            pOVar->GetVariant().m_variantStatus = a_status;
        else
            break;
    }
}












