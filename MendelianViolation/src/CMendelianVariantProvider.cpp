/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  CMendelianVariantProvider.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 1/31/17.
 *
 */

#include <algorithm>
#include "CMendelianVariantProvider.h"
#include "CSimplePEDParser.h"
#include "CSimpleBEDParser.h"
#include "Utils/CUtils.h"
#include <iostream>
#include <sstream>

using namespace mendelian;

bool CMendelianVariantProvider::OpenVcfFiles()
{
    bool bIsSuccessFather = true;
    bool bIsSuccessMother = true;
    bool bIsSuccessChild = true;
    
    //Open FATHER vcf file
    bIsSuccessFather = m_FatherVcf.Open(m_fatherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessFather)
        std::cerr << "Father VCF file is unable to open!: " << m_fatherChildConfig.m_pBaseVcfFileName << std::endl;
    
    //Open MOTHER vcf file
    bIsSuccessMother = m_MotherVcf.Open(m_motherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessMother)
        std::cerr << "Mother VCF file is unable to open!: " << m_motherChildConfig.m_pBaseVcfFileName << std::endl;
    
    //Open CHILD vcf file
    bIsSuccessChild = m_ChildVcf.Open(m_fatherChildConfig.m_pCalledVcfFileName);
    if(!bIsSuccessChild)
        std::cerr << "Child VCF file is unable to open!: " << m_fatherChildConfig.m_pCalledVcfFileName << std::endl;
    
    if(!bIsSuccessChild || !bIsSuccessFather || !bIsSuccessMother)
        return false;
    
    if(true == m_fatherChildConfig.m_bIsReadINFO)
    {
        m_FatherVcf.GetInfoNames(m_fatherChildConfig.m_infotags);
        m_MotherVcf.GetInfoNames(m_motherChildConfig.m_infotags);
        m_ChildVcf.GetInfoNames(m_motherChildConfig.m_infotags);
    }
    
    //Set Sample name from PED file
    if(true == m_motherChildConfig.m_bInitializeFromPED)
    {
        //This vector will be used if the program initialized with a PED file
        CSimplePEDParser pedParser;
        std::vector<std::string> parentChildSampleIds;
        
        pedParser.ParsePedigree(m_motherChildConfig.m_pPedigreeFileName);
        std::vector<std::string> sampleNamesChild;
        std::vector<std::string> sampleNamesMother;
        std::vector<std::string> sampleNamesFather;
        
        m_FatherVcf.GetSampleNames(sampleNamesFather);
        m_MotherVcf.GetSampleNames(sampleNamesMother);
        m_ChildVcf.GetSampleNames(sampleNamesChild);
        
        //Fill the MFC ids
        parentChildSampleIds = pedParser.GetIdsMFC(sampleNamesMother, sampleNamesFather, sampleNamesChild);
        if(parentChildSampleIds.size() != 3)
            return false;
        
        m_MotherVcf.SelectSample(parentChildSampleIds[0]);
        m_FatherVcf.SelectSample(parentChildSampleIds[1]);
        m_ChildVcf.SelectSample(parentChildSampleIds[2]);
    }
    
    else
    {
        //Set sample name of FATHER
        if (true == m_fatherChildConfig.m_bBaseSampleEnabled)
            m_FatherVcf.SelectSample(m_fatherChildConfig.m_pBaseSample);
        else
        {
            std::vector<std::string> sampleNames;
            m_FatherVcf.GetSampleNames(sampleNames);
            bIsSuccessFather = m_FatherVcf.SelectSample(sampleNames[0]);
        }
        
        if(!bIsSuccessFather)
        {
            std::cerr << "Father Sample name is incorrect!" << std::endl;
            return false;
        }
        
        //Set sample name of MOTHER
        if (true == m_motherChildConfig.m_bBaseSampleEnabled)
            m_MotherVcf.SelectSample(m_motherChildConfig.m_pBaseSample);
        else
        {
            std::vector<std::string> sampleNames;
            m_MotherVcf.GetSampleNames(sampleNames);
            bIsSuccessMother = m_MotherVcf.SelectSample(sampleNames[0]);
        }
        
        if(!bIsSuccessMother)
        {
            std::cerr << "Mother Sample name is incorrect!" << std::endl;
            return false;
        }
        
        //Set sample name of CHILD
        if (true == m_fatherChildConfig.m_bCalledSampleEnabled)
            m_ChildVcf.SelectSample(m_fatherChildConfig.m_pCalledSample);
        else
        {
            std::vector<std::string> sampleNames;
            m_ChildVcf.GetSampleNames(sampleNames);
            bIsSuccessChild = m_ChildVcf.SelectSample(sampleNames[0]);
        }
        
        if(!bIsSuccessChild)
        {
            std::cerr << "Child Sample name is incorrect!" << std::endl;
            return false;
        }
    }

    return bIsSuccessMother && bIsSuccessFather && bIsSuccessChild;
}

bool CMendelianVariantProvider::InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig)
{
    bool bIsSuccessVCFs = true;
    bool bIsSuccessFasta = true;
    
    m_motherChildConfig = a_rMotherChildConfig;
    m_fatherChildConfig = a_rFatherChildConfig;
    
    // OPEN VCF FILES
    bIsSuccessVCFs = OpenVcfFiles();

    if(!bIsSuccessVCFs)
        std::cerr << "VCF file(s) has error!" << std::endl;

    // OPEN FASTA FILE
    bIsSuccessFasta = m_referenceFasta.OpenFastaFile(a_rFatherChildConfig.m_pFastaFileName);

    if(!bIsSuccessFasta)
        std::cerr << "FASTA file is unable to open!: " << a_rFatherChildConfig.m_pFastaFileName << std::endl;
    
    if(bIsSuccessVCFs && bIsSuccessFasta)
    {
        //Fill the variants of 3 vcf file
        FillVariants();
        
        //Get the common chromosome ids and clear the uncommon variants from provider
        SetCommonChromosomes();
        
        //Fill the oriented variants of 3 vcf for genotype matching
        FillGenotypeMatchOrientedVariants(m_aCommonChromosomes);
        
        //Fill the oriented variants of 3 vcf for allele matching
        FillAlleleMatchOrientedVariants(m_aCommonChromosomes);
    }

    return bIsSuccessVCFs && bIsSuccessFasta;
}

void CMendelianVariantProvider::FillVariantForSample(int a_nSampleId, SConfig& a_rConfig)
{
    EMendelianVcfName sampleName = static_cast<EMendelianVcfName>(a_nSampleId);

    CSimpleBEDParser bedParser;
    unsigned int remainingBedContigCount = 0;
    
    if(true == a_rConfig.m_bInitializeFromBed)
    {
        bedParser.InitBEDFile(a_rConfig.m_pBedFileName);
        remainingBedContigCount = bedParser.m_nTotalContigCount;
    }
    
    int* pNonAssessedVariantCount;
    int* pAsteriskVariantCount;
    std::vector<std::vector<CVariant>>* pVariants;
    CVcfReader* pReader;
    std::string sampleNameStr;
    
    switch (sampleName)
    {
        case eFATHER:
            pNonAssessedVariantCount = &m_nFatherNotAssessedVariantCount;
            pVariants = &m_aFatherVariantList;
            pReader = &m_FatherVcf;
            pAsteriskVariantCount = &m_nFatherAsteriskCount;
            sampleNameStr = "father";
            break;
        case eMOTHER:
            pNonAssessedVariantCount = &m_nMotherNotAssessedVariantCount;
            pVariants = &m_aMotherVariantList;
            pReader = &m_MotherVcf;
            pAsteriskVariantCount = &m_nMotherAsteriskCount;
            sampleNameStr = "mother";
            break;
        case eCHILD:
            pNonAssessedVariantCount = &m_nChildNotAssessedVariantCount;
            pVariants = &m_aChildVariantList;
            pReader = &m_ChildVcf;
            pAsteriskVariantCount = &m_nChildAsteriskCount;
            sampleNameStr = "child";
            break;
            
        default:
            std::cerr << "Wrong Sample Enumeration in Mendelian Variant Provider!" << std::endl;
            break;
    }
    
    CVariant variant;
    int id = 0;
    std::string preChrId = "";
    unsigned int regionIterator = 0;
    
    std::vector<CVariant> multiTrimmableVarList;
    
    
    while(pReader->GetNextRecord(&variant, id, a_rConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            //We update the remaining contig count in BED file
            if(bedParser.m_regionMap[preChrId].size() > 0)
                remainingBedContigCount--;
            
            preChrId = variant.m_chrName;
            std::cerr << "Reading chromosome " << preChrId << " of Parent[" << sampleNameStr <<"] vcf" << std::endl;
            id = 0;
            variant.m_nId = id;
            
            regionIterator = 0;
        }
        
        if(true == a_rConfig.m_bInitializeFromBed)
        {
            //All BED regions are finished
            if(remainingBedContigCount == 0)
                break;
            
            //No Region exist for this chromosome
            if(regionIterator == bedParser.m_regionMap[variant.m_chrName].size())
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
            if(bedParser.m_regionMap[variant.m_chrName][regionIterator].m_nStartPos >= (variant.m_nOriginalPos + (int)variant.m_refSequence.length()))
                continue;
        }
            
        if(!variant.m_bIsNoCall && CUtils::IsHomRef(variant))
            continue;
        
        //Eliminate variants rather than diploid
        if(variant.m_nZygotCount != 2)
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            pAsteriskVariantCount++;
            continue;
        }
        
        else if(a_rConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            pNonAssessedVariantCount++;
        
        else if(CUtils::IsStructuralVariant(variant, a_rConfig.m_nMaxVariantSize))
            pNonAssessedVariantCount++;
        
        else if(true == variant.m_bHaveMultipleTrimOption)
        {
            multiTrimmableVarList.push_back(variant);
            id++;
        }
        
        else
        {
            (*pVariants)[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    FindOptimalTrimmings(multiTrimmableVarList, sampleName);
    AppendTrimmedVariants(multiTrimmableVarList, sampleName);
    
    for(unsigned int k = 0; k < pReader->GetContigs().size(); k++)
    {
        std::sort((*pVariants)[k].begin(), (*pVariants)[k].end(), CUtils::CompareVariants);
        (*pVariants)[k].shrink_to_fit();
    }
    
    (*pVariants).shrink_to_fit();
}

void CMendelianVariantProvider::FillVariants()
{
    //initialize variant lists
    m_aFatherVariantList = std::vector<std::vector<CVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherVariantList = std::vector<std::vector<CVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildVariantList = std::vector<std::vector<CVariant>>(m_ChildVcf.GetContigs().size());
    m_nFatherNotAssessedVariantCount = 0;
    m_nMotherNotAssessedVariantCount = 0;
    m_nChildNotAssessedVariantCount = 0;
    
    //TODO : this will be replaced ??
    m_nMotherAsteriskCount = 0;
    m_nFatherAsteriskCount = 0;
    m_nChildAsteriskCount  = 0;
    
    FillVariantForSample(eMOTHER, m_motherChildConfig);
    FillVariantForSample(eFATHER, m_fatherChildConfig);
    FillVariantForSample(eCHILD, m_motherChildConfig);
}


void CMendelianVariantProvider::FillGenotypeMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes)
{
    //INITIALIZE ORIENTED VARIANT LISTS
    m_aFatherOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_ChildVcf.GetContigs().size());
    
    for(unsigned int i=0; i < a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], true);
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], false);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid] = std::vector<core::COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aChildVariantList[a_aCommonChromosomes[i].m_nCid].size(); j++,k+=2)
        {
            m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid][k] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid][j], true);
            m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid][k+1] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid][j], false);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherOrientedVariantList[a_aCommonChromosomes[i].m_nMid] = std::vector<core::COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid].size(); j++, k+=2)
        {
            m_aMotherOrientedVariantList[a_aCommonChromosomes[i].m_nMid][k] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid][j], true);
            m_aMotherOrientedVariantList[a_aCommonChromosomes[i].m_nMid][k+1] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid][j], false);
        }
    }
}

void CMendelianVariantProvider::FillAlleleMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes)
{
    //INITIALIZE ORIENTED VARIANT LISTS
    m_aFatherAlleleMatchOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherAlleleMatchOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildAlleleMatchOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_ChildVcf.GetContigs().size());

    for(unsigned int i=0; i < a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aChildVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++,k+=2)
        {
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(unsigned int j=0, k=0; j < m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
    }
}

bool IsAutosome(const std::string& a_rChrName)
{
    std::stringstream convertor;
    int chrNumber;
    
    if(a_rChrName.length() > 3 && a_rChrName.substr(0,3) == "chr")
    {
        convertor << a_rChrName.substr(3);
        convertor >> chrNumber;
        
        if(convertor.fail())
            return false;
        else if(chrNumber > 0 && chrNumber < 23)
            return true;
    }
    
    else
    {
        convertor << a_rChrName;
        convertor >> chrNumber;
        
        if(convertor.fail())
            return false;
        else if(chrNumber > 0 && chrNumber < 23)
            return true;
    }
    
    return false;
}

void CMendelianVariantProvider::SetCommonChromosomes()
{
    int tripleIndex = 0;
    
    for(auto fatherItr = m_FatherVcf.m_chrIndexMap.begin(); fatherItr != m_FatherVcf.m_chrIndexMap.end(); fatherItr++)
    {
        bool isFound = false;
        
        for(auto motherItr = m_MotherVcf.m_chrIndexMap.begin(); motherItr != m_MotherVcf.m_chrIndexMap.end(); motherItr++)
        {
            for(auto childItr = m_ChildVcf.m_chrIndexMap.begin(); childItr != m_ChildVcf.m_chrIndexMap.end(); childItr++)
            {
                if(childItr->first == motherItr->first && childItr->first == fatherItr->first)
                {
                    if(m_aChildVariantList[childItr->second].size() > LEAST_VARIANT_THRESHOLD
                       &&
                       m_aFatherVariantList[fatherItr->second].size() > LEAST_VARIANT_THRESHOLD
                       &&
                       m_aMotherVariantList[motherItr->second].size() > LEAST_VARIANT_THRESHOLD)
                    {
                        if(m_motherChildConfig.m_bAutosomeOnly && !IsAutosome(motherItr->first))
                            continue;
                        
                        m_aCommonChromosomes.push_back(SChrIdTriplet(motherItr->second, fatherItr->second, childItr->second, motherItr->first, tripleIndex++));
                        isFound = true;
                        break;
                    }
                }
            }
            
            if(isFound == true)
                break;
        }
    }
        
    std::sort(m_aCommonChromosomes.begin(), m_aCommonChromosomes.end(), [](const SChrIdTriplet& t1, const SChrIdTriplet& t2){ return t1.m_nCid < t2.m_nCid; });
}

const std::vector<SChrIdTriplet>& CMendelianVariantProvider::GetCommonChromosomes() const
{
    return m_aCommonChromosomes;
}


std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList;
    
    switch (a_uFrom) {
        case eCHILD:
            for(unsigned int k = 0; k < m_aChildVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][k]);
            break;
        case eFATHER:
            for(unsigned int k = 0; k < m_aFatherVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][k]);
            break;
        case eMOTHER:
            for(unsigned int k = 0; k < m_aMotherVariantList[a_nChrNo].size();k++)
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
            for(unsigned int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eFATHER:
            for(unsigned int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eMOTHER:
            for(unsigned int k = 0; k < a_nIndexList.size(); k++)
                varList.push_back(&m_aMotherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
            
        default:
            break;
    }
    
    return varList;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetSortedVariantListByID(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList = GetVariantList(a_uFrom, a_nChrNo);
    std::sort(varList.begin(), varList.end(), [](const CVariant* pVar1, const CVariant* pVar2){return pVar1->m_nId < pVar2->m_nId;});
    return varList;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetSortedVariantListByIDandStartPos(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList = GetVariantList(a_uFrom, a_nChrNo);
    std::sort(varList.begin(), varList.end(), [](const CVariant* pVar1, const CVariant* pVar2)
    {
        if(pVar1->m_nOriginalPos != pVar2->m_nOriginalPos)
            return pVar1->m_nOriginalPos < pVar2->m_nOriginalPos;
        else if(pVar1->m_nStartPos != pVar2->m_nStartPos)
            return pVar1->m_nStartPos < pVar2->m_nStartPos;
        else if(pVar1->m_nEndPos != pVar2->m_nEndPos)
            return pVar1->m_nEndPos < pVar2->m_nEndPos;
        else
            return pVar1->m_nId < pVar2->m_nId;
    });
    return varList;
}


std::vector<const core::COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch) const
{
    std::vector<const core::COrientedVariant*> ovarList;
    const std::vector<std::vector<core::COrientedVariant>>* pBaseVarList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            pBaseVarList = a_bIsAlleleMatch ? &m_aChildAlleleMatchOrientedVariantList : &m_aChildOrientedVariantList;
            break;
        case eFATHER:
            pBaseVarList = a_bIsAlleleMatch ? &m_aFatherAlleleMatchOrientedVariantList : &m_aFatherOrientedVariantList;
            break;
        case eMOTHER:
            pBaseVarList = a_bIsAlleleMatch ? &m_aMotherAlleleMatchOrientedVariantList : &m_aMotherOrientedVariantList;
            break;
        default:
            break;
    }
    
    for(unsigned int k = 0; k < (*pBaseVarList)[a_nChrNo].size();k++)
        ovarList.push_back(&((*pBaseVarList)[a_nChrNo][k]));
    
    return ovarList;
}

std::vector<const core::COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch, const std::vector<int>& a_nIndexList) const
{
    std::vector<const core::COrientedVariant*> ovarList;
    const std::vector<std::vector<core::COrientedVariant>>* pBaseVarList;

    switch (a_uFrom)
    {
        case eCHILD:
            pBaseVarList = a_bIsAlleleMatch ? &m_aChildAlleleMatchOrientedVariantList : &m_aChildOrientedVariantList;
            break;
        case eFATHER:
            pBaseVarList = a_bIsAlleleMatch ? &m_aFatherAlleleMatchOrientedVariantList : &m_aFatherOrientedVariantList;
            break;
        case eMOTHER:
            pBaseVarList = a_bIsAlleleMatch ? &m_aMotherAlleleMatchOrientedVariantList : &m_aMotherOrientedVariantList;
            break;
        default:
            break;
    }
    
    for(unsigned int k = 0; k < a_nIndexList.size();k++)
    {
        ovarList.push_back(&((*pBaseVarList)[a_nChrNo][a_nIndexList[k]*2]));
        ovarList.push_back(&((*pBaseVarList)[a_nChrNo][a_nIndexList[k]*2+1]));
    }

    return ovarList;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_nIndexList) const
{
    std::vector<const CVariant*> resultList(a_nIndexList.size());
    for(unsigned int k = 0; k < a_nIndexList.size(); k++)
        resultList[k] = a_rVariantList[a_nIndexList[k]];
    
    return resultList;
}



int CMendelianVariantProvider:: GetVariantCount(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    if(a_uFrom == eMOTHER)
        return static_cast<int>(m_aMotherVariantList[a_nChrNo].size());
    if(a_uFrom == eFATHER)
        return static_cast<int>(m_aFatherVariantList[a_nChrNo].size());
    if(a_uFrom == eCHILD)
        return static_cast<int>(m_aChildVariantList[a_nChrNo].size());
    else
        return -1;
}


int CMendelianVariantProvider::GetSkippedVariantCount(EMendelianVcfName a_uFrom) const
{
    int totalCount = 0;
    
    const std::vector<std::vector<CVariant>>* pVariantList;
    switch (a_uFrom)
    {
        case eFATHER:
            pVariantList = &m_aFatherVariantList;
            break;
        case eMOTHER:
            pVariantList = &m_aMotherVariantList;
            break;
        case eCHILD:
            pVariantList = &m_aChildVariantList;
            break;
        default:
            break;
    }
    
    for(unsigned int k = 0; k < pVariantList->size(); k++)
    {
        for(CVariant var : (*pVariantList)[k])
        {
            if(var.m_variantStatus == eCOMPLEX_SKIPPED)
                totalCount++;
        }
    }
    
    return totalCount;
}


const std::vector<SVcfContig>& CMendelianVariantProvider::GetContigs() const
{
    return m_ChildVcf.GetContigs();
}

SVcfContig CMendelianVariantProvider::GetContig(const std::string& a_rChrName) const
{
    int contigId = m_ChildVcf.GetContigId(a_rChrName);
    
    if(contigId < 0)
    {
        std::cerr << "Unknown Chromosome : " << a_rChrName << std::endl;
        return SVcfContig();
    }
    
    return m_ChildVcf.GetContigs()[contigId];
}


int CMendelianVariantProvider::GetContigCount(EMendelianVcfName a_uFrom)
{
    switch (a_uFrom)
    {
        case eCHILD:
            return (int)m_aChildVariantList.size();
        case eFATHER:
            return (int)m_aFatherVariantList.size();
        case eMOTHER:
            return (int)m_aMotherVariantList.size();
        default:
            return -1;
    }
}


int CMendelianVariantProvider::GetNotAssessedVariantCount(EMendelianVcfName a_uFrom)
{
    unsigned int skippedVariantCount = 0;
    
    switch (a_uFrom) {
        case eMOTHER:
            skippedVariantCount = m_nMotherNotAssessedVariantCount;
            break;
        case eFATHER:
            skippedVariantCount = m_nFatherNotAssessedVariantCount;
            break;
        case eCHILD:
            skippedVariantCount = m_nChildNotAssessedVariantCount;
            break;
        default:
            break;
    }
    
    return static_cast<int>(skippedVariantCount);
}

void CMendelianVariantProvider::FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EMendelianVcfName a_uFrom)
{
    if(a_rVariantList.size() == 0)
        return;
    
    std::vector<std::vector<CVariant>>* allVarList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            allVarList = &m_aChildVariantList;
            break;
        case eFATHER:
            allVarList = &m_aFatherVariantList;
            break;
        case eMOTHER:
            allVarList = &m_aMotherVariantList;
            break;
        default:
            allVarList = 0;
            break;
    }
    
    CBaseVariantProvider::FindOptimalTrimmings(a_rVariantList, allVarList, m_motherChildConfig);
}


void CMendelianVariantProvider::AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EMendelianVcfName a_uFrom)
{
    std::vector<std::vector<CVariant>>* variantList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            variantList = &m_aChildVariantList;
            break;
        case eFATHER:
            variantList = &m_aFatherVariantList;
            break;
        case eMOTHER:
            variantList = &m_aMotherVariantList;
            break;
        default:
            variantList = 0;
            break;
    }
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
        (*variantList)[a_rVariantList[k].m_nChrId].push_back(a_rVariantList[k]);
}
