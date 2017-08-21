//
//  CMendelianVariantProvider.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include <algorithm>
#include "CMendelianVariantProvider.h"
#include "CSimplePEDParser.h"
#include "CSimpleBEDParser.h"
#include <iostream>
#include <sstream>

using namespace mendelian;

bool CMendelianVariantProvider::InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig)
{
    bool bIsSuccessFather = true;
    bool bIsSuccessMother = true;
    bool bIsSuccessChild = true;
    bool bIsSuccessFasta = true;
    
    m_motherChildConfig = a_rMotherChildConfig;
    m_fatherChildConfig = a_rFatherChildConfig;
    
    //This vector will be used if the program initialized with a PED file
    CSimplePEDParser pedParser;
    std::vector<std::string> parentChildSampleIds;
    
    //Open FATHER vcf file
    bIsSuccessFather = m_FatherVcf.Open(a_rFatherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessFather)
        std::cerr << "Father VCF file is unable to open!: " << a_rFatherChildConfig.m_pBaseVcfFileName << std::endl;
   
    //Open MOTHER vcf file
    bIsSuccessMother = m_MotherVcf.Open(a_rMotherChildConfig.m_pBaseVcfFileName);
    if(!bIsSuccessMother)
        std::cerr << "Mother VCF file is unable to open!: " << a_rMotherChildConfig.m_pBaseVcfFileName << std::endl;
    
    //Open CHILD vcf file
    bIsSuccessChild = m_ChildVcf.Open(a_rFatherChildConfig.m_pCalledVcfFileName);
    if(!bIsSuccessChild)
        std::cerr << "Child VCF file is unable to open!: " << a_rFatherChildConfig.m_pCalledVcfFileName << std::endl;
    
    
    if(!bIsSuccessChild || !bIsSuccessFather || !bIsSuccessMother)
        return false;
    
    
    //Set Sample name from PED file
    if(true == m_motherChildConfig.m_bInitializeFromPED)
    {
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

    // OPEN FASTA FILE
    bIsSuccessFasta = m_fastaParser.OpenFastaFile(a_rFatherChildConfig.m_pFastaFileName);
    if(!bIsSuccessFasta)
        std::cerr << "FASTA file is unable to open!: " << a_rFatherChildConfig.m_pFastaFileName << std::endl;
    
    if(bIsSuccessChild && bIsSuccessMother && bIsSuccessFather && bIsSuccessFasta)
    {
        //Fill the variants of 3 vcf file
        if(m_motherChildConfig.m_bInitializeFromBed)
            FillVariantsFromBED();
        else
            FillVariants();
        
        //Get the common chromosome ids and clear the uncommon variants from provider
        SetCommonChromosomes();
        
        //Filter NonAssessed variants coming from dependent sites. Update the variant ids
        //RemoveNonAssessedSites();

        //Fill the oriented variants of 3 vcf for genotype matching
        FillGenotypeMatchOrientedVariants(m_aCommonChromosomes);
        
        //Fill the oriented variants of 3 vcf for allele matching
        FillAlleleMatchOrientedVariants(m_aCommonChromosomes);
    }

    return bIsSuccessChild && bIsSuccessMother && bIsSuccessFather && bIsSuccessFasta;
}


bool CMendelianVariantProvider::CompareVariants(const CVariant& var1, const CVariant& var2)
{
    if(var1.m_nStartPos != var2.m_nStartPos)
        return var1.m_nStartPos < var2.m_nStartPos;
    else
        return var1.m_nEndPos < var2.m_nEndPos;
}

void CMendelianVariantProvider::FillVariantsFromBED()
{
    CSimpleBEDParser bedParser;
    bedParser.InitBEDFile(m_motherChildConfig.m_pBedFileName);

    SBedRegion bedRegion;
    bool hasNextRegion = bedParser.GetNextRegion(bedRegion);
    
    //There is no region in BED file
    if(false == hasNextRegion)
        return;
    
    CVariant variant;
    int id = 0;
    std::string preChrId = "";
    
    //initialize variant lists
    m_aFatherVariantList = std::vector<std::vector<CVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherVariantList = std::vector<std::vector<CVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildVariantList = std::vector<std::vector<CVariant>>(m_ChildVcf.GetContigs().size());
    m_aFatherNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_ChildVcf.GetContigs().size());
    
    //TODO : this will be replaced ??
    m_nMotherAsteriskCount = 0;
    m_nFatherAsteriskCount = 0;
    m_nChildAsteriskCount  = 0;
    
    //READ VARIANTS OF FATHER
    while(m_FatherVcf.GetNextRecord(&variant, id, m_fatherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Reading chromosome " << preChrId << " of Parent[FATHER] vcf" << std::endl;
        }

        //Pass to the next region
        if((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
            ||
            m_FatherVcf.m_chrIndexMap[variant.m_chrName] > m_FatherVcf.m_chrIndexMap[bedRegion.m_chrName])
        {
            hasNextRegion = bedParser.GetNextRegion(bedRegion);
            if(false == hasNextRegion)
                break;
        }
        
        //Variant Could not pass from BED region
        if(bedRegion.m_chrName != variant.m_chrName || bedRegion.m_nStartPos > variant.m_nStartPos || bedRegion.m_nEndPos < variant.m_nEndPos)
            continue;
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nFatherAsteriskCount++;
            continue;
        }
        
        else if(m_fatherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aFatherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_fatherChildConfig.m_nMaxVariantSize))
            m_aFatherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else
        {
            m_aFatherVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    preChrId = "";
    id = 0;
    bedParser.ResetIterator();
    bedParser.GetNextRegion(bedRegion);
    
    
    for(int k = 0; k < (int)m_FatherVcf.GetContigs().size(); k++)
    {
        std::sort(m_aFatherNotAssessedVariantList[k].begin(), m_aFatherNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aFatherVariantList[k].begin(), m_aFatherVariantList[k].end(), CompareVariants);
    }
    
    //READ VARIANTS OF MOTHER
    while(m_MotherVcf.GetNextRecord(&variant, id, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Reading chromosome " << preChrId << " of Parent[MOTHER] vcf" << std::endl;
        }
        
        //Pass to the next region
        if((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
           ||
           m_MotherVcf.m_chrIndexMap[variant.m_chrName] > m_MotherVcf.m_chrIndexMap[bedRegion.m_chrName])
        {
            hasNextRegion = bedParser.GetNextRegion(bedRegion);
            if(false == hasNextRegion)
                break;
        }
        
        //Variant Could not pass from BED region
        if(bedRegion.m_chrName != variant.m_chrName || bedRegion.m_nStartPos > variant.m_nStartPos || bedRegion.m_nEndPos < variant.m_nEndPos)
            continue;
        
        if(!variant.m_bIsNoCall &&  IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nMotherAsteriskCount++;
            continue;
        }
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aMotherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            m_aMotherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else
        {
            m_aMotherVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    preChrId = "";
    id = 0;
    bedParser.ResetIterator();
    bedParser.GetNextRegion(bedRegion);
    
    
    for(int k = 0; k < (int)m_MotherVcf.GetContigs().size(); k++)
    {
        std::sort(m_aMotherNotAssessedVariantList[k].begin(), m_aMotherNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aMotherVariantList[k].begin(), m_aMotherVariantList[k].end(), CompareVariants);
    }
    
    //READ VARIANTS OF CHILD
    while(m_ChildVcf.GetNextRecord(&variant, id, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Reading chromosome " << preChrId << " of child vcf" << std::endl;
        }

        //Pass to the next region
        if((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
           ||
           m_ChildVcf.m_chrIndexMap[variant.m_chrName] > m_ChildVcf.m_chrIndexMap[bedRegion.m_chrName])
        {
            hasNextRegion = bedParser.GetNextRegion(bedRegion);
            if(false == hasNextRegion)
                break;
        }
        
        //Variant Could not pass from BED region
        if(bedRegion.m_chrName != variant.m_chrName || bedRegion.m_nStartPos > variant.m_nStartPos || bedRegion.m_nEndPos < variant.m_nEndPos)
            continue;

        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nChildAsteriskCount++;
            continue;
        }
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aChildNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            m_aChildNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else
        {
            m_aChildVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    for(int k = 0; k < (int)m_ChildVcf.GetContigs().size(); k++)
    {
        std::sort(m_aChildNotAssessedVariantList[k].begin(), m_aChildNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aChildVariantList[k].begin(), m_aChildVariantList[k].end(), CompareVariants);
    }
    
   
    

}


void CMendelianVariantProvider::FillVariants()
{
    
    CVariant variant;
    int id = 0;
    std::string preChrId = "";

    //initialize variant lists
    m_aFatherVariantList = std::vector<std::vector<CVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherVariantList = std::vector<std::vector<CVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildVariantList = std::vector<std::vector<CVariant>>(m_ChildVcf.GetContigs().size());
    m_aFatherNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildNotAssessedVariantList = std::vector<std::vector<CVariant>>(m_ChildVcf.GetContigs().size());
    
    //TODO : this will be replaced ??
    m_nMotherAsteriskCount = 0;
    m_nFatherAsteriskCount = 0;
    m_nChildAsteriskCount  = 0;
    
    //READ VARIANTS OF FATHER
    while(m_FatherVcf.GetNextRecord(&variant, id, m_fatherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of Parent[FATHER] vcf" << std::endl;
        }
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nFatherAsteriskCount++;
            continue;
        }
        
        else if(m_fatherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aFatherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_fatherChildConfig.m_nMaxVariantSize))
            m_aFatherNotAssessedVariantList[variant.m_nChrId].push_back(variant);

        else
        {
            m_aFatherVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    preChrId = "";
    id = 0;

    
    for(int k = 0; k < (int)m_FatherVcf.GetContigs().size(); k++)
    {
        std::sort(m_aFatherNotAssessedVariantList[k].begin(), m_aFatherNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aFatherVariantList[k].begin(), m_aFatherVariantList[k].end(), CompareVariants);
    }

    //READ VARIANTS OF MOTHER
    while(m_MotherVcf.GetNextRecord(&variant, id, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of Parent[MOTHER] vcf" << std::endl;
        }
        
        if(!variant.m_bIsNoCall &&  IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nMotherAsteriskCount++;
            continue;
        }
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aMotherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            m_aMotherNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else
        {
            m_aMotherVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    preChrId = "";
    id = 0;

    
    for(int k = 0; k < (int)m_MotherVcf.GetContigs().size(); k++)
    {
        std::sort(m_aMotherNotAssessedVariantList[k].begin(), m_aMotherNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aMotherVariantList[k].begin(), m_aMotherVariantList[k].end(), CompareVariants);
    }
    
    //READ VARIANTS OF CHILD
    while(m_ChildVcf.GetNextRecord(&variant, id, m_motherChildConfig))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cout << "Processing chromosome " << preChrId << " of child vcf" << std::endl;
        }
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        std::size_t found = variant.m_allelesStr.find('*');
        if (found!=std::string::npos)
        {
            m_nChildAsteriskCount++;
            continue;
        }
        
        else if(m_motherChildConfig.m_bIsFilterEnabled && variant.m_bIsFilterPASS == false)
            m_aChildNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else if(IsStructuralVariant(variant, m_motherChildConfig.m_nMaxVariantSize))
            m_aChildNotAssessedVariantList[variant.m_nChrId].push_back(variant);
        
        else
        {
            m_aChildVariantList[variant.m_nChrId].push_back(variant);
            id++;
        }
    }
    
    for(int k = 0; k < (int)m_ChildVcf.GetContigs().size(); k++)
    {
        std::sort(m_aChildNotAssessedVariantList[k].begin(), m_aChildNotAssessedVariantList[k].end(), CompareVariants);
        std::sort(m_aChildVariantList[k].begin(), m_aChildVariantList[k].end(), CompareVariants);
    }

}


void CMendelianVariantProvider::FillGenotypeMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes)
{
    //INITIALIZE ORIENTED VARIANT LISTS
    m_aFatherOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_FatherVcf.GetContigs().size());
    m_aMotherOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_MotherVcf.GetContigs().size());
    m_aChildOrientedVariantList = std::vector<std::vector<core::COrientedVariant>>(m_ChildVcf.GetContigs().size());
    
    for(int i=0; i < (int)a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(int j=0, k=0; j < (int)m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], true);
            m_aFatherOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], false);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid] = std::vector<core::COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid].size() * 2);
        for(int j=0, k=0; j < (int)m_aChildVariantList[a_aCommonChromosomes[i].m_nCid].size(); j++,k+=2)
        {
            m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid][k] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid][j], true);
            m_aChildOrientedVariantList[a_aCommonChromosomes[i].m_nCid][k+1] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nCid][j], false);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherOrientedVariantList[a_aCommonChromosomes[i].m_nMid] = std::vector<core::COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid].size() * 2);
        for(int j=0, k=0; j < (int)m_aMotherVariantList[a_aCommonChromosomes[i].m_nMid].size(); j++, k+=2)
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

    for(int i=0; i < (int)a_aCommonChromosomes.size(); i++)
    {
        //GENERATE FATHER ORIENTED VARS
        m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(int j=0, k=0; j < (int)m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aFatherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aFatherVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
        
        //GENERATE CHILD ORIENTED VARS
        m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(int j=0, k=0; j < (int)m_aChildVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++,k+=2)
        {
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aChildAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aChildVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
        
        //GENERATE MOTHER ORIENTED VARS
        m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid] = std::vector<core::COrientedVariant>(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid].size() * 2);
        for(int j=0, k=0; j < (int)m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid].size(); j++, k+=2)
        {
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid][j], 0);
            m_aMotherAlleleMatchOrientedVariantList[a_aCommonChromosomes[i].m_nFid][k+1] = core::COrientedVariant(m_aMotherVariantList[a_aCommonChromosomes[i].m_nFid][j], 1);
        }
    }
}

bool CMendelianVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
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

bool CMendelianVariantProvider::IsHomRef(const CVariant& a_rVariant) const
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
    
    //Clear redundant variants TODO: We can remove redundant variants. Following part of code should be arranged to eliminate unique chromosomes
/*    for(int k = 0; k < m_aCommonChromosomes.size(); k++)
    {
        std::cerr << "Warning! Chromosome " << m_aCommonChromosomes[k].m_chrName << " is not contained by all three vcf files. Variants will be filtered out from comparison" << std::endl;
        m_aChildVariantList[m_aCommonChromosomes[k].m_nCid].clear();
        m_aFatherVariantList[m_aCommonChromosomes[k].m_nFid].clear();
        m_aMotherVariantList[m_aCommonChromosomes[k].m_nMid].clear();
    }
 */
    
    
    
}

std::vector<SChrIdTriplet>& CMendelianVariantProvider::GetCommonChromosomes()
{
    return m_aCommonChromosomes;
}


std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList;
    
    switch (a_uFrom) {
        case eCHILD:
            for(int k = 0; k < (int)m_aChildVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][k]);
            break;
        case eFATHER:
            for(int k = 0; k < (int)m_aFatherVariantList[a_nChrNo].size();k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][k]);
            break;
        case eMOTHER:
            for(int k = 0; k < (int)m_aMotherVariantList[a_nChrNo].size();k++)
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
            for(int k = 0; k < (int)a_nIndexList.size(); k++)
                varList.push_back(&m_aChildVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eFATHER:
            for(int k = 0; k < (int)a_nIndexList.size(); k++)
                varList.push_back(&m_aFatherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
        case eMOTHER:
            for(int k = 0; k < (int)a_nIndexList.size(); k++)
                varList.push_back(&m_aMotherVariantList[a_nChrNo][a_nIndexList[k]]);
            break;
            
        default:
            break;
    }
    
    return varList;
}

std::vector<const CVariant*> CMendelianVariantProvider::GetSortedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    std::vector<const CVariant*> varList = GetVariantList(a_uFrom, a_nChrNo);
    std::sort(varList.begin(), varList.end(), [](const CVariant* pVar1, const CVariant* pVar2){return pVar1->m_nId < pVar2->m_nId;});
    return varList;
}

std::vector<const core::COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch) const
{
    std::vector<const core::COrientedVariant*> ovarList;
    
    if(false == a_bIsAlleleMatch)
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < (int)m_aChildOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][k]);
                break;
            case eFATHER:
                for(int k = 0; k < (int)m_aFatherOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][k]);
                break;
            case eMOTHER:
                for(int k = 0; k < (int)m_aMotherOrientedVariantList[a_nChrNo].size();k++)
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
                for(int k = 0; k < (int)m_aChildAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
            case eFATHER:
                for(int k = 0; k < (int)m_aFatherAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
            case eMOTHER:
                for(int k = 0; k < (int)m_aMotherAlleleMatchOrientedVariantList[a_nChrNo].size();k++)
                    ovarList.push_back(&m_aMotherAlleleMatchOrientedVariantList[a_nChrNo][k]);
                break;
                
            default:
                break;
        }
    }
    
    return ovarList;
}

std::vector<const core::COrientedVariant*> CMendelianVariantProvider::GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch, const std::vector<int>& a_nIndexList) const
{
    std::vector<const core::COrientedVariant*> ovarList;

    if(false == a_bIsAlleleMatch)
    {
        switch (a_uFrom)
        {
            case eCHILD:
                for(int k = 0; k < (int)a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aChildOrientedVariantList[a_nChrNo][a_nIndexList[k]*2+1]);
                }

                break;
            case eFATHER:
                for(int k = 0; k < (int)a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aFatherOrientedVariantList[a_nChrNo][a_nIndexList[k]*2+1]);
                }
                break;
            case eMOTHER:
                for(int k = 0; k < (int)a_nIndexList.size();k++)
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
                for(int k = 0; k < (int)a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aChildAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2 +1]);
                }
                break;
            case eFATHER:
                for(int k = 0; k < (int)a_nIndexList.size();k++)
                {
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2]);
                    ovarList.push_back(&m_aFatherAlleleMatchOrientedVariantList[a_nChrNo][a_nIndexList[k]*2 +1]);
                }
                break;
            case eMOTHER:
                for(int k = 0; k < (int)a_nIndexList.size();k++)
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

void CMendelianVariantProvider::ReadContig(std::string a_chrId, SContig& a_rContig)
{
    m_fastaParser.FetchNewChromosome(a_chrId, a_rContig);
}

void CMendelianVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
    {
        if(pVar->m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pVar->m_variantStatus != eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pVar-> m_variantStatus == eNOT_ASSESSED)
            pVar->m_variantStatus = a_status;
        else
            continue;
    }
}


void CMendelianVariantProvider::SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const core::COrientedVariant* pOVar : a_rVariantList)
    {
        if(pOVar->GetVariant().m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pOVar->GetVariant().m_variantStatus != eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pOVar->GetVariant().m_variantStatus == eNOT_ASSESSED)
            pOVar->GetVariant().m_variantStatus = a_status;
        else
            continue;
    }
}


std::vector<const CVariant*> CMendelianVariantProvider::GetVariantList(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_nIndexList) const
{
    std::vector<const CVariant*> resultList(a_nIndexList.size());
    for(int k = 0; k < (int)a_nIndexList.size(); k++)
        resultList[k] = a_rVariantList[a_nIndexList[k]];
    
    return resultList;
}



int CMendelianVariantProvider:: GetVariantCount(EMendelianVcfName a_uFrom, int a_nChrNo) const
{
    if(a_uFrom == eMOTHER)
        return static_cast<int>(m_aMotherVariantList[a_nChrNo].size());
    else if(a_uFrom == eFATHER)
        return static_cast<int>(m_aFatherVariantList[a_nChrNo].size());
    else if(a_uFrom == eCHILD)
        return static_cast<int>(m_aChildVariantList[a_nChrNo].size());
    else
        return -1;
}

int CMendelianVariantProvider::Get0BasedVariantIndex(EMendelianVcfName a_uFrom, int a_nChr, int a_nVariantId) const
{
    int variantCountSoFar = 0;
    
    for(int k = 0; k < a_nChr; k++)
        variantCountSoFar += GetVariantCount(a_uFrom, k);
    
    
    return a_nVariantId - variantCountSoFar;
}


int CMendelianVariantProvider::GetSkippedVariantCount(EMendelianVcfName a_uFrom) const
{
    int totalCount = 0;
    
    
    if(a_uFrom == eFATHER)
    {
        for(int k = 0; k < (int)m_aFatherVariantList.size(); k++)
        {
            for(CVariant var : m_aFatherVariantList[k])
            {
                if(var.m_variantStatus == eCOMPLEX_SKIPPED)
                    totalCount++;
            }
        }
    }

    else if(a_uFrom == eMOTHER)
    {
        for(int k = 0; k < (int)m_aMotherVariantList.size(); k++)
        {
            for(CVariant var : m_aMotherVariantList[k])
            {
                if(var.m_variantStatus == eCOMPLEX_SKIPPED)
                    totalCount++;
            }
        }
    }
    
    
    else if(a_uFrom == eCHILD)
    {
        for(int k = 0; k < (int)m_aChildVariantList.size(); k++)
        {
            for(CVariant var : m_aChildVariantList[k])
            {
                if(var.m_variantStatus == eCOMPLEX_SKIPPED)
                    totalCount++;
            }
        }
    }
    
    return totalCount;
}


const std::vector<SVcfContig>& CMendelianVariantProvider::GetContigs()
{
    return m_ChildVcf.GetContigs();
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
            for(int k = 0; k < (int)m_aMotherNotAssessedVariantList.size(); k++)
                skippedVariantCount += m_aMotherNotAssessedVariantList[k].size();
            skippedVariantCount += m_nMotherAsteriskCount;
            break;
        case eFATHER:
            for(int k = 0; k < (int)m_aFatherNotAssessedVariantList.size(); k++)
                skippedVariantCount += m_aFatherNotAssessedVariantList[k].size();
            skippedVariantCount += m_nFatherAsteriskCount;
            break;
        case eCHILD:
            for(int k = 0; k < (int)m_aChildNotAssessedVariantList.size(); k++)
                skippedVariantCount += m_aChildNotAssessedVariantList[k].size();
            skippedVariantCount += m_nChildAsteriskCount;
            break;
        default:
            break;
    }
    
    return static_cast<int>(skippedVariantCount);
}


void CMendelianVariantProvider::RemoveNonAssessedSites()
{
    int motherId = 0;
    int fatherId = 0;
    int childId = 0;
    
    for(int k = 0; k < static_cast<int>(m_aCommonChromosomes.size()); k++)
    {
        SChrIdTriplet triplet = m_aCommonChromosomes[k];
        
        std::vector<int> mergedSorted_NA_VariantPositions;
        
        for(CVariant var : m_aFatherNotAssessedVariantList[triplet.m_nFid])
            mergedSorted_NA_VariantPositions.push_back(var.m_nOriginalPos);
        for(CVariant var : m_aMotherNotAssessedVariantList[triplet.m_nMid])
            mergedSorted_NA_VariantPositions.push_back(var.m_nOriginalPos);
        for(CVariant var : m_aChildNotAssessedVariantList[triplet.m_nCid])
            mergedSorted_NA_VariantPositions.push_back(var.m_nOriginalPos);

        //Sort positions and remove duplicates
        std::sort(mergedSorted_NA_VariantPositions.begin(), mergedSorted_NA_VariantPositions.end());
        std::vector<int>::iterator it;
        it = std::unique (mergedSorted_NA_VariantPositions.begin(), mergedSorted_NA_VariantPositions.end());
        mergedSorted_NA_VariantPositions.resize(std::distance(mergedSorted_NA_VariantPositions.begin(),it));
        
        int naItr = 0;
        
        //Remove dependent non-assessed sites from father side
        for(auto varItr = m_aFatherVariantList[triplet.m_nFid].begin(); varItr != m_aFatherVariantList[triplet.m_nFid].end();)
        {
            while(naItr != static_cast<int>(mergedSorted_NA_VariantPositions.size()) && mergedSorted_NA_VariantPositions[naItr] < varItr->m_nOriginalPos)
                naItr++;
            
            if(varItr->m_nOriginalPos == mergedSorted_NA_VariantPositions[naItr])
                m_aFatherVariantList[triplet.m_nFid].erase(varItr);
            
            else
            {
                varItr->m_nId = fatherId++;
                varItr++;
            }
        }
        
        naItr = 0;
        //Remove dependent non-assessed sites from mother side
        for(auto varItr = m_aMotherVariantList[triplet.m_nMid].begin(); varItr != m_aMotherVariantList[triplet.m_nMid].end();)
        {
            while(naItr != static_cast<int>(mergedSorted_NA_VariantPositions.size()) && mergedSorted_NA_VariantPositions[naItr] < varItr->m_nOriginalPos)
                naItr++;
            
            if(varItr->m_nOriginalPos == mergedSorted_NA_VariantPositions[naItr])
                m_aMotherVariantList[triplet.m_nMid].erase(varItr);
            
            else
            {
                varItr->m_nId = motherId++;
                varItr++;
            }
        }
        
        naItr = 0;
        //Remove dependent non-assessed sites from child side
        for(auto varItr = m_aChildVariantList[triplet.m_nCid].begin(); varItr != m_aChildVariantList[triplet.m_nCid].end();)
        {
            while(naItr != static_cast<int>(mergedSorted_NA_VariantPositions.size()) && mergedSorted_NA_VariantPositions[naItr] < varItr->m_nOriginalPos)
                naItr++;
            
            if(varItr->m_nOriginalPos == mergedSorted_NA_VariantPositions[naItr])
                m_aChildVariantList[triplet.m_nCid].erase(varItr);
            
            else
            {
                varItr->m_nId = childId++;
                varItr++;
            }
        }
        
    }
    
    

}

























