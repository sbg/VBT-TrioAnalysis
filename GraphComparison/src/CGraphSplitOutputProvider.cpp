//
//  CGraphSplitOutputProvider.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/15/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CGraphSplitOutputProvider.h"
#include "CGraphVariantProvider.h"
#include "CSimpleBEDParser.h"
#include "CGraphVcfReader.h"
#include "CVcfWriter.h"

#include <algorithm>
#include <iostream>

using namespace graphcomparison;

void CGraphSplitOutputProvider::SetVcfNames(const std::string& a_rBaseVcfPath, const std::string& a_rCalledVcfPath)
{
    m_baseVcfPath = a_rBaseVcfPath;
    m_calledVcfPath = a_rCalledVcfPath;

}

void CGraphSplitOutputProvider::SetExcludeIndexes(std::unordered_map<std::string, SGraphVarContainer>& a_rExcludedIndexContainers)
{
    m_excludedIndexContainers = a_rExcludedIndexContainers;
}

void CGraphSplitOutputProvider::SetOutputVcfFolder(const std::string& a_rVcfPath)
{
    m_vcfsFolder = a_rVcfPath;
}

void CGraphSplitOutputProvider::SetParameters(const std::string& a_rBedPath, bool a_bPassFilterEnabled, int a_nMaxBasePairLength)
{
    m_bedFilePath = a_rBedPath;
    if(m_bedFilePath != "")
        m_bIsBEDFileEnabled = true;
    m_bPassFilterEnabled = a_bPassFilterEnabled;
    m_nMaxBasePairLength = a_nMaxBasePairLength;
}

void CGraphSplitOutputProvider::GenerateVcfs(const std::string& a_rIncludeName, const std::string& a_rExcludeName, bool a_bIsBase)
{
    CVcfWriter includeWriter;
    CVcfWriter excludeWriter;
    CGraphVcfReader inputReader;
    
    //Iterator of excluded variant container
    int excludeIndexItr = 0;
    int includeIndexItr = 0;
    
    //Open vcf input and output vcf files
    includeWriter.CreateVcf(std::string(m_vcfsFolder + a_rIncludeName).c_str());
    excludeWriter.CreateVcf(std::string(m_vcfsFolder + a_rExcludeName).c_str());
    inputReader.Open(a_bIsBase ? m_baseVcfPath.c_str() : m_calledVcfPath.c_str());
    
    //Write Header of output vcf files
    bcf_hdr_t* inputHeader = inputReader.GetHeaderPointer();
    includeWriter.SetRawHeader(inputHeader);
    includeWriter.WriteHeaderToVcf();
    excludeWriter.SetRawHeader(inputHeader);
    excludeWriter.WriteHeaderToVcf();
    
    //BED region parameters if there is a defined bed file
    CSimpleBEDParser bedParser;
    unsigned int regionIterator;
    
    if(m_bIsBEDFileEnabled)
    {
        bedParser.InitBEDFile(m_bedFilePath);
    }
    
    CVariant variant;
    bcf1_t* pRecord = inputReader.GetRecordPointer();
    
    int id = 0;
    std::string preChrId = "";
    
    while(inputReader.GetNextRecord(&variant, id, true))
    {
        if(preChrId != variant.m_chrName)
        {
            preChrId = variant.m_chrName;
            std::cerr << "Writing chromosome " << preChrId << " of " << (a_bIsBase ? "base" : "called") << " vcf..." << std::endl;
            variant.m_nId = 0;
            id = 0;
            excludeIndexItr = 0;
            regionIterator++;
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
            if(a_bIsBase)
            {
                if(m_excludedIndexContainers[variant.m_chrName].baseExcludedIndexes[excludeIndexItr] == id)
                {
                    excludeWriter.AddRawRecord(pRecord);
                    excludeIndexItr++;
                }
                else if(m_excludedIndexContainers[variant.m_chrName].baseIncludedIndexes[includeIndexItr] == id)
                {
                    includeWriter.AddRawRecord(pRecord);
                    includeIndexItr++;
                }
            }
            else
            {
                if(m_excludedIndexContainers[variant.m_chrName].calledExcludedIndexes[excludeIndexItr] == id)
                {
                    excludeWriter.AddRawRecord(pRecord);
                    excludeIndexItr++;
                }
                else if(m_excludedIndexContainers[variant.m_chrName].calledIncludedIndexes[includeIndexItr] == id)
                {
                    includeWriter.AddRawRecord(pRecord);
                    includeIndexItr++;
                }
            }
            id++;
        }
    }
   
    
    
    
    
}



