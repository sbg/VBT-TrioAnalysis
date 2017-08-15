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
        }
        
        //Pass to the next region
        if(m_bIsBEDFileEnabled &&
           ((bedRegion.m_chrName == variant.m_chrName && variant.m_nStartPos > bedRegion.m_nEndPos)
            ||
            inputReader.m_chrIndexMap[variant.m_chrName] > inputReader.m_chrIndexMap[bedRegion.m_chrName]))
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
            if(a_bIsBase)
            {
                if(m_excludedIndexContainers[variant.m_chrName].baseExcludedIndexes[excludeIndexItr] == id)
                {
                    excludeWriter.AddRawRecord(pRecord);
                    excludeIndexItr++;
                }
                else
                    includeWriter.AddRawRecord(pRecord);
            }
            else
            {
                if(m_excludedIndexContainers[variant.m_chrName].calledExcludedIndexes[excludeIndexItr] == id)
                {
                    excludeWriter.AddRawRecord(pRecord);
                    excludeIndexItr++;
                }
                else
                    includeWriter.AddRawRecord(pRecord);
            }
            id++;
        }
    }
   
    
    
    
    
}



