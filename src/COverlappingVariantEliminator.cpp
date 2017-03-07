//
//  COverlapingVariantEliminator.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 3/2/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "COverlappingVariantEliminator.h"
#include <iostream>


void COverlappingVariantEliminator::FilterOverlaps(const std::string& a_rFileName, bool a_bIsFilterOverlap, bool a_bIsFilter00)
{

    std::string modeSubPath = "";
    if(a_bIsFilterOverlap)
        modeSubPath = modeSubPath + "_Overlap";
    if(a_bIsFilter00)
        modeSubPath = modeSubPath + "_00";
    
    modeSubPath = modeSubPath + ".vcf";

    //Set output filename
    m_fileName =  a_rFileName.substr(0,a_rFileName.size() - 4) + modeSubPath;

    
    //OPEN VCF FILE TO WRITE
    m_pHtsFileFiltered = hts_open(m_fileName.c_str() ,"w");
    
    //Open original vcf file to read
    m_vcfReader.Open(a_rFileName.c_str());
    
    //Copy the header file
    bcf_hdr_write(m_pHtsFileFiltered, m_vcfReader.GetHeaderPointer());

    int lastProcessedChrId = 1;
    
    int SkippedLineCount = 0;
    
    //Ending of the last variant
    int lastBound = 0;
    int curLast = 0;
    
    int sampleCount = m_vcfReader.GetNumberOfSamples();
    CVariant* variant = new CVariant[sampleCount];
    
    while(m_vcfReader.GetNextRecordMultiSample(variant))
    {
        //Eliminate 0/0 child variants
        if(variant[2].m_genotype[0] == 0 && variant[2].m_genotype[1] == 0)
            continue;
        
        //If chromosome changes, reset the last bound
        if(lastProcessedChrId != variant[0].m_nChrId)
        {
            lastBound = 0;
            lastProcessedChrId++;
            
            if(lastProcessedChrId == -1)
                break;
        }
        
        //Decide if that variant should be skipped
        bool shouldSkip = false;
        curLast = 0;
        for(int k = 0; k < sampleCount; k++)
        {
            if(variant[k].m_nStartPos <= lastBound)
            {
                shouldSkip = true;
                break;
            }
            
            curLast = std::max(curLast, variant[k].m_nEndPos);
        }
        
        //Write the variant into the new vcf file
        if(!shouldSkip)
        {
            bcf_write1(m_pHtsFileFiltered, m_vcfReader.GetHeaderPointer(), m_vcfReader.GetRecordPointer());
            lastBound = curLast;
        }
        else
            SkippedLineCount++;
    }
    
    //CLOSE VCF FILE
    int ret= hts_close(m_pHtsFileFiltered);
    if(ret != 0)
    {
        std::cerr << "A problem occured with saving the VCF file." << std::endl;
    }

    std::cout << "Skipped Line Count:" << SkippedLineCount << std::endl;
    
}
