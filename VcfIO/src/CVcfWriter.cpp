//
//  CVcfWriter.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/16/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVcfWriter.h"
#include <ctime>
#include <iostream>


CVcfWriter::CVcfWriter()
{
    m_HEADER_GUARD = 0;
    m_nSampleCount = 0;
}

void CVcfWriter::CreateVcf(const char* a_pFileName)
{
    m_pHtsFile = hts_open(a_pFileName,"w");
    m_pRecord  = bcf_init1();
}

void CVcfWriter::CloseVcf()
{
    bcf_destroy1(m_pRecord);
    bcf_hdr_destroy(m_pHeader);
    int ret= hts_close(m_pHtsFile);
    if(ret != 0)
    {
        std::cerr << "A problem occured with saving the VCF file." << std::endl;
    }
}

void CVcfWriter::AddHeaderLine(const std::string& a_rLine)
{
    bcf_hdr_append(m_pHeader, a_rLine.c_str());
}

void CVcfWriter::InitHeader()
{
    if(m_HEADER_GUARD == 0)
    {
        //Initialize header for write
        m_pHeader = bcf_hdr_init("w");
        //Add file date
        bcf_hdr_append(m_pHeader, ("##fileDate=" + GetTime()).c_str());
        m_HEADER_GUARD = 1;
    }
    else
        std::cerr << "Invalid Operation. Header cannot be initialized" << std::endl;

}

void CVcfWriter::WriteHeaderToVcf()
{
    if(m_HEADER_GUARD == 1)
    {
        //to update internal structures
        bcf_hdr_add_sample(m_pHeader, NULL);
        
        bcf_hdr_write(m_pHtsFile, m_pHeader);
        m_HEADER_GUARD = -1;
    }
    else
        std::cerr << "Invalid Operation. Header cannot be writtten into file" << std::endl;
    
}

void CVcfWriter::AddSampleName(const std::string &a_rSampleName)
{
    if(m_HEADER_GUARD == 1)
    {
        bcf_hdr_add_sample(m_pHeader, a_rSampleName.c_str());
        m_nSampleCount++;
    }
    else
        std::cerr << "Invalid Operation. Cannot add sample name" << std::endl;
}

void CVcfWriter::AddRecord(const SVcfRecord& a_rVcfRecord)
{
    if(m_HEADER_GUARD != -1)
    {
        std::cerr << "Invalid Operation. Cannot add record before submitting header" << std::endl;
        return;
    }
    
    //Clear previous record
    bcf_clear1(m_pRecord);
    
    //Set chromosome name
    m_pRecord->rid = bcf_hdr_name2id(m_pHeader, GetChrName(a_rVcfRecord.m_nChrId).c_str());
    
    //Set start position
    m_pRecord->pos = a_rVcfRecord.m_nPosition;
    
    //Set Quality
    if(a_rVcfRecord.m_nQuality != -1)
        m_pRecord->qual = a_rVcfRecord.m_nQuality;
    
    //Set alleles
    bcf_update_alleles_str(m_pHeader, m_pRecord, a_rVcfRecord.m_alleles.c_str());
    
    //Set filter
    if(a_rVcfRecord.m_aFilterString.size() != 0)
    {
        int32_t* tmpi = new int[a_rVcfRecord.m_aFilterString.size()];
        for(int k = 0; k < a_rVcfRecord.m_aFilterString.size(); k++)
            tmpi[k] = bcf_hdr_id2int(m_pHeader, BCF_DT_ID, a_rVcfRecord.m_aFilterString[k].c_str());
        bcf_update_filter(m_pHeader, m_pRecord, tmpi, static_cast<int>(a_rVcfRecord.m_aFilterString.size()));
    }
    
    //==Set Per Sample Data==

    //1.Genotype Set (GT)
    std::vector<int> genotypes;
    for(int k = 0; k < a_rVcfRecord.m_aSampleData.size(); k++)
    {
        for(int p = 0; p < a_rVcfRecord.m_aSampleData[k].m_nHaplotypeCount; p++)
        {
            if(a_rVcfRecord.m_aSampleData[k].m_bIsPhased)
                genotypes.push_back(bcf_gt_phased(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p]));
            else if(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p] == -1)
                genotypes.push_back(bcf_gt_missing);
            else
                genotypes.push_back(bcf_gt_unphased(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p]));
        }
        
        if(a_rVcfRecord.m_aSampleData[k].m_nHaplotypeCount == 1)
            genotypes.push_back(bcf_int32_vector_end);
    }

    
    if(m_nSampleCount != a_rVcfRecord.m_aSampleData.size())
    {
        genotypes.push_back(bcf_gt_missing);
        genotypes.push_back(bcf_gt_missing);
    }
    
    
    bcf_update_genotypes(m_pHeader, m_pRecord, static_cast<int*>(&genotypes[0]), bcf_hdr_nsamples(m_pHeader)*2);

    //2.Decision Set (BD)
    char* tmpstr[m_nSampleCount];
    int k;
    for(k = 0; k < a_rVcfRecord.m_aSampleData.size(); k++)
    {
        tmpstr[k] = new char[a_rVcfRecord.m_aSampleData[k].m_decisionBD.size()];
        strcpy(tmpstr[k], a_rVcfRecord.m_aSampleData[k].m_decisionBD.c_str());
    }
    
    if(a_rVcfRecord.m_aSampleData.size() != m_nSampleCount)
    {
        tmpstr[k] = new char[1];
        tmpstr[k][0] = bcf_str_missing;
    }
    bcf_update_format_string(m_pHeader, m_pRecord, "BD", (const char**)tmpstr, m_nSampleCount);
    

    //3.Match Type Set (BK)
    char* tmpstr2[m_nSampleCount];
    for(k = 0; k < a_rVcfRecord.m_aSampleData.size(); k++)
    {
        tmpstr2[k] = new char[a_rVcfRecord.m_aSampleData[k].m_matchTypeBK.size()];
        strcpy(tmpstr2[k], a_rVcfRecord.m_aSampleData[k].m_matchTypeBK.c_str());
    }
    
    if(a_rVcfRecord.m_aSampleData.size() != m_nSampleCount)
    {
        tmpstr2[k] = new char[1];
        tmpstr2[k][0] = bcf_str_missing;
    }
    bcf_update_format_string(m_pHeader, m_pRecord, "BK", (const char**)tmpstr2, m_nSampleCount);
    
 
    //Write record to created VCF File
    bcf_write1(m_pHtsFile, m_pHeader, m_pRecord);
    
    
    //Clean Temporary strings we used
    for(int k = 0; k < a_rVcfRecord.m_aSampleData.size(); k++)
    {
        delete[] tmpstr[k];
        delete[] tmpstr2[k];
    }
    
    
}

void CVcfWriter::AddMendelianRecord(const SVcfRecord& a_rVcfRecord)
{
    int success = 0;
    
    
    if(m_HEADER_GUARD != -1)
    {
        std::cerr << "Invalid Operation. Cannot add record before submitting header" << std::endl;
        return;
    }
    
    //Clear previous record
    bcf_clear1(m_pRecord);
    
    //Set chromosome name
    m_pRecord->rid = bcf_hdr_name2id(m_pHeader, GetChrName(a_rVcfRecord.m_nChrId).c_str());
    
    //Set start position
    m_pRecord->pos = a_rVcfRecord.m_nPosition;
    
    //Set Quality
    if(a_rVcfRecord.m_nQuality != -1)
        m_pRecord->qual = a_rVcfRecord.m_nQuality;
    
    //Set alleles
    success = bcf_update_alleles_str(m_pHeader, m_pRecord, a_rVcfRecord.m_alleles.c_str());
    if(success < 0)
        std::cerr << "Failed to update Alleles string for Record: " << "Chr" << a_rVcfRecord.m_nChrId << " Position: " << a_rVcfRecord.m_nPosition << std::endl;
    
    //Set Decision
    int decArray = atoi(a_rVcfRecord.m_mendelianDecision.c_str());
    success = bcf_update_info_int32(m_pHeader, m_pRecord, "MD", &decArray, 1);
                                    
    if(success < 0)
        std::cerr << "Failed to update MD INFO for Record: " << "Chr" << a_rVcfRecord.m_nChrId << " Position: " << a_rVcfRecord.m_nPosition << std::endl;
    
    
    //Set filter
    if(a_rVcfRecord.m_aFilterString.size() != 0)
    {
        int32_t* tmpi = new int[a_rVcfRecord.m_aFilterString.size()];
        for(int k = 0; k < a_rVcfRecord.m_aFilterString.size(); k++)
            tmpi[k] = bcf_hdr_id2int(m_pHeader, BCF_DT_ID, a_rVcfRecord.m_aFilterString[k].c_str());
        bcf_update_filter(m_pHeader, m_pRecord, tmpi, static_cast<int>(a_rVcfRecord.m_aFilterString.size()));
    }
    
    //==Set Per Sample Data==
    
    //1.Genotype Set (GT)
    int* genotypes = new int[bcf_hdr_nsamples(m_pHeader)*2];
    int genotypeItr = 0;
    for(int k = 0; k < a_rVcfRecord.m_aSampleData.size(); k++)
    {
        for(int p = 0; p < a_rVcfRecord.m_aSampleData[k].m_nHaplotypeCount; p++)
        {
            if(a_rVcfRecord.m_aSampleData[k].m_bIsPhased)
                genotypes[genotypeItr++] = bcf_gt_phased(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p]);
            else if(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p] == -1)
                genotypes[genotypeItr++] = bcf_gt_missing;
            else
                genotypes[genotypeItr++] = bcf_gt_unphased(a_rVcfRecord.m_aSampleData[k].m_aGenotype[p]);
        }
        
        if(a_rVcfRecord.m_aSampleData[k].m_nHaplotypeCount == 1)
            genotypes[genotypeItr++] = bcf_int32_vector_end;
    }

    if(m_nSampleCount != a_rVcfRecord.m_aSampleData.size())
    {
        for(int k = static_cast<int>(a_rVcfRecord.m_aSampleData.size()); k < m_nSampleCount; k++)
        {
            genotypes[genotypeItr++] = bcf_gt_missing;
            genotypes[genotypeItr++] = bcf_gt_missing;
        }
    }
    
    success = bcf_update_genotypes(m_pHeader, m_pRecord, genotypes, genotypeItr);
    if(success < 0)
        std::cerr << "Failed to update Genotypes for Record: " << "Chr" << a_rVcfRecord.m_nChrId << " Position: " << a_rVcfRecord.m_nPosition << std::endl;

    
    //Write record to created VCF File
    success = bcf_write1(m_pHtsFile, m_pHeader, m_pRecord);
    if(success < 0)
        std::cerr << "Failed to write Record to the file: " << "Chr" << a_rVcfRecord.m_nChrId << " Position: " << a_rVcfRecord.m_nPosition << std::endl;
    
    
    //Garbage Collection
    delete[] genotypes;
}



std::string CVcfWriter::GetChrName(int a_nChrId)
{
    if(a_nChrId < 23)
        return "chr" + std::to_string(a_nChrId);
    else if(a_nChrId == 23)
        return "chrX";
    else if(a_nChrId == 24)
        return "chrY";
    else
        return "chrMT";
}


std::string CVcfWriter::GetTime()
{
    std::tm* timeInfo;
    char buffer[16];
    
    std::time_t currentTime;
    
    time(&currentTime);
    timeInfo = localtime(&currentTime);

    strftime(buffer,16,"%Y%m%d",timeInfo);
    
    return std::string(buffer);
}

