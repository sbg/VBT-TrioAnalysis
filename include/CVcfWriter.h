//
//  CVcfWriter.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/16/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_VCF_WRITER_H_
#define _C_VCF_WRITER_H_

#include "htslib/vcf.h"
#include "CVariant.h"
#include <string>

struct SPerSampleData
{
    SPerSampleData()
    {
        m_decisionBD = bcf_str_missing;
        m_matchTypeBK = bcf_str_missing;
        m_nHaplotypeCount = 2;
        m_aGenotype[0] = -1;
        m_aGenotype[1] = -1;
    }
    
    //FORMAT :  Decison of variant (TP/FP/FN/N)
    std::string m_decisionBD;
    //FORMAT : Match type of variant (gt/allele match)
    std::string m_matchTypeBK;
    //Haplotype count of the variant
    int m_nHaplotypeCount;
    //Genotype of the variant
    int m_aGenotype[2];
    //Is the genotype is phased
    bool m_bIsPhased = false;
};


struct SVcfRecord
{
    //Position of the variant (0 based)
    int m_nPosition;
    //Quality of the variant
    int m_nQuality = -1;
    //Filter string of the variant (eg. "PASS")
    std::vector<std::string> m_aFilterString;
    //Alleles string separated by comma of the variant (eg. m_alleles = "AT,G")
    std::string m_alleles;
    //Chromosome id of the variant
    int m_nChrId;
    //Sample Data (Data to store for each sample)
    std::vector<SPerSampleData> m_aSampleData;
};


class CVcfWriter
{

public:
    
    CVcfWriter();
    
    //Creates VCF file at the given path
    void CreateVcf(const char* a_pFileName);
    
    //Creates the vcf header according to ga4gh and writes it to the opened vcf file
    void FillHeader();
    
    //Append the given variant to the opened vcf file
    void AddRecord(const SVcfRecord& a_rVcfRecord);
    
    //Close the VCF file
    void CloseVcf();
  
    //Add line to the header of vcf file
    void AddHeaderLine(const std::string& a_rLine);
    
    //Add new sample to the header
    void AddSampleName(const std::string& a_rSampleName);
    
    //Open header for writing
    void InitHeader();
    
    //Write the header to the vcf - No change should be made to the m_pHeader after this function is called.
    void WriteHeaderToVcf();
    
private:
    
    //Return the current time in YYYYMMDD format
    std::string GetTime();
    //Return chromosome name with given chromosome id
    std::string GetChrName(int a_nChrId);
    
    htsFile *   m_pHtsFile;
    bcf_hdr_t * m_pHeader;
    bcf1_t *    m_pRecord;
    
    //-1 : Close Header / 1: Inside Header / 0: Header Not Opened
    int m_HEADER_GUARD;
    
    //Stores sample count added
    int m_nSampleCount;
};




#endif //_C_VCF_WRITER_H_
