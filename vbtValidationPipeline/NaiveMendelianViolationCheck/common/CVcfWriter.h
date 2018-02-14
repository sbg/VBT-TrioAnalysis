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
 *  Created by Berke Cagkan Toptas
 *
 */

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
    //Is the genotype Nocall(./.)
    bool m_bIsNoCallVariant = false;
};


struct SVcfRecord
{
    //Position of the variant (0 based)
    int m_nPosition;
    //Quality of the variant
    int m_nQuality = -1;
    //Mendelian Decision INFO (Used for mendelian comparison feature)
    std::string m_mendelianDecision = "";
    //Filter string of the variant (eg. "PASS")
    std::vector<std::string> m_aFilterString;
    //Alleles string separated by comma of the variant (eg. m_alleles = "AT,G")
    std::string m_alleles;
    //Chromosome Name of the variant
    std::string m_chrName;
    //Sample Data (Data to store for each sample)
    std::vector<SPerSampleData> m_aSampleData;
};


class CVcfWriter
{

public:
    
    CVcfWriter();
    
    //Creates VCF file at the given path
    void CreateVcf(const char* a_pFileName);
        
    //Append the given variant to the opened vcf file
    void AddRecord(const SVcfRecord& a_rVcfRecord);
    
    //Append the given variant to the opned vcf file for Mendelian Trio Mode
    void AddMendelianRecord(const SVcfRecord& a_rVcfRecord);
    
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
    
    ///Set raw header to the output vcf
    void SetRawHeader(bcf_hdr_t* a_pHeader);

    ///Append the given raw variant to the opened vcf file
    void AddRawRecord(bcf1_t* a_rRecord);
    
private:
    
    //Return the current time in YYYYMMDD format
    std::string GetTime();
    
    htsFile *   m_pHtsFile;
    bcf_hdr_t * m_pHeader;
    bcf1_t *    m_pRecord;
    
    //-1 : Close Header / 1: Inside Header / 0: Header Not Opened
    int m_HEADER_GUARD;
    
    //Stores sample count added
    int m_nSampleCount;
};




#endif //_C_VCF_WRITER_H_
