//
//  CVcfReader.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include <stdio.h>
#include "CVcfReader.h"
#include <iostream>

CVcfReader::CVcfReader()
{
    m_bIsOpen = false;
}

CVcfReader::CVcfReader(const char * a_pFilename)
{
    m_bIsOpen = false;
    Open(a_pFilename);
}

CVcfReader::~CVcfReader()
{
    Close();
}

bool CVcfReader::Open(const char * a_pFilename)
{
    if (true == m_bIsOpen)
    {
        fprintf(stderr, "ERROR: The vcf has been opened.\n");
        return false;
    }
    
    assert(a_pFilename);
    m_pHtsFile = bcf_open(a_pFilename, "r");
    
    if (m_pHtsFile == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open the vcf file: %s.\n", a_pFilename);
        return false;
    }
    
    // Read header
    m_pHeader = bcf_hdr_read(m_pHtsFile);
    if (m_pHeader == NULL)
    {
        fprintf(stderr, "ERROR: Cannot load the header of vcf.\n");
        return false;
    }
    else
    {
        for (int i = 0; i < m_pHeader->n[BCF_DT_CTG]; ++i)
        {
            SVcfContig contig;
            contig.name   = m_pHeader->id[BCF_DT_CTG][i].key;
            contig.length = m_pHeader->id[BCF_DT_CTG][i].val->info[0];
            contigs_.push_back(contig);
        }
    }
        
    m_pRecord  = bcf_init();
    assert(m_pRecord);
    
    m_bIsOpen = true;
    return true;
}

bool CVcfReader::Close()
{
    if (m_bIsOpen)
    {
        bcf_hdr_destroy(m_pHeader);
        bcf_destroy(m_pRecord);
        bcf_close(m_pHtsFile);
        m_bIsOpen = false;
    }
    m_pHeader = NULL;
    m_pRecord = NULL;
    m_pHtsFile      = NULL;
    

    return true;
}

bool CVcfReader::GetNextRecord(CVariant * a_pVariant, int a_nId, const SConfig& a_rConfig)
{
    a_pVariant->Clear();
    a_pVariant->m_nVcfId = m_nVcfId;
    int* gt_arr = NULL;
    int ngt_arr = 0;

    int samplenumber = GetNumberOfSamples();
    int zygotCount;
    
    bcf_clear(m_pRecord);
    m_pRecord->d.m_allele = 0;
    const int ok = bcf_read(m_pHtsFile, m_pHeader, m_pRecord);
    bcf_unpack(m_pRecord, BCF_UN_ALL);
    
    if (ok == 0)
    {
        a_pVariant->m_nId = a_nId;
        a_pVariant->m_chrName = m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key;
        
        //READ CHROMOSOME ID: (TODO: this should be renewed. there should be sth that reads the chromosome id)
        if(a_pVariant->m_chrName == "x" || a_pVariant->m_chrName == "X" || a_pVariant->m_chrName == "chrX" || a_pVariant->m_chrName == "chrx")
            a_pVariant->m_nChrId = 23;
        else if(a_pVariant->m_chrName == "y" || a_pVariant->m_chrName == "Y" || a_pVariant->m_chrName == "chrY" || a_pVariant->m_chrName == "chry")
            a_pVariant->m_nChrId = 24;
        else if(a_pVariant->m_chrName == "MT" || a_pVariant->m_chrName == "mt" || a_pVariant->m_chrName == "chrMT" || a_pVariant->m_chrName == "chrmt")
            a_pVariant->m_nChrId = 25;
        else if(a_pVariant->m_chrName.length() == 5)
            a_pVariant->m_nChrId = atoi(a_pVariant->m_chrName.substr(3,2).c_str());
        else if(a_pVariant->m_chrName.length() == 4)
            a_pVariant->m_nChrId = atoi(a_pVariant->m_chrName.substr(3,1).c_str());
        else
            a_pVariant->m_nChrId = atoi(m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key);
    
        
        //READ FILTER DATA
        if(a_rConfig.m_bIsFilterEnabled)
        {
            bool isPassed = false;
            for(int k=0; k< m_pRecord->d.n_flt; k++)
            {
                if(k != 0)
                    a_pVariant->m_filterString += ",";
                a_pVariant->m_filterString += m_pHeader->id[BCF_DT_ID][m_pRecord->d.flt[k]].key;
                
                if(0 == strcmp(m_pHeader->id[BCF_DT_ID][m_pRecord->d.flt[k]].key, a_rConfig.m_pFilterName))
                    isPassed = true;
            }
            a_pVariant->m_bIsFilterPASS = isPassed;
        }
        
        //READ QUALITY DATA
        //TODO: read the quality
        
        
        //READ GENOTYPE DATA
        bcf_get_genotypes(m_pHeader, m_pRecord, &gt_arr, &ngt_arr);
        zygotCount = ngt_arr / samplenumber;
        a_pVariant->m_nAlleleCount = zygotCount;
        a_pVariant->m_bIsPhased = bcf_gt_is_phased(bcf_gt_allele(gt_arr[0]));
 
        
        //READ SEQUENCE DATA AND FILL ALLELES
        a_pVariant->m_refSequence = std::string(m_pRecord->d.allele[0]);
        
        for (int i = 0; i < zygotCount; ++i)
        {
            int index = bcf_gt_allele(gt_arr[i]) == -1 ? 0 : bcf_gt_allele(gt_arr[i]);
            a_pVariant->m_alleles[i].m_sequence = m_pRecord->d.allele[index];
            a_pVariant->m_alleles[i].m_nStartPos = m_pRecord->pos;
            a_pVariant->m_alleles[i].m_nEndPos = static_cast<int>(m_pRecord->pos + a_pVariant->m_refSequence.length());
            a_pVariant->m_nMaxLength = std::max(a_pVariant->m_nMaxLength,static_cast<int>(a_pVariant->m_alleles[i].m_sequence.length()));
        }
        
        //SET ZYGOSITY OF THE VARIANT (HOMOZYGOUS or HETEROZYGOUS)
        if(zygotCount == 2)
        {
            if(a_pVariant->m_alleles[0].m_sequence == a_pVariant->m_alleles[1].m_sequence)
            {
                a_pVariant->m_nAlleleCount = 1;
                a_pVariant->m_bIsHeterozygous = false;
            }
            else
            {
                a_pVariant->m_nAlleleCount = 2;
                a_pVariant->m_bIsHeterozygous = true;
            }
        }
        else
        {
            a_pVariant->m_nAlleleCount = 1;
            a_pVariant->m_bIsHeterozygous = false;
        }
        
        //TRIM FIRST NUCLEOTIDES IF THEY EXIST IN BOTH ALLELE AND IN REFERENCE
        if(HasRedundantFirstNucleotide())
        {
            a_pVariant->m_bIsFirstNucleotideTrimmed = true;
            for (int i = 0; i < ngt_arr; ++i)
            {
                TrimAllele(a_pVariant->m_alleles[i]);
            }
        }
        
        //SET START AND END POSITION OF VARIANT
        int maxEnd = -1;
        int minStart = INT_MAX;
        a_pVariant->m_nStartPos = m_pRecord->pos;
        for(int k=0; k < ngt_arr; k++)
        {
            maxEnd = std::max(maxEnd, static_cast<int>(a_pVariant->m_alleles[k].m_nEndPos));
            minStart = std::min(minStart, static_cast<int>(a_pVariant->m_alleles[k].m_nStartPos));
        }
        a_pVariant->m_nEndPos = maxEnd;
        a_pVariant->m_nStartPos = minStart;
 
        
        //FILL ORIGINAL ALLELE STR AND GENOTYPE FOR LATER ACCESS
        a_pVariant->m_nZygotCount = zygotCount;
        for(int k = 0; k < m_pRecord->n_allele; k++)
        {
            if(k != 0)
                a_pVariant->m_allelesStr += ",";
            a_pVariant->m_allelesStr += std::string(m_pRecord->d.allele[k]);
        }
       for(int k = 0; k < ngt_arr; k++)
           a_pVariant->m_genotype[k] = bcf_gt_allele(gt_arr[k]);
        
        //FREE BUFFERS
        free(gt_arr);
        return true;
    }
    else 
    {
        a_pVariant = 0;
        return false;
    }
}


bool CVcfReader::GetNextRecordMultiSample(CVariant* a_pVariant)
{
    int* gt_arr = NULL;
    int ngt_arr = 0;
    //int* fi_arr = NULL;
    
    
    if (!m_bIsOpen)
    {
        fprintf(stderr, "ERROR The vcf has not been opened.\n");
        return false;
    }

    //CLEAR THE OLD RECORD AND READ THE NEW RECORD
    bcf_clear(m_pRecord);
    m_pRecord->d.m_allele = 0;
    const int ok = bcf_read(m_pHtsFile, m_pHeader, m_pRecord);
    bcf_unpack(m_pRecord, BCF_UN_ALL);

    //UNABLE TO READ RECORD
    if(ok != 0)
        return false;

    
    //READ GENOTYPE DATA
    bcf_get_genotypes(m_pHeader, m_pRecord, &gt_arr, &ngt_arr);
    
    
    int samplenumber = GetNumberOfSamples();
    for(int k = 0; k < samplenumber; k++)
    {
        a_pVariant[k].Clear();
        a_pVariant[k].m_chrName = m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key;
        
        //READ CHROMOSOME ID: (TODO: this should be renewed. there should be sth that reads the chromosome id)
        if(a_pVariant[k].m_chrName.length() == 5)
            a_pVariant[k].m_nChrId = atoi(a_pVariant->m_chrName.substr(3,2).c_str());
        else if(a_pVariant[k].m_chrName.length() == 4)
            a_pVariant[k].m_nChrId = atoi(a_pVariant->m_chrName.substr(3,1).c_str());
        else
            a_pVariant[k].m_nChrId = atoi(m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key);
        
        //READ FILTER DATA
        //TODO: read filter data
        
        //READ QUALITY DATA
        //TODO: read the quality
        
        a_pVariant[k].m_nAlleleCount = ngt_arr/samplenumber;
        a_pVariant[k].m_bIsPhased = bcf_gt_is_phased(bcf_gt_allele(gt_arr[k * a_pVariant[k].m_nAlleleCount]));
        
        //READ SEQUENCE DATA AND FILL ALLELES
        a_pVariant->m_refSequence = std::string(m_pRecord->d.allele[0]);
        
        
        for (int i = 0; i < a_pVariant[k].m_nAlleleCount; ++i)
        {
            int index = bcf_gt_allele(gt_arr[2*k + i]) == -1 ? 1 : bcf_gt_allele(gt_arr[2*k + i]);
            a_pVariant[k].m_alleles[i].m_sequence = m_pRecord->d.allele[index];
            a_pVariant[k].m_alleles[i].m_nStartPos = m_pRecord->pos;
            a_pVariant[k].m_alleles[i].m_nEndPos = static_cast<int>(m_pRecord->pos + a_pVariant->m_refSequence.length());
            a_pVariant[k].m_nMaxLength = std::max(a_pVariant->m_nMaxLength,static_cast<int>(a_pVariant->m_alleles[i].m_sequence.length()));
        }
        
        if(a_pVariant[k].m_nAlleleCount == 2)
        {
            if(a_pVariant[k].m_alleles[0].m_sequence == a_pVariant[k].m_alleles[1].m_sequence)
            {
                a_pVariant[k].m_nAlleleCount = 1;
                a_pVariant[k].m_bIsHeterozygous = false;
            }
            else
            {
                a_pVariant[k].m_nAlleleCount = 2;
                a_pVariant[k].m_bIsHeterozygous = true;
            }
        }
        else
        {
            a_pVariant[k].m_nAlleleCount = 1;
            a_pVariant[k].m_bIsHeterozygous = false;
        }
        
        //CTODO TRIM ALLELE!!
        
        //SET START AND END POSITION OF VARIANT
        int maxEnd = -1;
        int minStart = INT_MAX;
        a_pVariant[k].m_nStartPos = m_pRecord->pos;
        for(int p=0; p < ngt_arr/samplenumber; p++)
        {
            maxEnd = std::max(maxEnd, static_cast<int>(a_pVariant[k].m_alleles[p].m_nEndPos));
            minStart = std::min(minStart, static_cast<int>(a_pVariant[k].m_alleles[p].m_nStartPos));
        }
        a_pVariant[k].m_nEndPos = maxEnd;
        a_pVariant[k].m_nStartPos = minStart;
        
        
    }

    //FREE BUFFERS
    free(gt_arr);
    //free(fi_arr);
    return true;
}

bool CVcfReader::SelectSample(std::string a_sampleName)
{
    int res = bcf_hdr_set_samples(m_pHeader, a_sampleName.c_str(), 0);
    
    if(res == 0)
        return true;
    else
        return false;
    
}

int CVcfReader::GetNumberOfSamples() const
{

    if (!m_bIsOpen)
    {
        fprintf(stderr, "ERROR: The vcf has not been opened.\n");
        return -1;
    }
    
    return bcf_hdr_nsamples(m_pHeader);
}

void CVcfReader::GetSampleNames(std::vector<std::string>& a_pSampleNameList)
{
    int samplecount = bcf_hdr_nsamples(m_pHeader);
    
    for(int k =0; k < samplecount; k++)
        a_pSampleNameList.push_back(m_pHeader->id[BCF_DT_SAMPLE][k].key);
}

int CVcfReader::GetContigId(std::string a_name) const
{
    for (unsigned int i = 0; i < contigs_.size(); ++i)
    {
        if (contigs_[i].name.compare(a_name) == 0)
            return i;
    }
    return -1;
}

const char* CVcfReader::getFilterString(int a_nFilterKey)
{
    return m_pHeader->id[BCF_DT_ID][a_nFilterKey].key;
}

int CVcfReader::getFilterKey(const char* a_pFilterValue)
{
    for(int k = 0; k < m_pHeader->n[BCF_DT_ID]; k++)
    {
        if(0 == strcmp(m_pHeader->id[BCF_DT_ID][k].key, a_pFilterValue))
            return k;
    }

    return -1;
}

void CVcfReader::TrimAllele(SAllele& a_rAllele)
{
    a_rAllele.m_sequence = a_rAllele.m_sequence.substr(1, a_rAllele.m_sequence.length() - 1);
    a_rAllele.m_nStartPos += 1;
}

bool CVcfReader::HasRedundantFirstNucleotide() const
{
    for(int k= 0; k < m_pRecord->n_allele; k++)
    {
        if(m_pRecord->d.allele[0][0] != m_pRecord->d.allele[k][0])
            return false;
    }
    return true;
}

const std::vector<SVcfContig>&  CVcfReader::GetContigs() const
{
    return contigs_;
}

void CVcfReader::PrintVCF(const SConfig& a_rConfig)
{
    CVariant variant;

    // READ SAMPLE VCF
    //fprintf(stderr,"#CHR\tPOS\tREF\tALTs\n");
    int k=0;
    while(GetNextRecord(&variant,0, a_rConfig))
    {
        std::cout << k++ << ": " <<  variant.ToString() << std::endl;
    }
}
