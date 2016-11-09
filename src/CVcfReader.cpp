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

bool CVcfReader::GetNextRecord(CVariant * a_pVariant)
{
    a_pVariant->Clear();
    a_pVariant->m_nVcfId = m_nVcfId;
    int* gt_arr = NULL;
    int ngt_arr = 0;

    if (!m_bIsOpen)
    {
        fprintf(stderr, "ERROR The vcf has not been opened.\n");
        return false;
    }
    
    bcf_clear(m_pRecord);
    m_pRecord->d.m_allele = 0;
    const int ok = bcf_read(m_pHtsFile, m_pHeader, m_pRecord);
    bcf_unpack(m_pRecord, BCF_UN_ALL);

    if (ok == 0)
    {
        a_pVariant->m_nPosition = m_pRecord->pos;
        a_pVariant->m_chrName = m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key;
        a_pVariant->m_aSequences.push_back(m_pRecord->d.als);
        
        //Chrid is 2 digit
        if(a_pVariant->m_chrName.length() == 5)
            a_pVariant->m_nChrId = atoi(a_pVariant->m_chrName.substr(3,2).c_str());
        else if(a_pVariant->m_chrName.length() == 4)
            a_pVariant->m_nChrId = atoi(a_pVariant->m_chrName.substr(3,1).c_str());
        else
            a_pVariant->m_nChrId = atoi(m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key);
       
        int ngt = bcf_get_genotypes(m_pHeader, m_pRecord, &gt_arr, &ngt_arr);
        a_pVariant->ngt_arr = ngt_arr;
        for(int k = 0; k < ngt_arr ; k++)
        {
            //std::cout << gt_arr[k] << std::endl;
            //std::cout << bcf_gt_allele(gt_arr[k]) << std::endl;
            a_pVariant->gt_arr[k] = bcf_gt_allele(gt_arr[k]);
        }
        
        a_pVariant->m_bIsPhased = bcf_gt_is_phased(a_pVariant->gt_arr[0]);

        for (int i = 1; i < m_pRecord->n_allele; ++i)
        {
            (*a_pVariant).m_aSequences.push_back(m_pRecord->d.allele[i]);
            a_pVariant->SetType(i);
        }

        //for(int k = 0; k < ngt_arr ; k++)
        //    std::cout << bcf_gt_allele(gt_arr[k]) << " ";
        //std::cout<< std::endl;
        return true;
    }
    else 
    {
        a_pVariant = 0;
        return false;
    }
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

int CVcfReader::GetContigId(const char* name) const
{
    for (unsigned int i = 0; i < contigs_.size(); ++i) {
        if (strcmp(contigs_[i].name.c_str(), name) == 0) return i;
    }
    return -1;
}

void CVcfReader::PrintVCF()
{
    CVariant variant;

    // READ SAMPLE VCF
    fprintf(stderr,"#CHR\tPOS\tREF\tALTs\n");
    while(GetNextRecord(&variant))
    {
        fprintf(stderr,"%s\t%d\t%s", variant.m_chrName.c_str(), variant.m_nPosition, variant.m_aSequences[0].c_str());
        if (variant.m_aSequences.size() > 1)
        {
            fprintf(stderr,"\t");
            for (unsigned int i = 1; i < variant.m_aSequences.size(); ++i) 
            {
                fprintf(stderr,"%s", variant.m_aSequences[i].c_str());
                if (i != variant.m_aSequences.size() - 1) fprintf(stderr,";");
                else fprintf(stderr,"\n");
            }
        }
    }
}
