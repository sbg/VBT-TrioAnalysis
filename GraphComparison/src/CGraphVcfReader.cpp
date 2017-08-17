//
//  CGraphVcfReader.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/11/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CGraphVcfReader.h"

using namespace graphcomparison;

CGraphVcfReader::CGraphVcfReader()
{
    m_bIsOpen = false;
}

CGraphVcfReader::~CGraphVcfReader()
{
    Close();
}

bool CGraphVcfReader::Open(const char * a_pFilename)
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
            m_contigs.push_back(contig);
            m_chrIndexMap[contig.name] = i;
        }
    }
    
    m_pRecord  = bcf_init();
    assert(m_pRecord);
    
    m_bIsOpen = true;
    return true;
}

bool CGraphVcfReader::Close()
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
    m_pHtsFile = NULL;
    
    
    return true;
}

bool CGraphVcfReader::GetNextRecord(CVariant * a_pVariant, int a_nId, bool a_bIsReadAll)
{
    a_pVariant->Clear();
    a_pVariant->m_nVcfId = 1;
    bcf_clear(m_pRecord);
    m_pRecord->d.m_allele = 0;
    const int ok = bcf_read(m_pHtsFile, m_pHeader, m_pRecord);
    
    if(!a_bIsReadAll)
        bcf_unpack(m_pRecord, BCF_UN_FLT);
    else
        bcf_unpack(m_pRecord, BCF_UN_ALL);
    
    if (ok == 0)
    {
        a_pVariant->m_nId = a_nId;
        a_pVariant->m_chrName = m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key;
        a_pVariant->m_nChrId = m_chrIndexMap[m_pHeader->id[BCF_DT_CTG][m_pRecord->rid].key];
        
        //READ FILTER DATA
        bool isPassed = false;
        
        if(m_pRecord->d.n_flt == 0)
            isPassed = true;
        
        for(int k=0; k< m_pRecord->d.n_flt; k++)
        {
            a_pVariant->m_filterString.push_back(std::string(m_pHeader->id[BCF_DT_ID][m_pRecord->d.flt[k]].key));
            
            if(0 == strcmp(m_pHeader->id[BCF_DT_ID][m_pRecord->d.flt[k]].key, "PASS"))
                isPassed = true;
            else if(0 == strcmp(m_pHeader->id[BCF_DT_ID][m_pRecord->d.flt[k]].key, "."))
                isPassed = true;
        }
        
        a_pVariant->m_bIsFilterPASS = isPassed;
        
        //READ SEQUENCE DATA AND FILL ALLELES
        a_pVariant->m_refSequence = std::string(m_pRecord->d.allele[0]);
        
        for (int i = 0; i < 2; ++i)
        {
            a_pVariant->m_alleles[i].m_sequence = m_pRecord->d.allele[1];
            a_pVariant->m_alleles[i].m_nStartPos = m_pRecord->pos;
            a_pVariant->m_alleles[i].m_nEndPos = static_cast<int>(m_pRecord->pos + a_pVariant->m_refSequence.length());
        }
        
        //SET ZYGOSITY OF THE VARIANT (HOMOZYGOUS or HETEROZYGOUS)
        a_pVariant->m_nAlleleCount = 1;
        a_pVariant->m_nZygotCount = 2;
        a_pVariant->m_bIsHeterozygous = false;
        
        //TRIM ALL REDUNDANT NUCLEOTIDES FROM BEGINNING AND END TO ALLOW MORE REFERENCE OVERLAPING
        a_pVariant->m_alleles[0].m_bIsIgnored = false;
        TrimRefOverlap(a_pVariant->m_alleles[0]);
        a_pVariant->m_alleles[1].m_bIsIgnored = false;
        TrimRefOverlap(a_pVariant->m_alleles[1]);

        
        //SET START AND END POSITION OF VARIANT
        int maxEnd = -1;
        int minStart = INT_MAX;
        for(int k=0; k < a_pVariant->m_nAlleleCount; k++)
        {
            if(!a_pVariant->m_alleles[k].m_bIsIgnored)
            {
                maxEnd = std::max(maxEnd, static_cast<int>(a_pVariant->m_alleles[k].m_nEndPos));
                minStart = std::min(minStart, static_cast<int>(a_pVariant->m_alleles[k].m_nStartPos));
            }
        }
        a_pVariant->m_nEndPos = maxEnd == -1 ? m_pRecord->pos : maxEnd;
        a_pVariant->m_nStartPos = minStart == INT_MAX ? m_pRecord->pos : minStart;
        
        //FILL ORIGINAL ALLELE STR AND GENOTYPE FOR LATER ACCESS
        for(int k = 0; k < m_pRecord->n_allele; k++)
        {
            if(k != 0)
                a_pVariant->m_allelesStr += ",";
            a_pVariant->m_allelesStr += std::string(m_pRecord->d.allele[k]);
        }
        
        //Set original genotypes
        a_pVariant->m_genotype[0] = 1;
        a_pVariant->m_genotype[1] = 1;
        a_pVariant->m_bIsPhased = false;
        
        //Set original position
        a_pVariant->m_nOriginalPos = m_pRecord->pos;
        
        return true;
    }
    else
    {
        a_pVariant = 0;
        return false;
    }
}


void CGraphVcfReader::TrimRefOverlap(SAllele& a_rAllele)
{
    //Overlapping deletion
    if(a_rAllele.m_sequence == "*")
    {
        //a_rAllele.m_sequence = "";
        return;
    }
    
    //Ref string
    std::string refString = m_pRecord->d.allele[0];
    
    int trimLengthFromBeginning = 0;
    int trimLengthFromEnd = 0;
    
    //Trim from the beginning
    int compSize = static_cast<int>(std::min(refString.size(), a_rAllele.m_sequence.size()));
    for(int k = 0; k < compSize; k++)
    {
        if(a_rAllele.m_sequence[k] == refString[k] && k == compSize -1)
        {
            trimLengthFromBeginning++;
            a_rAllele.m_nStartPos += trimLengthFromBeginning;
            break;
        }
        
        else if(a_rAllele.m_sequence[k] != refString[k])
        {
            a_rAllele.m_nStartPos += trimLengthFromBeginning;
            break;
        }
        
        else
            trimLengthFromBeginning++;
    }
    
    //Cut the beginning of the string
    a_rAllele.m_sequence = a_rAllele.m_sequence.substr(trimLengthFromBeginning);
    
    //Trim from the end
    for(int k = static_cast<int>(refString.size() - 1), p = static_cast<int>(a_rAllele.m_sequence.size() - 1); k >= trimLengthFromBeginning && p >= 0;  k--, p--)
    {
        if(a_rAllele.m_sequence[p] == refString[k] && p == 0)
        {
            trimLengthFromEnd = static_cast<int>(a_rAllele.m_sequence.size());
            a_rAllele.m_nEndPos -= trimLengthFromEnd;
            break;
        }
        
        else if(a_rAllele.m_sequence[p] == refString[k] && k == trimLengthFromBeginning)
        {
            trimLengthFromEnd = static_cast<int>(refString.size()) - trimLengthFromBeginning;
            a_rAllele.m_nEndPos -= trimLengthFromEnd;
            break;
        }
        
        
        else if(a_rAllele.m_sequence[p] != refString[k])
        {
            a_rAllele.m_nEndPos -= trimLengthFromEnd;
            break;
        }
        else
            trimLengthFromEnd++;
    }
    
    //Cut the end of the string
    a_rAllele.m_sequence = a_rAllele.m_sequence.substr(0, a_rAllele.m_sequence.length() - trimLengthFromEnd);
    
}

bcf_hdr_t* CGraphVcfReader::GetHeaderPointer()
{
    return m_pHeader;
}


bcf1_t* CGraphVcfReader::GetRecordPointer()
{
    return m_pRecord;
}






