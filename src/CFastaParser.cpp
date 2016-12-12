//
//  CFastaParser.cpp
//
//  Created by John Browning on 4/10/14.
//  Copyright (c) 2014 John Browning. All rights reserved.
//
//
//  Modified by: Berke Cagkan Toptas
//
//  This file opens a FASTA file and stores it in a string. It also contains a subroutine
//  to fetch a base given the index.
//

#include "CFastaParser.h"
#include <iostream>
#include <fstream>

bool CFastaParser::OpenFastaFile(const char *fn, std::string chrom)
{
    bool bIsSuccess = true;
    
    bIsSuccess = GenerateFastaIndex(fn);
    
    if(true == bIsSuccess)
    {
        char chrom_num[80] = {0};
        strcpy(chrom_num, chrom.c_str());
        fai = fai_load(fn);
        ref = faidx_fetch_seq(fai, chrom_num, 0, 0x7fffffff, &ref_len);
        
        if(ref_len == -1)
        {
            std::cout << "An Error is occured while reading FASTA file" << std::endl;
            bIsSuccess = false;
        }
        else if(ref_len == -2)
        {
            std::cout << "Specified chromosome " << chrom << " could not found in FASTA file" << std::endl;
            bIsSuccess = false;
        }
    }
    
    return bIsSuccess;
}

bool CFastaParser::GenerateFastaIndex(const char *fn)
{
    bool bIsSuccess = true;
    
    //Generate FASTA index file if it does not exists
    std::string faiName(fn);
    faiName = faiName + ".fai";
    
    //Try to open the file
    std::ifstream f(faiName.c_str());
    if(!f.good())
    {
        std::cout << "FASTA index file for " << fn << " does not exists. Generating the index file..." << std::endl;
        bIsSuccess = true;
        int res = fai_build(fn);
        if(res == 0)
            bIsSuccess = true;
        else
        {
            bIsSuccess = false;
            std::cout << "An error occured while generating FASTA index file for " << fn << std::endl;
        }
    }

    return bIsSuccess;
}

bool CFastaParser::FetchNewChromosome(std::string chrom)
{
    bool bIsSuccess = true;
    
    char chrom_num[80] = {0};
    strcpy( chrom_num, chrom.c_str());
    delete ref;
    ref = faidx_fetch_seq(fai, chrom_num, 0, 0x7fffffff, &ref_len);
    
    if(ref_len == -1)
    {
        std::cout << "An Error is occured while reading FASTA file" << std::endl;
        bIsSuccess = false;
    }
    else if(ref_len == -2)
    {
        std::cout << "Specified chromosome" << chrom  << " could not found in FASTA file" << std::endl;
        bIsSuccess = false;
    }
    
    return bIsSuccess;
    
}

void CFastaParser::DeleteReferenceMemory()
{
    free(ref);
    fai_destroy(fai);
}

std::string CFastaParser::GetRefBase(int a_nRefPos, int a_nLength)
{
    std::string temp; // = ref_sequence.substr(ref_pos-1, length);
    for (int i(0); i < a_nLength; ++i)
    {
        char tc = ref[a_nRefPos - 1 + i];
        if (tc > 'T')
        {
            tc = ref[a_nRefPos - 1 + i] - ' ';
            if ((tc != 'A') && (tc != 'C') && (tc != 'G') && (tc != 'T'))
                tc = 'N';
            temp += tc;
        }
        else
            temp += tc;
    }
    return temp;
}

char CFastaParser::GetRefBase(int a_nPosition)
{
    char tc = ref[a_nPosition - 1];
    
    if (tc > 'T')
        tc -= ' ';
    else
        tc = ref[a_nPosition - 1];
    
    if ((tc != 'A') && (tc != 'C') && (tc != 'G') && (tc != 'T'))
        tc = 'N';
    
    return tc;
}

int CFastaParser::GetRefLength()
{
    return ref_len;
}
