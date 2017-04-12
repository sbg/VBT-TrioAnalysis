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
#include <cstring>
#include <fstream>
#include <string>

bool CFastaParser::OpenFastaFile(const char *fn)
{
    bool bIsSuccess = true;
    
    bIsSuccess = GenerateFastaIndex(fn);

    if(true == bIsSuccess)
        fai = fai_load(fn);
    
    if(fai == 0)
        bIsSuccess = false;
    
    return bIsSuccess;
}

CFastaParser::~CFastaParser()
{
    if(fai != 0)
        fai_destroy(fai);
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

bool CFastaParser::FetchNewChromosome(std::string chromosome, SContig& a_rContig)
{
    bool bIsSuccess = true;
    
    char chrom_num[80] = {0};
    strcpy(chrom_num, chromosome.c_str());
    a_rContig.m_chromosome = chromosome;
    
    a_rContig.m_pRefSeq = faidx_fetch_seq(fai, chrom_num, 0, 0x7fffffff, &a_rContig.m_nRefLength);
    
    if(a_rContig.m_nRefLength == -1)
    {
        std::cout << "An Error is occured while reading FASTA file" << std::endl;
        bIsSuccess = false;
    }
    else if(a_rContig.m_nRefLength == -2)
    {
        std::cout << "Specified chromosome" << chromosome  << " could not found in FASTA file" << std::endl;
        bIsSuccess = false;
    }
    
    return bIsSuccess;
}














