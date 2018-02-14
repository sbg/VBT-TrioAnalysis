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
 *  CFastaParser.cpp
 *
 *  Created by John Browning on 4/10/14.
 *  Modified by: Berke Cagkan Toptas
 *
 */

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

CFastaParser::CFastaParser()
{
    fai = NULL;
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
    a_rContig.m_chromosomeName = chromosome;
    
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

bool SContig::Clean()
{
    bool bIsSuccess = true;
    if(m_pRefSeq != NULL)
        free(m_pRefSeq);
    else
        bIsSuccess = false;
    
    m_pRefSeq = NULL;
    return bIsSuccess;
}














