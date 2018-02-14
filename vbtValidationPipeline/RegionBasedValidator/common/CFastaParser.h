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
 */

#ifndef _C_FASTA_PARSER_H_
#define _C_FASTA_PARSER_H_

#include "htslib/faidx.h"
#include <string>

struct SContig
{
    bool Clean();
    
    //Chromosome name written on fasta/vcf
    std::string m_chromosome;
    //Chromosome index in variant provider
    int m_nChrId;
    char* m_pRefSeq = 0;
    int m_nRefLength;
};

class CFastaParser
{
    
public:
    
    //Destructor
    ~CFastaParser();
    
    //Open Given FASTA file with the given filename and creates FASTA index file if it does not exists
    bool OpenFastaFile(const char *fn);
    
    //Read contig from FASTA file name with the given chromosome name
    bool FetchNewChromosome(std::string chromosome, SContig& a_rContig);
    
    //Generate FASTA index file from given fasta file if it does not already exists
    bool GenerateFastaIndex(const char *fn);
    
    //Get the Subsequence from given chromosome and interval
    std::string GetSubSequence(const SContig& a_rContig, int a_nStartPos, int a_nEndPos) const;
    
    void PrintSubSequence(int a_nChrId, int a_nStartPos, int a_nLength);
    
private:
    
    faidx_t *fai;
};

#endif //_C_FASTA_PARSER_H_
