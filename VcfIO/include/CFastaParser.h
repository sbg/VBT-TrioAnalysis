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
 *  CFastaParser.h
 *  VariantBenchmarkingTools
 *
 *  Created by John Browning on 4/10/14.
 *
 */

#ifndef _C_FASTA_PARSER_H_
#define _C_FASTA_PARSER_H_

#include "faidx.h"
#include <string>

/**
 * @brief Container to store each FASTA contig
 *
 * SContig is used to store the reference sequence belongs to the given contig identified by the chromosome name
 */
struct SContig
{
    //Free the reference sequence
    bool Clean();
    
    //Chromosome name written on fasta/vcf
    std::string  m_chromosomeName;
    char* m_pRefSeq = 0;
    int m_nRefLength;
};


/**
 * @brief Read the provided Reference (FASTA) file
 *
 * CFastaParser is used to parse given reference files. Each contig present in FASTA file is accesible via FetchNewChromosome function
 */
class CFastaParser
{
    
public:
    
    ///Default Constructor
    CFastaParser();
    
    ///Destructor
    ~CFastaParser();
    
    ///Open Given FASTA file with the given filename and creates FASTA index file if it does not exists
    bool OpenFastaFile(const char *fn);
    
    ///Read contig from FASTA file name with the given chromosome name
    bool FetchNewChromosome(std::string chromosome, SContig& a_rContig);
    
    ///Generate FASTA index file from given fasta file if it does not already exists
    bool GenerateFastaIndex(const char *fn);
    
private:
    
    faidx_t *fai;
};

#endif //_C_FASTA_PARSER_H_
