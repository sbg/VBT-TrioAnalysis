//
//  CFastaParser.h
//
//  Created by John Browning on 4/10/14.
//  Copyright (c) 2014 John Browning. All rights reserved.
//
//  Modified by: Berke Cagkan Toptas
//
//  Header file, opens a FASTA file and stores it in a string. It also contains a subroutine
//  to fetch a base given the index.
//

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
