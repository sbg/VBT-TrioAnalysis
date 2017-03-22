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

struct SContig
{
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
    
private:
    
    faidx_t *fai;
};

#endif //_C_FASTA_PARSER_H_
