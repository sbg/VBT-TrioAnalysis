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

class CFastaParser
{
    
public:
    
    bool OpenFastaFile(const char *fn, std::string a_chrom);
    
    bool FetchNewChromosome(std::string chrom);
    
    void DeleteReferenceMemory();
    
    char GetRefBase(int a_nPosition);
    
    std::string GetRefBase(int a_nRefPos, int a_nLength);
    
    int GetRefLength();
    
    bool GenerateFastaIndex(const char *fn);
    
private:
    
    faidx_t *fai;
    int ref_len;
    char *ref;
    //    std::string ref_sequence;
};

#endif //_C_FASTA_PARSER_H_
