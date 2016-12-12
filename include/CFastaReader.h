/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
*/

#ifndef FASTA_READER_H_
#define FASTA_READER_H_

#include <stdlib.h>
#include <string>
#include <vector>
#include <zlib.h>
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread);

class CFastaReader
{
    public:
    CFastaReader();
    CFastaReader(const char * a_pFileName);
    ~CFastaReader();
        
    // Open a FASTA file and init Contig
    bool Open(const char * a_pFileName);
    // Close a FASTA file and free allocated Contig
    bool Close();            
    // Read the next Contig in FASTA file         
    void ReadContig();
    //Print the last read Contig in FASTA file
    void PrintContig();
    //Returns a pointer to the sequence array of last read Contig //TODO:Should take Chromosome id as input
    char* GetRefSeq() const;
    //Returns the size of ref sequence array of last read Contig //TODO:Should take Chromosome id as input
    int GetRefSeqSize() const;
    
    //Fills the m_contigList for the provided chromosomes from FASTA file
    void FillRefSeq(std::vector<std::string> a_chromosomeList);
    
    private:

    std::string m_fileName;
    //If FASTA file is open
    bool m_bIsOpen;
    gzFile m_fastaFile;
    kseq_t *m_pContig;
    
    std::vector<kseq_t*> m_contigList;
    
    int m_LastContigSequenceSize;
    
};





#endif // FASTA_READER_H_
