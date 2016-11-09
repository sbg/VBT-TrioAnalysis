#ifndef _C_REF_SEQ_PROVIDER_H_
#define _C_REF_SEQ_PROVIDER_H_

#include "CFastaReader.h"

class CRefSeqProvider
{
    public:
    CRefSeqProvider(const char* a_pFastaName);
    //Gets the reference nucleotide sequence from the FASTA file specified with the chromosome id
    const char* GetRefSeq(int a_nChromosomeId);

    private:
    //Fasta Reader instance
    CFastaReader m_fastaReader;
};

#endif //_C_REF_SEQ_PROVIDER_H_