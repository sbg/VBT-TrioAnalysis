/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _CVariant_H
#define _CVariant_H

#include <string>
#include "htslib/vcf.h"

enum EVariantType
{
    eSNP,
    eINDEL_ADD,
    eINDEL_DEL,
    eNO_OP
};

struct SAllele
{
    std::string m_sequence;
    int m_nStartPos = -1;
    int m_nEndPos = -1;
};

class CVariant
{
    public:
    //Default constructor
    CVariant();
    
    //Copy constructor
    CVariant(const CVariant& a_rObj);
    
    //Clear the Variant
    bool Clear();
    
    //Compare the given variant with this
    int CompareTo(const CVariant& a_rObj) const;
    
    //Return if genotype is phased
    bool IsPhased() const;
    
    //Gets the start position of the Variant
    int GetStart() const;
    
    //Gets the end position of the ref allele
    int GetEnd() const;
    
    // Return is the variant is heterozygous
    bool IsHeterozygous() const;
    
    // Print the Variant
    void Print() const;
    
    //Check if variant is initialized
    bool IsNull() const;
    
    //Return the unique id of the variant
    int GetId() const;
    
    //Return true if the variant filter column is PASS
    bool IsFilterPASS() const;
    
    //Return the reference sequences
    std::string GetRefSeq() const;
    
    //Return the allele sequence specified with the id (0 is first allele, 1 is second allele)
    SAllele GetAllele(int a_nAlleleId) const;

    // Print the variant [For Test Purpose]
    std::string ToString() const;
    
    //ID of which vcf file that the variant belongs to
    int m_nVcfId;
    //Id of the chromosome that variant belogs to
    int m_nChrId;
    //Chromosome name
    std::string m_chrName;
    //Unique Id of variant
    int m_nId;
 
    //Filter Data
    bool m_bIsFilterPASS;
    
    //True if the variant genotype is phased
    bool m_bIsPhased;
    bool m_bIsHeterozygous;

    //Allele array of the variant
    SAllele m_alleles[2];
    //Allele count of the variant (2 for diploid and 1 for haploid)
    int m_nAlleleCount;
    //Reference sequence
    std::string m_refSequence;
    
    //Start Position of the variant - min start pos of all alleles
    int m_nStartPos;
    //End Position of the variant - max end pos of all alleles
    int m_nEndPos;
    
    //Maximum length of the allele
    int m_nMaxLength;
};


#endif // _CVariant_H
