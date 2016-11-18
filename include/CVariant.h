/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _CVariant_H
#define _CVariant_H

#include <string>
#include <vector>
#include "htslib/vcf.h"
#include <iostream>


enum EVariantType
{
    eSNP,
    eINDEL_ADD,
    eINDEL_DEL,
    eNO_OP
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
    
    //Returns the maximum sequance size
    int GetMaxLength() const;

    //Return the reference sequences
    std::string GetRefSeq() const;
    //Return the allele sequence specified with the id (0 is first allele, 1 is second allele)
    std::string GetAllele(int a_nAlleleId) const;

    // Detects the type of of each alt with index 
    void SetType(int a_nAltIndex);    

    // Print the variant [For Test Purpose]
    std::string ToString() const;
    
    //ID of which vcf file that the variant belongs to
    int m_nVcfId;
    int m_nChrId;
    int m_nPosition;
    //Chromosome name
    std::string m_chrName;
    //Sequence array. m_aSequences[0] is the ref string
    std::vector<std::string> m_aSequences;
    
    //Genotype Data
    int gt_arr[2]; //Haplotype array
    int ngt_arr; //Haplotype count
    
    //Filter Data
    bool m_bIsFilterPASS;
    
    //Type of the variant array for each alt (SNP or INDEL)
    EVariantType m_aVarTypes[2];
    //Return true if the variant genotype is phased
    bool m_bIsPhased;
    //Unique Id of variant
    int m_nId;
    
    //Start Position of the variant
    int m_nStartPos;
    //End Position of the variant
    int m_nEndPos;

};


#endif // _CVariant_H
