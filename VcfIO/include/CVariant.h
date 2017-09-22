/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _CVariant_H
#define _CVariant_H

#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "EVariantMatch.h"
#include "EVariantCategory.h"


enum EVariantType
{
    eSNP,
    eINDEL,
    eSV
};

/**
 * @brief A Container that stores the necessary information of variant
 *
 * CVariant is a container for variants that stores necessary vcfrecord information read with htslib
 */
struct SAllele
{
    std::string m_sequence;
    int m_nStartPos = -1;
    int m_nEndPos = -1;
    bool m_bIsIgnored = false;
};

/**
 * @brief A Container that stores the necessary information of variant
 *
 * CVariant is a container for variants that stores necessary vcfrecord information read with htslib
 */
class CVariant
{
    public:
    ///Default constructor
    CVariant();
    
    ///Copy constructor
    CVariant(const CVariant& a_rObj);
    
    ///Clear the Variant
    bool Clear();
    
    ///Compare the given variant with this
    int CompareTo(const CVariant& a_rObj) const;
    
    ///Return if genotype is phased
    bool IsPhased() const;
    
    ///Gets the start position of the Variant
    int GetStart() const;
    
    ///Gets the end position of the ref allele
    int GetEnd() const;
    
    ///Gets the original position of the variant ignoring trimming and first nucleotide reduction
    int GetOriginalPos() const;
    
    ///Return is the variant is heterozygous
    bool IsHeterozygous() const;
    
    ///Print the Variant
    void Print() const;
    
    ///Check if variant is initialized
    bool IsNull() const;
    
    ///Return the unique id of the variant
    int GetId() const;
    
    ///Return true if the variant filter column is PASS
    bool IsFilterPASS() const;
    
    ///Return the original Alt string with the given index (Take back trim operation)
    std::string GetOriginalAlleleStr(int a_nAlleleIndex) const;
    
    ///Return the reference sequences
    std::string GetRefSeq() const;
    
    ///Return the allele sequence specified with the id (0 is first allele, 1 is second allele)
    SAllele GetAllele(int a_nAlleleId) const;

    ///Return the type of variant
    EVariantType GetVariantType() const;
    
    ///Return the category of variant to use in mendelian violation reporting
    EVariantCategory GetVariantCategory() const;
    
    ///Fill the Genotype list with the original genotype indexes and the genotype count
    void GetGenotypeArr(int* a_pGenotypeList, int& a_rGenotypeCount);
    
    ///Print the variant [For Test Purpose]
    std::string ToString() const;
    
    ///ID of which vcf file that the variant belongs to
    int m_nVcfId;

    ///Id of the chromosome that variant belogs to
    int m_nChrId;

    ///Unique Id of variant
    int m_nId;

    ///Allele count of the variant (2 for diploid and 1 for haploid)
    int m_nAlleleCount;
    
    ///Start Position of the variant - min start pos of all alleles
    int m_nStartPos;

    ///End Position of the variant - max end pos of all alleles
    int m_nEndPos;
    
    ///Original Genotype list read from vcf file
    int m_genotype[2];
    
    ///Zygot Count
    int m_nZygotCount;
    
    ///Original variant position
    int m_nOriginalPos;
    
    mutable EVariantMatch m_variantStatus;
    
    ///True if the variant genotype is phased
    bool m_bIsPhased;
    
    ///If the variant is heterozygous (ie. GT is 0/1 1/2 etc.)
    bool m_bIsHeterozygous;
    
    ///True if the first nucleotide is trimmed
    bool m_bIsFirstNucleotideTrimmed;
    
    ///True if genotype is ./.
    bool m_bIsNoCall;
    
    ///Allele array of the variant
    SAllele m_alleles[2];

    ///Reference sequence
    std::string m_refSequence;
    
    ///Filter Data
    bool m_bIsFilterPASS;
    std::vector<std::string> m_filterString;
    
    ///Chromosome name
    std::string m_chrName;
    
    ///Original Alleles string read from vcf file
    std::string m_allelesStr;
    
};


#endif // _CVariant_H
