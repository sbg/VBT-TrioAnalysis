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
 *  Created by Berke Cagkan Toptas
 *
 */

#ifndef _CVariant_H
#define _CVariant_H

#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "EVariantMatch.h"
#include "EMendelianDecision.h"

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
    bool m_bIsTrimmed = false;
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
    std::string GetOriginalAlleleStr(unsigned int a_nAlleleIndex) const;
    
    ///Return the reference sequences
    std::string GetRefSeq() const;
    
    ///Return the allele sequence specified with the id (0 is first allele, 1 is second allele)
    SAllele GetAllele(int a_nAlleleId) const;
    
    ///Return the type of variant
    EVariantType GetVariantType() const;
    
    ///Fill the Genotype list with the original genotype indexes and the genotype count
    void GetGenotypeArr(int* a_pGenotypeList, int& a_rGenotypeCount);
    
    ///Print the variant [For Test Purpose]
    std::string ToString() const;
    
    ///Gets the maximum number of nucleotides can be trimmed from beginning and ending of selected allele
    void GetMaxTrimStartEnd(int a_nAlleleIndex, unsigned int& trimLengthFromBeginning, unsigned int& trimLengthFromEnd);
    
    ///Trim the redundant nucleotides from beginning and ending of given allele (If allele can be trim multiple way, use the second parameter for order)
    void TrimVariant(int a_nAlleleIndex, bool a_bIsBeginFirst);
    
    ///Trim the redundant nucleotides from beginning and ending of given allele. Nucleotides to be clipped from beginning and ending of the allele are specified as parameter
    void TrimVariant(int a_nAlleleIndex, unsigned int trimLengthFromBeginning, unsigned int trimLengthFromEnd);
    
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
    
    mutable EMendelianDecision m_mendelianDecision;
    
    ///True if the variant genotype is phased
    bool m_bIsPhased;
    
    ///If the variant is heterozygous (ie. GT is 0/1 1/2 etc.)
    bool m_bIsHeterozygous;
    
    ///True if the first nucleotide is trimmed
    bool m_bIsFirstNucleotideTrimmed;
    
    ///True if genotype is ./.
    bool m_bIsNoCall;
    
    ///If the variant can be trimmed more than 1 way (for -ref-overlap mode)
    bool m_bHaveMultipleTrimOption;
    
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
 
private:
    
    ///Trim the redundant nucleotides from beginning and ending of given allele (Beginning nucleotides will be trimmed first)
    void TrimVariantBeginFirst(int a_nAlleleIndex);
    
    ///Trim the redundant nucleotides from beginning and ending of given allele (Ending nucleotides will be trimmed first)
    void TrimVariantEndFirst(int a_nAlleleIndex);

    
};


#endif // _CVariant_H
