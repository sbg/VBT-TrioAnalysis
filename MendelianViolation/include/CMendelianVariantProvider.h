//
//  CMendelianVariantProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#ifndef _C_MENDELIAN_VARIANT_PROVIDER_H_
#define _C_MENDELIAN_VARIANT_PROVIDER_H_

#include "CVcfReader.h"
#include "COrientedVariant.h"
#include "Constants.h"
#include "EMendelianVcfName.h"
#include "SChrIdTriplet.h"
#include "CBaseVariantProvider.h"

namespace mendelian
{

/**
 * @brief Reads and stores variants of family members
 *
 * CMendelianVariantProvider contains functions to parse vcf files of family members. All variants are stored in this object and
 * other classes can get access variant lists via this class.
 */
class CMendelianVariantProvider : public CBaseVariantProvider
{

public:
    
    ///Initialize the vcf and fasta files for mendelian violation mode
    bool InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig);
    
    ///Return all the variants belongs to given chromosome
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const;

    ///Return all the variants belongs to given chromosome sorted by variant ids
    std::vector<const CVariant*> GetSortedVariantListByID(EMendelianVcfName a_uFrom, int a_nChrNo) const;
    
    ///Return all the variants belongs to given chromosome sorted by variant ids
    std::vector<const CVariant*> GetSortedVariantListByIDandStartPos(EMendelianVcfName a_uFrom, int a_nChrNo) const;

    ///Return all the variants belongs to given chromosome according to given index list
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_nIndexList) const;

    ///Return all the variants belongs to given List according to given index list
    std::vector<const CVariant*> GetVariantList(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_nIndexList) const;

    ///Return the count of not assessed variants in vcf file - variants that contains any of *, <, >, [, ], {, } symbols at their allele string
    int GetNotAssessedVariantCount(EMendelianVcfName a_uFrom);
    
    ///Return contig information from header of child vcf
    const std::vector<SVcfContig>& GetContigs() const;
    
    ///Return the total contig count of requested vcf
    int GetContigCount(EMendelianVcfName a_uFrom);
    
    ///Get the total variant count for given chromosome
    int GetVariantCount(EMendelianVcfName a_uFrom, int a_nChrNo) const;
    
    ///Return all the oriented variants belongs to given chromosome
    std::vector<const core::COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch = false) const;
    
    ///Return all the oriented variants belongs to given chromosome with provided index list
    std::vector<const core::COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom,
                                                                      int a_nChrNo,
                                                                      bool a_bIsAlleleMatch,
                                                                      const std::vector<int>& a_nIndexList) const;
        
    //Return a list of common chromosome id triplets found in all 3 vcf file
    std::vector<SChrIdTriplet>& GetCommonChromosomes();

    //Set the status of each variant in the given list
    void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
    void SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
        
    //Returns the count of eNOT_ASSESSED variants for all chromosome belong to input VCF sample
    int GetSkippedVariantCount(EMendelianVcfName a_uFrom) const;
        
private:
    
    //Fill the common chromosome list
    void SetCommonChromosomes();
    
    //Fill Variants for given sample Id
    void FillVariantForSample(int a_nSampleId, SConfig& a_rConfig);
    
    //Fill Variant sets for parent and child
    void FillVariants();
        
    //Fill Oriented variant sets for parent and child
    void FillGenotypeMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes);

    //Fill Oriented variant sets for parent and child
    void FillAlleleMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes);
        
    //Find the optimal trimmings for given variant list
    void FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EMendelianVcfName a_uFrom);
    
    //Merge trimmed variants with the original variant list
    void AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EMendelianVcfName a_uFrom);

    
    //VCF FILES FOR COMPARISON
    CVcfReader m_FatherVcf;
    CVcfReader m_MotherVcf;
    CVcfReader m_ChildVcf;
    
    //Config objects for variant provider
    SConfig m_motherChildConfig;
    SConfig m_fatherChildConfig;
    
    //List that stores Father Variants in order
    std::vector<std::vector<CVariant>> m_aFatherVariantList;
    //List that stores Mother Variants in order
    std::vector<std::vector<CVariant>> m_aMotherVariantList;
    //List that stores Child Variants in order
    std::vector<std::vector<CVariant>> m_aChildVariantList;
    
    //List that store the genotype match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aMotherOrientedVariantList;
    //List that store the genotype match called Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aFatherOrientedVariantList;
    //List that store the genotype match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aChildOrientedVariantList;

    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aMotherAlleleMatchOrientedVariantList;
    //List that store the allele match called Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aFatherAlleleMatchOrientedVariantList;
    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aChildAlleleMatchOrientedVariantList;
    
    //Father variants which are filtered out from comparison
    int m_nFatherNotAssessedVariantCount;
    //Mother variants which are filtered out from comparison
    int m_nMotherNotAssessedVariantCount;
    //Child variants which are filtered out from comparison
    int m_nChildNotAssessedVariantCount;
    
    //Chromosomes that all three samples shares variant
    std::vector<SChrIdTriplet> m_aCommonChromosomes;
    
    //Variant counts contains asterisk which will are eliminated from comparison
    int m_nMotherAsteriskCount;
    int m_nFatherAsteriskCount;
    int m_nChildAsteriskCount;
};

}

#endif //_C_MENDELIAN_VARIANT_PROVIDER_H_

