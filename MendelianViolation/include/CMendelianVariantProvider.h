//
//  CMendelianVariantProvider.h
//  VCFComparison
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_MENDELIAN_VARIANT_PROVIDER_H_
#define _C_MENDELIAN_VARIANT_PROVIDER_H_

#include "CVcfReader.h"
#include "CFastaParser.h"
#include "COrientedVariant.h"
#include "Constants.h"
#include "EMendelianVcfName.h"
#include "SChrIdTriplet.h"

class CMendelianVariantProvider
{

public:
    
    //Initialize the vcf and fasta files for mendelian violation mode
    bool InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig);
    
    //Return all the variants belongs to given chromosome
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const;

    //Return all the variants belongs to given chromosome sorted by variant ids
    std::vector<const CVariant*> GetSortedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const;

    //Return all the variants belongs to given chromosome according to given index list
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_nIndexList) const;

    //Return all the variants belongs to given List according to given index list
    std::vector<const CVariant*> GetVariantList(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_nIndexList) const;

    //Return the count of not assessed variants in vcf file - variants that contains any of *, <, >, [, ], {, } symbols at their allele string
    int GetNotAssessedVariantCount(EMendelianVcfName a_uFrom);
    
    //Return contig information from header of child vcf
    const std::vector<SVcfContig>& GetContigs();
    
    //Return the total contig count of requested vcf
    int GetContigCount(EMendelianVcfName a_uFrom);
    
    //Get the total variant count for given chromosome
    int GetVariantCount(EMendelianVcfName a_uFrom, int a_nChrNo) const;
    
    //Return all the oriented variants belongs to given chromosome
    std::vector<const COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch = false) const;
    
    //Return all the oriented variants belongs to given chromosome with provided index list
    std::vector<const COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch, const std::vector<int>& a_nIndexList) const;
    
    //Return contig object given by the chromosome Id
    void GetContig(int a_nChrId, SContig& a_rContig) const;
    
    //Read contig given by the chromosome id
    void ReadContig(std::string a_chrId, SContig& a_rContig);
    
    //Return a list of common chromosome id triplets found in all 3 vcf file
    std::vector<SChrIdTriplet>& GetCommonChromosomes();

    //Set the status of each variant in the given lust
    void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
    void SetVariantStatus(const std::vector<const COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
    
    //Returns the 0 based index of variants starting from the first variant for that chromosome
    int Get0BasedVariantIndex(EMendelianVcfName a_uFrom, int a_nChr, int a_nVariantId) const;
    
    //Returns the count of eNOT_ASSESSED variants for all chromosome belong to input VCF sample
    int GetSkippedVariantCount(EMendelianVcfName a_uFrom) const;
        
private:
    
    //Fill the common chromosome list
    void SetCommonChromosomes();
    
    //Checks whether given variant is a structural variant type (A complex type)
    bool IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const;
    
    //Checks if the variant genotype is 0/0 (homref)
    bool IsHomRef(const CVariant& a_rVariant) const;
    
    //Fill Variant sets for parent and child
    void FillVariants();
    
    //Fill Oriented variant sets for parent and child
    void FillGenotypeMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes);

    //Fill Oriented variant sets for parent and child
    void FillAlleleMatchOrientedVariants(std::vector<SChrIdTriplet>& a_aCommonChromosomes);
    
    static bool CompareVariants(const CVariant& var1, const CVariant& var2);
    
    //VCF FILES FOR COMPARISON
    CVcfReader m_FatherVcf;
    CVcfReader m_MotherVcf;
    CVcfReader m_ChildVcf;
    
    //REFERENCE FASTA
    CFastaParser m_referenceFasta;
    
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
    std::vector<std::vector<COrientedVariant>> m_aMotherOrientedVariantList;
    //List that store the genotype match called Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<COrientedVariant>> m_aFatherOrientedVariantList;
    //List that store the genotype match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<COrientedVariant>> m_aChildOrientedVariantList;

    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<COrientedVariant>> m_aMotherAlleleMatchOrientedVariantList;
    //List that store the allele match called Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<COrientedVariant>> m_aFatherAlleleMatchOrientedVariantList;
    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<COrientedVariant>> m_aChildAlleleMatchOrientedVariantList;
    
    //List that stores Father variants which are filtered out from comparison
    std::vector<std::vector<CVariant>> m_aFatherNotAssessedVariantList;
    //List that stores Mother variants which are filtered out from comparison
    std::vector<std::vector<CVariant>> m_aMotherNotAssessedVariantList;
    //List that stores Child variants which are filtered out from comparison
    std::vector<std::vector<CVariant>> m_aChildNotAssessedVariantList;
    
    //Reference to the fasta reader object
    CFastaParser m_fastaParser;
    std::vector<SChrIdTriplet> m_aCommonChromosomes;
    
};

#endif //_C_MENDELIAN_VARIANT_PROVIDER_H_

