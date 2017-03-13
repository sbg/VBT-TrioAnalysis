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


class CMendelianVariantProvider
{

public:
    
    //Initialize the vcf and fasta files for mendelian violation mode
    bool InitializeReaders(const SConfig &a_rFatherChildConfig, const SConfig& a_rMotherChildConfig);
    
    //Return all the variants belongs to given chromosome
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo) const;
    
    //Return all the variants belongs to given chromosome according to given index list
    std::vector<const CVariant*> GetVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_nIndexList) const;

    //Return all the variants belongs to given List according to given index list
    std::vector<const CVariant*> GetVariantList(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_nIndexList) const;

    //Return all the oriented variants belongs to given chromosome
    std::vector<const COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch = false) const;
    
    //Return all the oriented variants belongs to given chromosome with provided index list
    std::vector<const COrientedVariant*> GetOrientedVariantList(EMendelianVcfName a_uFrom, int a_nChrNo, bool a_bIsAlleleMatch, const std::vector<int>& a_nIndexList) const;
    
    //Return contig object given by the chromosome Id
    void GetContig(int a_nChrId, SContig& a_rContig) const;
    
    //Return a list of common chromosome ids found in all 3 vcf file
    std::vector<int> GetCommonChromosomes(bool a_bIsCalledInProviderInitialization = false);

    //Set the status of each variant in the given lust
    void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
    void SetVariantStatus(const std::vector<const COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
        
private:
    
    //Checks whether given variant is a structural variant type (A complex type)
    bool IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const;
    
    //Fill Variant sets for parent and child
    void FillVariants();
    
    //Fill Oriented variant sets for parent and child
    void FillGenotypeMatchOrientedVariants(std::vector<int>& a_aCommonChromosomes);

    //Fill Oriented variant sets for parent and child
    void FillAlleleMatchOrientedVariants(std::vector<int>& a_aCommonChromosomes);
    
    //Push the variant to the variantlist in the order of starting point (ascending order)
    void PushVariant(CVariant& a_rVariant, std::vector<CVariant>& a_rVecToPush);
    
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
    std::vector<CVariant> m_aFatherVariantList[CHROMOSOME_COUNT];
    //List that stores Mother Variants in order
    std::vector<CVariant> m_aMotherVariantList[CHROMOSOME_COUNT];
    //List that stores Child Variants in order
    std::vector<CVariant> m_aChildVariantList[CHROMOSOME_COUNT];
    
    
    //List that store the genotype match base Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aMotherOrientedVariantList[CHROMOSOME_COUNT];
    //List that store the genotype match called Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aFatherOrientedVariantList[CHROMOSOME_COUNT];
    //List that store the genotype match base Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aChildOrientedVariantList[CHROMOSOME_COUNT];

    
    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aMotherAlleleMatchOrientedVariantList[CHROMOSOME_COUNT];
    //List that store the allele match called Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aFatherAlleleMatchOrientedVariantList[CHROMOSOME_COUNT];
    //List that store the allele match base Oriented variant tuples (In the order of genotype)
    std::vector<COrientedVariant> m_aChildAlleleMatchOrientedVariantList[CHROMOSOME_COUNT];
    
    //List that stores Father variants which are filtered out from comparison
    std::vector<CVariant> m_aFatherNotAssessedVariantList[CHROMOSOME_COUNT];
    //List that stores Mother variants which are filtered out from comparison
    std::vector<CVariant> m_aMotherNotAssessedVariantList[CHROMOSOME_COUNT];
    //List that stores Child variants which are filtered out from comparison
    std::vector<CVariant> m_aChildNotAssessedVariantList[CHROMOSOME_COUNT];
    
    //Reference to the fasta reader object
    CFastaParser m_fastaParser;
    //Reference contig list
    SContig m_aContigList[CHROMOSOME_COUNT];
    
};

#endif //_C_MENDELIAN_VARIANT_PROVIDER_H_

