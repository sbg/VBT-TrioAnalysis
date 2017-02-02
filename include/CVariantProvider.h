/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_VARIANT_PROVIDER_H_
#define _C_VARIANT_PROVIDER_H_

#include "EVcfName.h"
#include "CVcfReader.h"
#include "SConfig.h"
#include "CFastaParser.h"
#include "Constants.h"

class COrientedVariant;

class CVariantProvider
{
    public:
        //Destructor
        ~CVariantProvider();
    
        //Initialize the VCF readers for base and called vcf file
        bool InitializeReaders(const SConfig& a_rConfig);
    
        //Return the variant list with the given index list
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes);
    
        //Return the all the variants belongs to given chromosome
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo);
    
        //Return the list of variants from varlist according to given variant indexes
        std::vector<const CVariant*> GetVariantList(std::vector<const CVariant*>& a_varList, const std::vector<int>& a_VariantIndexes);
    
        //Return all the oriented variants belongs to given chromosome
        std::vector<const COrientedVariant*> GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, bool a_bIsGenotypeMatch);
    
        //Return contig object given by the chromosome Id
        void GetContig(int a_nChrId, SContig& a_rContig) const;
    
        //Return the count of chromosome which both contained by baseline and called VCF
        void GetUniqueChromosomeIds(std::vector<int>& a_rChrIds);
    
        //Intersect base and called VCF files to find common chromosomes
        std::vector<SVcfContig> GetCommonChromosomes();
    
        //Read the header of Called vcf and return the filter names and descriptions
        void GetFilterInfo(EVcfName a_vcfType, std::vector<std::string>& a_rFilterNames, std::vector<std::string>& a_rFilterDescriptions);

        //Return the access of not-asessed variants
        std::vector<CVariant>& GetNotAssessedVariantList(EVcfName a_uFrom, int a_nChrNo);
    
        //Initialize Homozygous Oriented Variant Lists with given base and called variant set
        void FillAlleleMatchVariantList(int a_nChrId, std::vector<const CVariant*>& a_rBaseVariants, std::vector<const CVariant*>& a_rCalledVariants);
    
        //Set the status of each variant in the given lust
        void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
        void SetVariantStatus(const std::vector<const COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
    
    
    private:
    
        //Read through the variant file and fill the variant lists. It assumes that positions are sorted.
        void FillVariantLists();
    
        //Read through the variant lists and generate oriented variant list for call and base
        void FillOrientedVariantLists();

        //Checks whether given variant is a structural variant type (A complex type)
        bool IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const;
    
        //Push the variant to the variantlist in the order of starting point (ascending order)
        void PushVariant(CVariant& a_rVariant, std::vector<CVariant>& a_rVecToPush);
    
        //VCF Readers
        CVcfReader m_baseVCF;
        CVcfReader m_calledVCF;

        //List that stores base Variants in order
        std::vector<CVariant> m_aBaseVariantList[CHROMOSOME_COUNT];
        //List that stores called Variants in order
        std::vector<CVariant> m_aCalledVariantList[CHROMOSOME_COUNT];
    
        //THESE TWO LIST STORES VARIANT ORIENTATIONS WHILE PERFORMING GT MATCHING
        //List that store the base Oriented variant tuples (In the order of genotype)
        std::vector<COrientedVariant> m_aBaseOrientedVariantList[CHROMOSOME_COUNT];
        //List that store the called Oriented variant tuples (In the order of genotype)
        std::vector<COrientedVariant> m_aCalledOrientedVariantList[CHROMOSOME_COUNT];
    
        //THESE TWO LIST STORES WHILE VARIANT ORIENTATIONS PERFORMING ALLELE MATCH OPERATION AFTER GT MATCHING
        //List that store the base Oriented variant tuples (Allele match homozygous variants)
        std::vector<COrientedVariant> m_aBaseHomozygousOrientedVariantList[CHROMOSOME_COUNT];
        //List that store the called Oriented variant tuples (Allele match homozygous variants)
        std::vector<COrientedVariant> m_aCalledHomozygousOrientedVariantList[CHROMOSOME_COUNT];
    
        //List that stores baseline variants which are filtered out from comparison
        std::vector<CVariant> m_aBaseNotAssessedVariantList[CHROMOSOME_COUNT];
        //List that stores called variants which are filtered out from comparison
        std::vector<CVariant> m_aCalledNotAssessedVariantList[CHROMOSOME_COUNT];
    
        //Parameters come from command line arguments
        SConfig m_config;
    
        //Reference to the fasta reader object
        CFastaParser m_fastaParser;
        //Reference contig list
        SContig m_aContigList[CHROMOSOME_COUNT];
    
};


#endif // _C_VARIANT_PROVIDER_H_
