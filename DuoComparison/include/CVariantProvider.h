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
#include "SChrIdTuple.h"
#include "CBaseVariantProvider.h"

namespace core
{
    class COrientedVariant;
}

namespace duocomparison
{

/**
 * @brief Reads and stores variants of base/called vcf files
 *
 * CVariantProvider contains functions to parse vcf files of truth and query. All variants are stored in this object and
 * other classes can get access variant lists via this class.
 */
class CVariantProvider : public CBaseVariantProvider
{
public:
    ///Destructor
    ~CVariantProvider();
    ///Constructor
    CVariantProvider();

    ///Initialize the VCF readers for base and called vcf file
    bool InitializeReaders(const SConfig& a_rConfig);

    ///Return the variant list with the given index list
    std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes);

    ///Return the all the variants belongs to given chromosome
    std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo);

    ///Return the list of variants from varlist according to given variant indexes
    std::vector<const CVariant*> GetVariantList(std::vector<const CVariant*>& a_varList, const std::vector<int>& a_VariantIndexes);

    ///Return the variants that are not assessed due to being in a complex regions during variant comparison operation
    std::vector<const CVariant*> GetSkippedComplexVariantList(EVcfName a_uFrom, int a_nChrNo);
    
    ///Return all the oriented variants belongs to given chromosome
    std::vector<const core::COrientedVariant*> GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, bool a_bIsGenotypeMatch);

    ///Return the index tuples of chromosomes which both contained by baseline and called VCF
    std::vector<SChrIdTuple>& GetChromosomeIdTuples();

    ///Return contig information from header of query vcf
    const std::vector<SVcfContig>& GetContigs() const;

    ///Read the header of Called vcf and return the filter names and descriptions
    void GetFilterInfo(EVcfName a_vcfType, std::vector<std::string>& a_rFilterNames, std::vector<std::string>& a_rFilterDescriptions);

    ///Return the access of not-asessed variants
    std::vector<const CVariant*> GetNotAssessedVariantList(EVcfName a_uFrom, int a_nChrNo);

    ///Initialize Homozygous Oriented Variant Lists with given base and called variant set
    void FillAlleleMatchVariantList(SChrIdTuple& a_rTuple,
                                    std::vector<const CVariant*>& a_rBaseVariants,
                                    std::vector<const CVariant*>& a_rCalledVariants);

private:

    //Finds the tuple index list of chromosome which is contained by both baseline and called vcf
    void SetChromosomeIdTuples();

    void FillVariantForSample(int a_nSampleId, SConfig& a_rConfig);
    
    //Read through the variant file and fill the variant lists. It assumes that positions are sorted.
    void FillVariantLists();
    
    //Read through the variant lists and generate oriented variant list for call and base
    void FillOrientedVariantLists();
    
    //Find the optimal trimmings for given variant list
    void FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom);
    
    //Merge trimmed variants with the original variant list
    void AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom);

    //VCF Readers
    CVcfReader m_baseVCF;
    CVcfReader m_calledVCF;

    //List that stores base Variants in order
    std::vector<std::vector<CVariant>> m_aBaseVariantList;
    //List that stores called Variants in order
    std::vector<std::vector<CVariant>> m_aCalledVariantList;

    //THESE TWO LIST STORES VARIANT ORIENTATIONS WHILE PERFORMING GT MATCHING
    //List that store the base Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aBaseOrientedVariantList;
    //List that store the called Oriented variant tuples (In the order of genotype)
    std::vector<std::vector<core::COrientedVariant>> m_aCalledOrientedVariantList;

    //THESE TWO LIST STORES WHILE VARIANT ORIENTATIONS PERFORMING ALLELE MATCH OPERATION AFTER GT MATCHING
    //List that store the base Oriented variant tuples (Allele match homozygous variants)
    bool m_bIsHomozygousOvarListInitialized;
    std::vector<std::vector<core::COrientedVariant>> m_aBaseHomozygousOrientedVariantList;
    //List that store the called Oriented variant tuples (Allele match homozygous variants)
    std::vector<std::vector<core::COrientedVariant>> m_aCalledHomozygousOrientedVariantList;

    //List that stores baseline variants which are filtered out from comparison
    std::vector<std::vector<CVariant>> m_aBaseNotAssessedVariantList;
    //List that stores called variants which are filtered out from comparison
    std::vector<std::vector<CVariant>> m_aCalledNotAssessedVariantList;

    //Parameters come from command line arguments
    SConfig m_config;

    //Chromosome id tuples for each common chromosome
    std::vector<SChrIdTuple> m_aCommonChrTupleList;
};

}

#endif // _C_VARIANT_PROVIDER_H_
