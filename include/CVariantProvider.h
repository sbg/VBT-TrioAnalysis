/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_VARIANT_PROVIDER_H_
#define _C_VARIANT_PROVIDER_H_

#include "EVcfName.h"
#include "CVcfReader.h"
#include "CPath.h"
#include "SConfig.h"
#include "CFastaReader.h"

const int CHROMOSOME_COUNT = 23;


class CVariantProvider
{
    public:
        
        //Initialize the VCF readers for base and called vcf file
        void InitializeReaders(const SConfig& a_rConfig);
    
        //Get the variant with given Chr number and variant id
        const CVariant* GetVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId) const;
    
        //Get the oriented variant with the given chr variant id and orientation
        COrientedVariant* GetOrientedVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId, bool a_bOrientation);
    
        // Returns the size of variant list for that chromosome
        int GetVariantListSize(EVcfName a_uFrom, int a_nChrNo) const;
    
        //Sets the fasta reader reference object
        void SetFastaReader(const CFastaReader& m_rFastaReader);
    
        //Return the variant list with the given index list
        std::vector<CVariant> GetVariantList(EVcfName a_uFrom, int a_nChrNo, std::vector<int> a_VariantIndexes);
    
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
    
        //List that store the base Oriented variant tuples (In the order of genotype)
        std::vector<COrientedVariant> m_aBaseOrientedVariantList[CHROMOSOME_COUNT];
        //List that store the called Oriented variant tuples (In the order of genotype)
        std::vector<COrientedVariant> m_aCalledOrientedVariantList[CHROMOSOME_COUNT];
    
        SConfig m_config;
    
        //Reference to the fasta reader object
        const CFastaReader* m_pFastaReader;
};


#endif // _C_VARIANT_PROVIDER_H_
