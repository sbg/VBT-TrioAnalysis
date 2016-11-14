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

const int CHROMOSOME_COUNT = 23;


class CVariantProvider
{
    public:
        
        //Initialize the VCF readers for base and called vcf file
        void InitializeReaders(const SConfig& a_rConfig);
    
        //Read through the variant file and fill the variant lists. It assumes that positions are sorted.
        void FillVariantLists();
    
        //Get the variant with given Chr number and variant id
        const CVariant* GetVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId) const;
    
        // Returns the size of variant list for that chromosome
        int GetVariantListSize(EVcfName a_uFrom, int a_nChrNo) const;
        
    private:
        //VCF Readers
        CVcfReader m_baseVCF;
        CVcfReader m_calledVCF;

        //List that stores base Variants in order
        std::vector<CVariant> m_aBaseVariantList[CHROMOSOME_COUNT];
        //List that stores called Variants in order
        std::vector<CVariant> m_aCalledVariantList[CHROMOSOME_COUNT];
    
        bool m_bIsFilterPASS;
};


#endif // _C_VARIANT_PROVIDER_H_
