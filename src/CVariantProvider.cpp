#include "CVariantProvider.h"
#include <iostream>

void CVariantProvider::InitializeReaders(const SConfig& a_rConfig)
{
    bool bIsSuccess;

    m_bIsFilterPASS = a_rConfig.m_bIsFilterPASS;
    
    bIsSuccess = m_baseVCF.Open(a_rConfig.m_pBaseVcfFileName);
    if(!bIsSuccess)
        std::cout << "Baseline VCF file is unable to open!: " << a_rConfig.m_pBaseVcfFileName << std::endl;
    m_baseVCF.setID(0);
    bIsSuccess = m_calledVCF.Open(a_rConfig.m_pCalledVcfFileName);
    if(!bIsSuccess)
        std::cout << "Called VCF file is unable to open!: " << a_rConfig.m_pCalledVcfFileName << std::endl;
    m_calledVCF.setID(1);
    

    FillVariantLists();
}


void CVariantProvider::FillVariantLists()
{
    CVariant variant;
    int id = 0;
    while(m_baseVCF.GetNextRecord(&variant, id++))
    {
       // if(m_bIsFilterPASS == true && !variant.IsFilterPASS())
       //    continue;
           
        m_aBaseVariantList[variant.m_nChrId].push_back(variant);
    }
    
    while(m_calledVCF.GetNextRecord(&variant, id++))
    {
      //  if(m_bIsFilterPASS == true && !variant.IsFilterPASS())
      //     continue;
           
        m_aCalledVariantList[variant.m_nChrId].push_back(variant);
    }
}

const CVariant* CVariantProvider::GetVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId) const
{
    switch(a_uFrom)
    {
        case eBASE:
            return &m_aBaseVariantList[a_nChrNo][a_nVariantId];
        case eCALLED:
            return &m_aCalledVariantList[a_nChrNo][a_nVariantId];
    }
}

int CVariantProvider::GetVariantListSize(EVcfName a_uFrom, int a_nChrNo) const
{
    switch(a_uFrom)
    {
        case eBASE:
            return (int)m_aBaseVariantList[a_nChrNo].size();
        case eCALLED:
            return (int)m_aCalledVariantList[a_nChrNo].size();
    }
}
