#include "CVariantProvider.h"
#include <iostream>

void CVariantProvider::InitializeReaders(const char* a_pBaseVcfFile, const char* a_pCalledVcfFile)
{
    bool bIsSuccess;

    bIsSuccess = m_baseVCF.Open(a_pBaseVcfFile);
    if(bIsSuccess)
        std::cout << "Successfully open Base: " << a_pBaseVcfFile << std::endl;
    m_baseVCF.setID(0);
    bIsSuccess = m_calledVCF.Open(a_pCalledVcfFile);
    if(bIsSuccess)
        std::cout << "Successfully open Caller: " << a_pCalledVcfFile << std::endl; 
    m_calledVCF.setID(1);
    
    FillVariantLists();
}


void CVariantProvider::FillVariantLists()
{
    CVariant variant;
    int id = 0;
    while(m_baseVCF.GetNextRecord(&variant, id++))
    {
        m_aBaseVariantList[variant.m_nChrId].push_back(variant);
    }
    
    while(m_calledVCF.GetNextRecord(&variant, id++))
    {
        m_aCalledVariantList[variant.m_nChrId].push_back(variant);
    }
}

const CVariant* CVariantProvider::GetVariant(EVcfName a_uFrom, int a_nChrNo, int a_nVariantId)
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
