#include "CVariantProvider.h"
#include <iostream>

void CVariantProvider::InitializeReaders(const SConfig& a_rConfig)
{
    bool bIsSuccess;

    m_config = a_rConfig;
    
    bIsSuccess = m_baseVCF.Open(a_rConfig.m_pBaseVcfFileName);
    if(!bIsSuccess)
        std::cout << "Baseline VCF file is unable to open!: " << a_rConfig.m_pBaseVcfFileName << std::endl;
    m_baseVCF.setID(0);
    bIsSuccess = m_calledVCF.Open(a_rConfig.m_pCalledVcfFileName);
    if(!bIsSuccess)
        std::cout << "Called VCF file is unable to open!: " << a_rConfig.m_pCalledVcfFileName << std::endl;
    m_calledVCF.setID(1);

    FillVariantLists();
//    int k = 0;
//    for(CVariant variant : m_aBaseVariantList[21])
//          std::cout << k++ << ": " <<  variant.ToString() << std::endl;
//    
    
}

void CVariantProvider::SetFastaReader(const CFastaReader& a_rFastaReader)
{
    m_pFastaReader = &a_rFastaReader;
}


void CVariantProvider::FillVariantLists()
{
    CVariant variant;
    int id = 0;
    while(m_baseVCF.GetNextRecord(&variant, id++, m_config))
    {
        //if(id == 10000)
        //    break;
        
        if(variant.m_bIsFilterPASS == false)
            continue;
        
        if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            continue;
        
        else if(m_pFastaReader->GetRefSeqSize() < variant.GetEnd())
            continue;
        
        m_aBaseVariantList[variant.m_nChrId].push_back(variant);
    }
    
    while(m_calledVCF.GetNextRecord(&variant, id++, m_config))
    {
        //if(id == 20000)
        //    break;

        if(variant.m_bIsFilterPASS == false)
            continue;
        
        if(IsStructuralVariant(variant, m_config.m_nMaxVariantSize))
            continue;
        
        else if(m_pFastaReader->GetRefSeqSize() < variant.GetEnd())
            continue;
        
        m_aCalledVariantList[variant.m_nChrId].push_back(variant);
    }
}

void CVariantProvider::FillOrientedVariantLists()
{
    for(int i=0; i < CHROMOSOME_COUNT; i++)
    {
        for(int j=0; j < m_aBaseVariantList[i].size(); j++)
        {
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], true));
            m_aBaseOrientedVariantList[i].push_back(COrientedVariant(m_aBaseVariantList[i][j], false));
        }
        
        for(int j=0; j < m_aCalledVariantList[i].size(); j++)
        {
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], true));
            m_aCalledOrientedVariantList[i].push_back(COrientedVariant(m_aCalledVariantList[i][j], false));
        }
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



std::vector<CVariant> CVariantProvider::GetVariantList(EVcfName a_uFrom, int a_nChrNo, std::vector<int> a_VariantIndexes)
{
    std::vector<CVariant> result;

    switch (a_uFrom)
    {
        case eBASE:
            for(int k=0; k< a_VariantIndexes.size(); k++)
                result.push_back(m_aBaseVariantList[a_nChrNo][a_VariantIndexes[k]]);
            break;
            
        case eCALLED:
            for(int k=0; k< a_VariantIndexes.size(); k++)
                result.push_back(m_aCalledVariantList[a_nChrNo][a_VariantIndexes[k]]);
            
        default:
            break;
    }

    return result;
}

bool CVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
{
    std::size_t found;
    
    for(int k = 0; k < a_rVariant.m_nAlleleCount; k++)
    {
        std::string allele = a_rVariant.m_alleles[k].m_sequence;
        
        if(allele.size() > a_nMaxLength)
            return true;
        found = allele.find('[');
        if(found != std::string::npos)
            return true;
        found = allele.find('<');
        if(found != std::string::npos)
            return true;
        found = allele.find('*');
        if(found != std::string::npos)
            return true;
    }
    return false;
}








































