#include "CVariant.h"
#include <iostream>


CVariant::CVariant(): m_nVcfId(-1),
            m_nChrId(-1), 
            m_nStartPos(-1),
            m_chrName(),
            m_bIsPhased(false) 
{
    m_nId = -1;
    m_nAlleleCount = 0;
    m_nMaxLength = 0;
    m_bIsDuplicate = false;
}

CVariant::CVariant(const CVariant& a_rObj)
{
    m_nVcfId = a_rObj.m_nVcfId;
    m_nChrId = a_rObj.m_nChrId;
    m_chrName = a_rObj.m_chrName;
    m_bIsPhased = a_rObj.m_bIsPhased;
    m_nAlleleCount = a_rObj.m_nAlleleCount;
    m_alleles[0].m_nEndPos = a_rObj.m_alleles[0].m_nEndPos;
    m_alleles[0].m_nStartPos = a_rObj.m_alleles[0].m_nStartPos;
    m_alleles[0].m_sequence = a_rObj.m_alleles[0].m_sequence;
    m_alleles[1].m_nEndPos = a_rObj.m_alleles[1].m_nEndPos;
    m_alleles[1].m_nStartPos = a_rObj.m_alleles[1].m_nStartPos;
    m_alleles[1].m_sequence = a_rObj.m_alleles[1].m_sequence;m_nStartPos = a_rObj.m_nStartPos;
    m_nEndPos = a_rObj.m_nEndPos;
    m_nId = a_rObj.m_nId;
    m_refSequence = a_rObj.m_refSequence;
    m_nMaxLength = a_rObj.m_nMaxLength;
    m_bIsDuplicate = a_rObj.m_bIsDuplicate;
}


int CVariant::CompareTo(const CVariant& a_rObj) const
{
    if(GetStart() < a_rObj.GetStart())
        return -1;
    else if(GetStart() > a_rObj.GetStart())
        return 1;
    else
        return GetEnd() - a_rObj.GetEnd();
}

bool CVariant::Clear()
{
    m_nVcfId = -1;
    m_nChrId = -1;
    m_nStartPos = -1;
    m_chrName.clear();
    m_nAlleleCount = 0;
    m_nMaxLength = 0;
    m_bIsDuplicate = false;
    return true;
}

bool CVariant::IsHeterozygous() const
{
    if(m_nAlleleCount != 2)
        return false;
    
    return m_alleles[0].m_sequence != m_alleles[1].m_sequence;
}

void CVariant::Print() const
{
    if(m_nVcfId == 0)
        std::cout << "Belongs to  : Baseline" << std::endl;
    else
        std::cout << "Belongs to  : Called" << std::endl;
    std::cout <<     "Ref : " << m_refSequence << std::endl;
    for(int k = 0; k < m_nAlleleCount; k++)
    {
        std::cout << "Alt" << k << ": " << m_alleles[k].m_sequence << std::endl;
    }
}

int CVariant::GetId() const
{
    return m_nId;
}

int CVariant::GetStart() const
{
    return m_nStartPos;
}

int CVariant::GetEnd() const
{
    return m_nEndPos;
}

bool CVariant::IsPhased() const
{
    return m_bIsPhased;
}

std::string CVariant::GetRefSeq() const
{
    return m_refSequence;
}

SAllele CVariant::GetAllele(int a_nAlleleId) const
{
    return m_alleles[a_nAlleleId];
}

bool CVariant::IsNull() const
{
    if(m_nVcfId == -1)
        return true;
    else
        return false;
}

std::string CVariant::ToString() const
{
    std::string toRet = "";
    
    toRet = m_chrName + ":" + std::to_string(GetStart() + 1) + "-" + std::to_string(GetEnd() + 1) + " (";
    
    for(int k=0; k < m_nAlleleCount; k++)
    {
        if(k > 0)
            toRet = toRet + ":";
        
        toRet = toRet + m_alleles[k].m_sequence;
    }
    toRet = toRet + ")";
    
    return toRet;
}

bool CVariant::IsFilterPASS() const
{
    return m_bIsFilterPASS;
}











