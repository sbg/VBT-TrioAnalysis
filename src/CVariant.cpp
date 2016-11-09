#include "CVariant.h"

CVariant::CVariant(): m_nVcfId(-1),
            m_nChrId(-1), 
            m_nPosition(-1),
            m_chrName(), 
            m_aSequences(),
            ngt_arr(0),
            m_bIsPhased(false) 
{
    gt_arr[0] = -1;
    gt_arr[1] = -1;
}

CVariant::CVariant(const CVariant& a_rObj)
: m_aSequences(a_rObj.m_aSequences)
{
    m_nVcfId = a_rObj.m_nVcfId;
    m_nChrId = a_rObj.m_nChrId;
    m_chrName = a_rObj.m_chrName;
    m_bIsPhased = a_rObj.m_bIsPhased;
    m_nPosition = a_rObj.m_nPosition;
    gt_arr[0] = a_rObj.gt_arr[0];
    gt_arr[1] = a_rObj.gt_arr[1];
    ngt_arr = a_rObj.ngt_arr;
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
        m_nPosition    = -1;
        m_chrName.clear();
        m_aSequences.clear();
        return true;
}

bool CVariant::IsHeterozygous() const
{
        if(ngt_arr != 2)
            return false;
    
        //std::cout << gt_arr[0] << " " << gt_arr[1] << std::endl;
        return gt_arr[0] != gt_arr[1];
}

void CVariant::SetType(int a_nAltIndex)
{
    EVariantType uType;
        if(gt_arr[a_nAltIndex -1] == 0)
        uType = eNO_OP;
    else if(m_aSequences[a_nAltIndex].length() == m_aSequences[0].length())
        uType = eSNP;
    else if(m_aSequences[a_nAltIndex].length() < m_aSequences[0].length())
        uType = eINDEL_DEL;
    else
        uType = eINDEL_ADD; 

    m_aVarTypes[a_nAltIndex -1] = uType;      
}

void CVariant::Print() const
{
    if(m_nVcfId == 0)
        std::cout << "Belongs to  : Baseline" << std::endl;
    else
        std::cout << "Belongs to  : Called" << std::endl;
    std::cout <<     "Ref : " << m_aSequences[0] << std::endl;
    for(int k = 1; k < m_aSequences.size(); k++)
    {
        std::cout << "Alt" << k << ": " << m_aSequences[k] << std::endl;
        if(m_aVarTypes[k] == eSNP)
            std::cout << "Type: SNP" << std::endl;
        else if(m_aVarTypes[k] == eINDEL_DEL)
            std::cout << "Type: DELETION" << std::endl;
        else if(m_aVarTypes[k] == eINDEL_ADD)
            std::cout << "Type: ADDITION" << std::endl;
        else
            std::cout << "Type: NO OP" << std::endl;
    }
    std::cout << "GenoType: " << gt_arr[0] << "/" << gt_arr[1] << std::endl;

}

int CVariant::GetStart() const
{
    return m_nPosition;
}

int CVariant::GetEnd() const
{
    return m_nPosition + (int)m_aSequences[0].length();
}

bool CVariant::IsPhased() const
{
    return m_bIsPhased;
}

std::string CVariant::GetRefSeq() const
{
    return m_aSequences[0];
}

std::string CVariant::GetAllele(int a_nAlleleId) const
{
    return m_aSequences[a_nAlleleId - 1];
}

bool CVariant::IsNull() const
{
    if(m_nVcfId == -1)
        return true;
    else
        return false;
}
