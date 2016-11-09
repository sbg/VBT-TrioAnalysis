#include "COrientedVariant.h"


COrientedVariant::COrientedVariant()
{
    m_nAlleleIndex = -1;
}

COrientedVariant::COrientedVariant(const CVariant& a_rObj, bool a_bIsOrder)
: m_variant(a_rObj)
{
    if(a_rObj.IsHeterozygous() && false == a_bIsOrder)
    {
        m_nAlleleIndex = a_rObj.gt_arr[1];
    }
    else
    {
        m_nAlleleIndex = a_rObj.gt_arr[0];
    }
    m_bIsOrderOfGenotype = a_bIsOrder;
}

COrientedVariant::COrientedVariant(const COrientedVariant& a_rObj)
: m_variant(a_rObj.m_variant)
{
    m_nAlleleIndex = a_rObj.m_nAlleleIndex;
    m_bIsOrderOfGenotype = a_rObj.m_bIsOrderOfGenotype;
}


int COrientedVariant::CompareTo(const COrientedVariant& a_rObj) const
{
    return m_variant.CompareTo(a_rObj.m_variant);
}

std::string COrientedVariant::GetAlleleString() const
{
    return m_variant.m_aSequences[m_nAlleleIndex];
}

int COrientedVariant::GetAlleleStartPos() const
{
    return m_variant.GetStart();
}

int COrientedVariant::GetAlleleEndPos() const
{
    return m_variant.GetStart() + (int)m_variant.m_aSequences[m_nAlleleIndex].length();
}

COrientedVariant COrientedVariant::Other() const
{
    return COrientedVariant(m_variant, !m_bIsOrderOfGenotype);
}

bool COrientedVariant::IsNull() const
{
    return m_variant.IsNull();
}

void COrientedVariant::SetToNull()
{
    m_variant.m_nVcfId = -1;
}

const CVariant& COrientedVariant::GetVariant() const
{
    return m_variant;
}
