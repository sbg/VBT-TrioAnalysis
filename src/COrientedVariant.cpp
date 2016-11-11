#include "COrientedVariant.h"


COrientedVariant::COrientedVariant()
{
    m_nAlleleIndex = -1;
    m_nOtherAlleleIndex = -1;
}

COrientedVariant::COrientedVariant(const CVariant& a_rObj, bool a_bIsOrder)
: m_variant(a_rObj)
{
    if(a_rObj.IsHeterozygous() && false == a_bIsOrder)
    {
        m_nAlleleIndex = a_rObj.gt_arr[1];
        m_nOtherAlleleIndex = a_rObj.gt_arr[0];
    }
    else
    {
        m_nAlleleIndex = a_rObj.gt_arr[0];
        m_nOtherAlleleIndex = a_rObj.gt_arr[1];
    }
    m_bIsOrderOfGenotype = a_bIsOrder;
}

COrientedVariant::COrientedVariant(const COrientedVariant& a_rObj)
: m_variant(a_rObj.m_variant)
{
    m_nAlleleIndex = a_rObj.m_nAlleleIndex;
    m_nOtherAlleleIndex = a_rObj.m_nOtherAlleleIndex;
    m_bIsOrderOfGenotype = a_rObj.m_bIsOrderOfGenotype;
}


int COrientedVariant::CompareTo(const COrientedVariant& a_rObj) const
{
    //decide by variant id
    int id = (m_variant.GetId() < a_rObj.m_variant.GetId()) ? -1 : ((m_variant.GetId() == a_rObj.m_variant.GetId()) ? 0 : 1);
    if(id != 0)
        return id;
    
    //decide by genotype order
    int genotype = (m_bIsOrderOfGenotype == a_rObj.m_bIsOrderOfGenotype) ? 0 : (m_bIsOrderOfGenotype ? 1 : -1);
    if(genotype != 0)
        return genotype;
    
    // Decide by allele index
    int alleleId = (m_nAlleleIndex < a_rObj.m_nAlleleIndex) ? -1 : ((m_nAlleleIndex == a_rObj.m_nAlleleIndex) ? 0 : 1);
    if (alleleId != 0)
        return alleleId;
    
    //Decide by other allele index
    return (m_nOtherAlleleIndex < a_rObj.m_nOtherAlleleIndex) ? -1 : ((m_nOtherAlleleIndex == a_rObj.m_nOtherAlleleIndex) ? 0 : 1);
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

int COrientedVariant::GetAlleleIndex() const
{
    return m_nAlleleIndex;
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
