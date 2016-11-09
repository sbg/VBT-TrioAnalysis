#include "CSemiPath.h"

inline int max(int a, int b)
{
    return a > b ? a : b;
}

CSemiPath::CSemiPath(const char* a_aRefSequence, int a_nRefSize, EVcfName a_uVcfName) 
: m_haplotypeA(a_aRefSequence, a_nRefSize),
  m_haplotypeB(a_aRefSequence, a_nRefSize)
{
    m_uVcfName = a_uVcfName;
    m_nVariantIndex = -1;
    m_nIncludedVariantEndPosition = 0;
    m_nVariantEndPosition = 0;
    m_bFinishedHapA = false;
    m_bFinishedHapB = false;   
}

CSemiPath::CSemiPath(const CSemiPath& a_rObj)
: m_haplotypeA(a_rObj.m_haplotypeA),
  m_haplotypeB(a_rObj.m_haplotypeB)
{
    m_uVcfName = a_rObj.m_uVcfName;
    m_nVariantIndex = a_rObj.m_nVariantIndex;
    m_nIncludedVariantEndPosition = a_rObj.m_nIncludedVariantEndPosition;  
    m_nVariantEndPosition = a_rObj.m_nVariantEndPosition;

    m_aIncludedVariants = std::vector<COrientedVariant>(a_rObj.m_aIncludedVariants);
    m_aExcludedVariants = std::vector<CVariant>(a_rObj.m_aExcludedVariants);

    m_bFinishedHapA = a_rObj.m_bFinishedHapA;
    m_bFinishedHapB = a_rObj.m_bFinishedHapB;
}

EVcfName CSemiPath::GetVcfName() const
{
    return m_uVcfName;
}

void CSemiPath::IncludeVariant(const COrientedVariant& a_rVariant, int a_nVariantIndex)
{
    if(a_nVariantIndex <= m_nVariantIndex)
    {
        std::cout << "ERROR INCLUDING FILE" << std::endl;
        return;
    }

    m_aIncludedVariants.push_back(a_rVariant);
    m_nVariantIndex = a_nVariantIndex;
    m_nVariantEndPosition = max(m_nVariantEndPosition, a_rVariant.GetVariant().GetEnd());
    m_nIncludedVariantEndPosition = max(m_nIncludedVariantEndPosition, a_rVariant.GetVariant().GetEnd());

    m_haplotypeA.AddVariant(a_rVariant);
    m_haplotypeB.AddVariant(a_rVariant.Other());
}

void CSemiPath::ExcludeVariant(const CVariant& a_rVariant, int a_nVariantIndex)
{
    if (a_nVariantIndex <= m_nVariantIndex)
    {
        std::cout << "ERROR EXCLUDING FILE" << std::endl;
        return;
    }

    m_aExcludedVariants.push_back(a_rVariant);
    m_nVariantEndPosition = max(m_nVariantEndPosition, a_rVariant.GetEnd());
    m_nVariantIndex = a_nVariantIndex;
}

int CSemiPath::CompareHaplotypePositions() const
{
    return m_haplotypeA.GetTemplatePosition() - m_haplotypeB.GetTemplatePosition();
}

int CSemiPath::GetPosition() const
{
    if (m_haplotypeA.GetTemplatePosition() > m_haplotypeB.GetTemplatePosition()) 
        return m_haplotypeA.GetTemplatePosition();
    else 
        return m_haplotypeB.GetTemplatePosition();
}

int CSemiPath::GetVariantIndex() const
{
    return m_nVariantIndex;
}

int CSemiPath::GetVariantEndPosition() const
{
    return m_nVariantEndPosition;    
}

bool CSemiPath::IsOnTemplate() const 
{
    return m_haplotypeA.IsOnTemplate() && m_haplotypeB.IsOnTemplate();
}

int CSemiPath::GetIncludedVariantEndPosition() const
{
    return m_nIncludedVariantEndPosition;
}

int CSemiPath::CompareTo(const CSemiPath& a_rObj) const
{
    int res = m_haplotypeA.CompareTo(a_rObj.m_haplotypeA);
    if(res != 0)
        return res;
    else
        return m_haplotypeB.IsEqual(a_rObj.m_haplotypeB);
}

bool CSemiPath::IsEqual(const CSemiPath& a_rObj) const
{
    return (CompareTo(a_rObj) == 0);
}

bool CSemiPath::HasFinished() const
{
    return m_bFinishedHapA && m_bFinishedHapB;
}

void CSemiPath::MoveForward(int a_nPosition)
{
    m_haplotypeA.MoveForward(a_nPosition);
    m_haplotypeB.MoveForward(a_nPosition);
}

bool CSemiPath::Matches(const CSemiPath& a_rOther)
{
    if (!FinishedHaplotypeA() && !a_rOther.FinishedHaplotypeA() && NextHaplotypeABase() != a_rOther.NextHaplotypeABase()) 
        return false;
    else if (!FinishedHaplotypeB() && !a_rOther.FinishedHaplotypeB() && NextHaplotypeBBase() != a_rOther.NextHaplotypeBBase())
        return false;
    else
        return true;
}

bool CSemiPath::WantsFutureVariantBases() const
{
    return m_haplotypeA.WantsFutureVariantBases() || m_haplotypeB.WantsFutureVariantBases();
}

bool CSemiPath::FinishedHaplotypeA() const
{
    return m_bFinishedHapA;
}

bool CSemiPath::FinishedHaplotypeB() const
{
    return m_bFinishedHapB;
}

char CSemiPath::NextHaplotypeABase() const
{
    return m_haplotypeA.NextBase();
}

char CSemiPath::NextHaplotypeBBase() const
{
    return m_haplotypeB.NextBase();
}

void CSemiPath::StepHaplotypeA()
{
    if (m_haplotypeA.HasNext())
        m_haplotypeA.Next();
    else
        m_bFinishedHapA = true;
}

void CSemiPath::StepHaplotypeB()
{
    if (m_haplotypeB.HasNext())
        m_haplotypeB.Next();
    else
        m_bFinishedHapB = true;
}
