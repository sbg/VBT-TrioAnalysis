/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/HalfPath.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 * Licensed under the Simplified BSD License: https://github.com/RealTimeGenomics/rtg-tools/blob/master/LICENSE.txt
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 *
 */

#include "CSemiPath.h"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace core;

inline int max(int a, int b)
{
    return a > b ? a : b;
}

CSemiPath::CSemiPath()
{}

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

    m_aIncludedVariants = std::vector<const COrientedVariant*>(a_rObj.m_aIncludedVariants);
    m_aExcludedVariants = std::vector<int>(a_rObj.m_aExcludedVariants);

    m_bFinishedHapA = a_rObj.m_bFinishedHapA;
    m_bFinishedHapB = a_rObj.m_bFinishedHapB;
}

EVcfName CSemiPath::GetVcfName() const
{
    return m_uVcfName;
}

void CSemiPath::IncludeVariant(const COrientedVariant& a_rVariant, int a_nVariantIndex)
{
    assert(a_nVariantIndex > m_nVariantIndex);

    m_aIncludedVariants.push_back(&a_rVariant);
    m_nVariantIndex = a_nVariantIndex;
    m_nVariantEndPosition = max(m_nVariantEndPosition, a_rVariant.GetVariant().GetEnd());
    m_nIncludedVariantEndPosition = std::max(m_nIncludedVariantEndPosition, a_rVariant.GetVariant().GetEnd());
    
    m_haplotypeA.AddVariant(a_rVariant);
    m_haplotypeB.AddVariant(a_rVariant.Other());
}

void CSemiPath::ExcludeVariant(const CVariant& a_rVariant, int a_nVariantIndex)
{
    assert(a_nVariantIndex > m_nVariantIndex);

    m_aExcludedVariants.push_back(a_nVariantIndex);
    m_nVariantEndPosition = max(m_nVariantEndPosition, a_rVariant.GetEnd());
    m_nVariantIndex = a_nVariantIndex;
}

int CSemiPath::CompareHaplotypePositions() const
{
    return m_haplotypeA.GetTemplatePosition() - m_haplotypeB.GetTemplatePosition();
}


const std::vector<int>& CSemiPath::GetExcluded() const
{
    return m_aExcludedVariants;
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

void CSemiPath::SetVariantIndex(int a_nVariantIndex)
{
    m_nVariantIndex = a_nVariantIndex;
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

const std::vector<const COrientedVariant*>& CSemiPath::GetIncludedVariants() const
{
    return m_aIncludedVariants;
}

int CSemiPath::CompareTo(const CSemiPath& a_rObj) const
{
    int res = m_haplotypeA.CompareTo(a_rObj.m_haplotypeA);
    if(res != 0)
        return res;
    else
        return m_haplotypeB.CompareTo(a_rObj.m_haplotypeB);
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


bool CSemiPath::IsNew(const COrientedVariant& a_rVar) const
{
    if(a_rVar.GetStartPos() >= m_nIncludedVariantEndPosition)
        return true;
    else
        return m_haplotypeA.IsNew(a_rVar) && m_haplotypeB.IsNew(a_rVar.Other());
}


bool CSemiPath::Matches(const CSemiPath& a_rOther)
{    
    if (!FinishedHaplotypeA() && !a_rOther.FinishedHaplotypeA() && toupper(NextHaplotypeABase()) != toupper(a_rOther.NextHaplotypeABase()))
        return false;
    else if (!FinishedHaplotypeB() && !a_rOther.FinishedHaplotypeB() && toupper(NextHaplotypeBBase()) != toupper(a_rOther.NextHaplotypeBBase()))
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

void CSemiPath::ClearIncludedVariants()
{
    m_aIncludedVariants.clear();
}

void CSemiPath::AddIncludedVariants(std::vector<const COrientedVariant*>& a_rIncludedVarList)
{
    m_aIncludedVariants = a_rIncludedVarList;
}

void CSemiPath::ClearExcludedVariants()
{
    m_aExcludedVariants.clear();
}

void CSemiPath::AddExcludedVariants(std::vector<int> &a_rExcludedVarList)
{
    m_aExcludedVariants = a_rExcludedVarList;
}

void CSemiPath::SortIncludedVariants()
{
    std::sort(m_aIncludedVariants.begin(), m_aIncludedVariants.end(), [](const COrientedVariant* pOvar1, const COrientedVariant* pOvar2){return pOvar1->GetVariant().m_nId < pOvar2->GetVariant().m_nId;});
}

void CSemiPath::Print() const
{
    std::cout<< "Pos:" << GetPosition() << " VarEnd Pos:" << GetVariantEndPosition() << " VarEnd Ind:" << GetVariantIndex() << std::endl;
    std::cout<< "Excluded Var Count:" << m_aExcludedVariants.size() << " Included Var Count:" << m_aIncludedVariants.size() << std::endl;
    
    if(m_aIncludedVariants.size() > 0 && !m_aIncludedVariants.back()->IsNull())
        m_aIncludedVariants.back()->Print();
        
    std::cout<< "Haplotype A:" << std::endl;
    m_haplotypeA.Print();
    std::cout<< "Haplotype B:" << std::endl;
    m_haplotypeB.Print();
}

