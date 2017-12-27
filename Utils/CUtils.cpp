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
 *  CUtils.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 12/4/17.
 *
 */

#include "Utils/CUtils.h"
#include "CVariant.h"
#include "COrientedVariant.h"
#include <fstream>

//Checks if the given two range is overlapping
bool CUtils::IsOverlap(int left1, int right1, int left2, int right2)
{
    //If the interval length is 0 (eg. 974791-974791) we need to check if the boundaries matches
    if(left1 == left2)
        return true;
    
    if(right1 == left1)
        return (right2 > left1 && left2 <= left1);
    
    if(right2 == left2)
        return (right1 > left2 && left1 <= left2);
    
    return std::min(right1, right2) - std::max(left1, left2) > 0;
}

bool CUtils::CompareVariants(const CVariant& var1, const CVariant& var2)
{
    if(var1.m_nStartPos != var2.m_nStartPos)
        return var1.m_nStartPos < var2.m_nStartPos;
    
    if(var1.m_nEndPos != var2.m_nEndPos)
        return var1.m_nEndPos < var2.m_nEndPos;

    return var1.m_nId < var2.m_nId;
}

bool CUtils::IsHomRef(const CVariant& a_rVariant)
{
    bool res = true;
    
    for(int k = 0; k < a_rVariant.m_nZygotCount; k++)
    {
        if(a_rVariant.m_genotype[k] != 0)
        {
            res = false;
            break;
        }
    }
    return res;
}

bool CUtils::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength)
{
    std::size_t found;
    
    for(int k = 0; k < a_rVariant.m_nAlleleCount; k++)
    {
        std::string allele = a_rVariant.m_alleles[k].m_sequence;
        
        if((int)allele.size() > a_nMaxLength)
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
        found = allele.find('.');
        if(found != std::string::npos)
            return true;
    }
    return false;
}

//Compare variants according to id for sort operation
bool CUtils::CompareVariantsById(const CVariant* v1, const CVariant* v2)
{
    return v1->m_nId < v2->m_nId;
}

bool CUtils::CompareOrientedVariantsById(const core::COrientedVariant* v1, const core::COrientedVariant* v2)
{
    return v1->GetVariant().m_nId < v2->GetVariant().m_nId;
}

bool CUtils::IsFileExists (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
