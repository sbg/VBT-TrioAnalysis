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
 *  CPathSet.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas
 *
 */

#include "CPathSet.h"
#include <iostream>

using namespace core;

CPathSet::CPathSet()
{
}

int CPathSet::Size() const
{
    return (int)(m_set.size());
}

void CPathSet::Clear()
{
    m_set.clear();
}

void CPathSet::Add(const CPathContainer& item)
{
    m_set.insert(item);
}

void CPathSet::Erase(const CPathContainer& item)
{
    m_set.erase(item);
}

// Get the least advanced path
void CPathSet::GetLeastAdvanced(CPathContainer& item)
{
    item = *m_set.begin();
    m_set.erase(m_set.begin());
}

bool CPathSet::Empty()
{
    return m_set.empty();
}

bool CPathSet::Contains(const CPathContainer& item) const
{
    std::set<CPathContainer>::iterator it;
    
    for(it = m_set.begin(); it != m_set.end(); ++it)
    {
        if(it->m_pPath->IsEqual(*item.m_pPath))
            return true;
    }
    
    return false;
}

CPathContainer CPathSet::floor(const CPathContainer& a_rObj)
{
    return *m_set.lower_bound(a_rObj);
}

std::set<CPathContainer>::iterator CPathSet::Find(const CPathContainer& item)
{
    return m_set.find(item);
}

std::set<CPathContainer>::iterator CPathSet::End() const
{
    return m_set.end();
}


void CPathSet::Print() const
{
    std::set<CPathContainer>::iterator it;
    std::cout << "Paths:";
    for(it = m_set.begin(); it != m_set.end(); ++it)
    {
        std::cout << it->m_pPath->m_nPathId << " ";
    }
    std::cout << std::endl;
}
