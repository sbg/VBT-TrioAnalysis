//
//  CThreadSafePathList.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 11/10/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CThreadSafePathList.h"

CThreadSafePathList::CThreadSafePathList()
{
}

int CThreadSafePathList::Size() const
{
    return (int)(m_set.size());
}

void CThreadSafePathList::Clear()
{
    m_set.clear();
}

void CThreadSafePathList::Add(const CPath& item)
{
    m_set.insert(CPath(item));
}

void CThreadSafePathList::Erase(const CPath& item)
{
    m_set.erase(item);
}

// Get the least advanced path
void CThreadSafePathList::GetLeastAdvanced(CPath& item)
{
    item = *m_set.begin();
    m_set.erase(m_set.begin());
}

bool CThreadSafePathList::Empty()
{
    return m_set.empty();
}

bool CThreadSafePathList::Contains(const CPath& item) const
{
    std::set<CPath>::iterator it;
    
    for(it = m_set.begin(); it != m_set.end(); ++it)
    {
        if(it->IsEqual(item))
            return true;
    }
    
    return false;
}

CPath CThreadSafePathList::floor(const CPath& a_rObj)
{
    return *m_set.lower_bound(a_rObj);
}

std::set<CPath>::iterator CThreadSafePathList::Find(const CPath& item)
{
    return m_set.find(item);
}

std::set<CPath>::iterator CThreadSafePathList::End() const
{
    return m_set.end();
}


void CThreadSafePathList::Print() const
{
    std::set<CPath>::iterator it;
    std::cout << "Paths:";
    for(it = m_set.begin(); it != m_set.end(); ++it)
    {
        std::cout << it->m_nPathId << " ";
    }
    std::cout << std::endl;
}
