//
//  CThreadSafePathList.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 11/10/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CThreadSafePathList.h"



bool compare2Path (const CPath& lhs, const CPath& rhs)
{
    bool res = !lhs.IsEqual(rhs) && !rhs.IsEqual(lhs);
    
    return lhs.IsEqual(rhs);
}

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
CPath CThreadSafePathList::GetLeastAdvanced()
{
    CPath toRet = *m_set.begin();
    m_set.erase(m_set.begin());
    return toRet;
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
