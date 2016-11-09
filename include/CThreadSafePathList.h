/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

//TODO: Read/Write lock for multithread support will be added

#ifndef _C_THREAD_SAFE_PATH_STACK_H_
#define _C_THREAD_SAFE_PATH_STACK_H_

#include "CPath.h"
#include <set>

struct classcomp 
{
  bool operator() (const CPath& lhs, const CPath& rhs) const
  {
      return lhs.IsEqual(rhs);
  }
};

class CThreadSafePathList
{
  public:

    int Size() const
    {
      return (int)(m_set.size());
    }

    void Clear()
    {
      m_set.clear();
    }

    // Add the path to the processed path list if there is no better path with the same sequence number
    void Add(const CPath& item) 
    {
      m_set.insert(item);
    }

    void Erase(const CPath& item)
    {
      m_set.erase(item);
    }

    // Get the least advanced path
    CPath GetLeastAdvanced() 
    { 
      CPath toRet = *m_set.begin();
      m_set.erase(m_set.begin());
      return toRet;
    }

    bool Empty()
    {
      return m_set.empty();
    }

    bool Contains(const CPath& item) const
    {
      std::set<CPath>::iterator it;

      for(it = m_set.begin(); it != m_set.end(); ++it)
      {
        if(true == it->IsEqual(item))
          return true;
      }

      return false;
    }

    CPath floor(const CPath& a_rObj)
    {
      return *m_set.lower_bound(a_rObj);  
    }

  private:  
    std::set<CPath,classcomp>  m_set;

};

#endif //_C_THREAD_SAFE_PATH_STACK_H_
