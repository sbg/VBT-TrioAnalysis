/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

//TODO: Read/Write lock for multithread support will be added

#ifndef _C_THREAD_SAFE_PATH_STACK_H_
#define _C_THREAD_SAFE_PATH_STACK_H_

#include "CPath.h"
#include <set>

class CThreadSafePathList
{
  public:

    CThreadSafePathList();
    
    int Size() const;
    
    std::set<CPathContainer>::iterator End() const;

    void Clear();

    void Add(const CPathContainer& item);

    void Erase(const CPathContainer& item);
    
    std::set<CPathContainer>::iterator Find(const CPathContainer& item);
    
    // Get the least advanced path
    void GetLeastAdvanced(CPathContainer& items);

    bool Empty();

    bool Contains(const CPathContainer& item) const;

    CPathContainer floor(const CPathContainer& a_rObj);
    
    void Print() const;

  private:
    
    std::set<CPathContainer> m_set;

};

#endif //_C_THREAD_SAFE_PATH_STACK_H_
