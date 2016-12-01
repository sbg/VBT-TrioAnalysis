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
    
    std::set<CPath>::iterator End() const;

    void Clear();

    void Add(const CPath& item);

    void Erase(const CPath& item);
    
    std::set<CPath>::iterator Find(const CPath& item);
    
    // Get the least advanced path
    void GetLeastAdvanced(CPath& items);

    bool Empty();

    bool Contains(const CPath& item) const;

    CPath floor(const CPath& a_rObj);
    
    void Print() const;

  private:
    
    std::set<CPath> m_set;

};

#endif //_C_THREAD_SAFE_PATH_STACK_H_
