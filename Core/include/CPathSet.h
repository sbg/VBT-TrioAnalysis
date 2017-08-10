/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_PATH_SET_H_
#define _C_PATH_SET_H_

#include "CPath.h"
#include <set>

namespace core
{

class CPathSet
{
  public:

    CPathSet();
    
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

}

#endif //_C_PATH_SET_H_
