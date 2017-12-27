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

#ifndef _C_PATH_SET_H_
#define _C_PATH_SET_H_

#include "CPath.h"
#include <set>

namespace core
{

/**
 * @brief Search Tree of CPaths - current countainer is std::set from STL
 *
 */
class CPathSet
{
  public:

    CPathSet();
    
    ///Return size of Tree
    int Size() const;
    
    std::set<CPathContainer>::iterator End() const;

    ///Clear the search tree
    void Clear();

    ///Add CPathContainer to search tree
    void Add(const CPathContainer& item);

    ///Delete CPathContainer from search tree
    void Erase(const CPathContainer& item);
    
    ///Find closest CPathContainer in search tree to the given path
    std::set<CPathContainer>::iterator Find(const CPathContainer& item);
    
    ///Pops the least advanced CPathContainer from the search tree
    void GetLeastAdvanced(CPathContainer& items);

    ///Checks if the search tree is empty
    bool Empty();

    ///Chekcs if the search tree contains given CPathContainer
    bool Contains(const CPathContainer& item) const;

    CPathContainer floor(const CPathContainer& a_rObj);
    
    ///Print the search tree [FOR TEST]
    void Print() const;

  private:
    
    std::set<CPathContainer> m_set;

};

}

#endif //_C_PATH_SET_H_
