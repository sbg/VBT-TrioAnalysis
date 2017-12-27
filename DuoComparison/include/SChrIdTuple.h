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
 *  SChrIdTuple.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 5/2/17.
 *
 */

#ifndef _S_CHR_ID_TUPLE_H_
#define _S_CHR_ID_TUPLE_H_

#include <string>

namespace duocomparison
{

/**
 * @brief Groups indexes of common chromosomes for truth and query vcfs
 *
 */
struct SChrIdTuple
{
    SChrIdTuple(int a_base, int a_called, std::string a_chrName, int a_nTupleIndex)
    {
        m_nBaseId = a_base;
        m_nCalledId = a_called;
        m_chrName = a_chrName;
        m_nTupleIndex = a_nTupleIndex;
    }
    
    SChrIdTuple()
    {
        m_nCalledId = -1;
        m_nBaseId = -1;
        m_chrName = "none";
        m_nTupleIndex = -1;
    }
    
    ///Index of chromosome in truth vcf
    int m_nBaseId;
    ///Index of chromosome in query vcf
    int m_nCalledId;
    ///Tuple Index
    int m_nTupleIndex;
    ///Chromosome name
    std::string m_chrName;
};

}

#endif /* _S_CHR_ID_TUPLE_H_ */
