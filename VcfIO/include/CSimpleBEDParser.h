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
 *  CSimpleBEDParser.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Toptas on 7/7/17.
 *
 */

#ifndef _C_SIMPLE_BED_PARSER_H_
#define _C_SIMPLE_BED_PARSER_H_

#include <string>
#include <vector>
#include <unordered_map>

/**
 * @brief Defines a BED region specified by the contig name and the interval
 *
 */
struct SBedRegion
{
    std::string m_chrName;
    int m_nStartPos;
    int m_nEndPos;
};

/**
 * @brief Parse BED region for input filtering
 *
 */
class CSimpleBEDParser
{
    
public:
    
    ///Reads the bed file and save it inside. If there is an unknown region format, it returns FALSE
    bool InitBEDFile(const std::string& a_rBEDFilePath);
        
    ///Bed Region map use chromosome name as key [This one is not accesible with GetNextRegion]
    std::unordered_map<std::string, std::vector<SBedRegion>> m_regionMap;
    
    //Total number of contigs in BED file
    unsigned int m_nTotalContigCount;
    
};

#endif // _C_SIMPLE_BED_PARSER_H_
