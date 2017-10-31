//
//  CSimpleBEDParser.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 7/7/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

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
    
private:
    
    ///List of regions from BED file
    std::vector<SBedRegion> m_regionArray;
    
    ///Points to the Next Region in m_regionArray
    unsigned int m_nIterator;
    
};

#endif // _C_SIMPLE_BED_PARSER_H_
