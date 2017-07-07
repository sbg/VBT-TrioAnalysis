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

struct SBedRegion
{
    std::string m_chrName;
    int m_nStartPos;
    int m_nEndPos;
};

class CSimpleBEDParser
{
    
public:
    
    //Reads the bed file and save it inside. If there is an unknown region format, it returns FALSE
    bool InitBEDFile(const std::string& a_rBEDFilePath);
    
    //Return the next region from bed file
    bool GetNextRegion(SBedRegion& a_rRegion);
    
    //Set the next Region as first region
    void ResetIterator();
    
private:
    
    //List of regions from BED file
    std::vector<SBedRegion> m_regionArray;
    
    //Points to the Next Region in m_regionArray
    int m_nIterator;
    
};

#endif // _C_SIMPLE_BED_PARSER_H_
