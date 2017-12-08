//
//  CSimpleBEDParser.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 7/7/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#include "CSimpleBEDParser.h"
#include <fstream>
#include <sstream>
#include <iostream>

bool CSimpleBEDParser::InitBEDFile(const std::string& a_rBEDFilePath)
{
    bool bIsSuccess = true;
    
    //Open bed file
    std::ifstream bedFile;
    bedFile.open(a_rBEDFilePath.c_str());
    
    try
    {
        std::string bedRecord;
        while (std::getline(bedFile, bedRecord))
        {
            //Ignore comment Line
            if(bedRecord[0] == '#')
                continue;
            
            //Ignore track line
            std::string lineTag = bedRecord.substr(0,5);
            if(lineTag == "track")
                continue;
            
            std::stringstream ss(bedRecord);
            
            SBedRegion region;
            ss >> region.m_chrName;
            
            if (ss.peek() == '\t')
                ss.ignore();
            
            ss >> region.m_nStartPos;
            
            if (ss.peek() == '\t')
                ss.ignore();
            
            ss >> region.m_nEndPos;
            
            if (ss.peek() == '\t')
                ss.ignore();
            
            m_regionMap[region.m_chrName].push_back(region);
        }
    }
    catch(std::exception ex)
    {
        std::cerr << ex.what() << std::endl;
        bIsSuccess = false;
    }
        
    bedFile.close();

    //Specify the total contig count
    m_nTotalContigCount = 0;
    for (auto it= m_regionMap.begin(); it!=m_regionMap.end(); ++it)
        m_nTotalContigCount++;
    
    return bIsSuccess;
}




