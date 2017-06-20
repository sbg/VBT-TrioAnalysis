//
//  CSimplePEDParser.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 6/20/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CSimplePEDParser.h"
#include <fstream>
#include <sstream>
#include <iostream>

void CSimplePEDParser::ParsePedigree(const std::string& a_rPedFilePath)
{
    std::ifstream pedFile;
    
    pedFile.open(a_rPedFilePath);
    
    std::string pedRecord;
    while (std::getline(pedFile, pedRecord))
    {
        std::stringstream ss(pedRecord);
        
        SPerson person;
        
        std::string familyID;
        ss >> familyID;
        
        if (ss.peek() == '\t')
            ss.ignore();
        
        ss >> person.m_id;
        
        if (ss.peek() == '\t')
            ss.ignore();
        
        ss >> person.m_fatherId;
        
        if (ss.peek() == '\t')
            ss.ignore();
        
        ss >> person.m_motherId;
        
        if (ss.peek() == '\t')
            ss.ignore();
        
        std::string gender;
        ss >> gender;
        person.m_bIsMale = std::stoi(gender) == 1;
        
        m_familyMap[familyID].push_back(person);
    }
}

std::vector<std::string> CSimplePEDParser::GetIdsMFC(const std::string& id1, const std::string& id2, const std::string& id3)
{
    for (auto it = m_familyMap.begin(); it != m_familyMap.end(); ++it )
    {
        for(int k = 0; k < it->second.size(); k++)
        {
            if(it->second[k].m_id == id1)
            {
                if((it->second[k].m_fatherId == id2 && it->second[k].m_motherId == id3)
                   ||
                   (it->second[k].m_fatherId == id3 && it->second[k].m_motherId == id2))
                {
                    std::vector<std::string> res;
                    res.push_back(it->second[k].m_motherId);
                    res.push_back(it->second[k].m_fatherId);
                    res.push_back(it->second[k].m_id);
                    return res;
                }
            }
            
            else if(it->second[k].m_id == id2)
            {
                if((it->second[k].m_fatherId == id1 && it->second[k].m_motherId == id3)
                   ||
                   (it->second[k].m_fatherId == id3 && it->second[k].m_motherId == id1))
                {
                    std::vector<std::string> res;
                    res.push_back(it->second[k].m_motherId);
                    res.push_back(it->second[k].m_fatherId);
                    res.push_back(it->second[k].m_id);
                    return res;
                }
            }
            
            
            else if(it->second[k].m_id == id3)
            {
                if((it->second[k].m_fatherId == id2 && it->second[k].m_motherId == id1)
                   ||
                   (it->second[k].m_fatherId == id1 && it->second[k].m_motherId == id2))
                {
                    std::vector<std::string> res;
                    res.push_back(it->second[k].m_motherId);
                    res.push_back(it->second[k].m_fatherId);
                    res.push_back(it->second[k].m_id);
                    return res;
                }
            }
        }
    }
    
    std::cerr << "Parents and Child could not be identified from given pedigree file." << std::endl;
    std::vector<std::string> res;
    return res;
}

