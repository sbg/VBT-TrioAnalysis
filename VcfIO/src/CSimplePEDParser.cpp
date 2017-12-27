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
 *  CSimplePEDParser.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 6/20/17.
 *
 */

#include "CSimplePEDParser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

void CSimplePEDParser::ParsePedigree(const std::string& a_rPedFilePath)
{
    std::ifstream pedFile;
    
    pedFile.open(a_rPedFilePath);
    
    std::string pedRecord;
    while (std::getline(pedFile, pedRecord))
    {
        //Ignore comment Line
        if(pedRecord[0] == '#')
            continue;
        
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


std::vector<std::string> CSimplePEDParser::GetIdsMFC(const std::vector<std::string> motherList,
                                                     const std::vector<std::string> fatherList,
                                                     const std::vector<std::string> childList)
{
    for (auto it = m_familyMap.begin(); it != m_familyMap.end(); ++it )
    {
        for(unsigned int k = 0; k < childList.size(); k++)
        {
            for(unsigned int m = 0; m < it->second.size(); m++)
            {
                if(it->second[m].m_id == childList[k])
                {
                    auto fatherItr = std::find(fatherList.begin(), fatherList.end(), it->second[m].m_fatherId);
                    auto motherItr = std::find(motherList.begin(), motherList.end(), it->second[m].m_motherId);
                    if(fatherItr != fatherList.end() && motherItr != motherList.end())
                    {
                        std::vector<std::string> res;
                        res.push_back(it->second[m].m_motherId);
                        res.push_back(it->second[m].m_fatherId);
                        res.push_back(it->second[m].m_id);
                        return res;
                    }
                }
            }
        }

    }

    std::cerr << "Parents and Child could not be identified from given pedigree file." << std::endl;
    std::vector<std::string> res;
    return res;
}

std::vector<std::string> CSimplePEDParser::GetIdsMFC(const std::string& id1, const std::string& id2, const std::string& id3)
{
    for (auto it = m_familyMap.begin(); it != m_familyMap.end(); ++it )
    {
        for(unsigned int k = 0; k < it->second.size(); k++)
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

