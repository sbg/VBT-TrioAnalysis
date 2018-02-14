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
 *  Created by Berke Cagkan Toptas
 *
 */

#ifndef _C_SIMPLE_PED_PARSER_H_
#define _C_SIMPLE_PED_PARSER_H_

#include <vector>
#include <string>
#include <unordered_map>

/**
 * @brief Defines a simplified pedigree entry
 *
 */
struct SPerson
{
    bool m_bIsMale;
    std::string m_id;
    std::string m_fatherId;
    std::string m_motherId;
};


/**
 * @brief Parse pedigree(PED) files for trio analysis where the input is a trio rather than 3 single individual
 *
 */
class CSimplePEDParser
{
    
public:
    
    ///Reads the pedigree file and fill the family map
    void ParsePedigree(const std::string& a_rPedFilePath);
    
    ///Searches family map for provided 3 person ids and identift mother, father and child ids in order. Returns an empty vector if it could not identify
    std::vector<std::string> GetIdsMFC(const std::string& id1, const std::string& id2, const std::string& id3);
    ///Searches family map for provided family id lists. Returns an empty vector if it could not identify
    std::vector<std::string> GetIdsMFC(const std::vector<std::string> motherList,
                                       const std::vector<std::string> fatherList,
                                       const std::vector<std::string> childList);
    
private:

///A Family pedigree tree constructed as the PED file is parsed
std::unordered_map<std::string, std::vector<SPerson>> m_familyMap;

};


#endif /* _C_SIMPLE_PED_PARSER_H_ */
