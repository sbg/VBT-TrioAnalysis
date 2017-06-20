//
//  CSimplePEDParser.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 6/20/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_SIMPLE_PED_PARSER_H_
#define _C_SIMPLE_PED_PARSER_H_

#include <vector>
#include <string>
#include <unordered_map>

struct SPerson
{
    bool m_bIsMale;
    std::string m_id;
    std::string m_fatherId;
    std::string m_motherId;
};

class CSimplePEDParser
{
    
public:
    
    //Reads the pedigree file and fill the family map
    void ParsePedigree(const std::string& a_rPedFilePath);
    
    //Searches family map for provided 3 person ids and identift mother, father and child ids in order. Returns an empty vector if it could not identify
    std::vector<std::string> GetIdsMFC(const std::string& id1, const std::string& id2, const std::string& id3);
    
private:
    
std::unordered_map<std::string, std::vector<SPerson>> m_familyMap;

};


#endif /* _C_SIMPLE_PED_PARSER_H_ */
