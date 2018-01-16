//
//  SInfo.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 1/11/18.
//  Copyright © 2018 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_INFO_H_
#define _S_INFO_H_

#include <vector>

struct SInfoEntry
{
    int n; //Len
    int type; //Info Type
    std::string key = ""; //Info Tag
    void* values = NULL; //Data
};

struct SInfo
{
    SInfo() {};
    
    SInfo(const SInfo& a_rInfo) {m_infoArray = std::vector<SInfoEntry>(a_rInfo.m_infoArray);};
    
    void Clear() {m_infoArray.clear();};
    
    std::vector<SInfoEntry> m_infoArray;
};


#endif // _S_INFO_H_
