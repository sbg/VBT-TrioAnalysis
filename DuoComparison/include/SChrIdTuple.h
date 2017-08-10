//
//  SChrIdTuple.h
//  VCFComparison
//
//  Created by Berke.Toptas on 5/2/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CHR_ID_TUPLE_H_
#define _S_CHR_ID_TUPLE_H_

#include <string>

namespace duocomparison
{

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
    
    int m_nBaseId;
    int m_nCalledId;
    int m_nTupleIndex;
    std::string m_chrName;
};

}

#endif /* _S_CHR_ID_TUPLE_H_ */
