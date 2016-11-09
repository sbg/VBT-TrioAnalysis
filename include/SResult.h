/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _S_RESULT_H_
#define _S_RESULT_H_

struct SResult
{
    SResult()
    {
        m_nFP = 0;
        m_nTP = 0;
        m_nFN = 0;
    }

    //False Positive Count
    int m_nFP;
    //True Positive Count
    int m_nTP;
    //False Negative Count
    int m_nFN;
};


#endif //_S_RESULT_H_
