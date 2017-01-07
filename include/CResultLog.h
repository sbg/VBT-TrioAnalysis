//
//  CResultLog.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/29/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_RESULT_LOG_H_
#define _C_RESULT_LOG_H_

#include <string>

struct SLogEntry
{
    bool m_bIsNull = true;
    int m_nTpCalled;
    int m_nTpBase;
    int m_nHalfTpCalled;
    int m_nHalfTpBase;
    int m_nFp;
    int m_nFn;
};


class CResultLog
{
    
public:
    
    //Sets the path of log file
    void SetLogPath(const std::string& a_rLogPath);
    
    //Records the result for given chromosome
    void LogStatistic(int a_nChromosomeId, int a_nTpCalled, int a_nTpBaseline, int a_nHalfTPCalled, int a_nHalfTPBaseline, int a_nFalsePositive, int a_nFalseNegative);
    
    //Write the results in log.txt file
    void WriteStatistics();
    
private:
    
    SLogEntry m_aResultEntries[25];
    
    std::string m_aLogPath;
    
};



#endif // _C_RESULT_LOG_H_
