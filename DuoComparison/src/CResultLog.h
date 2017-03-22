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
#include "Constants.h"

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

struct SMendelianLogEntry
{
    bool m_bIsNull = true;
    int m_nTpCalled;
    int m_nTpBase;
    int m_nFP;
    int m_nFn;
};

struct SMendelianMergeLogEntry
{
    bool m_bIsNull = true;
    int m_nCommonHomozygous;
    int m_nCommonHeterozygous;
    int m_nCommonFilteredHeterozygous;
    int m_nFCHomozygous;
    int m_nFCHeterozygous;
    int m_nMCHomozygous;
    int m_nMCHeterozygous;
    int m_nChildUnique;
};


class CResultLog
{
    
public:
    
    //Sets the path of log file
    void SetLogPath(const std::string& a_rLogPath);
    
    //Records the result for given chromosome
    void LogStatistic(int a_nChromosomeId, int a_nTpCalled, int a_nTpBaseline, int a_nHalfTPCalled, int a_nHalfTPBaseline, int a_nFalsePositive, int a_nFalseNegative);
    
    //Records log entry for mendelian comparison
    void LogMendelianStatistic(bool a_bIsFatherChild, int a_nChromosomeId, int a_nTpCalled, int a_nTpBaseline, int a_nFalsePositive, int a_nFalseNegative);
    
    void LogMendelianIntersectionStatistic(int a_nChromosomeId,
                                           int a_nCommonHomozygous,
                                           int a_nCommonHeterozygous,
                                           int a_nCommonFilteredHeterozygous,
                                           int a_nFCHomozygous,
                                           int a_nFCHeterozgous,
                                           int a_nMCHomozygous,
                                           int a_nMCHeterozygous,
                                           int a_nChildUnique);
    
    //Write the results in log.txt file
    void WriteStatistics();
    
    //Write the results of mendelian violation check in log.txt file
    void WriteMendelianStatistics();
    
private:
    
    //Log entry array for ga4gh output [THE STANTARD OUTPUT]
    SLogEntry m_aResultEntries[CHROMOSOME_COUNT];
    
    //Log entry array for Mendelian Violation mode output
    SMendelianLogEntry m_aFatherChildLogEntries[CHROMOSOME_COUNT];
    SMendelianLogEntry m_aMotherChildLogEntries[CHROMOSOME_COUNT];
    SMendelianMergeLogEntry m_aChildMergeLogEntries[CHROMOSOME_COUNT];
    
    std::string m_aLogPath;
    
};



#endif // _C_RESULT_LOG_H_
