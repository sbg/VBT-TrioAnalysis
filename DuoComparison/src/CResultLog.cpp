//
//  CResultLog.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/29/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CResultLog.h"
#include <fstream>


void CResultLog::SetLogPath(const std::string& a_rLogFolder)
{
    m_aLogPath = a_rLogFolder + "/log.txt";
}


//Records the result for given chromosome
void CResultLog::LogStatistic(std::string a_chromosomeName, int a_nTpCalled, int a_nTpBaseline, int a_nHalfTPCalled, int a_nHalfTPBaseline, int a_nFalsePositive, int a_nFalseNegative)
{
    SLogEntry entry;
    
    entry.m_nTpCalled = a_nTpCalled;
    entry.m_nTpBase = a_nTpBaseline;
    entry.m_nHalfTpCalled = a_nHalfTPCalled;
    entry.m_nHalfTpBase = a_nHalfTPBaseline;
    entry.m_nFp = a_nFalsePositive;
    entry.m_nFn = a_nFalseNegative;
    entry.m_bIsNull = false;
    entry.m_chrName = a_chromosomeName;
    
    m_aResultEntries.push_back(entry);
}

//Write the results in log.txt file
void CResultLog::WriteStatistics(int a_nMode)
{
    //a_nMode :  0- Genotype Match (SPLIT) 1-Allele Match (SPLIT)  2- Both (GA4GH)
    
    std::ofstream outputLog;
    outputLog.open(m_aLogPath);
    
    if(a_nMode != 1)
    {
        outputLog << "====== GENOTYPE MATCHING MODE (GA4GH Method 3) ======" << std::endl;
        outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
        outputLog << "\t" << "F-measure" << std::endl;

        for(int k = 0; k < (int)m_aResultEntries.size(); k++)
        {
            if(m_aResultEntries[k].m_bIsNull)
                continue;
            
            int TPbase = m_aResultEntries[k].m_nTpBase;
            int TPcalled = m_aResultEntries[k].m_nTpCalled;
            int FP = m_aResultEntries[k].m_nFp + m_aResultEntries[k].m_nHalfTpCalled;
            int FN = m_aResultEntries[k].m_nFn + m_aResultEntries[k].m_nHalfTpBase;
            
            outputLog << m_aResultEntries[k].m_chrName << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
            
            double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
            double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
            double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
        
            outputLog.precision(4);
            outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
        }
    }
    
    if(a_nMode == 2)
        outputLog << std::endl << std::endl;
    
    if(a_nMode != 0)
    {
        outputLog << "====== ALLELE MATCHING MODE (GA4GH Method 2) ======" << std::endl;
        outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
        outputLog << "\t" << "F-measure" << std::endl;
        
        for(int k = 0; k < (int)m_aResultEntries.size(); k++)
        {
            if(m_aResultEntries[k].m_bIsNull)
                continue;
            
            int TPbase = m_aResultEntries[k].m_nTpBase + m_aResultEntries[k].m_nHalfTpBase;
            int TPcalled = m_aResultEntries[k].m_nTpCalled + m_aResultEntries[k].m_nHalfTpCalled;
            int FP = m_aResultEntries[k].m_nFp;
            int FN = m_aResultEntries[k].m_nFn;
            
            outputLog << m_aResultEntries[k].m_chrName << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
            
            double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
            double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
            double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
            
            outputLog.precision(4);
            outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
        }
    }
    
    if(a_nMode == 2)
    {
        int totalBaseTP = 0;
        for(int k= 0; k < (int)m_aResultEntries.size(); k++)
            totalBaseTP += m_aResultEntries[k].m_nTpBase;
        
        int totalCalledTP = 0;
        for(int k= 0; k < (int)m_aResultEntries.size(); k++)
            totalCalledTP += m_aResultEntries[k].m_nTpCalled;
        
        outputLog << std::endl << std::endl;
        outputLog << "Total TP Base :" << totalBaseTP << std::endl;
        outputLog << "Total TP Called:" << totalCalledTP << std::endl;
    }
    
    outputLog.close();
}





