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
void CResultLog::LogStatistic(int a_nChromosomeId, int a_nTpCalled, int a_nTpBaseline, int a_nHalfTPCalled, int a_nHalfTPBaseline, int a_nFalsePositive, int a_nFalseNegative)
{
    m_aResultEntries[a_nChromosomeId].m_nTpCalled = a_nTpCalled;
    m_aResultEntries[a_nChromosomeId].m_nTpBase = a_nTpBaseline;
    m_aResultEntries[a_nChromosomeId].m_nHalfTpCalled = a_nHalfTPCalled;
    m_aResultEntries[a_nChromosomeId].m_nHalfTpBase = a_nHalfTPBaseline;
    m_aResultEntries[a_nChromosomeId].m_nFp = a_nFalsePositive;
    m_aResultEntries[a_nChromosomeId].m_nFn = a_nFalseNegative;
    m_aResultEntries[a_nChromosomeId].m_bIsNull = false;
}

void CResultLog::LogMendelianStatistic(bool a_bIsFatherChild, int a_nChromosomeId, int a_nTpCalled, int a_nTpBaseline, int a_nFalsePositive, int a_nFalseNegative)
{
    SMendelianLogEntry* entry;
    
    if(true == a_bIsFatherChild)
        entry = &m_aFatherChildLogEntries[a_nChromosomeId];
    else
        entry = &m_aMotherChildLogEntries[a_nChromosomeId];

    entry->m_nTpBase = a_nTpBaseline;
    entry->m_nTpCalled = a_nTpCalled;
    entry->m_nFn = a_nFalseNegative;
    entry->m_nFP = a_nFalsePositive;
    entry->m_bIsNull = false;
}


//Write the results in log.txt file
void CResultLog::WriteStatistics()
{
    std::ofstream outputLog;
    outputLog.open(m_aLogPath);
    
    outputLog << "====== GENOTYPE MATCHING MODE (GA4GH Method 3) ======" << std::endl;
    outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
    outputLog << "\t" << "F-measure" << std::endl;

    for(int k = 0; k < 25; k++)
    {
        if(m_aResultEntries[k].m_bIsNull)
            continue;
        
        if(k < 22)
            outputLog << "Chr" << k+1;
        else if(k == 23)
            outputLog << "ChrX";
        else if(k == 24)
            outputLog << "ChrY";
        else
            outputLog << "ChrM";
        
        int TPbase = m_aResultEntries[k].m_nTpBase;
        int TPcalled = m_aResultEntries[k].m_nTpCalled;
        int FP = m_aResultEntries[k].m_nFp + m_aResultEntries[k].m_nHalfTpCalled;
        int FN = m_aResultEntries[k].m_nFn + m_aResultEntries[k].m_nHalfTpBase;
        
        
        outputLog << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
        
        double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
        double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
        double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
    
        outputLog.precision(4);
        outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
    }
    
    outputLog << std::endl << std::endl;
    outputLog << "====== ALLELE MATCHING MODE (GA4GH Method 2) ======" << std::endl;
    outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
    outputLog << "\t" << "F-measure" << std::endl;
    
    for(int k = 0; k < 25; k++)
    {
        if(m_aResultEntries[k].m_bIsNull)
            continue;
        
        if(k < 22)
            outputLog << "Chr" << k+1;
        else if(k == 23)
            outputLog << "ChrX";
        else if(k == 24)
            outputLog << "ChrY";
        else
            outputLog << "ChrM";
        
        int TPbase = m_aResultEntries[k].m_nTpBase + m_aResultEntries[k].m_nHalfTpBase;
        int TPcalled = m_aResultEntries[k].m_nTpCalled + m_aResultEntries[k].m_nHalfTpCalled;
        int FP = m_aResultEntries[k].m_nFp;
        int FN = m_aResultEntries[k].m_nFn;
        
        outputLog << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
        
        double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
        double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
        double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
        
        outputLog.precision(4);
        outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
    }
    
    
    
    outputLog.close();
}

void CResultLog::WriteMendelianStatistics()
{
    std::ofstream outputLog;
    outputLog.open(m_aLogPath);
    
    outputLog << "====== FATHER CHILD COMPARISON ======" << std::endl;
    outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
    outputLog << "\t" << "F-measure" << std::endl;
    
    for(int k = 0; k < 25; k++)
    {
        if(m_aFatherChildLogEntries[k].m_bIsNull)
            continue;
        
        if(k < 22)
            outputLog << "Chr" << k+1;
        else if(k == 23)
            outputLog << "ChrX";
        else if(k == 24)
            outputLog << "ChrY";
        else
            outputLog << "ChrM";
        
        int TPbase = m_aFatherChildLogEntries[k].m_nTpBase;
        int TPcalled = m_aFatherChildLogEntries[k].m_nTpCalled;
        int FP = m_aFatherChildLogEntries[k].m_nFP;
        int FN = m_aFatherChildLogEntries[k].m_nFn;
        
        outputLog << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
        
        double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
        double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
        double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
        
        outputLog.precision(4);
        outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
    }
    
    outputLog << "====== MOTHER CHILD COMPARISON ======" << std::endl;
    outputLog << "ID" << "\t" << "True-Pos-Called" << "\t\t" << "True-Pos-Baseline" << "\t" << "False-Pos" << "\t" << "False-Neg" << "\t" << "Precision" << "\t" << "Recall";
    outputLog << "\t" << "F-measure" << std::endl;
    
    for(int k = 0; k < 25; k++)
    {
        if(m_aMotherChildLogEntries[k].m_bIsNull)
            continue;
        
        if(k < 22)
            outputLog << "Chr" << k+1;
        else if(k == 23)
            outputLog << "ChrX";
        else if(k == 24)
            outputLog << "ChrY";
        else
            outputLog << "ChrM";
        
        int TPbase = m_aMotherChildLogEntries[k].m_nTpBase;
        int TPcalled = m_aMotherChildLogEntries[k].m_nTpCalled;
        int FP = m_aMotherChildLogEntries[k].m_nFP;
        int FN = m_aMotherChildLogEntries[k].m_nFn;
        
        
        outputLog << "\t" << TPcalled << "\t" << TPbase << "\t" << FP << "\t" << FN;
        
        double Precision = static_cast<double>(TPbase) / static_cast<double>(TPbase + FP);
        double Recall = static_cast<double>(TPbase) / static_cast<double>(TPbase + FN);
        double Fmeasure = Precision + Recall == 0.0 ? 0.0 : (2.0 * Precision * Recall) / (Precision + Recall);
        
        outputLog.precision(4);
        outputLog << "\t" << Precision << "\t" << Recall << "\t" << Fmeasure << std::endl;
    }
}





