//
//  CMendelianResultLog.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/28/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_MENDELIAN_RESULT_LOG_H_
#define _C_MENDELIAN_RESULT_LOG_H_

#include "Constants.h"
#include <string>
#include <vector>
#include "SChrIdTriplet.h"

namespace mendelian
{

struct SMendelianBestPathLogEntry
{
    bool m_bIsNull = true;
    int m_nTpCalled;
    int m_nTpBase;
    int m_nFP;
    int m_nFn;
    std::string m_chrName;
};

struct SMendelianDetailedLogEntry
{
    void clear();
    
    //0 -> compliant 1 -> violation 2 -> no call parent 3-> no call child
    int m_nSNP[4];
    int m_nInsertSmall[4];
    int m_nInsertMedium[4];
    int m_nInsertLarge[4];
    int m_nDeleteSmall[4];
    int m_nDeleteMedium[4];
    int m_nDeleteLarge[4];
    int m_nComplexSmall[4];
    int m_nComplexMedium[4];
    int m_nComplexLarge[4];
};


struct SMendelianShortLogEntry
{
    int m_nSNPcompliant;
    int m_nSNPviolation;
    int m_nINDELcompliant;
    int m_nINDELviolation;
    std::string m_chrName;
};


struct SMendelianDetailedLogGenotypes
{
    void clear();
    
    //For 10 variant categories, 27 different genotype possibilities for compliant and violation(27 + MultiAllelic)
    int genotypesCompliant[10][27+1];
    
    int genotypesViolation[10][27+1];
};


class CMendelianResultLog
{
    
public:
    
    //Sets the result log file directory
    void SetLogDirectory(std::string a_rLogDirectory);
    
    
    //For given chromosome it reports the result of best path algorithm for parent-child
    void LogBestPathStatistic(bool a_bIsFatherChild, SChrIdTriplet a_triplet, int a_nTpCalled, int a_nTpBaseline, int a_nFalsePositive, int a_nFalseNegative);
    
    //For given chromosome it reports the short result table (SNP and INDEL counts for non 0/0 variants)
    void LogShortReport(std::string& a_rChrName, int a_nSNPcompliant, int a_nSNPviolation, int a_nINDELcompliant, int a_nINDELviolation);
    
    //For given chromosome logs the number of compliant and violations for each of 10 variant category
    void LogDetailedReport(SMendelianDetailedLogEntry& a_rLogEntry);
    
    //For ALL chromosomes, logs the genotype possibility counts for each 10 variant category
    void LogGenotypeMatrix(SMendelianDetailedLogGenotypes& a_rLogEntry);
    
    //For ALL chromosomes, logs the skipped variant counts for father/mother/child
    void LogSkippedVariantCounts(int a_nChildSkipped, int a_nFatherSkipped, int a_nMotherSkipped);
    
    //For ALL chromosomes, logs the filtered variant counts for father/mother/child
    void LogFilteredComplexVariantCounts(int a_nChildFiltered, int a_nFatherFiltered, int a_nMotherFiltered);
    
    //Write Statistic of parent-child comparison from Best Path Algorithm
    void WriteBestPathStatistics(const std::string& a_rPrefixName);
    
    //Write Detailed report (Extended statistic of each variant group) to the detailedLog.txt
    void WriteDetailedReportTable(const std::string& a_rPrefixName);
    
    //Write Detailed report in TSV format
    void WriteDetailedReportTabDelimited(const std::string& a_rPrefixName);
    
    //Write a summary report to the shortLog.txt
    void WriteShortReportTable(const std::string& a_rPrefixName);

private:
    
    //Write a row to the 4 column result log table in SQL console table format
    void WriteRow(std::ofstream* a_pOutput, const std::string& a_name, int* a_pColumns);
    
    //Write a row to the 4 column result log table in SQL console table format
    void WriteRowSum(std::ofstream* a_pOutput, const std::string& a_name);
    

    //Write genotypes row to the result log table in SQL console table format
    void WriteGenotypeRow(std::ofstream* a_pOutput,
                          const std::string& a_name,
                          SMendelianDetailedLogGenotypes& a_genotypes,
                          int a_nVariantCategoryIndex,
                          bool a_bIsCompliantTable);
    
    //Write genotypes row to the result log table in SQL console table format
    void WriteGenotypeRowSum(std::ofstream* a_pOutput, const std::string& a_name, bool a_bIsCompliantTable);


    //Stores best Path algorithm results (TP/FP/FN of comparisons)
    std::vector<SMendelianBestPathLogEntry> m_aMotherChildLogEntries;
    std::vector<SMendelianBestPathLogEntry> m_aFatherChildLogEntries;

    //Detailed log entry for each chromosome
    SMendelianDetailedLogEntry m_DetailedLogEntries;
    
    //Detailed log genotype TOTAL possibility array for each variant category
    SMendelianDetailedLogGenotypes m_DetailedLogGenotypes;
    
    //Short log entry for each chromosome
    std::vector<SMendelianShortLogEntry> m_aShortLogEntries;
    
    //File directory where the result logs will be stored
    std::string m_aLogDirectory;
    
    //Skipped variant counts [Complex Region]
    int m_nTotalSkippedCountFather;
    int m_nTotalSkippedCountMother;
    int m_nTotalSkippedCountChild;
    
    //Filtered  variant counts [Complex variant representation]
    int m_nTotalNonAssessedVarCountFather;
    int m_nTotalNonAssessedVarCountMother;
    int m_nTotalNonAssessedVarCountChild;
    
};

}
    
#endif // _C_MENDELIAN_RESULT_LOG_H_
