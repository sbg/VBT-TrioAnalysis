//
//  CVcfAnalyzer.h
//  VCFComparison
//
//  Created by Berke.Toptas on 12/8/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef C_VCF_ANALYZER_H
#define C_VCF_ANALYZER_H

#include <thread>
#include "CPathReplay.h"
#include "SConfig.h"
#include "CVariantProvider.h"
#include "CResultLog.h"

class CVcfAnalyzer
{

public:
    
    //Run the application
    void Run(int argc, char** argv);
    
    
private:

    //Read Parameters from command line. If the mandatory arguments are given, return true.
    bool ReadParameters(int argc, char** argv);
    
    //Prints the help menu at console
    void PrintHelp() const;
    
    //Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
    int AssignJobsToThreads(int a_nThreadCount);
    
    //Calculate the thread count that will be generated and initialize the threads
    void SetThreadsCustom(int a_nMemoryInMB);
    
    //Function that process chromosome in bulk for SPLIT mode (process either genotype or allele match)
    void ThreadFunctionSPLIT(std::vector<int> a_nChrArr, bool a_bIsGenotypeMatch);

    //Function that process chromosome in bulk for GA4GH mode (process both genotype and allele matches)
    void ThreadFunctionGA4GH(std::vector<int> a_nChrArr);

    
    void PrintVariants(std::string a_outputDirectory, std::string a_FileName, const std::vector<const COrientedVariant*>& a_rOvarList) const;
    void PrintVariants(std::string a_outputDirectory, std::string a_FileName, const std::vector<const CVariant*>& a_rVarList) const;
    
    
    //Configurations for Vcf comparison
    SConfig m_config;
    
    //Object to write variant statistic into a file
    CResultLog m_resultLogger;
    
    //Variant provider instance
    CVariantProvider m_provider;
    
    //Best Paths written by each thread for each unique chromosome exists
    CPath m_aBestPaths[CHROMOSOME_COUNT];
    
    //Best Paths written by each thread to find Allele matches for each unique chromosome exists
    CPath m_aBestPathsAllele[CHROMOSOME_COUNT];
    
    //Thread pool we have for multitasking by per chromosome
    std::thread *m_pThreadPool;

};

#endif //C_VCF_ANALYZER_H
