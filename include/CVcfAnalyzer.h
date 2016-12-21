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
#include "Sconfig.h"
#include "CVariantProvider.h"

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
    
    //Process all chromosomes in parallel. Called if there is enough memory and core
    void SetThreadsPlatform();
    
    //Calculate the thread count that will be generated and initialize the threads
    void SetThreadsCustom(int a_nMemoryInMB);
    
    //Function that will be send to each thread
    void ThreadFunc(int a_nChromosomeId);
    
    //Function that process chromosome in bulk
    void ThreadFunc2(std::vector<int> a_nChrArr);
    
    //Maximum # of thread that program executes
    int m_nMaxThreadCount;
    //Maximum memory will be used by the program
    int m_nMaxMemoryInByte;
    
    //Configurations for Vcf comparison
    SConfig m_config;
    
    //Variant provider instance
    CVariantProvider m_provider;
    
    //Best Paths written by each thread for each unique chromosome exists
    CPath m_aBestPaths[CHROMOSOME_COUNT];

};

#endif //C_VCF_ANALYZER_H
