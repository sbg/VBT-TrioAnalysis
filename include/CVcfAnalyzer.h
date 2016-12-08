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
    
    
    //Maximum # of thread that program executes
    int m_nMaxThreadCount;
    //Maximum memory will be used by the program
    int m_nMaxMemoryInByte;
    
    //Configurations for Vcf comparison
    SConfig m_config;
    
};

#endif //C_VCF_ANALYZER_H
