//
//  CVcfAnalyzer.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/8/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVcfAnalyzer.h"
#include "iostream"


void CVcfAnalyzer::Run(int argc, char** argv)
{


    
}



bool CVcfAnalyzer::ReadParameters(int argc, char** argv)
{

    const char* PARAM_BASE = "-b";
    const char* PARAM_CALLED = "-c";
    const char* PARAM_FILTER = "-f";
    const char* PARAM_REFERENCE = "-r";
    
    bool bBaselineSet = false;
    bool bCalledSet = false;
    bool bReferenceSet = false;
    
    int it = 1;
        
    while(it < argc)
    {
        if(0 == strcmp(argv[it], PARAM_BASE))
        {
            m_config.m_pBaseVcfFileName = argv[it+1];
            bBaselineSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CALLED))
        {
            m_config.m_pCalledVcfFileName = argv[it+1];
            bCalledSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REFERENCE))
        {
            m_config.m_pCalledVcfFileName = argv[it+1];
            bReferenceSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_FILTER))
        {
            if(0 == strcmp("none", argv[it+1]))
                m_config.m_bIsFilterEnabled = false;
            else
            {
                m_config.m_bIsFilterEnabled = true;
                m_config.m_pFilterName = argv[it+1];
            }
        }
        
        it += 2;
    }
    
    if(!bBaselineSet)
        std::cout << "Baseline vcf file is not set" << std::endl;
    else if(!bCalledSet)
        std::cout << "Called vcf file is not set" << std::endl;
    
    
    
    return true;
}

void CVcfAnalyzer::PrintHelp() const
{
    std::cout << "==== SBG VCF COMPARISON TOOL ==== " << std::endl;
    std::cout << "AUTHOR: Berke Cagkan Toptas (berke.toptas@sbgenomics.com)" << std::endl;
    std::cout << "VERSION: 1.0 " << "DATE: DEC 8 2016" << std::endl;
    std::cout << "COPYRIGHT 2016 SEVEN BRIDGES GENOMICS. ALL RIGHTS RESERVED." << std::endl;
    
    std::cout << std::endl << std::endl;
    std::cout << "<<< INSTRUCTIONS >>>" << std::endl;
    std::cout << "-b baseline_vcf_path (required*)" << std::endl;
    std::cout << "-c called_vcf_path (required*)" << std::endl;
    std::cout << "-r reference_fasta_path (required*)" << std::endl;
    std::cout << "-f filter_name (optional. Default is PASS. Use none for unfilter vcf)" << std::endl;
}
