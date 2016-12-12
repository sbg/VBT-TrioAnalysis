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
    //Read command line parameters
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return;
    
    //Initialize FASTA reader
    isSuccess = m_fastaReader.Open(m_config.m_pFastaFileName);
  
    if(!isSuccess)
        return;
    
    //Initialize Variant providers
    m_provider.SetFastaReader(m_fastaReader);
    m_provider.InitializeReaders(m_config);
    
}


bool CVcfAnalyzer::ReadParameters(int argc, char** argv)
{

    const char* PARAM_BASE = "-base";
    const char* PARAM_CALLED = "-called";
    const char* PARAM_FILTER = "-filter";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_HELP = "-help";
    
    bool bBaselineSet = false;
    bool bCalledSet = false;
    bool bReferenceSet = false;
    
    int it = 1;
        
    while(it < argc)
    {
        if(0 == strcmp(argv[it], PARAM_HELP))
        {
            PrintHelp();
            return false;
        }
        
        else if(0 == strcmp(argv[it], PARAM_BASE))
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
    else if(!bReferenceSet)
        std::cout << "Reference fasta file is not set" << std::endl;
    
    
    return bCalledSet & bBaselineSet & bReferenceSet;
}

void CVcfAnalyzer::PrintHelp() const
{
    std::cout << "==== SBG VCF COMPARISON TOOL VERSION 1.0 ==== " << std::endl;
    std::cout << "Based on paper: http://biorxiv.org/content/biorxiv/early/2015/08/02/023754.full.pdf" << std::endl;
    std::cout << "AUTHOR: Berke Cagkan Toptas (berke.toptas@sbgenomics.com)" << std::endl;
    std::cout << "COPYRIGHT (C) 2016 SEVEN BRIDGES GENOMICS." << std::endl;
    std::cout << std::endl;
    std::cout << " --- PARAMETERS --- " << std::endl;
    std::cout << "-help                      [Prints the parameter options]" << std::endl;
    std::cout << "-base baseline_vcf_path    [Required.Add baseline VCF file.]" << std::endl;
    std::cout << "-called called_vcf_path    [Required.Add called VCF file.]" << std::endl;
    std::cout << "-ref reference_fasta_path  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-filter filter_name        [Optional.Filter variants based on filter column. Default value is PASS. Use 'none' to unfilter]" << std::endl;
}
