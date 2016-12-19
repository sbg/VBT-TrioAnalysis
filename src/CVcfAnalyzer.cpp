//  CVcfAnalyzer.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/8/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVcfAnalyzer.h"
#include "iostream"
#include <math.h>


void CVcfAnalyzer::Run(int argc, char** argv)
{
    
    std::clock_t start;
    double duration;
    
    //Read command line parameters to m_config object
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return;

   // m_config.m_pBaseVcfFileName = "/Users/c1ms21p6h3qk/Desktop/BigTestData/HG00171-30x.vcf";
   // m_config.m_pCalledVcfFileName = "/Users/c1ms21p6h3qk/Desktop/BigTestData/HG00171-30x.vcf";
    m_config.m_pFastaFileName = "/Users/c1ms21p6h3qk/Desktop/BigTestData/human_g1k_v37_decoy.fasta";
    m_config.m_pBaseVcfFileName = "/Users/c1ms21p6h3qk/Desktop/BigTestData/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf.vcf";
    m_config.m_pCalledVcfFileName = "/Users/c1ms21p6h3qk/Desktop/BigTestData/gral0.9.sorted.and_more.concat.vcf";
    
    start = std::clock();
    
    //Initialize Variant providers which contains VCF and FASTA files
    isSuccess = m_provider.InitializeReaders(m_config);
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Vcf and fasta Parser read completed in " << duration << " secs" << std::endl;

    if(!isSuccess)
        return;
    
    start = std::clock();
    SetThreadsCustom(8 * 1024);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Program Completed in " << duration << " secs" << std::endl;
}

void CVcfAnalyzer::SetThreadsPlatform()
{
    std::vector<int> chromosomeIds;
 
    std::thread threadPool[CHROMOSOME_COUNT];
    
    m_provider.GetUniqueChromosomeIds(chromosomeIds);

    for(int k = 0; k < chromosomeIds.size(); k++)
    {
        threadPool[k] = std::thread(&CVcfAnalyzer::ThreadFunc, this, chromosomeIds[k]);
    }
    
    for(int k = 0; k < chromosomeIds.size(); k++)
    {
        threadPool[k].join();
    }
}

void CVcfAnalyzer::SetThreadsCustom(int a_nMemoryInMB)
{
    std::vector<int> chromosomeIds;
    m_provider.GetUniqueChromosomeIds(chromosomeIds);
    
    //IF WE HAVE TOO MUCH MEMORY GO RUN ALL THREADS CONSECUTIVELY
    if(a_nMemoryInMB > 32 * 1024)
    {
        SetThreadsPlatform();
        return;
    }

    //IF WE HAVE LESS MEMORY USE SINGLE THREADING
    if(a_nMemoryInMB < 2 * 1024 || chromosomeIds.size() < 4)
    {
        for(int k = 0; k < chromosomeIds.size(); k++)
            ThreadFunc(chromosomeIds[k]);
        return;
    }

    //ELSE ALLOCATE 1 GB MEMORY (ON AVG) FOR EACH THREAD
    const int maxThreadCount = ceil(a_nMemoryInMB / 1024);
    std::vector<int>* chrArrs = new std::vector<int>[maxThreadCount];

    for(int k = 0, p = 0; k < chromosomeIds.size(); k++)
    {
        chrArrs[p].push_back(chromosomeIds[k]);
        p++;
        p = p % maxThreadCount;
    }
    
    std::thread *threadPool = new std::thread[maxThreadCount];
    
    for(int k = 0; k < maxThreadCount; k++)
    {
        threadPool[k] = std::thread(&CVcfAnalyzer::ThreadFunc2, this, chrArrs[k]);
    }
        
    for(int k = 0; k < maxThreadCount; k++)
    {
        threadPool[k].join();
    }
}

void CVcfAnalyzer::ThreadFunc(int a_nChromosomeId)
{
    CPathReplay pathReplay;
    SContig ctg;
    
    pathReplay.SetVariantProvider(m_provider);
    m_provider.GetContig(a_nChromosomeId, ctg);
    
    m_aBestPaths[a_nChromosomeId] = pathReplay.FindBestPath(ctg);
    
    std::cout << "===CHROMOSOME " << a_nChromosomeId + 1 << "===" << std::endl;
    std::cout << "Called Included Count: " << m_aBestPaths[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants().size() << std::endl;
    std::cout << "Called Excluded Count: " << m_aBestPaths[a_nChromosomeId].m_calledSemiPath.GetExcluded().size() << std::endl;
    std::cout << "Baseline Included Count: " << m_aBestPaths[a_nChromosomeId].m_baseSemiPath.GetIncludedVariants().size() << std::endl;
    std::cout << "Baseline Excluded Count: " << m_aBestPaths[a_nChromosomeId].m_baseSemiPath.GetExcluded().size() << std::endl;
    return;
}

void CVcfAnalyzer::ThreadFunc2(std::vector<int> a_nChrArr)
{
    for(int k = 0; k < a_nChrArr.size(); k++)
    {
        CPathReplay pathReplay;
        SContig ctg;
        
        pathReplay.SetVariantProvider(m_provider);
        m_provider.GetContig(a_nChrArr[k], ctg);
        m_aBestPaths[a_nChrArr[k]] = pathReplay.FindBestPath(ctg);
        
        std::cout << "===CHROMOSOME " << a_nChrArr[k] + 1 << "===" << std::endl;
        std::cout << "Called Included Count: " << m_aBestPaths[a_nChrArr[k]].m_calledSemiPath.GetIncludedVariants().size() << std::endl;
        std::cout << "Called Excluded Count: " << m_aBestPaths[a_nChrArr[k]].m_calledSemiPath.GetExcluded().size() << std::endl;
        std::cout << "Baseline Included Count: " << m_aBestPaths[a_nChrArr[k]].m_baseSemiPath.GetIncludedVariants().size() << std::endl;
        std::cout << "Baseline Excluded Count: " << m_aBestPaths[a_nChrArr[k]].m_baseSemiPath.GetExcluded().size() << std::endl;
    }
}


bool CVcfAnalyzer::ReadParameters(int argc, char** argv)
{

    const char* PARAM_BASE = "-base";
    const char* PARAM_CALLED = "-called";
    const char* PARAM_FILTER = "-filter";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_HELP = "-help";
    const char* PARAM_SAMPLE_BASE = "-SampleBase";
    const char* PARAM_SAMPLE_CALLED = "-SampleCalled";
    const char* PARAM_SNP_ONLY = "-SNP_ONLY";
    const char* PARAM_INDEL_ONLY = "-INDEL_ONLY";
    
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
            m_config.m_pFastaFileName = argv[it+1];
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
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_BASE))
        {
            m_config.m_bBaseSampleEnabled = true;
            m_config.m_pBaseSample = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_CALLED))
        {
            m_config.m_bCalledSampleEnabled = true;
            m_config.m_pCalledSample = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_SNP_ONLY))
        {
            m_config.m_bSNPOnly = true;
            it--;
        }
        
        else if(0 == strcmp(argv[it], PARAM_INDEL_ONLY))
        {
            m_config.m_bINDELOnly = true;
            it--;
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
    std::cout << "==== SBG VCF COMPARISON TOOL VERSION 1.0 (Beta) ==== " << std::endl;
    std::cout << "Based on paper: http://biorxiv.org/content/biorxiv/early/2015/08/02/023754.full.pdf" << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgenomics.com)" << std::endl;
    std::cout << "Please notify me if program fails or return unexpected results" << std::endl;
    std::cout << "COPYRIGHT (C) 2016 SEVEN BRIDGES GENOMICS." << std::endl;
    std::cout << std::endl;
    std::cout << " --- PARAMETERS --- " << std::endl;
    std::cout << "-base <baseline_vcf_path>    [Required.Add baseline VCF file.]" << std::endl;
    std::cout << "-called <called_vcf_path>    [Required.Add called VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-filter <filter_name>        [Optional.Filter variants based on filter column. Default value is PASS. Use 'none' to unfilter]" << std::endl;
    std::cout << "-SNP_ONLY                    [Optional.Filter INDELs out from both base and called VCF file.]" << std::endl;
    std::cout << "-INDEL_ONLY                  [Optional.Filter SNPs out from both base and called VCF file.]" << std::endl;
    std::cout << "-SampleBase <sample_name>    [Optional.Read only the given sample in base VCF. Default is the first sample.]" << std::endl;
    std::cout << "-SampleCalled <sample_name>  [Optional.Read only the given sample in called VCF. Default is the first sample.]" << std::endl;
    std::cout << std::endl;
    std::cout << "Example Commands:" << std::endl;
    std::cout << "./sbgVcfComp -called called.vcf -base base.vcf -ref reference.fa -filter none" << std::endl;
    std::cout << "./sbgVcfComp -called called.vcf -base base.vcf -ref reference.fa -filter PASS -SNP_ONLY -SampleBase sample0 -SampleCalled sample07" << std::endl;
}
