//
//  CGraphVcfAnalyzer.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/10/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <thread>
#include "CGraphVcfAnalyzer.h"
#include "CGraphSplitOutputProvider.h"
#include "CPathReplay.h"

using namespace graphcomparison;

CGraphVcfAnalyzer::CGraphVcfAnalyzer()
{
    m_bIsPassFilterEnabled = false;
    m_bedFilePath = "";
    m_nThreadCount = DEFAULT_THREAD_COUNT;
    m_nMaxBasePairLength = DEFAULT_MAX_BP_LENGTH;
    m_nMaxPathSize = DEFAULT_MAX_PATH_SIZE;
    m_nMaxIterationCount = DEFAULT_MAX_ITERATION_SIZE;
}

int CGraphVcfAnalyzer::run(int argc, char** argv)
{
    std::time_t programStart;
    double duration;
    
    programStart = std::time(0);

    //Read command line parameters to m_config object
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return -1;
    
    //Initialize Variant providers which contains VCF and FASTA files
    m_provider.SetParameters(m_baseVcfPath, m_calledVcfPath, m_referencePath, m_bedFilePath, m_bIsPassFilterEnabled, m_nMaxBasePairLength);
    isSuccess = m_provider.InitializeReaders();
    
    if(!isSuccess)
        return -1;
    
    duration = std::difftime(std::time(0) ,programStart);
    std::cout << "Vcf read completed in " << duration << " secs" << std::endl;
    
    std::time_t bestpathAlgorithmStart = std::time(0);
    
    //Creates the threads according to given memory and process the data
    AssignJobsToThreads(m_nThreadCount);

    std::cout << "Generating Outputs..." << std::endl;
    CGraphSplitOutputProvider outputProvider;
    outputProvider.SetOutputVcfFolder(m_outputDirectoryPath);
    outputProvider.SetVcfNames(m_baseVcfPath, m_calledVcfPath);
    outputProvider.SetParameters(m_bedFilePath, m_bIsPassFilterEnabled, m_nMaxBasePairLength);
    outputProvider.SetExcludeIndexes(m_uniqueVariantsListPerChromosome);
    outputProvider.GenerateVcfs("/TPCalled.vcf", "/FP.vcf", false);
    outputProvider.GenerateVcfs("/TPBase.vcf", "/FN.vcf", true);

    duration = std::difftime(std::time(0), bestpathAlgorithmStart);
    std::cout << "Processing Chromosomes completed in " << duration << " secs" << std::endl;

    PrintLogs();
    
    duration = std::difftime(std::time(0), programStart);
    std::cout << std::endl << "Total execution time is " << duration << " secs" << std::endl;

    return 0;
}

int CGraphVcfAnalyzer::AssignJobsToThreads(int a_nThreadCount)
{
    //Get the list of chromosomes to be processed
    std::vector<duocomparison::SChrIdTuple> chromosomeListToProcess = m_provider.GetChromosomeIdTuples();
    
    int exactThreadCount = std::min(a_nThreadCount, (int)chromosomeListToProcess.size());
    
    //Allocate threads
    std::thread* pThreadPoolList = new std::thread[exactThreadCount];
    
    int threadPoolIt = 0;
    std::vector<duocomparison::SChrIdTuple> *chromosomeLists = new std::vector<duocomparison::SChrIdTuple>[exactThreadCount];
    
    //Divide tasks into threads
    for(int k = 0; k < (int)chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(chromosomeListToProcess[k]);
        threadPoolIt = (threadPoolIt+1) % exactThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < exactThreadCount; k++)
    {
        pThreadPoolList[k] = std::thread(&CGraphVcfAnalyzer::ThreadFunctionSPLIT, this, chromosomeLists[k]);
    }
    
    for(int k = 0; k < exactThreadCount; k++)
        pThreadPoolList[k].join();
    
    //Clear allocated memory
    delete[] chromosomeLists;
    delete[] pThreadPoolList;
    
    return exactThreadCount;
}

void CGraphVcfAnalyzer::ThreadFunctionSPLIT(std::vector<duocomparison::SChrIdTuple> a_aTuples)
{
    
    for(int k = 0; k < (int)a_aTuples.size(); k++)
    {
        std::vector<const CVariant*> varListBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId);
        std::vector<const CVariant*> varListCalled = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId);
        
        std::vector<const core::COrientedVariant*> ovarListBase = m_provider.GetOrientedVariantList(eBASE, a_aTuples[k].m_nBaseId);
        std::vector<const core::COrientedVariant*> ovarListCalled = m_provider.GetOrientedVariantList(eCALLED, a_aTuples[k].m_nCalledId);
        
        SContig ctg;
        m_provider.GetContig(a_aTuples[k].m_chrName, ctg);
        
        
        core::CPath path;
        int replayIterationCount = 0;
        
        std::vector<int> baseExcludedVariantIndexes;
        std::vector<int> calledExcludedVariantIndexes;
        
        while(replayIterationCount < GRAPH_REPLAY_MAX_ITERATION)
        {
            core::CPathReplay pathReplay(varListBase, varListCalled, ovarListBase, ovarListCalled);
            pathReplay.SetMaxPathAndIteration(m_nMaxPathSize, m_nMaxIterationCount);
            path = pathReplay.FindBestPath(ctg, false);
            
            const std::vector<const core::COrientedVariant*>& includedVarsBase = path.m_baseSemiPath.GetIncludedVariants();
            const std::vector<const core::COrientedVariant*>& includedVarsCall = path.m_calledSemiPath.GetIncludedVariants();
            
            baseExcludedVariantIndexes = CGraphVariantProvider::GetExcludedIndexes(varListBase, path.m_baseSemiPath.GetExcluded());
            calledExcludedVariantIndexes = CGraphVariantProvider::GetExcludedIndexes(varListCalled, path.m_calledSemiPath.GetExcluded());
            
            if(includedVarsBase.size() == 0 && includedVarsCall.size() == 0)
            {
                SGraphVarContainer container;
                container.chrName = a_aTuples[k].m_chrName;
                container.baseExcludedIndexes = baseExcludedVariantIndexes;
                container.calledExcludedIndexes = calledExcludedVariantIndexes;
                m_uniqueVariantsListPerChromosome[container.chrName] = container;
                std::cout << "Total Iteration : " << replayIterationCount + 1 << std::endl;
                break;
            }
            
            else
            {
                varListBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId, baseExcludedVariantIndexes);
                varListCalled = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId, calledExcludedVariantIndexes);
                
                ovarListBase   = m_provider.GetOrientedVariantList(eBASE,   a_aTuples[k].m_nBaseId,   baseExcludedVariantIndexes);
                ovarListCalled = m_provider.GetOrientedVariantList(eCALLED, a_aTuples[k].m_nCalledId, calledExcludedVariantIndexes);
                
                pathReplay.Clear();
                replayIterationCount++;
            }
            
        }
    }
}

bool CGraphVcfAnalyzer::ReadParameters(int argc, char** argv)
{
    
    const char* PARAM_BASE = "-base";
    const char* PARAM_CALLED = "-called";
    const char* PARAM_FILTER = "--pass-filter";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_BED = "-bed";
    const char* PARAM_HELP = "--help";
    const char* PARAM_OUTPUT_DIR = "-outDir";
    const char* PARAM_THREAD_COUNT = "-thread-count";
    const char* PARAM_MAX_PATH_SIZE = "-max-path-size";
    const char* PARAM_MAX_ITERATION_COUNT = "-max-iteration-count";
    const char* PARAM_MAX_BP_LENGTH = "-max-bp-length";
    
    bool bBaselineSet = false;
    bool bCalledSet = false;
    bool bReferenceSet = false;
    bool bOutputDirSet = false;
    
    //Start from index 2 since first parameter will be variant comparison mode indicator
    int it = 2;
    
    if(argc < 2)
    {
        PrintHelp();
        return false;
    }
    
    while(it < argc)
    {
        if(0 == strcmp(argv[it], PARAM_HELP))
        {
            PrintHelp();
            return false;
        }
        
        else if(0 == strcmp(argv[it], PARAM_BASE))
        {
            m_baseVcfPath = argv[it+1];
            bBaselineSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CALLED))
        {
            m_calledVcfPath = argv[it+1];
            bCalledSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REFERENCE))
        {
            m_referencePath = argv[it+1];
            bReferenceSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_OUTPUT_DIR))
        {
            m_outputDirectoryPath = argv[it+1];
            bOutputDirSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_FILTER))
        {
            m_bIsPassFilterEnabled = true;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_BED))
        {
            m_bedFilePath = argv[it+1];
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_THREAD_COUNT))
        {
            m_nThreadCount = atoi(argv[it+1]);
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MAX_PATH_SIZE))
        {
            m_nMaxPathSize = atoi(argv[it+1]);
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MAX_ITERATION_COUNT))
        {
            m_nMaxIterationCount = atoi(argv[it+1]);
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MAX_BP_LENGTH))
        {
            m_nMaxBasePairLength = atoi(argv[it+1]);
            it+=2;
        }
        
        else
        {
            std::cerr << "Unkonown parameter is detected: " << argv[it] << std::endl;
            return false;
        }
    }
    
    if(!bBaselineSet)
        std::cerr << "Baseline vcf file is not set" << std::endl;
    else if(!bCalledSet)
        std::cerr << "Called vcf file is not set" << std::endl;
    else if(!bReferenceSet)
        std::cerr << "Reference fasta file is not set" << std::endl;
    else if(!bOutputDirSet)
        std::cerr << "Output Directory is not set" << std::endl;
    
    
    return bCalledSet & bBaselineSet & bReferenceSet & bOutputDirSet;

}

void CGraphVcfAnalyzer::PrintLogs()
{
    const char separator = ' ';

    std::cout << std::endl;
    std::cout << "=====================================================================" << std::endl;
    auto chrIdTuples = m_provider.GetChromosomeIdTuples();

    std::cout << std::left << std::setw(20)  << std::setfill(separator) << "Chromosome";
    std::cout << std::left << std::setw(15)  << std::setfill(separator) << "TP Base";
    std::cout << std::left << std::setw(15)  << std::setfill(separator) << "TP Called";
    std::cout << std::left << std::setw(19)  << std::setfill(separator) << "False Positive";
    std::cout << std::left << std::setw(15)  << std::setfill(separator) << "False Negative";
    std::cout << std::endl;
    //std::cout << "Chr" << "\t" << "TP Base" << "\t" <<  "TP Called" << "\t" << "False Positive" << "\t" << "False Negative" << std::endl;
    
    for(int k = 0; k < chrIdTuples.size(); k++)
    {
        unsigned long excludedSizeBase = m_uniqueVariantsListPerChromosome[chrIdTuples[k].m_chrName].baseExcludedIndexes.size();
        unsigned long includedSizeBase = m_provider.GetVariantList(eBASE, chrIdTuples[k].m_nBaseId).size() - excludedSizeBase;
        unsigned long excludedSizeCalled = m_uniqueVariantsListPerChromosome[chrIdTuples[k].m_chrName].baseExcludedIndexes.size();
        unsigned long includedSizeCalled = m_provider.GetVariantList(eCALLED, chrIdTuples[k].m_nCalledId).size() - excludedSizeCalled;
        
        //std::cout << chrIdTuples[k].m_chrName << "\t" << includedSizeBase << "\t"  << includedSizeCalled << "\t" << excludedSizeCalled << "\t" << excludedSizeBase << std::endl;
        std::cout << std::left << std::setw(20)  << std::setfill(separator) << chrIdTuples[k].m_chrName;
        std::cout << std::left << std::setw(15)  << std::setfill(separator) << includedSizeBase;
        std::cout << std::left << std::setw(15)  << std::setfill(separator) << includedSizeCalled;
        std::cout << std::left << std::setw(19)  << std::setfill(separator) << excludedSizeCalled;
        std::cout << std::left << std::setw(15)  << std::setfill(separator) << excludedSizeBase;
        std::cout << std::endl;
    }
}

void CGraphVcfAnalyzer::PrintHelp() const
{
    std::cout << "==== VARIANT BENCHMARKING TOOL VERSION 1.0 (Beta) ==== " << std::endl;
    std::cout << "Based on paper: http://biorxiv.org/content/biorxiv/early/2015/08/02/023754.full.pdf" << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgdinc.com)" << std::endl;
    std::cout << "Please notify me if program fails or return unexpected results" << std::endl;
    std::cout << "COPYRIGHT (C) 2017 SBGD INC" << std::endl;
    std::cout << std::endl;
    std::cout << " --- GRAPH COMPARISON PARAMETERS --- " << std::endl;
    std::cout << "-base <baseline_vcf_path>    [Required.Add baseline VCF file.]" << std::endl;
    std::cout << "-called <called_vcf_path>    [Required.Add called VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-outDir <output_directory>   [Required.Add output directory]" << std::endl;
    std::cout << "--pass-filter                [Optional.Process only 'PASS' or '.' variants. Default value is false]" << std::endl;
    std::cout << "-thread-count                [Optional.Specify the number of threads that program will use. Default value is 2]" << std::endl;
    std::cout << "-max-bp-length               [*Optional.Specify the maximum base pair length of variant to process. Default value is 1000]" << std::endl;
    std::cout << "-max-path-size <size>        [*Optional.Specify the maximum size of path that core algorithm can store inside. Default value is 150,000]" << std::endl;
    std::cout << "-max-iteration-count <count> [*Optional.Specify the maximum iteration count that core algorithm can decide to include/exclude variant. Default value is 10,000,000]" << std::endl;
    std::cout << "(*) - for advanced usage" << std::endl;
    std::cout << std::endl;
}





