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
#include "CSyncPoint.h"

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
        //List that store the base Oriented variant tuples (In the order of genotype)
        std::deque<core::COrientedVariant> BaseOrientedVariantList;
        //List that store the called Oriented variant tuples (In the order of genotype)
        std::deque<core::COrientedVariant> CalledOrientedVariantList;

        m_provider.FillOrientedVariantList(a_aTuples[k], BaseOrientedVariantList, CalledOrientedVariantList);
        
        std::vector<const CVariant*> varListBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId);
        std::vector<const CVariant*> varListCalled = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId);
        
        std::vector<const core::COrientedVariant*> ovarListBase = m_provider.GetOrientedVariantList(BaseOrientedVariantList);
        std::vector<const core::COrientedVariant*> ovarListCalled = m_provider.GetOrientedVariantList(CalledOrientedVariantList);
        
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
            
            if(includedVarsBase.size() == 0 && includedVarsCall.size() == 0)
            {
                SGraphVarContainer container;
                container.chrName = a_aTuples[k].m_chrName;
                
                container.baseIncludedIndexes = m_provider.GetVariantIndexesByStatus(eBASE, a_aTuples[k].m_nBaseId, eALLELE_MATCH);
                container.baseExcludedIndexes = m_provider.GetVariantIndexesByStatus(eBASE, a_aTuples[k].m_nBaseId, eNO_MATCH);
                
                container.calledIncludedIndexes = m_provider.GetVariantIndexesByStatus(eCALLED, a_aTuples[k].m_nCalledId, eALLELE_MATCH);
                container.calledExcludedIndexes = m_provider.GetVariantIndexesByStatus(eCALLED, a_aTuples[k].m_nCalledId, eNO_MATCH);
                
                m_uniqueVariantsListPerChromosome[container.chrName] = container;
                std::cout << "Total Iteration : " << replayIterationCount + 1 << std::endl;
                break;
            }
            
            else
            {
                //Get the excluded variants for base and called graph
                baseExcludedVariantIndexes = CGraphVariantProvider::GetExcludedIndexes(varListBase, path.m_baseSemiPath.GetExcluded());
                calledExcludedVariantIndexes = CGraphVariantProvider::GetExcludedIndexes(varListCalled, path.m_calledSemiPath.GetExcluded());
                
                std::vector<const CVariant*> excludedVarsBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId, baseExcludedVariantIndexes);
                std::vector<const CVariant*> excludedVarsCall = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId, calledExcludedVariantIndexes);
                
                //Set Variant status
                m_provider.SetVariantStatus(excludedVarsBase, EVariantMatch::eNO_MATCH);
                m_provider.SetVariantStatus(excludedVarsCall, EVariantMatch::eNO_MATCH);
                m_provider.SetVariantStatus(includedVarsBase, EVariantMatch::eALLELE_MATCH);
                m_provider.SetVariantStatus(includedVarsCall, EVariantMatch::eALLELE_MATCH);

                std::vector<core::CSyncPoint> syncPointList;
                GetSyncPointList(includedVarsBase, includedVarsCall, excludedVarsBase, excludedVarsCall, path.m_aSyncPointList, syncPointList);
                
                std::vector<const CVariant*> tmpBase;
                std::vector<const CVariant*> tmpCalled;
                FilterGuaranteedUniqueVariants(syncPointList, tmpBase, tmpCalled);
                
                varListBase = tmpBase;
                varListCalled = tmpCalled;
                
                ovarListBase   = m_provider.GetOrientedVariantList(BaseOrientedVariantList,   varListBase);
                ovarListCalled = m_provider.GetOrientedVariantList(CalledOrientedVariantList, varListCalled);
                
                pathReplay.Clear();
                replayIterationCount++;
            }
            
        }
    }
}

void CGraphVcfAnalyzer::GetSyncPointList(const std::vector<const core::COrientedVariant*> a_rBaseIncluded,
                                         const std::vector<const core::COrientedVariant*> a_rCalledIncluded,
                                         const std::vector<const CVariant*>& a_rBaseExcluded,
                                         const std::vector<const CVariant*>& a_rCalledExcluded,
                                         const std::vector<int>& a_rSyncPointCoordinates,
                                         std::vector<core::CSyncPoint>& a_rSyncPointList)
{
    
    int baseIncludedItr = 0;
    int baseExcludedItr = 0;
    int calledIncludedItr = 0;
    int calledExcludedItr = 0;
    
    for(int k = 0; k < (int)a_rSyncPointCoordinates.size(); k++)
    {
        core::CSyncPoint ssPoint;
        ssPoint.m_nStartPosition = k > 0 ? a_rSyncPointCoordinates[k-1] : 0;
        ssPoint.m_nEndPosition = a_rSyncPointCoordinates[k];
        ssPoint.m_nIndex = k;
        int bound = ssPoint.m_nStartPosition == ssPoint.m_nEndPosition ? 1 : 0;
        
        while(baseIncludedItr < (int)a_rBaseIncluded.size() && a_rBaseIncluded[baseIncludedItr]->GetStartPos() < (a_rSyncPointCoordinates[k] + bound))
        {
            const core::COrientedVariant* pOvar = a_rBaseIncluded[baseIncludedItr];
            ssPoint.m_baseVariantsIncluded.push_back(pOvar);
            baseIncludedItr++;
        }
        
        while(calledIncludedItr < (int)a_rCalledIncluded.size() && a_rCalledIncluded[calledIncludedItr]->GetStartPos() < (a_rSyncPointCoordinates[k] + bound))
        {
            const core::COrientedVariant* pOvar = a_rCalledIncluded[calledIncludedItr];
            ssPoint.m_calledVariantsIncluded.push_back(pOvar);
            calledIncludedItr++;
        }
        
        while(baseExcludedItr < (int)a_rBaseExcluded.size() && a_rBaseExcluded[baseExcludedItr]->m_nStartPos < (a_rSyncPointCoordinates[k] + bound))
        {
            ssPoint.m_baseVariantsExcluded.push_back(a_rBaseExcluded[baseExcludedItr]);
            baseExcludedItr++;
        }
        
        while(calledExcludedItr < (int)a_rCalledExcluded.size() && a_rCalledExcluded[calledExcludedItr]->m_nStartPos < (a_rSyncPointCoordinates[k] + bound))
        {
            ssPoint.m_calledVariantsExcluded.push_back(a_rCalledExcluded[calledExcludedItr]);
            calledExcludedItr++;
        }
        
        a_rSyncPointList.push_back(ssPoint);
    }
    
    //Add Remaining variants to the last syncPoint
    core::CSyncPoint sPoint;
    sPoint.m_nStartPosition = a_rSyncPointCoordinates[a_rSyncPointCoordinates.size()-1];
    sPoint.m_nEndPosition = INT_MAX;
    sPoint.m_nIndex = static_cast<int>(a_rSyncPointCoordinates.size()-1);
  
    while(baseIncludedItr < (int)a_rBaseIncluded.size() && a_rBaseIncluded[baseIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsIncluded.push_back(a_rBaseIncluded[baseIncludedItr]);
        baseIncludedItr++;
    }
    while(calledIncludedItr < (int)a_rCalledIncluded.size() && a_rCalledIncluded[calledIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsIncluded.push_back(a_rCalledIncluded[calledIncludedItr]);
        calledIncludedItr++;
    }
    
    while(baseExcludedItr < (int)a_rBaseExcluded.size() && a_rBaseExcluded[baseExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsExcluded.push_back(a_rBaseExcluded[baseExcludedItr]);
        baseExcludedItr++;
    }
    
    while(calledExcludedItr < (int)a_rCalledExcluded.size() && a_rCalledExcluded[calledExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsExcluded.push_back(a_rCalledExcluded[calledExcludedItr]);
        calledExcludedItr++;
    }
    a_rSyncPointList.push_back(sPoint);

}

void CGraphVcfAnalyzer::FilterGuaranteedUniqueVariants(std::vector<core::CSyncPoint>& a_rSyncPointList,
                                                       std::vector<const CVariant*>& a_rBaseExcluded,
                                                       std::vector<const CVariant*>& a_rCalledExcluded)
{
    for(core::CSyncPoint syncPoint : a_rSyncPointList)
    {
        if(syncPoint.m_baseVariantsIncluded.size() == 0 && syncPoint.m_calledVariantsIncluded.size() == 0)
            continue;
        else
        {
            for(const CVariant* var : syncPoint.m_baseVariantsExcluded)
                a_rBaseExcluded.push_back(var);
            for(const CVariant* var : syncPoint.m_calledVariantsExcluded)
                a_rCalledExcluded.push_back(var);
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
        unsigned long includedSizeBase = m_uniqueVariantsListPerChromosome[chrIdTuples[k].m_chrName].baseIncludedIndexes.size();
        unsigned long excludedSizeCalled = m_uniqueVariantsListPerChromosome[chrIdTuples[k].m_chrName].calledExcludedIndexes.size();
        unsigned long includedSizeCalled = m_uniqueVariantsListPerChromosome[chrIdTuples[k].m_chrName].calledIncludedIndexes.size();
        
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
    std::cout << std::endl;
    std::cout << " --- GRAPH COMPARISON PARAMETERS --- " << std::endl;
    std::cout << "-base <baseline_vcf_path>    [Required.Add baseline VCF file.]" << std::endl;
    std::cout << "-called <called_vcf_path>    [Required.Add called VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-outDir <output_directory>   [Required.Add output directory]" << std::endl;
    std::cout << "-bed <bed_file_path>         [Optional.Process only given regions in BED file]" << std::endl;
    std::cout << "--pass-filter                [Optional.Process only 'PASS' or '.' variants. Default value is false]" << std::endl;
    std::cout << "-thread-count                [Optional.Specify the number of threads that program will use. Default value is 2]" << std::endl;
    std::cout << "-max-bp-length               [*Optional.Specify the maximum base pair length of variant to process. Default value is 1000]" << std::endl;
    std::cout << "-max-path-size <size>        [*Optional.Specify the maximum size of path that core algorithm can store inside. Default value is 150,000]" << std::endl;
    std::cout << "-max-iteration-count <count> [*Optional.Specify the maximum iteration count that core algorithm can decide to include/exclude variant. Default value is 10,000,000]" << std::endl;
    std::cout << "(*) - for advanced usage" << std::endl;
    std::cout << std::endl;
}





