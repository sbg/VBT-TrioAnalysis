//  CVcfAnalyzer.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 12/8/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVcfAnalyzer.h"
#include "iostream"
#include <math.h>
#include <algorithm>
#include "CGa4ghOutputProvider.h"
#include "CSplitOutputProvider.h"
#include <fstream>

using namespace duocomparison;

void CVcfAnalyzer::Run(int argc, char** argv)
{
    
    std::time_t start;
    double duration;
    
    //Read command line parameters to m_config object
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return;

    std::cout << "MaxPath: " << m_config.m_nMaxPathSize << std::endl;
    std::cout << "MaxIteration: " << m_config.m_nMaxIterationCount << std::endl;

    start = std::time(0);
    
    //Initialize Variant providers which contains VCF and FASTA files
    isSuccess = m_provider.InitializeReaders(m_config);

    if(!isSuccess)
        return;
    
    duration = std::difftime(std::time(0) ,start);
    std::cout << "Vcf and fasta Parser read completed in " << duration << " secs" << std::endl;
    
    std::time_t start1 = std::time(0);
    
    //Creates the threads according to given memory and process the data
    AssignJobsToThreads(m_config.m_nThreadCount);
    
    
    if(0 == strcmp(m_config.m_pOutputMode, "SPLIT"))
    {
        std::cout << "Generating Outputs [SPLIT MODE]..." << std::endl;
        CSplitOutputProvider outputprovider;
        outputprovider.SetVcfPath(m_config.m_pOutputDirectory);
        outputprovider.SetVariantProvider(&m_provider);
        outputprovider.SetBestPaths(m_aBestPaths);
        outputprovider.SetContigList(m_provider.GetContigs());
        outputprovider.GenerateSplitVcfs(m_provider.GetChromosomeIdTuples());
    }
    
    else
    {
        std::cout << "Generating Outputs [GA4GH MODE]..." << std::endl;
        CGa4ghOutputProvider outputprovider;
        outputprovider.SetVcfPath(m_config.m_pOutputDirectory);
        outputprovider.SetVariantProvider(&m_provider);
        outputprovider.SetBestPaths(m_aBestPaths, m_aBestPathsAllele);
        outputprovider.SetContigList(m_provider.GetContigs());
        outputprovider.GenerateGa4ghVcf(m_provider.GetChromosomeIdTuples());
    }
    
    if(true == m_config.m_bGenerateSyncPoints)
    {
        std::vector<SChrIdTuple> chromosomeListToProcess = m_provider.GetChromosomeIdTuples();
        m_resultLogger.OpenSyncPointFile(std::string(m_config.m_pOutputDirectory) + "/SyncPointList.txt");
        for(unsigned int k = 0; k < chromosomeListToProcess.size(); k++)
        {
            std::vector<core::CSyncPoint> syncPointList;
            CalculateSyncPointList(chromosomeListToProcess[k], syncPointList);
            m_resultLogger.WriteSyncPointList(chromosomeListToProcess[k].m_chrName, syncPointList);
        }
        m_resultLogger.CloseSyncPointFile();
    }

    m_resultLogger.SetLogPath(m_config.m_pOutputDirectory);
    int logMode = (0 == strcmp(m_config.m_pOutputMode, "SPLIT") ? 0 : 2) + (m_config.m_bIsGenotypeMatch ? 0 : 1);
    m_resultLogger.WriteStatistics(logMode);
    
    duration = std::difftime(std::time(0), start1);
    std::cout << "Processing Chromosomes completed in " << duration << " secs" << std::endl;
    duration = std::difftime(std::time(0), start);
    std::cout << "Total execution time is " << duration << " secs" << std::endl;
}

int CVcfAnalyzer::AssignJobsToThreads(int a_nThreadCount)
{
    //Get the list of chromosomes to be processed
    std::vector<SChrIdTuple> chromosomeListToProcess = m_provider.GetChromosomeIdTuples();
    
    //Initialize best Path vectors
    m_aBestPaths = std::vector<core::CPath>(chromosomeListToProcess.size());
    m_aBestPathsAllele = std::vector<core::CPath>(chromosomeListToProcess.size());
    
    int exactThreadCount = std::min(a_nThreadCount, (int)chromosomeListToProcess.size());
        
    //Allocate threads
    m_pThreadPool = new std::thread[exactThreadCount];
    
    int threadPoolIt = 0;
    std::vector<SChrIdTuple> *chromosomeLists = new std::vector<SChrIdTuple>[exactThreadCount];
    
    //Divide tasks into threads
    for(unsigned int k = 0; k < chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(chromosomeListToProcess[k]);
        threadPoolIt = (threadPoolIt+1) % exactThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < exactThreadCount; k++)
    {
        if(0 == strcmp("SPLIT", m_config.m_pOutputMode))
            m_pThreadPool[k] = std::thread(&CVcfAnalyzer::ThreadFunctionSPLIT, this, chromosomeLists[k], m_config.m_bIsGenotypeMatch);
        else
            m_pThreadPool[k] = std::thread(&CVcfAnalyzer::ThreadFunctionGA4GH, this, chromosomeLists[k]);
    }
    
    for(int k = 0; k < exactThreadCount; k++)
        m_pThreadPool[k].join();

    //Clear allocated memory
    delete[] chromosomeLists;
    delete[] m_pThreadPool;

    return exactThreadCount;
}


void CVcfAnalyzer::ThreadFunctionGA4GH(std::vector<SChrIdTuple> a_aTuples)
{
    for(unsigned int k = 0; k < a_aTuples.size(); k++)
    {
        std::vector<const CVariant*> varListBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId);
        std::vector<const CVariant*> varListCalled = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId);
        std::vector<const core::COrientedVariant*> ovarListBase = m_provider.GetOrientedVariantList(eBASE, a_aTuples[k].m_nBaseId, true);
        std::vector<const core::COrientedVariant*> ovarListCalled = m_provider.GetOrientedVariantList(eCALLED, a_aTuples[k].m_nCalledId, true);
        
        core::CPathReplay pathReplay(varListBase, varListCalled, ovarListBase, ovarListCalled);
        pathReplay.SetMaxPathAndIteration(m_config.m_nMaxPathSize, m_config.m_nMaxIterationCount);
        SContig ctg;
        m_provider.GetContig(a_aTuples[k].m_chrName, ctg);
        
        //Find Best Path [GENOTYPE MATCH]
        m_aBestPaths[a_aTuples[k].m_nTupleIndex] = pathReplay.FindBestPath(ctg,true);
        
        //Genotype Match variants
        const std::vector<const core::COrientedVariant*>& includedVarsBase = m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsCall = m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetIncludedVariants();
        
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId, m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsCall = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId, m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetExcluded());
        
        //Fill oriented variants for allele match
        m_provider.FillAlleleMatchVariantList(a_aTuples[k], excludedVarsBase, excludedVarsCall);
        
        //Clear old variant pointers
        ovarListBase.clear();
        ovarListCalled.clear();
        varListBase.clear();
        varListCalled.clear();
        
        //Set new variant pointers for allele match comparison
        varListBase = excludedVarsBase;
        varListCalled = excludedVarsCall;
        ovarListBase = m_provider.GetOrientedVariantList(eBASE, a_aTuples[k].m_nBaseId, false);
        ovarListCalled = m_provider.GetOrientedVariantList(eCALLED, a_aTuples[k].m_nCalledId, false);
        pathReplay.Clear();
        
        //Find Best Path [ALLELE MATCH]
        m_aBestPathsAllele[a_aTuples[k].m_nTupleIndex] = pathReplay.FindBestPath(ctg, false);
        
        //No Match variants
        std::vector<const CVariant*> excludedVarsBase2 = m_provider.GetVariantList(excludedVarsBase,
                                                                                   m_aBestPathsAllele[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsCall2 = m_provider.GetVariantList(excludedVarsCall,
                                                                                   m_aBestPathsAllele[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetExcluded());
        
        //Allele Match variants
        const std::vector<const core::COrientedVariant*>& includedVarsBase2 = m_aBestPathsAllele[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsCall2 = m_aBestPathsAllele[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetIncludedVariants();
        
        mtx.lock();
        m_resultLogger.LogStatistic(a_aTuples[k].m_chrName,
                                    a_aTuples[k].m_nBaseId,
                                    static_cast<int>(includedVarsCall.size()),
                                    static_cast<int>(includedVarsBase.size()),
                                    static_cast<int>(includedVarsCall2.size()),
                                    static_cast<int>(includedVarsBase2.size()),
                                    static_cast<int>(excludedVarsCall2.size()),
                                    static_cast<int>(excludedVarsBase2.size()));
        mtx.unlock();
        
        m_provider.SetVariantStatus(excludedVarsBase2, eNO_MATCH);
        m_provider.SetVariantStatus(excludedVarsCall2, eNO_MATCH);
        m_provider.SetVariantStatus(includedVarsCall,  eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(includedVarsBase,  eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(includedVarsCall2, eALLELE_MATCH);
        m_provider.SetVariantStatus(includedVarsBase2, eALLELE_MATCH);
        
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("FP_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , excludedVarsCall);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("TP_BASE_") + std::to_string(a_nChrArr[k] +1) + std::string(".txt")  , includedVarsBase);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("TP_CALLED_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , includedVarsCall);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("FN_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , excludedVarsBase);
    }
}

void CVcfAnalyzer::ThreadFunctionSPLIT(std::vector<SChrIdTuple> a_aTuples, bool a_bIsGenotypeMatch)
{
    for(unsigned int k = 0; k < a_aTuples.size(); k++)
    {
        std::vector<const CVariant*> varListBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId);
        std::vector<const CVariant*> varListCalled = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId);
        
        std::vector<const core::COrientedVariant*> ovarListBase;
        std::vector<const core::COrientedVariant*> ovarListCalled;
        
        if(!a_bIsGenotypeMatch)
        {
            //Fill oriented variants for allele match
            m_provider.FillAlleleMatchVariantList(a_aTuples[k], varListBase, varListCalled);
        }
        
        ovarListBase = m_provider.GetOrientedVariantList(eBASE, a_aTuples[k].m_nBaseId, a_bIsGenotypeMatch);
        ovarListCalled = m_provider.GetOrientedVariantList(eCALLED, a_aTuples[k].m_nCalledId, a_bIsGenotypeMatch);
        
        core::CPathReplay pathReplay(varListBase, varListCalled, ovarListBase, ovarListCalled);
        pathReplay.SetMaxPathAndIteration(m_config.m_nMaxPathSize, m_config.m_nMaxIterationCount);
        SContig ctg;
        m_provider.GetContig(a_aTuples[k].m_chrName, ctg);
        
        m_aBestPaths[a_aTuples[k].m_nTupleIndex] = pathReplay.FindBestPath(ctg,a_bIsGenotypeMatch);
        
        //Genotype Match variants
        const std::vector<const core::COrientedVariant*>& includedVarsBase = m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsCall = m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetIncludedVariants();
        
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsBase = m_provider.GetVariantList(eBASE, a_aTuples[k].m_nBaseId, m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsCall = m_provider.GetVariantList(eCALLED, a_aTuples[k].m_nCalledId, m_aBestPaths[a_aTuples[k].m_nTupleIndex].m_calledSemiPath.GetExcluded());
        
        mtx.lock();
        m_resultLogger.LogStatistic(a_aTuples[k].m_chrName,
                                    a_aTuples[k].m_nBaseId,
                                    static_cast<int>(includedVarsCall.size()),
                                    static_cast<int>(includedVarsBase.size()),
                                    static_cast<int>(0),
                                    static_cast<int>(0),
                                    static_cast<int>(excludedVarsCall.size()),
                                    static_cast<int>(excludedVarsBase.size()));
        mtx.unlock();
        
        EVariantMatch match = a_bIsGenotypeMatch ? eGENOTYPE_MATCH : eALLELE_MATCH;
        m_provider.SetVariantStatus(excludedVarsBase, eNO_MATCH);
        m_provider.SetVariantStatus(excludedVarsCall, eNO_MATCH);
        m_provider.SetVariantStatus(includedVarsCall,  match);
        m_provider.SetVariantStatus(includedVarsBase,  match);
        
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("FP_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , excludedVarsCall);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("TP_BASE_") + std::to_string(a_nChrArr[k] +1) + std::string(".txt")  , includedVarsBase);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("TP_CALLED_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , includedVarsCall);
        //PrintVariants(std::string(m_config.m_pOutputDirectory), std::string("FN_") + std::to_string(a_nChrArr[k] + 1) + std::string(".txt")  , excludedVarsBase);
    }
}

void CVcfAnalyzer::CalculateSyncPointList(const SChrIdTuple& a_rTuple, std::vector<core::CSyncPoint>& a_rSyncPointList)
{
    std::vector<const core::COrientedVariant*> pBaseIncluded = m_aBestPaths[a_rTuple.m_nTupleIndex].m_baseSemiPath.GetIncludedVariants();
    std::vector<const core::COrientedVariant*> pCalledIncluded = m_aBestPaths[a_rTuple.m_nTupleIndex].m_calledSemiPath.GetIncludedVariants();;
    
    std::vector<const CVariant*> pBaseExcluded = m_provider.GetVariantList(eBASE, a_rTuple.m_nBaseId, m_aBestPaths[a_rTuple.m_nTupleIndex].m_baseSemiPath.GetExcluded());
    std::vector<const CVariant*> pCalledExcluded = m_provider.GetVariantList(eCALLED, a_rTuple.m_nCalledId, m_aBestPaths[a_rTuple.m_nTupleIndex].m_calledSemiPath.GetExcluded());
    
    core::CPath *pPath = &m_aBestPaths[a_rTuple.m_nTupleIndex];
    
    unsigned int baseIncludedItr = 0;
    unsigned int baseExcludedItr = 0;
    unsigned int calledIncludedItr = 0;
    unsigned int calledExcludedItr = 0;

    for(unsigned int k = 0; k < pPath->m_aSyncPointList.size(); k++)
    {
        core::CSyncPoint ssPoint;
        ssPoint.m_nStartPosition = k > 0 ? pPath->m_aSyncPointList[k-1] : 0;
        ssPoint.m_nEndPosition = pPath->m_aSyncPointList[k];
        ssPoint.m_nIndex = (int)k;
        
        while(baseIncludedItr < pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            const core::COrientedVariant* pOvar = pBaseIncluded[baseIncludedItr];
            ssPoint.m_baseVariantsIncluded.push_back(pOvar);
            baseIncludedItr++;
        }
        
        while(calledIncludedItr < pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            const core::COrientedVariant* pOvar = pCalledIncluded[calledIncludedItr];
            ssPoint.m_calledVariantsIncluded.push_back(pOvar);
            calledIncludedItr++;
        }
        
        while(baseExcludedItr < pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            ssPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
            baseExcludedItr++;
        }
        
        while(calledExcludedItr < pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            ssPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
            calledExcludedItr++;
        }
        
        a_rSyncPointList.push_back(ssPoint);
    }
    
    //Add Remaining variants to the last syncPoint
    core::CSyncPoint sPoint;
    sPoint.m_nStartPosition = pPath->m_aSyncPointList[pPath->m_aSyncPointList.size()-1];
    sPoint.m_nEndPosition = INT_MAX;
    sPoint.m_nIndex = static_cast<int>(pPath->m_aSyncPointList.size()-1);
    
    while(baseIncludedItr < pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsIncluded.push_back(pBaseIncluded[baseIncludedItr]);
        baseIncludedItr++;
    }
    while(calledIncludedItr < pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsIncluded.push_back(pCalledIncluded[calledIncludedItr]);
        calledIncludedItr++;
    }
    
    while(baseExcludedItr < pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
        baseExcludedItr++;
    }
    
    while(calledExcludedItr < pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
        calledExcludedItr++;
    }
    a_rSyncPointList.push_back(sPoint);
    
}


void CVcfAnalyzer::PrintVariants(std::string a_outputDirectory, std::string a_FileName, const std::vector<const core::COrientedVariant*>& a_rOvarList) const
{
    std::ofstream outputFile;
    
    std::string fileName = a_outputDirectory + std::string("/") + a_FileName;
    
    outputFile.open(fileName.c_str());
    
    for(int k = (int)a_rOvarList.size()-1; k >= 0; k--)
        outputFile << a_rOvarList[k]->GetVariant().ToString() << std::endl;
    
    outputFile.close();
}

void CVcfAnalyzer::PrintVariants(std::string a_outputDirectory, std::string a_FileName, const std::vector<const CVariant*>& a_rVarList) const
{
    std::ofstream outputFile;
    
    std::string fileName = a_outputDirectory + std::string("/") + a_FileName;
    
    outputFile.open(fileName.c_str());
    
    for(const CVariant* pVar : a_rVarList)
        outputFile << pVar->ToString() << std::endl;

    outputFile.close();
}


bool CVcfAnalyzer::ReadParameters(int argc, char** argv)
{

    const char* PARAM_BASE = "-base";
    const char* PARAM_CALLED = "-called";
    const char* PARAM_FILTER = "-filter";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_BED = "-bed";
    const char* PARAM_HELP = "--help";
    const char* PARAM_SAMPLE_BASE = "-sample-base";
    const char* PARAM_SAMPLE_CALLED = "-sample-called";
    const char* PARAM_SNP_ONLY = "--SNP_ONLY";
    const char* PARAM_INDEL_ONLY = "--INDEL_ONLY";
    const char* PARAM_OUTPUT_DIR = "-outDir";
    const char* PARAM_REF_OVERLAP = "--ref-overlap";
    const char* PARAM_CLIP_FROM_END = "--trim-endings-first";
    const char* PARAM_PLATFORM = "--platform-mode";
    const char* PARAM_THREAD_COUNT = "-thread-count";
    const char* PARAM_OUTPUT_MODE = "-output-mode";
    const char* PARAM_ALLELE_MATCH = "--allele-match";
    const char* PARAM_GENERATE_SYNC_POINT = "--generate-sync-point";
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
            m_config.m_pBaseVcfFileName = argv[it+1];
            bBaselineSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CALLED))
        {
            m_config.m_pCalledVcfFileName = argv[it+1];
            bCalledSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REFERENCE))
        {
            m_config.m_pFastaFileName = argv[it+1];
            bReferenceSet = true;
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_OUTPUT_DIR))
        {
            m_config.m_pOutputDirectory = argv[it+1];
            bOutputDirSet = true;
            it += 2;
        }

        else if(0 == strcmp(argv[it], PARAM_BED))
        {
            m_config.m_bInitializeFromBed = true;
            m_config.m_pBedFileName = argv[it+1];
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_OUTPUT_MODE))
        {
            if(0 == strcmp("SPLIT", argv[it+1]))
                m_config.m_pOutputMode = "SPLIT";
            else if(0 == strcmp("GA4GH", argv[it+1]))
                m_config.m_pOutputMode = "GA4GH";
            else
            {
                std::cout << "Invalid Output Mode. SPLIT mode is selected instead!" << std::endl;
                m_config.m_pOutputMode = "SPLIT";
            }
            
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REF_OVERLAP))
        {
            m_config.m_bIsRefOverlap = true;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CLIP_FROM_END))
        {
            m_config.m_bTrimBeginningFirst = false;
            it++;
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
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_BASE))
        {
            m_config.m_bBaseSampleEnabled = true;
            m_config.m_pBaseSample = argv[it+1];
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_CALLED))
        {
            m_config.m_bCalledSampleEnabled = true;
            m_config.m_pCalledSample = argv[it+1];
            it += 2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_SNP_ONLY))
        {
            m_config.m_bSNPOnly = true;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_INDEL_ONLY))
        {
            m_config.m_bINDELOnly = true;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_PLATFORM))
        {
            m_config.m_nThreadCount = MAX_THREAD_COUNT;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_THREAD_COUNT))
        {
            m_config.m_nThreadCount = std::min(std::max(1, atoi(argv[it+1])), MAX_THREAD_COUNT);
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_ALLELE_MATCH))
        {
            m_config.m_bIsGenotypeMatch = false;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_GENERATE_SYNC_POINT))
        {
            m_config.m_bGenerateSyncPoints = true;
            it++;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MAX_PATH_SIZE))
        {
            m_config.m_nMaxPathSize = atoi(argv[it+1]);
            it+=2;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MAX_ITERATION_COUNT))
        {
            m_config.m_nMaxIterationCount = atoi(argv[it+1]);
            it+=2;
        }

        else if(0 == strcmp(argv[it], PARAM_MAX_BP_LENGTH))
        {
            m_config.m_nMaxVariantSize = atoi(argv[it+1]);
            it+=2;
        }
        
        else
            it++; //break;
    }
    
    if(!bBaselineSet)
        std::cout << "Baseline vcf file is not set" << std::endl;
    else if(!bCalledSet)
        std::cout << "Called vcf file is not set" << std::endl;
    else if(!bReferenceSet)
        std::cout << "Reference fasta file is not set" << std::endl;
    else if(!bOutputDirSet)
        std::cout << "Output Directory is not set" << std::endl;
    
    
    return bCalledSet & bBaselineSet & bReferenceSet & bOutputDirSet;
}

void CVcfAnalyzer::PrintHelp() const
{
    std::cout << std::endl;
    std::cout << " --- VARCOMP PARAMETERS --- " << std::endl;
    std::cout << "-base <baseline_vcf_path>    [Required.Add baseline VCF file.]" << std::endl;
    std::cout << "-called <called_vcf_path>    [Required.Add called VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-outDir <output_directory>   [Required.Add output directory]" << std::endl;
    std::cout << "-bed <bed_file_path>         [Optional.Filter variants out of comparison for the given regions]" << std::endl;
    std::cout << "-output-mode <output_mode>   [Optional.Choose the output mode. SPLIT creates 4 vcf files. GA4GH creates a single merged vcf. Default value is SPLIT]" << std::endl;
    std::cout << "-filter <filter_name>        [Optional.Filter variants based on filter column. Default value is PASS. Use 'none' to unfilter]" << std::endl;
    std::cout << "--allele-match               [Optional.Execute the variant comparison engine in allele matching mode]" << std::endl;
    std::cout << "--SNP_ONLY                   [Optional.Filter INDELs out from both base and called VCF file.]" << std::endl;
    std::cout << "--INDEL_ONLY                 [Optional.Filter SNPs out from both base and called VCF file.]" << std::endl;
    std::cout << "-sample-base <sample_name>   [Optional.Read only the given sample in base VCF. Default value is the first sample.]" << std::endl;
    std::cout << "-sample-called <sample_name> [Optional.Read only the given sample in called VCF. Default value is the first sample.]" << std::endl;
    std::cout << "--ref-overlap                [Optional.Allow reference overlapping by trimming nucleotides and ignoring 0 genotype.]" << std::endl;
    std::cout << "--generate-sync-point        [Optional.Prints the sync point list of two vcf file. Default value is false.]" << std::endl;
    std::cout << "--trim-endings-first         [Optional.If set, starts trimming variants from ending base pairs. Default is from beginning]" << std::endl;
    std::cout << "--platform-mode              [Optional.Allow to run program with the thread number of different chromosome count.]" << std::endl;
    std::cout << "-thread-count                [Optional.Specify the number of threads that program will use. Default value is 2]" << std::endl;
    std::cout << "-max-bp-length               [*Optional.Specify the maximum base pair length of variant to process. Default value is 1000]" << std::endl;
    std::cout << "-max-path-size <size>        [*Optional.Specify the maximum size of path that core algorithm can store inside. Default value is 150,000]" << std::endl;
    std::cout << "-max-iteration-count <count> [*Optional.Specify the maximum iteration count that core algorithm can decide to include/exclude variant. Default value is 10,000,000]" << std::endl;
    std::cout << "(*) - advanced usage" << std::endl;
    std::cout << std::endl;
    std::cout << "Example Commands:" << std::endl;
    std::cout << "./vbt varcomp -called called.vcf -base base.vcf -ref reference.fa -outDir SampleResultDir -filter none" << std::endl;
    std::cout << "./vbt varcomp -called duo.vcf -base duo.vcf -ref reference.fa -outDir SampleResultDir -filter PASS -SNP_ONLY -sample-base sample00 -sample-called sample01" << std::endl;
}
