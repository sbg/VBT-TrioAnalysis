//
//  CMendelianViolationAnalyzer.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 1/31/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianAnalyzer.h"
#include "SConfig.h"
#include <string>
#include "CPathReplay.h"
#include <iostream>


void CMendelianAnalyzer::run(int argc, char **argv)
{
    std::clock_t start;
    double duration;
    
    //Reads the command line parameters
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return;
    
    start = std::clock();
    
    //Initialize variant provider
    isSuccess = m_provider.InitializeReaders(m_fatherChildConfig, m_motherChildConfig);
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Vcf and fasta Parser read completed in " << duration << " secs" << std::endl;
    
    if(!isSuccess)
        return;
    
    //Decide Thread count and start parent child comparison
    if(m_fatherChildConfig.m_bIsPlatformMode)
    {
        AssignJobsToThreads(CHROMOSOME_COUNT);
        for(int k = 0; k < CHROMOSOME_COUNT; k++)
            m_pThreadPool[k].join();
    }
    else
    {
        AssignJobsToThreads(MAC_THREAD_COUNT);
        for(int k = 0; k < MAC_THREAD_COUNT; k++)
            m_pThreadPool[k].join();
    }
    
    std::vector<int> chrIds = m_provider.GetCommonChromosomes();
    for(int k = 0; k < chrIds.size(); k++)
        MergeFunc(chrIds[k]);
    
    //Write results to log file
    m_resultLog.SetLogPath(m_fatherChildConfig.m_pOutputDirectory);
    m_resultLog.WriteMendelianStatistics();
    
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Program Completed in " << duration << " secs" << std::endl;
}

void CMendelianAnalyzer::ThreadFunc(const std::vector<SMendelianThreadParam>& a_rThreadParams)
{
    for(int k = 0; k < a_rThreadParams.size(); k++)
    {
        int a_nChromosomeId = a_rThreadParams[k].m_nChromosomeIndex;
        bool a_bIsFatherChild = a_rThreadParams[k].m_bIsFatherChildComparison;
        
        //Perform operations on mother or father according to selection
        EMendelianVcfName baseName = a_bIsFatherChild ? eFATHER : eMOTHER;
        
        std::vector<const CVariant*> varListBase = m_provider.GetVariantList(baseName, a_nChromosomeId);
        std::vector<const CVariant*> varListCalled = m_provider.GetVariantList(eCHILD, a_nChromosomeId);

        std::vector<const COrientedVariant*> ovarListBase = m_provider.GetOrientedVariantList(baseName, a_nChromosomeId);
        std::vector<const COrientedVariant*> ovarListCalled = m_provider.GetOrientedVariantList(eCHILD, a_nChromosomeId);

        //Values to Log
        int TPBase, TPCalled, Fp, Fn;
        
        //Get the chromosome ref seq
        SContig ctg;
        m_provider.GetContig(a_nChromosomeId, ctg);
        
        //Create path replay for father child;
        CPathReplay replay(varListBase, varListCalled, ovarListBase, ovarListCalled);
        
        //Find Best Path Parent Child
        if(true == a_bIsFatherChild)
        {
            m_aBestPathsFatherChild[a_nChromosomeId] = replay.FindBestPath(ctg, false);
            TPBase = static_cast<int>(m_aBestPathsFatherChild[a_nChromosomeId].m_baseSemiPath.GetIncludedVariants().size());
            TPCalled = static_cast<int>(m_aBestPathsFatherChild[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants().size());
            Fp = static_cast<int>(m_aBestPathsFatherChild[a_nChromosomeId].m_calledSemiPath.GetExcluded().size());
            Fn = static_cast<int>(m_aBestPathsFatherChild[a_nChromosomeId].m_baseSemiPath.GetExcluded().size());
        }
        else
        {
            m_aBestPathsMotherChild[a_nChromosomeId] = replay.FindBestPath(ctg, false);
            TPBase = static_cast<int>(m_aBestPathsMotherChild[a_nChromosomeId].m_baseSemiPath.GetIncludedVariants().size());
            TPCalled = static_cast<int>(m_aBestPathsMotherChild[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants().size());
            Fp = static_cast<int>(m_aBestPathsMotherChild[a_nChromosomeId].m_calledSemiPath.GetExcluded().size());
            Fn = static_cast<int>(m_aBestPathsMotherChild[a_nChromosomeId].m_baseSemiPath.GetExcluded().size());
        }
        
        //Log statistics of current process
        m_resultLog.LogMendelianStatistic(a_bIsFatherChild, a_nChromosomeId, TPCalled, TPBase, Fp, Fn);
        
//        std::cout << "====CHROMOSOME " << a_nChromosomeId << (a_bIsFatherChild ? "  (father-child)" : "   (mother-child)") << "====" << std::endl;
//        std::cout << "TP Base:" << TPBase << std::endl;
//        std::cout << "TP Called:" << TPCalled << std::endl;
//        std::cout << "FP:" << Fp << std::endl;
//        std::cout << "FN:" << Fn << std::endl;
        
    }
}

bool CMendelianAnalyzer::ReadParameters(int argc, char **argv)
{
    const char* PARAM_HELP = "-help";
    
    const char* PARAM_FATHER = "-father";
    const char* PARAM_MOTHER = "-mother";
    const char* PARAM_CHILD = "-child";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_FILTER = "-filter";
    
    const char* PARAM_SAMPLE_FATHER = "-SampleFather";
    const char* PARAM_SAMPLE_MOTHER = "-SampleMother";
    const char* PARAM_SAMPLE_CHILD = "-SampleChild";
    
    const char* PARAM_OUTPUT_DIR = "-outDir";
    const char* PARAM_REF_OVERLAP = "-ref-overlap";
    const char* PARAM_PLATFORM = "-platform-mode";
    
    bool bFatherSet = false;
    bool bMotherSet = false;
    bool bChildSet = false;
    bool bReferenceSet = false;
    bool bOutputDirSet = false;
    
    //Start from index 2 since first parameter will be mendelian mode indicator
    int it = 2;
    
    while(it < argc)
    {
        if(0 == strcmp(argv[it], PARAM_HELP))
        {
            //PrintHelp();
            return false;
        }
        
        else if(0 == strcmp(argv[it], PARAM_FATHER))
        {
            m_fatherChildConfig.m_pBaseVcfFileName = argv[it+1];
            bFatherSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_MOTHER))
        {
            m_motherChildConfig.m_pBaseVcfFileName = argv[it+1];
            bMotherSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CHILD))
        {
            m_motherChildConfig.m_pCalledVcfFileName = argv[it+1];
            m_fatherChildConfig.m_pCalledVcfFileName = argv[it+1];
            bChildSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REFERENCE))
        {
            m_motherChildConfig.m_pFastaFileName = argv[it+1];
            m_fatherChildConfig.m_pFastaFileName = argv[it+1];
            bReferenceSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_OUTPUT_DIR))
        {
            m_motherChildConfig.m_pOutputDirectory = argv[it+1];
            m_fatherChildConfig.m_pOutputDirectory = argv[it+1];
            bOutputDirSet = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_REF_OVERLAP))
        {
            m_motherChildConfig.m_bIsRefOverlap = true;
            m_fatherChildConfig.m_bIsRefOverlap = true;
            it--;
        }
        
        else if(0 == strcmp(argv[it], PARAM_FILTER))
        {
            if(0 == strcmp("none", argv[it+1]))
            {
                m_motherChildConfig.m_bIsFilterEnabled = false;
                m_fatherChildConfig.m_bIsFilterEnabled = false;
            }
            else
            {
                m_motherChildConfig.m_bIsFilterEnabled = true;
                m_fatherChildConfig.m_bIsFilterEnabled = true;
                
                m_motherChildConfig.m_pFilterName = argv[it+1];
                m_fatherChildConfig.m_pFilterName = argv[it+1];
            }
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_FATHER))
        {
            m_fatherChildConfig.m_bBaseSampleEnabled = true;
            m_fatherChildConfig.m_pBaseSample = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_MOTHER))
        {
            m_motherChildConfig.m_bBaseSampleEnabled = true;
            m_motherChildConfig.m_pBaseSample = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_SAMPLE_CHILD))
        {
            m_motherChildConfig.m_bCalledSampleEnabled = true;
            m_motherChildConfig.m_pCalledSample = argv[it+1];
            
            m_fatherChildConfig.m_bCalledSampleEnabled = true;
            m_fatherChildConfig.m_pCalledSample = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_PLATFORM))
        {
            m_motherChildConfig.m_bIsPlatformMode = true;
            m_fatherChildConfig.m_bIsPlatformMode = true;
            it--;
        }
        
        it += 2;
    }
    
    if(!bChildSet)
        std::cout << "Child vcf file is not set" << std::endl;
    else if(!bFatherSet)
        std::cout << "Father vcf file is not set" << std::endl;
    else if(!bMotherSet)
        std::cout << "Mother vcf file is not set" << std::endl;
    else if(!bReferenceSet)
        std::cout << "Reference fasta file is not set" << std::endl;
    else if(!bOutputDirSet)
        std::cout << "Output Directory is not set" << std::endl;
    
    
    return bFatherSet && bMotherSet && bChildSet && bReferenceSet && bOutputDirSet;
    
}

void CMendelianAnalyzer::AssignJobsToThreads(int a_nThreadCount)
{
    //Get the list of chromosomes to be processed
    std::vector<int> chromosomeListToProcess = m_provider.GetCommonChromosomes();

    //Allocate threads
    m_pThreadPool = new std::thread[std::min(a_nThreadCount, static_cast<int>(chromosomeListToProcess.size() *2))];
    
    int threadPoolIt = 0;
    std::vector<SMendelianThreadParam> *chromosomeLists = new std::vector<SMendelianThreadParam>[a_nThreadCount];
    
    //Divide tasks into threads
    for(int k = 0; k < chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(SMendelianThreadParam(chromosomeListToProcess[k], true));
        threadPoolIt = (threadPoolIt+1) % a_nThreadCount;
        chromosomeLists[threadPoolIt].push_back(SMendelianThreadParam(chromosomeListToProcess[k], false));
        threadPoolIt = (threadPoolIt+1) % a_nThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < a_nThreadCount; k++)
    {
        m_pThreadPool[k] = std::thread(&CMendelianAnalyzer::ThreadFunc, this, chromosomeLists[k]);
    }
    
}

void CMendelianAnalyzer::MergeFunc(int a_nChromosomeId)
{
    //Included lists of child
    const std::vector<const COrientedVariant*>& FCchildIncluded = m_aBestPathsFatherChild[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants();
    const std::vector<const COrientedVariant*>& MCchildIncluded = m_aBestPathsMotherChild[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants();
    
    //Excluded list of child (Intersection of the 2 is the unique child set)
    const std::vector<int>& FCexcludedIndexes = m_aBestPathsFatherChild[a_nChromosomeId].m_calledSemiPath.GetExcluded();
    const std::vector<int>& MCexcludedIndexes = m_aBestPathsMotherChild[a_nChromosomeId].m_calledSemiPath.GetExcluded();
    std::vector<int> childUniqueVarIndexList;

    for(int FC = 0, MC = 0; FC < FCexcludedIndexes.size() && MC < MCexcludedIndexes.size();)
    {
        if(FCexcludedIndexes[FC] == MCexcludedIndexes[MC])
        {
            childUniqueVarIndexList.push_back(FCexcludedIndexes[FC]);
            FC++;
            MC++;
        }
        else if(FCexcludedIndexes[FC] > MCexcludedIndexes[MC])
            MC++;
        else
            FC++;
    }
    
    //Unique Excluded variant list of child
    //const std::vector<const COrientedVariant*>& childOnlyExcluded = m_provider.GetOrientedVariantList(eCHILD, a_nChromosomeId, common);
    
    
    std::vector<const COrientedVariant*> childCommonHeterozygous;
    std::vector<const COrientedVariant*> childFilteredHeterozygous;
    std::vector<const COrientedVariant*> childCommonHomozygous;
    std::vector<const COrientedVariant*> motherChildHeterozygous;
    std::vector<const COrientedVariant*> motherChildHomozygous;
    std::vector<const COrientedVariant*> fatherChildHeterozygous;
    std::vector<const COrientedVariant*> fatherChildHomozygous;

    std::vector<const COrientedVariant*> MendelianViolationVars;
    std::vector<const COrientedVariant*> MendelianCompliantVars;
    
    //Assure that given chromosome is processed already;
    assert(FCchildIncluded.size() > 0);
    assert(MCchildIncluded.size() > 0);
    
    int uniqueHomozygousVariantInFC = 0;
    int uniqueHomozygousVariantInMC = 0;

    int uniqueHeterozygousVariantInFC = 0;
    int uniqueHeterozygousVariantInMC = 0;

    
    //Get the Intersection of the two list
    int mcIndex, fcIndex;
    for(mcIndex = 0, fcIndex = 0; mcIndex < MCchildIncluded.size() && fcIndex < FCchildIncluded.size();)
    {
        if(MCchildIncluded[mcIndex]->GetVariant().m_nStartPos == FCchildIncluded[fcIndex]->GetVariant().m_nStartPos)
        {
            if(MCchildIncluded[mcIndex]->GetVariant().m_bIsHeterozygous)
            {
                if(MCchildIncluded[mcIndex]->GetAlleleIndex() == FCchildIncluded[fcIndex]->GetAlleleIndex())
                    childFilteredHeterozygous.push_back(MCchildIncluded[mcIndex]);
                else
                    childCommonHeterozygous.push_back(MCchildIncluded[mcIndex]);
            }
            else
                childCommonHomozygous.push_back(MCchildIncluded[mcIndex]);
            
            mcIndex++;
            fcIndex++;
        }
        
        else if(MCchildIncluded[mcIndex]->GetVariant().m_nStartPos > FCchildIncluded[fcIndex]->GetVariant().m_nStartPos)
        {
            if(FCchildIncluded[fcIndex]->GetVariant().m_bIsHeterozygous)
                uniqueHeterozygousVariantInFC++;
            else
                uniqueHomozygousVariantInFC++;
            
            fcIndex++;
        }
        
        else if(MCchildIncluded[mcIndex]->GetVariant().m_nStartPos < FCchildIncluded[fcIndex]->GetVariant().m_nStartPos)
        {
            if(MCchildIncluded[mcIndex]->GetVariant().m_bIsHeterozygous)
                uniqueHeterozygousVariantInMC++;
            else
                uniqueHomozygousVariantInMC++;
            
            mcIndex++;
        }
        
        else
        {
            if(MCchildIncluded[mcIndex]->GetVariant().m_bIsHeterozygous)
                uniqueHeterozygousVariantInMC++;
            else
                uniqueHomozygousVariantInMC++;
            
            mcIndex++;
        }
        
    }
    
    //Add the remaining MC child variants if exists
    if(mcIndex < MCchildIncluded.size())
    {
        for(int k = mcIndex; k < MCchildIncluded.size(); k++)
        {
            if(MCchildIncluded[k]->GetVariant().m_bIsHeterozygous)
                uniqueHeterozygousVariantInMC++;
            else
                uniqueHomozygousVariantInMC++;
        }
    }

    //Add the remaining FC child variatns if exists
    if(fcIndex < FCchildIncluded.size())
    {
        for(int k = fcIndex; k < FCchildIncluded.size(); k++)
        {
            if(FCchildIncluded[k]->GetVariant().m_bIsHeterozygous)
                uniqueHeterozygousVariantInFC++;
            else
                uniqueHomozygousVariantInFC++;
        }
    }
    
    //Send result to the log file
    m_resultLog.LogMendelianIntersectionStatistic(a_nChromosomeId,
                                                  static_cast<int>(childCommonHomozygous.size()),
                                                  static_cast<int>(childCommonHeterozygous.size()),
                                                  static_cast<int>(childFilteredHeterozygous.size()),
                                                  uniqueHomozygousVariantInFC,
                                                  uniqueHeterozygousVariantInFC,
                                                  uniqueHomozygousVariantInMC,
                                                  uniqueHeterozygousVariantInMC,
                                                  static_cast<int>(childUniqueVarIndexList.size()));
}















