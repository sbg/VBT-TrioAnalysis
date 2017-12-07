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
#include "CVariantIterator.h"
#include "CSyncPoint.h"
#include "Utils/CUtils.h"
#include <algorithm>

using namespace mendelian;

CMendelianAnalyzer::CMendelianAnalyzer() :
m_mendelianDecider(m_aBestPathsFatherChildGT, m_aBestPathsFatherChildAM, m_aBestPathsMotherChildGT, m_aBestPathsMotherChildAM, m_provider, m_resultLog)
{
    m_noCallMode = ENoCallMode::eExplicitNoCall;
}

int CMendelianAnalyzer::run(int argc, char **argv)
{
    std::time_t start, start1;
    double duration;
    
    //Reads the command line parameters
    bool isSuccess = ReadParameters(argc, argv);
    
    //Parameter Read fails, return
    if(!isSuccess)
        return -1;
    
    start = std::time(0);
    
    //Initialize variant provider
    isSuccess = m_provider.InitializeReaders(m_fatherChildConfig, m_motherChildConfig);
    
    //Variant Initialization fails, return
    if(!isSuccess)
        return -1;

    duration = std::difftime(std::time(0), start);
    std::cerr << "[stderr] Vcf and fasta Parser read completed in " << duration << " secs" << std::endl;
    start1 = std::time(0);
    
    std::cerr << "[stderr] initializing output writer" << std::endl;
    
    //Initialize output writer
    std::string directory = std::string(m_fatherChildConfig.m_pOutputDirectory);
    std::string trioPath = directory + (directory[directory.length()-1] != '/' ? "/" + std::string(m_fatherChildConfig.m_output_prefix) + "_trio.vcf" : std::string(m_fatherChildConfig.m_output_prefix) + "_trio.vcf");
    
    m_trioWriter.SetTrioPath(trioPath);
    m_trioWriter.SetNoCallMode(m_noCallMode);
    m_trioWriter.SetResultLogPointer(&m_resultLog);
    m_trioWriter.SetContigList(m_provider.GetContigs(),
                               static_cast<int>(m_provider.GetCommonChromosomes().size()),
                               m_provider.GetContigCount(eCHILD),
                               m_provider.GetContigCount(eFATHER),
                               m_provider.GetContigCount(eMOTHER));
    
    std::cerr << "[stderr] Running best path algorithm pipeline for each chromosome..." << std::endl;
    
    //Run core comparison engine on parallel
    AssignJobsToThreads(m_fatherChildConfig.m_nThreadCount);
    
    std::cerr << "[stderr] Evaluating mendelian consistency of variants..." << std::endl;
    
    //Perform merge process
    m_mendelianDecider.SetNocallMode(m_noCallMode);
    std::vector<SChrIdTriplet> chrIds = m_provider.GetCommonChromosomes();
    for(unsigned int k = 0; k < chrIds.size(); k++)
    {
        //Initialize the decision arrays
        std::vector<EMendelianDecision> childDecisions  = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eCHILD,  chrIds[k].m_nCid));
        std::vector<EMendelianDecision> motherDecisions = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eMOTHER, chrIds[k].m_nMid));
        std::vector<EMendelianDecision> fatherDecisions = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eFATHER, chrIds[k].m_nFid));
        
        //Set all decisions to unknown at the beginning
        for(unsigned int m = 0; m < childDecisions.size(); m++)
            childDecisions[m] = eUnknown;
        for(unsigned int m = 0; m < motherDecisions.size(); m++)
            motherDecisions[m] = eUnknown;
        for(unsigned int m = 0; m < fatherDecisions.size(); m++)
            fatherDecisions[m] = eUnknown;
        
        //Merge the chromosome and fill the decisions arrays
        m_mendelianDecider.MergeFunc(chrIds[k], motherDecisions, fatherDecisions, childDecisions);
        
        //Set decision arrays and variants to the output Trio Merger
        m_trioWriter.SetDecisionsAndVariants(chrIds[k], eCHILD,  childDecisions, m_provider.GetSortedVariantListByIDandStartPos(eCHILD, chrIds[k].m_nCid));
        m_trioWriter.SetDecisionsAndVariants(chrIds[k], eMOTHER, motherDecisions, m_provider.GetSortedVariantListByIDandStartPos(eMOTHER, chrIds[k].m_nMid));
        m_trioWriter.SetDecisionsAndVariants(chrIds[k], eFATHER, fatherDecisions, m_provider.GetSortedVariantListByIDandStartPos(eFATHER, chrIds[k].m_nFid));
    }
    
    std::cerr << "[stderr] Generating the output trio vcf..." << std::endl;
    //Generate trio output vcf from common chromosomes
    m_trioWriter.GenerateTrioVcf(chrIds);
    
    std::cerr << "[stderr] Generating detailed output logs.." << std::endl;
    
    m_resultLog.LogSkippedVariantCounts(m_provider.GetSkippedVariantCount(eCHILD),
                                        m_provider.GetSkippedVariantCount(eFATHER),
                                        m_provider.GetSkippedVariantCount(eMOTHER));
    
    m_resultLog.LogFilteredComplexVariantCounts(m_provider.GetNotAssessedVariantCount(eCHILD),
                                                m_provider.GetNotAssessedVariantCount(eFATHER),
                                                m_provider.GetNotAssessedVariantCount(eMOTHER));
    
    //Write results to log file
    m_resultLog.SetLogDirectory(m_fatherChildConfig.m_pOutputDirectory);
    m_resultLog.WriteBestPathStatistics(m_fatherChildConfig.m_output_prefix);
    m_resultLog.WriteDetailedReportTable(m_fatherChildConfig.m_output_prefix);
    m_resultLog.WriteDetailedReportTabDelimited(m_fatherChildConfig.m_output_prefix);
    m_resultLog.WriteShortReportTable(m_fatherChildConfig.m_output_prefix);
    
    duration = std::difftime(std::time(0), start1);
    std::cerr << "[stderr] Processing Chromosomes completed in " << duration << " secs" << std::endl;
    duration = std::difftime(std::time(0), start);
    std::cerr << "[stderr] Total execution time is " << duration << " secs" << std::endl;
    
    return 0;
}

bool CMendelianAnalyzer::ReadParameters(int argc, char **argv)
{
    const char* PARAM_HELP = "--help";
    
    const char* PARAM_FATHER = "-father";
    const char* PARAM_MOTHER = "-mother";
    const char* PARAM_CHILD = "-child";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_FILTER = "-filter";
    const char* PARAM_PEDIGREE = "-pedigree";
    const char* PARAM_BED = "-bed";
    
    const char* PARAM_SAMPLE_FATHER = "-sample-father";
    const char* PARAM_SAMPLE_MOTHER = "-sample-mother";
    const char* PARAM_SAMPLE_CHILD = "-sample-child";
    
    const char* PARAM_OUTPUT_DIR = "-outDir";
    const char* PARAM_REF_OVERLAP = "--ref-overlap";
    const char* PARAM_CLIP_FROM_END = "--trim-endings-first";
    const char* PARAM_THREAD_COUNT = "-thread-count";
    const char* PARAM_NO_CALL = "-no-call";
    
    const char* PARAM_OUTPUT_PREFIX = "-out-prefix";
    
    const char* PARAM_AUTOSOME_ONLY = "--autosome-only";
    
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
            PrintHelp();
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

        else if(0 == strcmp(argv[it], PARAM_PEDIGREE))
        {
            m_motherChildConfig.m_pPedigreeFileName = argv[it+1];
            m_motherChildConfig.m_bInitializeFromPED = true;
            m_fatherChildConfig.m_pPedigreeFileName = argv[it+1];
            m_fatherChildConfig.m_bInitializeFromPED = true;
        }
        
        else if(0 == strcmp(argv[it], PARAM_BED))
        {
            m_motherChildConfig.m_bInitializeFromBed = true;
            m_motherChildConfig.m_pBedFileName = argv[it+1];
            m_fatherChildConfig.m_bInitializeFromBed = true;
            m_fatherChildConfig.m_pBedFileName = argv[it+1];
        }
        
        else if(0 == strcmp(argv[it], PARAM_REF_OVERLAP))
        {
            m_motherChildConfig.m_bIsRefOverlap = true;
            m_fatherChildConfig.m_bIsRefOverlap = true;
            it--;
        }
        
        else if(0 == strcmp(argv[it], PARAM_CLIP_FROM_END))
        {
            m_motherChildConfig.m_bTrimBeginningFirst = false;
            m_fatherChildConfig.m_bTrimBeginningFirst = false;
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
        
        else if(0 == strcmp(argv[it], PARAM_NO_CALL))
        {
            if(0 == strcmp("none", argv[it+1]))
                m_noCallMode = ENoCallMode::eNone;
            else if(0 == strcmp("explicit", argv[it+1]))
                m_noCallMode = ENoCallMode::eExplicitNoCall;
            else if(0 == strcmp("implicit", argv[it+1]))
                m_noCallMode = ENoCallMode::eImplicitNoCall;
            else
                m_noCallMode = ENoCallMode::eNone;
        }
        
        else if(0 == strcmp(argv[it], PARAM_AUTOSOME_ONLY))
        {
            m_motherChildConfig.m_bAutosomeOnly = true;
            m_fatherChildConfig.m_bAutosomeOnly = true;
            it--;
        }
        
        else if(0 == strcmp(argv[it], PARAM_OUTPUT_PREFIX))
        {
            m_motherChildConfig.m_output_prefix = argv[it+1];
            m_fatherChildConfig.m_output_prefix = argv[it+1];
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
        
        else if(0 == strcmp(argv[it], PARAM_THREAD_COUNT))
        {
            m_motherChildConfig.m_nThreadCount = std::min(std::max(1, atoi(argv[it+1])), MAX_THREAD_COUNT);
            m_fatherChildConfig.m_nThreadCount = std::min(std::max(1, atoi(argv[it+1])), MAX_THREAD_COUNT);
        }
        
        else
        {
            std::cerr << "Unknown Command or Argument: " << argv[it] << std::endl;
            std::cerr << "Use the following command for parameter usage:" << std::endl;
            std::cerr << "./vbt mendelian --help" << std::endl;
            return false;
        }

        it += 2;
    }
    
    bool bIsVcfFilesAccesible = true;
    
    if(!bChildSet)
        std::cerr << "Child vcf file is not set" << std::endl;
    else if(!bFatherSet)
        std::cerr << "Father vcf file is not set" << std::endl;
    else if(!bMotherSet)
        std::cerr << "Mother vcf file is not set" << std::endl;
    else if(!bReferenceSet)
        std::cerr << "Reference fasta file is not set" << std::endl;
    else if(!bOutputDirSet)
        std::cerr << "Output Directory is not set" << std::endl;
    else if(!CUtils::IsFileExists(m_motherChildConfig.m_pBaseVcfFileName) || !CUtils::IsFileExists(m_motherChildConfig.m_pCalledVcfFileName) || !CUtils::IsFileExists(m_fatherChildConfig.m_pBaseVcfFileName))
    {
        bIsVcfFilesAccesible = false;
        std::cerr << "One of vcf file paths is wrong" << std::endl;
    }
    
    return bFatherSet && bMotherSet && bChildSet && bReferenceSet && bOutputDirSet && bIsVcfFilesAccesible;
    
}

int CMendelianAnalyzer::AssignJobsToThreads(int a_nThreadCount)
{
    //Thread pool we have for multitasking by per chromosome
    std::thread *pThreadPool;
    
    //Get the list of chromosomes to be processed
    std::vector<SChrIdTriplet> chromosomeListToProcess = m_provider.GetCommonChromosomes();

    //Initialize best path vectors
    m_aBestPathsFatherChildGT = std::vector<core::CPath>(m_provider.GetCommonChromosomes().size());
    m_aBestPathsMotherChildGT = std::vector<core::CPath>(m_provider.GetCommonChromosomes().size());
    m_aBestPathsFatherChildAM = std::vector<core::CPath>(m_provider.GetCommonChromosomes().size());
    m_aBestPathsMotherChildAM = std::vector<core::CPath>(m_provider.GetCommonChromosomes().size());
    
    int exactThreadCount = std::min(a_nThreadCount, (int)chromosomeListToProcess.size());
    
    //Allocate threads
    pThreadPool = new std::thread[exactThreadCount];
    
    int threadPoolIt = 0;
    std::vector<SChrIdTriplet> *chromosomeLists = new std::vector<SChrIdTriplet>[exactThreadCount];
    
    //Divide tasks into threads
    for(unsigned int k = 0; k < chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(chromosomeListToProcess[k]);
        threadPoolIt = (threadPoolIt+1) % exactThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < exactThreadCount; k++)
        pThreadPool[k] = std::thread(&CMendelianAnalyzer::ProcessChromosome, this, chromosomeLists[k]);
    
    for(int k = 0; k < exactThreadCount; k++)
        pThreadPool[k].join();

    //Clean allocation
    delete[] chromosomeLists;
    
    //Clean threads
    delete[] pThreadPool;
    
    return exactThreadCount;
    
}

void CMendelianAnalyzer::ProcessChromosome(const std::vector<SChrIdTriplet>& a_nChromosomeIds)
{    
    for(SChrIdTriplet triplet : a_nChromosomeIds)
    {
        //Lock data read
        mtx.lock();
        //Get variant list of parent-child for given chromosome
        std::vector<const CVariant*> varListFather = m_provider.GetVariantList(eFATHER, triplet.m_nFid);
        std::vector<const CVariant*> varListMother = m_provider.GetVariantList(eMOTHER, triplet.m_nMid);
        std::vector<const CVariant*> varListChild = m_provider.GetVariantList(eCHILD, triplet.m_nCid);
        
        //Get oriented variant list of parent-child for given chromosome
        std::vector<const core::COrientedVariant*> ovarListGTFather = m_provider.GetOrientedVariantList(eFATHER, triplet.m_nFid);
        std::vector<const core::COrientedVariant*> ovarListGTMother = m_provider.GetOrientedVariantList(eMOTHER, triplet.m_nMid);
        std::vector<const core::COrientedVariant*> ovarListGTChild = m_provider.GetOrientedVariantList(eCHILD, triplet.m_nCid);

        //Get the chromosome ref seq
        SContig ctg;
        m_provider.ReadContig(triplet.m_chrName, ctg);
        //Unlock data read
        mtx.unlock();

        // === PROCESS FATHER-CHILD ===
        
        //Create path replay for parent child;
        core::CPathReplay replayFatherChildGT(varListFather, varListChild, ovarListGTFather, ovarListGTChild);
        
        //Find Best Path Father-Child GT Match
        m_aBestPathsFatherChildGT[triplet.m_nTripleIndex] = replayFatherChildGT.FindBestPath(ctg, true);
        
        //Genotype Match variants
        const std::vector<const core::COrientedVariant*>& includedVarsChildGT = m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsFatherGT = m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants();
        
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsFather = m_provider.GetVariantList(eFATHER, triplet.m_nFid, m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsChild = m_provider.GetVariantList(eCHILD, triplet.m_nCid,m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        
        //Allele Match oriented variants
        std::vector<const core::COrientedVariant*> ovarListAMFather = m_provider.GetOrientedVariantList(eFATHER, triplet.m_nFid, true, m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const core::COrientedVariant*> ovarListAMChildFC = m_provider.GetOrientedVariantList(eCHILD, triplet.m_nCid, true, m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        
        //Clear Father child replay object
        replayFatherChildGT.Clear();
        
        //Change the variant list to process
        core::CPathReplay replayFatherChildAM(excludedVarsFather, excludedVarsChild, ovarListAMFather, ovarListAMChildFC);
        
        //Find Best Path Father-Child AM Match
        m_aBestPathsFatherChildAM[triplet.m_nTripleIndex] = replayFatherChildAM.FindBestPath(ctg, false);
        const std::vector<const core::COrientedVariant*>& includedVarsChildAM = m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsFatherAM = m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants();

        const std::vector<const CVariant*> excludedVarsChildFC = m_provider.GetVariantList(excludedVarsChild, m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        const std::vector<const CVariant*> excludedVarsFatherFC = m_provider.GetVariantList(excludedVarsFather, m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
        
        //Set Variant status of child variants
        m_provider.SetVariantStatus(includedVarsChildAM, eALLELE_MATCH);
        m_provider.SetVariantStatus(includedVarsChildGT, eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(excludedVarsChildFC, eNO_MATCH);
        
        //Set Variant status of father variants
        m_provider.SetVariantStatus(includedVarsFatherGT, eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(includedVarsFatherAM, eALLELE_MATCH);
        m_provider.SetVariantStatus(excludedVarsFatherFC, eNO_MATCH);

        //Clear Father child replay object
        replayFatherChildAM.Clear();
        
        // === PROCESS MOTHER-CHILD ===
     
        //Create path replay for parent child;
        core::CPathReplay replayMotherChildGT(varListMother, varListChild, ovarListGTMother, ovarListGTChild);
        
        //Find Best Path Father-Child GT Match
        m_aBestPathsMotherChildGT[triplet.m_nTripleIndex] = replayMotherChildGT.FindBestPath(ctg, true);
        
        //Genotype Match variants
        const std::vector<const core::COrientedVariant*>& includedVarsChildGTMC = m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsMotherGT = m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants();
        
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsMother = m_provider.GetVariantList(eMOTHER, triplet.m_nMid, m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsChild2 = m_provider.GetVariantList(eCHILD, triplet.m_nCid,m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        
        std::vector<const core::COrientedVariant*> ovarListAMMother = m_provider.GetOrientedVariantList(eMOTHER, triplet.m_nMid, true, m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());
        std::vector<const core::COrientedVariant*> ovarListAMChildMC = m_provider.GetOrientedVariantList(eCHILD, triplet.m_nCid, true, m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        
        //Clear Mother child replay object
        replayMotherChildGT.Clear();
        //Change the variant list to process
        core::CPathReplay replayMotherChildAM(excludedVarsMother, excludedVarsChild2, ovarListAMMother, ovarListAMChildMC);
        
        //Find Best Path Mother-Child AM Match
        m_aBestPathsMotherChildAM[triplet.m_nTripleIndex] = replayMotherChildAM.FindBestPath(ctg, false);
        const std::vector<const core::COrientedVariant*>& includedVarsChildAMMC = m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetIncludedVariants();
        const std::vector<const core::COrientedVariant*>& includedVarsMotherAM = m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants();

        const std::vector<const CVariant*> excludedVarsChildMC = m_provider.GetVariantList(excludedVarsChild2, m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded());
        const std::vector<const CVariant*> excludedVarsMotherMC = m_provider.GetVariantList(excludedVarsMother, m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded());

        
        //Set Variant status of child variants
        m_provider.SetVariantStatus(includedVarsChildAMMC, eALLELE_MATCH);
        m_provider.SetVariantStatus(includedVarsChildGTMC, eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(excludedVarsChildMC, eNO_MATCH);
        
        //Set Variant status of mother variants
        m_provider.SetVariantStatus(includedVarsMotherGT, eGENOTYPE_MATCH);
        m_provider.SetVariantStatus(includedVarsMotherAM, eALLELE_MATCH);
        m_provider.SetVariantStatus(excludedVarsMotherMC, eNO_MATCH);
        
    
        //Clear Father child replay object
        replayMotherChildAM.Clear();
        
        //Lock the logging mechanism
        mtx.lock();
        
        //Send TP/FP/FN values to the log file
        m_resultLog.LogBestPathStatistic(true,
                                         triplet,
                                         static_cast<int>(includedVarsChildGT.size() + includedVarsChildAM.size()),
                                         static_cast<int>(m_aBestPathsFatherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants().size() +
                                                          m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants().size()),
                                         static_cast<int>(m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded().size()),
                                         static_cast<int>(m_aBestPathsFatherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded().size()));
        
        m_resultLog.LogBestPathStatistic(false,
                                         triplet,
                                         static_cast<int>(includedVarsChildGTMC.size() + includedVarsChildAMMC.size()),
                                         static_cast<int>(m_aBestPathsMotherChildGT[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants().size() +
                                                          m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetIncludedVariants().size()),
                                         static_cast<int>(m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_calledSemiPath.GetExcluded().size()),
                                         static_cast<int>(m_aBestPathsMotherChildAM[triplet.m_nTripleIndex].m_baseSemiPath.GetExcluded().size()));
        
        //Unlock the logging mechanism
        mtx.unlock();
        
        if(!ctg.Clean())
        {
            std::cerr << ctg.m_chromosomeName << " not cleaned.." << std::endl;
        }
    }
}

void CMendelianAnalyzer::PrintHelp() const
{
    std::cout << std::endl;
    std::cout << " --- MENDELIAN PARAMETERS --- " << std::endl;
    std::cout << "-father <father_vcf_path>    [Required.Add father VCF file.]" << std::endl;
    std::cout << "-mother <mother_vcf_path>    [Required.Add mother VCF file.]" << std::endl;
    std::cout << "-child <child_vcf_path>      [Required.Add child VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-outDir <output_directory>   [Required.Add output directory]" << std::endl;
    std::cout << "-pedigree <PED_file_path>    [Optional.Indentifies parent-child indexes from given PED file" << std::endl;
    std::cout << "-no-call <no_call_mode>      [Optional. Decides what to do with no call variants. There are 3 modes:" << std::endl;
    std::cout << "\t" << "implicit : mark boths implicit and explicit no call variant as NoCall" << std::endl;
    std::cout << "\t" << "explicit : mark explicit no call variants only as NoCall. Implicit no call variants will be treated as 0/0" << std::endl;
    std::cout << "\t" << "none : [Default Value] Treat all of no call variants as 0/0" << std::endl;
    std::cout << "-filter <filter_name>        [Optional.Filter variants based on filter column. Default value is PASS. Use 'none' to unfilter]" << std::endl;
    std::cout << "--ref-overlap                [Optional.Allow reference overlapping by trimming nucleotides and ignoring 0 genotype.]" << std::endl;
    std::cout << "--trim-endings-first         [Optional.If set, starts trimming variants from ending base pairs. Default is from beginning]" << std::endl;
    std::cout << "-sample-father <sample_name>  [Optional.Read only the given sample in father VCF. Default value is the first sample.]" << std::endl;
    std::cout << "-sample-mother <sample_name>  [Optional.Read only the given sample in mother VCF. Default value is the first sample.]" << std::endl;
    std::cout << "-sample-child <sample_name>   [Optional.Read only the given sample in child VCF. Default value is the first sample.]" << std::endl;
    std::cout << "-thread-count                 [Optional.Specify the number of threads that program will use. Default value is 2]" << std::endl;
    std::cout << std::endl;
    std::cout << "Example Commands:" << std::endl;
    std::cout << "./vbt mendelian -mother mother.vcf -father father.vcf -child child.vcf -ref reference.fasta -outDir SampleResultDir -filter none -no-call explicit" << std::endl;
    std::cout << "./vbt mendelian -mother trio.vcf -father trio.vcf -child trio.vcf -ref reference.fasta -outDir SampleResultDir -filter PASS -sample-mother mother -sample-father father -sample-child child" << std::endl;
}




