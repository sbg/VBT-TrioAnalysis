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
#include <algorithm>


//#include <fstream>


//Compare variants according to id for sort operation
bool variantCompare(const CVariant* v1, const CVariant* v2) {return v1->m_nId < v2->m_nId;}

//Checks if the given two range is overlapping
bool isOverlap(int left1, int right1, int left2, int right2)
{
    return std::min(right1, right2) - std::max(left1, left2) >= 0;
}


CMendelianAnalyzer::CMendelianAnalyzer()
{
    m_noCallMode = ENoCallMode::eNone;
}

void CMendelianAnalyzer::run(int argc, char **argv)
{
    std::time_t start, start1;
    double duration;
    
    //Reads the command line parameters
    bool isSuccess = ReadParameters(argc, argv);
    
    if(!isSuccess)
        return;
    
    start = std::time(0);
    
    //Initialize variant provider
    isSuccess = m_provider.InitializeReaders(m_fatherChildConfig, m_motherChildConfig);
    
    if(!isSuccess)
        return;

    duration = std::difftime(std::time(0), start);
    std::cout << "Vcf and fasta Parser read completed in " << duration << " secs" << std::endl;
    start1 = std::time(0);
    
    std::cout << "Running best path algorithm pipeline for each chromosome..." << std::endl;
    
    //Run duo comparison engine on parallel
    int threadCount = AssignJobsToThreads(m_fatherChildConfig.m_nThreadCount);
    for(int k = 0; k < threadCount; k++)
        m_pThreadPool[k].join();
    
    std::cout << "Evaluating mendelian consistency of variants..." << std::endl;
    
    //Perform merge process
    std::vector<int> chrIds = m_provider.GetCommonChromosomes();
    for(int k = 0; k < (int)chrIds.size(); k++)
        MergeFunc(chrIds[k]);
    
    std::cout << "Generating the output trio vcf..." << std::endl;
    
    //Initialize output writer
    std::string trioPath = std::string(m_fatherChildConfig.m_pOutputDirectory) + "/trio.vcf";
    m_trioWriter.SetTrioPath(trioPath);
    m_trioWriter.SetNoCallMode(m_noCallMode);
    m_trioWriter.SetResultLogPointer(&m_resultLog);
    
    for(int k = 0; k < (int)chrIds.size(); k++)
    {
        m_trioWriter.SetVariants(chrIds[k], eFATHER, m_provider.GetVariantList(eFATHER, chrIds[k]));
        m_trioWriter.SetVariants(chrIds[k], eMOTHER, m_provider.GetVariantList(eMOTHER, chrIds[k]));
        m_trioWriter.SetVariants(chrIds[k], eCHILD,  m_provider.GetVariantList(eCHILD,  chrIds[k]));
        
        m_trioWriter.SetDecisions(chrIds[k], eCHILD,  m_aChildDecisions[chrIds[k]]);
        m_trioWriter.SetDecisions(chrIds[k], eMOTHER, m_aMotherDecisions[chrIds[k]]);
        m_trioWriter.SetDecisions(chrIds[k], eFATHER, m_aFatherDecisions[chrIds[k]]);
    }
    m_trioWriter.GenerateTrioVcf();
    
    
    //Write results to log file
    m_resultLog.SetLogDirectory(m_fatherChildConfig.m_pOutputDirectory);
    m_resultLog.WriteBestPathStatistics();
    m_resultLog.WriteDetailedReportTable();
    m_resultLog.WriteShortReportTable();
    
    
    duration = std::difftime(std::time(0), start1);
    std::cout << "Processing Chromosomes completed in " << duration << " secs" << std::endl;
    duration = std::difftime(std::time(0), start);
    std::cout << "Total execution time is " << duration << " secs" << std::endl;

}

bool CMendelianAnalyzer::ReadParameters(int argc, char **argv)
{
    const char* PARAM_HELP = "--help";
    
    const char* PARAM_FATHER = "-father";
    const char* PARAM_MOTHER = "-mother";
    const char* PARAM_CHILD = "-child";
    const char* PARAM_REFERENCE = "-ref";
    const char* PARAM_FILTER = "-filter";
    
    const char* PARAM_SAMPLE_FATHER = "-sampleFather";
    const char* PARAM_SAMPLE_MOTHER = "-sampleMother";
    const char* PARAM_SAMPLE_CHILD = "-sampleChild";
    
    const char* PARAM_OUTPUT_DIR = "-outDir";
    const char* PARAM_REF_OVERLAP = "--ref-overlap";
    const char* PARAM_PLATFORM = "--platform-mode";
    const char* PARAM_THREAD_COUNT = "-thread-count";
    const char* PARAM_NO_CALL = "-no-call";
    
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
            m_motherChildConfig.m_nThreadCount = CHROMOSOME_COUNT;
            m_fatherChildConfig.m_nThreadCount = CHROMOSOME_COUNT;
            it--;
        }
        
        else if(0 == strcmp(argv[it], PARAM_THREAD_COUNT))
        {
            m_motherChildConfig.m_nThreadCount = std::min(std::max(1, atoi(argv[it+1])), CHROMOSOME_COUNT);
            m_fatherChildConfig.m_nThreadCount = std::min(std::max(1, atoi(argv[it+1])), CHROMOSOME_COUNT);
            it+=2;
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

void CMendelianAnalyzer::GetSyncPointList(int a_nChrId, bool a_bIsFatherChild, std::vector<CSyncPoint>& a_rSyncPointList, bool a_bIsGT)
{
    std::vector<const COrientedVariant*> pBaseIncluded;
    std::vector<const COrientedVariant*> pCalledIncluded;
    
    std::vector<const CVariant*> pBaseExcluded;
    std::vector<const CVariant*> pCalledExcluded;
    
    CPath *pPath;

    
    if(a_bIsGT == false)
    {
        pPath = a_bIsFatherChild ? &m_aBestPathsFatherChildAM[a_nChrId] : &m_aBestPathsMotherChildAM[a_nChrId];
        CPath *pPathGT = a_bIsFatherChild ? &m_aBestPathsFatherChildGT[a_nChrId] : &m_aBestPathsMotherChildGT[a_nChrId];
        
        std::vector<const CVariant*> excludedVarsBase = m_provider.GetVariantList(a_bIsFatherChild ? eFATHER : eMOTHER, a_nChrId, pPathGT->m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsCalled = m_provider.GetVariantList(eCHILD, a_nChrId, pPathGT->m_calledSemiPath.GetExcluded());

        pBaseIncluded = pPath->m_baseSemiPath.GetIncludedVariants();
        pCalledIncluded = pPath->m_calledSemiPath.GetIncludedVariants();
        
        pBaseExcluded = m_provider.GetVariantList(excludedVarsBase, pPath->m_baseSemiPath.GetExcluded());
        pCalledExcluded = m_provider.GetVariantList(excludedVarsCalled, pPath->m_calledSemiPath.GetExcluded());
    }
    else
    {
        pPath = a_bIsFatherChild ? &m_aBestPathsFatherChildGT[a_nChrId] : &m_aBestPathsMotherChildGT[a_nChrId];
        
        pBaseIncluded = pPath->m_baseSemiPath.GetIncludedVariants();
        pCalledIncluded = pPath->m_calledSemiPath.GetIncludedVariants();
        
        pBaseExcluded = m_provider.GetVariantList(a_bIsFatherChild ? eFATHER : eMOTHER, a_nChrId, pPath->m_baseSemiPath.GetExcluded());
        pCalledExcluded = m_provider.GetVariantList(eCHILD, a_nChrId, pPath->m_calledSemiPath.GetExcluded());
    }
    
    
    int baseIncludedItr = 0;
    int baseExcludedItr = 0;
    int calledIncludedItr = 0;
    int calledExcludedItr = 0;

    for(int k = 0; k < (int)pPath->m_aSyncPointList.size(); k++)
    {
        CSyncPoint sPoint;
        sPoint.m_nStartPosition = k > 0 ? pPath->m_aSyncPointList[k-1] : 0;
        sPoint.m_nEndPosition = pPath->m_aSyncPointList[k];
        sPoint.m_nIndex = k;
        
        while(baseIncludedItr < (int)pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_baseVariantsIncluded.push_back(pBaseIncluded[baseIncludedItr]);
            baseIncludedItr++;
        }

        while(calledIncludedItr < (int)pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_calledVariantsIncluded.push_back(pCalledIncluded[calledIncludedItr]);
            calledIncludedItr++;
        }
        
        while(baseExcludedItr < (int)pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
            baseExcludedItr++;
        }
        
        while(calledExcludedItr < (int)pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
            calledExcludedItr++;
        }
        
        a_rSyncPointList.push_back(sPoint);
    }
    
    //Add Remaining variants to the last syncPoint
    CSyncPoint sPoint;
    sPoint.m_nStartPosition = pPath->m_aSyncPointList[pPath->m_aSyncPointList.size()-1];
    sPoint.m_nEndPosition = INT_MAX;
    sPoint.m_nIndex = static_cast<int>(pPath->m_aSyncPointList.size()-1);
    
    while(baseIncludedItr < (int)pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsIncluded.push_back(pBaseIncluded[baseIncludedItr]);
        baseIncludedItr++;
    }
    while(calledIncludedItr < (int)pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsIncluded.push_back(pCalledIncluded[calledIncludedItr]);
        calledIncludedItr++;
    }
    
    while(baseExcludedItr < (int)pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
        baseExcludedItr++;
    }
    
    while(calledExcludedItr < (int)pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= sPoint.m_nEndPosition)
    {
        sPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
        calledExcludedItr++;
    }
    a_rSyncPointList.push_back(sPoint);
    
}

void CMendelianAnalyzer::CheckFor0PathFor00(int a_nChrId,
                                            bool a_bIsFatherChild,
                                            std::vector<const CVariant*>& a_rVarList,
                                            std::vector<const CVariant*>& a_rViolationList,
                                            std::vector<const CVariant*>& a_rCompliantList)
{
    //Get sync point list
    std::vector<CSyncPoint> a_rSyncPointList;
    GetSyncPointList(a_nChrId, a_bIsFatherChild, a_rSyncPointList, true);
    
    int varlistItr = 0;
    for(int k = 0; k < (int)a_rSyncPointList.size() && varlistItr < (int)a_rVarList.size(); k++)
    {
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < (int)a_rVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_rVarList[varlistItr]->m_nStartPos)
        {
            a_rViolationList.push_back(a_rVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(int m = 0; m < (int)a_rSyncPointList[k].m_calledVariantsIncluded.size(); m++)
        {
            if(a_rSyncPointList[k].m_calledVariantsIncluded[m]->GetVariant().m_nId == a_rVarList[varlistItr]->m_nId)
            {
                tmpVarList.push_back(a_rVarList[varlistItr]);
                varlistItr++;
                bDoCheck = true;
            }
        }
        
        if(true == bDoCheck)
        {
            bool bIsCompliant = true;
            
            for(int m = 0; m < (int)a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                {
                    if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                    {
                        bIsCompliant = false;
                        break;
                    }
                }
            }
            
            if(bIsCompliant == true)
            {
                for(int aa = 0; aa < (int)tmpVarList.size(); aa++)
                    a_rCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(int aa = 0; aa < (int)tmpVarList.size(); aa++)
                    a_rViolationList.push_back(tmpVarList[aa]);
            }
        }
    }
    
}


void CMendelianAnalyzer::CheckFor0Path(int a_nChrId,
                                       bool a_bIsFatherChild,
                                       std::vector<const CVariant *> &a_pVarList,
                                       std::vector<const CVariant *> &a_pViolantList,
                                       std::vector<const CVariant *> &a_pCompliantList,
                                       bool a_bIsUpdateDecisionList)
{
    
    //Get sync point list
    std::vector<CSyncPoint> a_rSyncPointList;
    GetSyncPointList(a_nChrId, a_bIsFatherChild, a_rSyncPointList);
    
    
    int varlistItr = 0;
    for(int k = 0; k < (int)a_rSyncPointList.size() && varlistItr < (int)a_pVarList.size(); k++)
    {
        
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < (int)a_pVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_pVarList[varlistItr]->m_nStartPos)
        {
            a_pViolantList.push_back(a_pVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(int m = 0; m < (int)a_rSyncPointList[k].m_calledVariantsExcluded.size(); m++)
        {
            if(a_rSyncPointList[k].m_calledVariantsExcluded[m]->m_nId == a_pVarList[varlistItr]->m_nId)
            {
                tmpVarList.push_back(a_pVarList[varlistItr]);
                varlistItr++;
                bDoCheck = true;
            }
        }
        
        if(true == bDoCheck)
        {
            bool bIsCompliant = true;
            
            for(int m = 0; m < (int)a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                {
                    
                    if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                    {
                        if(a_bIsUpdateDecisionList)
                        {
                            //We are marking decision of mother/father variant as violation
                            if(true == a_bIsFatherChild)
                                m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, pVar->m_nId)] = eViolation;
                            else
                                m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, pVar->m_nId)] = eViolation;
                        }
                        
                        bIsCompliant = false;
                        break;
                    }
                    
                    else
                    {
                        if(a_bIsUpdateDecisionList)
                        {
                            //We are marking decision of mother/father variant as compliant
                            if(true == a_bIsFatherChild)
                                m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, pVar->m_nId)] = eCompliant;
                            else
                                m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, pVar->m_nId)] = eCompliant;
                        }
                    }
                }
            }
            
            if(bIsCompliant == true)
            {
                for(int aa = 0; aa < (int)tmpVarList.size(); aa++)
                    a_pCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(int aa = 0; aa < (int)tmpVarList.size(); aa++)
                    a_pViolantList.push_back(tmpVarList[aa]);
            }
        }
    }
    
}

void CMendelianAnalyzer::CheckFor00Child(int a_nChrId,
                                         std::vector<const CVariant*>& a_rOvarList,
                                         std::vector<const CVariant*>& a_rViolationList,
                                         std::vector<const CVariant*>& a_rCompliantList,
                                         bool a_bIsGTMatch)
{
    std::vector<const CVariant*> fatherCompliants;
    std::vector<const CVariant*> motherCompliants;
    std::vector<const CVariant*> fatherViolants;
    std::vector<const CVariant*> motherViolants;
    
    std::sort(a_rOvarList.begin(), a_rOvarList.end(), variantCompare);
    
    if(false == a_bIsGTMatch)
    {
        CheckFor0Path(a_nChrId, true,  a_rOvarList, fatherViolants, fatherCompliants, false);
        CheckFor0Path(a_nChrId, false, a_rOvarList, motherViolants, motherCompliants, false);
    }
    else
    {
        CheckFor0PathFor00(a_nChrId, true,  a_rOvarList, fatherViolants, fatherCompliants);
        CheckFor0PathFor00(a_nChrId, false, a_rOvarList, motherViolants, motherCompliants);
    }
    
    //Intersect father and mother compliant variants and write to a_rCompliantList
    for (std::vector<const CVariant*>::iterator i = fatherCompliants.begin(); i != fatherCompliants.end(); ++i)
    {
        if (std::find(motherCompliants.begin(), motherCompliants.end(), *i) != motherCompliants.end())
        {
            a_rCompliantList.push_back(*i);
        }
    }
    
    //Write all other variants to a_rViolationList
    for (std::vector<const CVariant*>::iterator i = a_rOvarList.begin(); i != a_rOvarList.end(); ++i)
    {
        if (std::find(a_rCompliantList.begin(), a_rCompliantList.end(), *i) == a_rCompliantList.end())
        {
            a_rViolationList.push_back(*i);
        }
    }

}


int CMendelianAnalyzer::AssignJobsToThreads(int a_nThreadCount)
{
    //Get the list of chromosomes to be processed
    std::vector<int> chromosomeListToProcess = m_provider.GetCommonChromosomes();

    int exactThreadCount = std::min(a_nThreadCount, (int)chromosomeListToProcess.size());
    
    //Allocate threads
    m_pThreadPool = new std::thread[exactThreadCount];
    
    int threadPoolIt = 0;
    std::vector<int> *chromosomeLists = new std::vector<int>[exactThreadCount];
    
    //Divide tasks into threads
    for(int k = 0; k < (int)chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(chromosomeListToProcess[k]);
        threadPoolIt = (threadPoolIt+1) % exactThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < exactThreadCount; k++)
    {
        m_pThreadPool[k] = std::thread(&CMendelianAnalyzer::ProcessChromosome, this, chromosomeLists[k]);
    }
    
    //Clean allocation
    delete[] chromosomeLists;
    
    return exactThreadCount;
    
}

void CMendelianAnalyzer::ProcessChromosome(const std::vector<int>& a_nChromosomeIds)
{    
    for(int chromosomeId : a_nChromosomeIds)
    {
    
        //Get variant list of parent-child for given chromosome
        std::vector<const CVariant*> varListFather = m_provider.GetVariantList(eFATHER, chromosomeId);
        std::vector<const CVariant*> varListMother = m_provider.GetVariantList(eMOTHER, chromosomeId);
        std::vector<const CVariant*> varListChild = m_provider.GetVariantList(eCHILD, chromosomeId);
        
        //Get oriented variant list of parent-child for given chromosome
        std::vector<const COrientedVariant*> ovarListGTFather = m_provider.GetOrientedVariantList(eFATHER, chromosomeId);
        std::vector<const COrientedVariant*> ovarListGTMother = m_provider.GetOrientedVariantList(eMOTHER, chromosomeId);
        std::vector<const COrientedVariant*> ovarListGTChild = m_provider.GetOrientedVariantList(eCHILD, chromosomeId);

        //Get the chromosome ref seq
        SContig ctg;
        m_provider.GetContig(chromosomeId, ctg);
        
        
        // === PROCESS FATHER-CHILD ===
        
        //Create path replay for parent child;
        CPathReplay replayFatherChildGT(varListFather, varListChild, ovarListGTFather, ovarListGTChild);
        
        //Find Best Path Father-Child GT Match
        m_aBestPathsFatherChildGT[chromosomeId] = replayFatherChildGT.FindBestPath(ctg, true);
        
        //Genotype Match variants
        const std::vector<const COrientedVariant*>& includedVarsChildGT = m_aBestPathsFatherChildGT[chromosomeId].m_calledSemiPath.GetIncludedVariants();
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsFather = m_provider.GetVariantList(eFATHER, chromosomeId, m_aBestPathsFatherChildGT[chromosomeId].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsChild = m_provider.GetVariantList(eCHILD, chromosomeId,m_aBestPathsFatherChildGT[chromosomeId].m_calledSemiPath.GetExcluded());
        
        //Allele Match oriented variants
        std::vector<const COrientedVariant*> ovarListAMFather = m_provider.GetOrientedVariantList(eFATHER, chromosomeId, true, m_aBestPathsFatherChildGT[chromosomeId].m_baseSemiPath.GetExcluded());
        std::vector<const COrientedVariant*> ovarListAMChildFC = m_provider.GetOrientedVariantList(eCHILD, chromosomeId, true, m_aBestPathsFatherChildGT[chromosomeId].m_calledSemiPath.GetExcluded());
        
        //Clear Father child replay object
        replayFatherChildGT.Clear();
        
        //Change the variant list to process
        CPathReplay replayFatherChildAM(excludedVarsFather, excludedVarsChild, ovarListAMFather, ovarListAMChildFC);
        
        //Find Best Path Father-Child AM Match
        m_aBestPathsFatherChildAM[chromosomeId] = replayFatherChildAM.FindBestPath(ctg, false);
        const std::vector<const COrientedVariant*>& includedVarsChildAM = m_aBestPathsFatherChildAM[chromosomeId].m_calledSemiPath.GetIncludedVariants();

        //Set Variant status of child variants
        m_provider.SetVariantStatus(includedVarsChildAM, eALLELE_MATCH);
        m_provider.SetVariantStatus(includedVarsChildGT, eGENOTYPE_MATCH);
        
        //Clear Father child replay object
        replayFatherChildAM.Clear();
        
        
        // === PROCESS MOTHER-CHILD ===
     
        //Create path replay for parent child;
        CPathReplay replayMotherChildGT(varListMother, varListChild, ovarListGTMother, ovarListGTChild);
        
        //Find Best Path Father-Child GT Match
        m_aBestPathsMotherChildGT[chromosomeId] = replayMotherChildGT.FindBestPath(ctg, true);
        
        //Genotype Match variants
        const std::vector<const COrientedVariant*>& includedVarsChildGTMC = m_aBestPathsMotherChildGT[chromosomeId].m_calledSemiPath.GetIncludedVariants();
        //Variants that will be passed for allele match check
        std::vector<const CVariant*> excludedVarsMother = m_provider.GetVariantList(eMOTHER, chromosomeId, m_aBestPathsMotherChildGT[chromosomeId].m_baseSemiPath.GetExcluded());
        std::vector<const CVariant*> excludedVarsChild2 = m_provider.GetVariantList(eCHILD, chromosomeId,m_aBestPathsMotherChildGT[chromosomeId].m_calledSemiPath.GetExcluded());
        
        std::vector<const COrientedVariant*> ovarListAMMother = m_provider.GetOrientedVariantList(eMOTHER, chromosomeId, true, m_aBestPathsMotherChildGT[chromosomeId].m_baseSemiPath.GetExcluded());
        std::vector<const COrientedVariant*> ovarListAMChildMC = m_provider.GetOrientedVariantList(eCHILD, chromosomeId, true, m_aBestPathsMotherChildGT[chromosomeId].m_calledSemiPath.GetExcluded());
        
        //Clear Mother child replay object
        replayMotherChildGT.Clear();
        //Change the variant list to process
        CPathReplay replayMotherChildAM(excludedVarsMother, excludedVarsChild2, ovarListAMMother, ovarListAMChildMC);
        
        //Find Best Path Mother-Child AM Match
        m_aBestPathsMotherChildAM[chromosomeId] = replayMotherChildAM.FindBestPath(ctg, false);
        const std::vector<const COrientedVariant*>& includedVarsChildAMMC = m_aBestPathsMotherChildAM[chromosomeId].m_calledSemiPath.GetIncludedVariants();
        
        //Set Variant status of child variants
        m_provider.SetVariantStatus(includedVarsChildAMMC, eALLELE_MATCH);
        m_provider.SetVariantStatus(includedVarsChildGTMC, eGENOTYPE_MATCH);
        
        
        //Clear Father child replay object
        replayMotherChildAM.Clear();
        
        //Initialize the decision arrays
        m_aMotherDecisions[chromosomeId] = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eMOTHER, chromosomeId));
        m_aFatherDecisions[chromosomeId] = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eFATHER, chromosomeId));
        m_aChildDecisions[chromosomeId]  = std::vector<EMendelianDecision>(m_provider.GetVariantCount(eCHILD,  chromosomeId));
        
        //Send TP/FP/FN values to the log file
        m_resultLog.LogBestPathStatistic(true,
                                         chromosomeId,
                                         static_cast<int>(includedVarsChildGT.size() + includedVarsChildAM.size()),
                                         static_cast<int>(m_aBestPathsFatherChildGT[chromosomeId].m_baseSemiPath.GetIncludedVariants().size() + m_aBestPathsFatherChildAM[chromosomeId].m_baseSemiPath.GetIncludedVariants().size()),
                                         static_cast<int>(m_aBestPathsFatherChildAM[chromosomeId].m_calledSemiPath.GetExcluded().size()),
                                         static_cast<int>(m_aBestPathsFatherChildAM[chromosomeId].m_baseSemiPath.GetExcluded().size()));
        
        m_resultLog.LogBestPathStatistic(false,
                                         chromosomeId,
                                         static_cast<int>(includedVarsChildGTMC.size() + includedVarsChildAMMC.size()),
                                         static_cast<int>(m_aBestPathsMotherChildGT[chromosomeId].m_baseSemiPath.GetIncludedVariants().size() + m_aBestPathsMotherChildAM[chromosomeId].m_baseSemiPath.GetIncludedVariants().size()),
                                         static_cast<int>(m_aBestPathsMotherChildAM[chromosomeId].m_calledSemiPath.GetExcluded().size()),
                                         static_cast<int>(m_aBestPathsMotherChildAM[chromosomeId].m_baseSemiPath.GetExcluded().size()));
        
    }

}

void CMendelianAnalyzer::CheckUniqueVars(EMendelianVcfName a_checkSide, int a_nChrId, const std::vector<const CVariant*>& a_rVariantList, std::vector<bool>& a_rSideDecisions)
{
    
    std::vector<const CVariant*> varListToCheckChild = m_provider.GetVariantList(eCHILD, a_nChrId);;
    std::vector<const CVariant*> varListToCheckParent = a_checkSide == eMOTHER ? m_provider.GetVariantList(eFATHER, a_nChrId) : m_provider.GetVariantList(eMOTHER, a_nChrId);
    
    std::vector<bool> decChild(a_rVariantList.size());
    std::vector<bool> decParent(a_rVariantList.size());
    
    //initialize all variants as true at first. We will then eliminate those with overlapping non-0 variants.
    for(int k = 0; k < (int)a_rVariantList.size(); k++)
    {
        decChild[k] = true;
        decParent[k] = true;
    }
    
    
    //Check overlaps for child
    int varItr = 0;
    for(int k = 0; k < (int)a_rVariantList.size(); k++)
    {
        //Check if the variant is already marked
        if(a_checkSide == eMOTHER && m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, a_rVariantList[k]->m_nId)] != eUnknown)
        {
            decChild[k] = m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, a_rVariantList[k]->m_nId)] == eCompliant ? true : false;
            continue;
        }
        //Check if the variant is already marked
        else if(a_checkSide == eFATHER && m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, a_rVariantList[k]->m_nId)] != eUnknown)
        {
            decChild[k] = m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, a_rVariantList[k]->m_nId)] == eCompliant ? true : false;
            continue;
        }

        
        if(a_rVariantList[k]->m_genotype[0] != 0 && a_rVariantList[k]->m_genotype[1] != 0)
        {
            decChild[k] = false;
            continue;
        }
        
        while(varItr <  (int)varListToCheckChild.size() && varListToCheckChild[varItr]->m_nEndPos < a_rVariantList[k]->m_nStartPos)
            varItr++;
    
        if(varItr == (int)varListToCheckChild.size())
            break;
        
        else if(isOverlap(a_rVariantList[k]->m_nStartPos, a_rVariantList[k]->m_nEndPos, varListToCheckChild[varItr]->m_nStartPos, varListToCheckChild[varItr]->m_nEndPos))
        {
            if(m_aChildDecisions[a_nChrId][varItr] == eCompliant)
                decChild[k] = true;
            
            else if(m_aChildDecisions[a_nChrId][varItr] == eViolation)
                decChild[k] = false;
            
            else if(varListToCheckChild[varItr]->m_genotype[0] != 0 && varListToCheckChild[varItr]->m_genotype[1] != 0)
                decChild[k] = false;
            
        }
    }
    
    
    //Check overlaps for other parent
    varItr = 0;
    for(int k = 0; k < (int)a_rVariantList.size(); k++)
    {
        //Check if the variant is already marked
        if(a_checkSide == eMOTHER && m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, a_rVariantList[k]->m_nId)] != eUnknown)
        {
            decParent[k] = m_aMotherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChrId, a_rVariantList[k]->m_nId)] == eCompliant ? true : false;
            continue;
        }
        //Check if the variant is already marked
        else if(a_checkSide == eFATHER && m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, a_rVariantList[k]->m_nId)] != eUnknown)
        {
            decParent[k] = m_aFatherDecisions[a_nChrId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChrId, a_rVariantList[k]->m_nId)] == eCompliant ? true : false;
            continue;
        }
        
        
        if(a_rVariantList[k]->m_genotype[0] != 0 && a_rVariantList[k]->m_genotype[1] != 0)
        {
            decParent[k] = false;
            continue;
        }
        
        while(varItr <  (int)varListToCheckParent.size() && varListToCheckParent[varItr]->m_nEndPos < a_rVariantList[k]->m_nStartPos)
            varItr++;
        
        if(varItr == (int)a_rVariantList.size())
            break;
        
        else if(isOverlap(a_rVariantList[k]->m_nStartPos, a_rVariantList[k]->m_nEndPos, varListToCheckParent[varItr]->m_nStartPos, varListToCheckParent[varItr]->m_nEndPos))
        {
            if(varListToCheckParent[varItr]->m_genotype[0] != 0 && varListToCheckParent[varItr]->m_genotype[1] != 0)
            {
                decParent[k] = false;
            }
        }
    }
    
    //Output the final decision
    for(int k = 0; k < (int)a_rVariantList.size(); k++)
    {
        a_rSideDecisions[k] = decChild[k] && decParent[k];
    }
}


void CMendelianAnalyzer::MergeFunc(int a_nChromosomeId)
{
    std::vector<const CVariant*> compliants;
    std::vector<const CVariant*> violations;

    std::vector<const CVariant*> SameAlleleMatchViolationVars;
    std::vector<const CVariant*> MendelianCompliantVars;
    
    std::vector<const CVariant*> motherChildOnly;
    std::vector<const CVariant*> fatherChildOnly;
    
    std::vector<const CVariant*> check0atMotherSide;
    std::vector<const CVariant*> check0atFatherSide;
    std::vector<const CVariant*> check00Child;
    std::vector<const CVariant*> childUniqueList;
    std::vector<const CVariant*> check00ChildGTMatched;

    //Included lists of child
    CVariantIterator FatherChildVariants(m_aBestPathsFatherChildGT[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants(),
                                         m_aBestPathsFatherChildAM[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants());
    
    CVariantIterator MotherChildVariants(m_aBestPathsMotherChildGT[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants(),
                                         m_aBestPathsMotherChildAM[a_nChromosomeId].m_calledSemiPath.GetIncludedVariants());

    //Check if the two list have common variants
    if(FatherChildVariants.hasNext() == false && MotherChildVariants.hasNext() == false)
        return;
    
    const COrientedVariant* varMC = MotherChildVariants.Next();
    const COrientedVariant* varFC = FatherChildVariants.Next();
    
    while(true)
    {
        //If Mother variant and Father variant is Common check for same allele match condition
        if(varMC->GetVariant().m_nId == varFC->GetVariant().m_nId)
        {
            if(varMC->GetVariant().m_genotype[0] == 0 && varMC->GetVariant().m_genotype[1] == 0)
                check00ChildGTMatched.push_back(&varMC->GetVariant());
            
            else if(varMC->GetVariant().m_bIsHeterozygous)
            {
                //ELIMINATE SAME ALLELE MATCHING(EXCEPTION 1)
                if(varMC->GetAlleleIndex() == varFC->GetAlleleIndex())
                {
                    if(varMC->GetVariant().m_variantStatus == eGENOTYPE_MATCH
                       ||
                       varFC->GetVariant().m_variantStatus == eGENOTYPE_MATCH)
                        MendelianCompliantVars.push_back(&varMC->GetVariant());
                    else
                        SameAlleleMatchViolationVars.push_back(&varMC->GetVariant());
                }
                else
                    MendelianCompliantVars.push_back(&varMC->GetVariant());
            }
            else
                MendelianCompliantVars.push_back(&varMC->GetVariant());
            
            if(MotherChildVariants.hasNext() && FatherChildVariants.hasNext())
            {
                varMC = MotherChildVariants.Next();
                varFC = FatherChildVariants.Next();
            }
            
            else
            {
                varMC = NULL;
                varFC = NULL;
                break;
            }
        }
        
        //If we have variant match with Father side only, we filter 0/x variants and rest of them are marked as violation
        else if(varMC->GetVariant().m_nId > varFC->GetVariant().m_nId)
        {
            if(varFC->GetVariant().m_bIsHeterozygous)
            {
                if(varFC->GetVariant().m_genotype[0] == 0 || varFC->GetVariant().m_genotype[1] == 0)
                    check0atMotherSide.push_back(&varFC->GetVariant());
                else
                    fatherChildOnly.push_back(&varFC->GetVariant());
            }
            else
            {
                if(varFC->GetVariant().m_genotype[0] == 0)
                    check0atMotherSide.push_back(&varFC->GetVariant());
                else
                    fatherChildOnly.push_back(&varFC->GetVariant());
            }
            
            if(FatherChildVariants.hasNext())
                varFC = FatherChildVariants.Next();
            else
            {
                varFC = NULL;
                break;
            }
        }
        
        //If we have variant match with Mother side only, we filter 0/x variants and rest of them are marked as violation
        else
        {
            if(varMC->GetVariant().m_bIsHeterozygous)
            {
                if(varMC->GetVariant().m_genotype[0] == 0 || varMC->GetVariant().m_genotype[1] == 0)
                    check0atFatherSide.push_back(&varMC->GetVariant());
                else
                    motherChildOnly.push_back(&varMC->GetVariant());
            }
            else
            {
                if(varMC->GetVariant().m_genotype[0] == 0)
                    check0atFatherSide.push_back(&varMC->GetVariant());
                else
                    motherChildOnly.push_back(&varMC->GetVariant());
            }
            
            if(MotherChildVariants.hasNext())
                varMC = MotherChildVariants.Next();
            else
            {
                varMC = NULL;
                break;
            }
        }
        
    }
    
    
    //Process remaining vars in FatherChild explained as above
    while(true)
    {
        if(varFC == NULL)
        {
            if(FatherChildVariants.hasNext())
                varFC = FatherChildVariants.Next();
            else
                break;
        }
        
        if(varFC->GetVariant().m_bIsHeterozygous)
        {
            if(varFC->GetVariant().m_genotype[0] == 0 || varFC->GetVariant().m_genotype[1] == 0)
                check0atMotherSide.push_back(&varFC->GetVariant());
            else
                fatherChildOnly.push_back(&varFC->GetVariant());
        }
        else
        {
            if(varFC->GetVariant().m_genotype[0] == 0)
                check0atMotherSide.push_back(&varFC->GetVariant());
            else
                fatherChildOnly.push_back(&varFC->GetVariant());
        }
        
        if(!FatherChildVariants.hasNext())
            break;
        else
            varFC = FatherChildVariants.Next();
    }
    
    //Process remaining vars in MotherChild explined as above
    while(true)
    {
        if(varMC == NULL)
        {
            if(MotherChildVariants.hasNext())
                varMC = MotherChildVariants.Next();
            else
                break;
        }
        
        if(varMC->GetVariant().m_bIsHeterozygous)
        {
            if(varMC->GetVariant().m_genotype[0] == 0 || varMC->GetVariant().m_genotype[1] == 0)
                check0atFatherSide.push_back(&varMC->GetVariant());
            else
                motherChildOnly.push_back(&varMC->GetVariant());
        }
        else
        {
            if(varMC->GetVariant().m_genotype[0] == 0)
                check0atFatherSide.push_back(&varMC->GetVariant());
            else
                motherChildOnly.push_back(&varMC->GetVariant());
        }
        
        
        if(!MotherChildVariants.hasNext())
            break;
        else
            varMC = MotherChildVariants.Next();
    }
    
    std::vector<const CVariant*> compliantVarsFrom0CheckMother;
    std::vector<const CVariant*> compliantVarsFrom0CheckFather;
    std::vector<const CVariant*> violationVarsFrom0CheckMother;
    std::vector<const CVariant*> violationVarsFrom0CheckFather;
    
    //Check for 0/x child variant set at father side
    CheckFor0Path(a_nChromosomeId, true, check0atFatherSide, violationVarsFrom0CheckFather, compliantVarsFrom0CheckFather);
    //Check for 0/x child variant set at the mother side
    CheckFor0Path(a_nChromosomeId, false, check0atMotherSide, violationVarsFrom0CheckMother, compliantVarsFrom0CheckMother);
    
    std::vector<const CVariant*> compliantVarsFrom00CheckGT;
    std::vector<const CVariant*> violationVarsFrom00CheckGT;
    
    //Check for 0/0 child variant set for both parent
    CheckFor00Child(a_nChromosomeId, check00ChildGTMatched, violationVarsFrom00CheckGT, compliantVarsFrom00CheckGT, true);


    //Gather all compliant variants of child we found so far
    compliants.insert(std::end(compliants), std::begin(MendelianCompliantVars), std::end(MendelianCompliantVars));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckFather), std::end(compliantVarsFrom0CheckFather));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckMother), std::end(compliantVarsFrom0CheckMother));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom00CheckGT), std::end(compliantVarsFrom00CheckGT));
    std::sort(compliants.begin(), compliants.end(), variantCompare);

    
    //Gather all violation variants of child we found so far
    violations.insert(std::end(violations), std::begin(SameAlleleMatchViolationVars), std::end(SameAlleleMatchViolationVars));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckFather), std::end(violationVarsFrom0CheckFather));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckMother), std::end(violationVarsFrom0CheckMother));
    violations.insert(std::end(violations), std::begin(fatherChildOnly), std::end(fatherChildOnly));
    violations.insert(std::end(violations), std::begin(motherChildOnly), std::end(motherChildOnly));
    violations.insert(std::end(violations), std::begin(violationVarsFrom00CheckGT), std::end(violationVarsFrom00CheckGT));
    std::sort(violations.begin(), violations.end(), variantCompare);

    
    //Find Child Unique variants
    std::vector<const CVariant*> childVariants = m_provider.GetVariantList(eCHILD, a_nChromosomeId);
    std::vector<int>childProcessedArray(childVariants.size());
    for(int k = 0; k <  (int)childProcessedArray.size(); k++)
        childProcessedArray[k] = 0;
    
    //Mark mendelian compliant vars as processed
    for(int k = 0, m = 0; k < (int)childVariants.size() && m < (int)compliants.size(); k++)
    {
        if(childVariants[k]->m_nId == compliants[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    //Mark mendelian violation vars as processed
    for(int k = 0, m = 0; k < (int)childVariants.size() && m < (int)violations.size(); k++)
    {
        if(childVariants[k]->m_nId == violations[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    for(int childItr = 0; childItr < (int)childProcessedArray.size(); childItr++)
    {
        if(childProcessedArray[childItr] == 0)
        {
            if(childVariants[childItr]->m_genotype[0] == 0 && childVariants[childItr]->m_genotype[1] == 0)
            {
                check00Child.push_back(childVariants[childItr]);
            }
            else
            {
                const CVariant* pVar = childVariants[childItr];
                childUniqueList.push_back(pVar);
            }
        }
    }
    
    std::vector<const CVariant*> compliantVarsFrom00Check;
    std::vector<const CVariant*> violationVarsFrom00Check;
    
    //Check for 0/0 unique child variants for both parent
    CheckFor00Child(a_nChromosomeId, check00Child, violationVarsFrom00Check, compliantVarsFrom00Check, false);
    
    //Add the new compliants we found to compliants list
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom00Check), std::end(compliantVarsFrom00Check));
    std::sort(compliants.begin(), compliants.end(), variantCompare);

    //Add the new violations we found to violation list
    violations.insert(std::end(violations), std::begin(violationVarsFrom00Check), std::end(violationVarsFrom00Check));
    violations.insert(std::end(violations), std::begin(childUniqueList), std::end(childUniqueList));
    std::sort(violations.begin(), violations.end(), variantCompare);

    
    //We looked up all child variants. Now, we will look at parent variants where there is no corresponding child variant exist in the child.vcf (check for hidden 0/0 child variants)

    //Fill the child decision array
    std::vector<const CVariant*>::iterator compliantsIterator = compliants.begin();
    std::vector<const CVariant*>::iterator violationsIterator = violations.begin();
    for(int k = 0; k < (int)childVariants.size(); k++)
    {
        if(compliantsIterator != compliants.end() && childVariants[k]->m_nId == (*compliantsIterator)->m_nId)
        {
            m_aChildDecisions[a_nChromosomeId][k] = EMendelianDecision::eCompliant;
            compliantsIterator++;
        }
        else if(violationsIterator != violations.end() && childVariants[k]->m_nId == (*violationsIterator)->m_nId)
        {
            m_aChildDecisions[a_nChromosomeId][k] = EMendelianDecision::eViolation;
            violationsIterator++;
        }
        else
            m_aChildDecisions[a_nChromosomeId][k] = EMendelianDecision::eUnknown;
    }
    
    //Mother Checks
    std::vector<const CVariant*> uniqueMotherVars = m_provider.GetVariantList(m_provider.GetVariantList(eMOTHER, a_nChromosomeId,  m_aBestPathsMotherChildGT[a_nChromosomeId].m_baseSemiPath.GetExcluded()),
                                                                              m_aBestPathsMotherChildAM[a_nChromosomeId].m_baseSemiPath.GetExcluded());
    std::vector<bool> motherDecisions(uniqueMotherVars.size());
    CheckUniqueVars(eMOTHER, a_nChromosomeId, uniqueMotherVars, motherDecisions);
    
    //Father Checks
    std::vector<const CVariant*> uniqueFatherVars = m_provider.GetVariantList(m_provider.GetVariantList(eFATHER, a_nChromosomeId,  m_aBestPathsFatherChildGT[a_nChromosomeId].m_baseSemiPath.GetExcluded()),
                                                                              m_aBestPathsFatherChildAM[a_nChromosomeId].m_baseSemiPath.GetExcluded());
    std::vector<bool> fatherDecisions(uniqueFatherVars.size());
    CheckUniqueVars(eFATHER, a_nChromosomeId, uniqueFatherVars, fatherDecisions);

    //Fill the mother decision array
    for(int k = 0; k < (int)motherDecisions.size(); k++)
        m_aMotherDecisions[a_nChromosomeId][m_provider.Get0BasedVariantIndex(eMOTHER, a_nChromosomeId, uniqueMotherVars[k]->m_nId)] = (motherDecisions[k] ? eCompliant : eViolation);
    
    //Fill the father decision array
    for(int k = 0; k < (int)fatherDecisions.size(); k++)
        m_aFatherDecisions[a_nChromosomeId][m_provider.Get0BasedVariantIndex(eFATHER, a_nChromosomeId, uniqueFatherVars[k]->m_nId)] = (fatherDecisions[k] ? eCompliant : eViolation);

    
    
    //If NoCall Mode Is Enabled, mark all decisions of nocall childs as NoCallChild and all nocall parents as NoCallParent
    if(m_noCallMode != eNone)
    {
        std::vector<const CVariant*> motherVariants = m_provider.GetVariantList(eMOTHER, a_nChromosomeId);
        std::vector<const CVariant*> fatherVariants = m_provider.GetVariantList(eFATHER, a_nChromosomeId);

        for(int k = 0; k < (int)motherVariants.size(); k ++)
        {
            if(motherVariants[k]->m_bIsNoCall)
                m_aMotherDecisions[a_nChromosomeId][k] = eNoCallParent;
        }
        
        for(int k = 0; k < (int)fatherVariants.size(); k ++)
        {
            if(fatherVariants[k]->m_bIsNoCall)
                m_aFatherDecisions[a_nChromosomeId][k] = eNoCallParent;
        }
        
        for(int k = 0; k < (int)childVariants.size(); k ++)
        {
            if(childVariants[k]->m_bIsNoCall)
                m_aChildDecisions[a_nChromosomeId][k] = eNoCallChild;
        }
    }
    
    ReportChildChromosomeData(a_nChromosomeId, compliants, violations);
    
    std::cout << "===================== STATISTICS " << a_nChromosomeId + 1 << " ===================" << std::endl;
    std::cout << "Total Compliants:" << compliants.size() << std::endl;
    std::cout << "Total Violations:" << violations.size() << std::endl;
    std::cout << "Child Var Size:" << childVariants.size()<< std::endl;
    std::cout << "=====================================================" << std::endl << std::endl;

    /*
    std::string commonPath = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/CHR1/chr" + std::to_string(a_nChromosomeId + 1);

    std::ofstream compliantsAll;
    compliantsAll.open(commonPath + "_CompliantsALL.txt");
    for(const CVariant* k : compliants)
        compliantsAll << k->m_nOriginalPos << std::endl;
    compliantsAll.close();

    std::ofstream violationsAll;
    violationsAll.open(commonPath + "_ViolationsALL.txt");
    for(const CVariant* k : violations)
        violationsAll << k->m_nOriginalPos << std::endl;
    violationsAll.close();
    
    std::ofstream outputFileAlleleMatchVio;
    outputFileAlleleMatchVio.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21SameAlleleMatchViolation.txt");
    for(const CVariant* k : MendelianViolationVars)
        outputFileAlleleMatchVio << k->m_nOriginalPos << std::endl;
    outputFileAlleleMatchVio.close();

    std::ofstream outputFileMotherOnly;
    outputFileMotherOnly.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21MotherChildOnly.txt");
    for(const CVariant* k : motherChildOnly)
        outputFileMotherOnly << k->m_nOriginalPos << std::endl;
    outputFileMotherOnly.close();

    std::ofstream outputFileFatherOnly;
    outputFileFatherOnly.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21FatherChildOnly.txt");
    for(const CVariant* k : fatherChildOnly)
        outputFileFatherOnly << k->m_nOriginalPos << std::endl;
    outputFileFatherOnly.close();

    std::ofstream outputFile;
    outputFile.open(commonPath + "_InitialCompliants.txt");
    for(const CVariant* k : MendelianCompliantVars)
        outputFile << k->m_nOriginalPos << std::endl;
    outputFile.close();

    std::ofstream outputFileF;
    outputFileF.open(commonPath + "_Father0sCompliant.txt");
    for(const CVariant* k : compliantVarsFrom0CheckFather)
        outputFileF << k->m_nOriginalPos << std::endl;
    outputFileF.close();

    std::ofstream outputFileFvio;
    outputFileFvio.open(commonPath + "_Father0sViolation.txt");
    for(const CVariant* k : violationVarsFrom0CheckFather)
        outputFileFvio << k->m_nOriginalPos << std::endl;
    outputFileFvio.close();
    

    std::ofstream outputFileM;
    outputFileM.open(commonPath + "_Mother0sCompliant.txt");
    for(const CVariant* k : compliantVarsFrom0CheckMother)
        outputFileM << k->m_nOriginalPos << std::endl;
    outputFileM.close();

    std::ofstream outputFileMvio;
    outputFileMvio.open(commonPath + "_Mother0sViolation.txt");
    for(const CVariant* k : violationVarsFrom0CheckMother)
        outputFileMvio << k->m_nOriginalPos << std::endl;
    outputFileMvio.close();

    std::ofstream outputFileCC1;
    outputFileCC1.open(commonPath + "_00Compliants.txt");
    for(const CVariant* k : compliantVarsFrom00Check)
        outputFileCC1 << k->m_nOriginalPos << std::endl;
    outputFileCC1.close();

    std::ofstream outputFileCC1v;
    outputFileCC1v.open(commonPath + "_00Violations.txt");
    for(const CVariant* k : violationVarsFrom00Check)
        outputFileCC1v << k->m_nOriginalPos << std::endl;
    outputFileCC1v.close();

    std::ofstream outputFileCC11;
    outputFileCC11.open(commonPath + "_00CompliantsGT.txt");
    for(const CVariant* k : compliantVarsFrom00CheckGT)
        outputFileCC11 << k->m_nOriginalPos << std::endl;
    outputFileCC11.close();

    std::ofstream outputFileCC11v;
    outputFileCC11v.open(commonPath + "_00ViolationsGT.txt");
    for(const CVariant* k : violationVarsFrom00CheckGT)
        outputFileCC11v << k->m_nOriginalPos << std::endl;
    outputFileCC11v.close();
     
    */
}


void CMendelianAnalyzer::ReportChildChromosomeData(int a_nChromosomeId, std::vector<const CVariant*>& a_rCompliants, std::vector<const CVariant*>& a_rViolations)
{
    int compliantSNPcount = 0;
    int compliantINDELcount = 0;
    int violationSNPcount = 0;
    int violationINDELcount = 0;

    for(const CVariant* pVar : a_rCompliants)
    {
        if(pVar->GetVariantType() == eSNP)
            compliantSNPcount++;
        else
            compliantINDELcount++;
    }
    
    for(const CVariant* pVar : a_rViolations)
    {
        if(pVar->GetVariantType() == eSNP)
            violationSNPcount++;
        else
            violationINDELcount++;
    }
    
    m_resultLog.LogShortReport(a_nChromosomeId, compliantSNPcount, violationSNPcount, compliantINDELcount, violationINDELcount);
}


void CMendelianAnalyzer::PrintHelp() const
{
    std::cout << "==== SBG VCF COMPARISON TOOL VERSION 1.0 (Beta) ==== " << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgenomics.com || berke.toptas@sbgdinc.com)" << std::endl;
    std::cout << "Please notify me if program fails or return unexpected results" << std::endl;
    std::cout << "COPYRIGHT (C) 2016 SEVEN BRIDGES GENOMICS." << std::endl;
    std::cout << "COPYRIGHT (C) 2017 SBGD INC" << std::endl;
    std::cout << std::endl;
    std::cout << " --- MENDELIAN PARAMETERS --- " << std::endl;
    std::cout << "-father <father_vcf_path>    [Required.Add father VCF file.]" << std::endl;
    std::cout << "-mother <mother_vcf_path>    [Required.Add mother VCF file.]" << std::endl;
    std::cout << "-child <child_vcf_path>      [Required.Add child VCF file.]" << std::endl;
    std::cout << "-ref <reference_fasta_path>  [Required.Add reference FASTA file]" << std::endl;
    std::cout << "-outDir <output_directory>   [Required.Add output directory]" << std::endl;
    std::cout << "-no-call <no_call_mode>      [Optional. Decides what to do with no call variants. There are 3 modes:" << std::endl;
    std::cout << "\t" << "implicit : mark boths implicit and explicit no call variant as NoCall" << std::endl;
    std::cout << "\t" << "explicit : mark explicit no call variants only as NoCall. Implicit no call variants will be treated as 0/0" << std::endl;
    std::cout << "\t" << "none : [DEFAULT] Treat all of no call variants as 0/0" << std::endl;
    std::cout << "-filter <filter_name>        [Optional.Filter variants based on filter column. Default value is PASS. Use 'none' to unfilter]" << std::endl;
    std::cout << "-SampleFather <sample_name>  [Optional.Read only the given sample in father VCF. Default is the first sample.]" << std::endl;
    std::cout << "-SampleMother <sample_name>  [Optional.Read only the given sample in mother VCF. Default is the first sample.]" << std::endl;
    std::cout << "-SampleChild <sample_name>   [Optional.Read only the given sample in child VCF. Default is the first sample.]" << std::endl;
    std::cout << "-threadcount                 [Optional.Specify the number of threads that program will use. Default value is 2]" << std::endl;
    std::cout << std::endl;
    std::cout << "Example Commands:" << std::endl;
    std::cout << "./vbt mendelian -mother mother.vcf -father father.vcf -child child.vcf -ref reference.fa -outDir SampleResultDir -filter none -no-call explicit" << std::endl;
    std::cout << "./vbt mendelian -mother trio.vcf -father trio.vcf -child trio.vcf -ref reference.fa -outDir SampleResultDir -filter PASS -sampleMother mother -sampleFather father -sampleChild child" << std::endl;
}




