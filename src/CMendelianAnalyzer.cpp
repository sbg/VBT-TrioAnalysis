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

#include <fstream>


//Compare variants according to id for sort operation
bool variantCompare(const CVariant* v1, const CVariant* v2) {return v1->m_nId < v2->m_nId;}

//Checks if the given two range is overlapping
bool isOverlap(int left1, int right1, int left2, int right2)
{
    return std::min(right1, right2) - std::max(left1, left2) >= 0;
}



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
        int threadCount = AssignJobsToThreads(CHROMOSOME_COUNT);
        for(int k = 0; k < threadCount; k++)
            m_pThreadPool[k].join();
    }
    else
    {
        int threadCount = AssignJobsToThreads(MAC_THREAD_COUNT);
        for(int k = 0; k < threadCount; k++)
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

    for(int k = 0; k < pPath->m_aSyncPointList.size(); k++)
    {
        CSyncPoint sPoint;
        sPoint.m_nStartPosition = k > 0 ? pPath->m_aSyncPointList[k-1] : 0;
        sPoint.m_nEndPosition = pPath->m_aSyncPointList[k];
        sPoint.m_nIndex = k;
        
        while(baseIncludedItr < pBaseIncluded.size() && pBaseIncluded[baseIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_baseVariantsIncluded.push_back(pBaseIncluded[baseIncludedItr]);
            baseIncludedItr++;
        }

        while(calledIncludedItr < pCalledIncluded.size() && pCalledIncluded[calledIncludedItr]->GetStartPos() <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_calledVariantsIncluded.push_back(pCalledIncluded[calledIncludedItr]);
            calledIncludedItr++;
        }
        
        while(baseExcludedItr < pBaseExcluded.size() && pBaseExcluded[baseExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_baseVariantsExcluded.push_back(pBaseExcluded[baseExcludedItr]);
            baseExcludedItr++;
        }
        
        while(calledExcludedItr < pCalledExcluded.size() && pCalledExcluded[calledExcludedItr]->m_nStartPos <= pPath->m_aSyncPointList[k])
        {
            sPoint.m_calledVariantsExcluded.push_back(pCalledExcluded[calledExcludedItr]);
            calledExcludedItr++;
        }
        
        a_rSyncPointList.push_back(sPoint);
    }
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
    for(int k = 0; k < a_rSyncPointList.size() && varlistItr < a_rVarList.size(); k++)
    {
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < a_rVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_rVarList[varlistItr]->m_nStartPos)
        {
            a_rViolationList.push_back(a_rVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(int m = 0; m < a_rSyncPointList[k].m_calledVariantsIncluded.size(); m++)
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
            
            for(int m = 0; m < a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(eSNP == pVar->GetVariantType())
                {
                    if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                    {
                        if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                        {
                            bIsCompliant = false;
                            break;
                        }
                    }
                }
                
                else if(eINDEL == pVar->GetVariantType())
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
                for(int aa = 0; aa < tmpVarList.size(); aa++)
                    a_rCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(int aa = 0; aa < tmpVarList.size(); aa++)
                    a_rViolationList.push_back(tmpVarList[aa]);
            }
        }
    }
    
}


void CMendelianAnalyzer::CheckFor0Path(int a_nChrId,
                                       bool a_bIsFatherChild,
                                       std::vector<const CVariant *> &a_pVarList,
                                       std::vector<const CVariant *> &a_pViolantList,
                                       std::vector<const CVariant *> &a_pCompliantList)
{
    
    //Get sync point list
    std::vector<CSyncPoint> a_rSyncPointList;
    GetSyncPointList(a_nChrId, a_bIsFatherChild, a_rSyncPointList);
    
    
    int varlistItr = 0;
    for(int k = 0; k < a_rSyncPointList.size() && varlistItr < a_pVarList.size(); k++)
    {
        //If we check that syncpoint
        bool bDoCheck = false;
        std::vector<const CVariant*> tmpVarList;
        
        //Exclude variant if we somehow skip the syncpoint intervals
        while(varlistItr < a_pVarList.size() -1 && a_rSyncPointList[k].m_nStartPosition > a_pVarList[varlistItr]->m_nStartPos)
        {
            a_pViolantList.push_back(a_pVarList[varlistItr]);
            varlistItr++;
        }
        
        //Check if the sync interval contains 0/x child variants
        for(int m = 0; m < a_rSyncPointList[k].m_calledVariantsExcluded.size(); m++)
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
            
            for(int m = 0; m < a_rSyncPointList[k].m_baseVariantsExcluded.size(); m++)
            {
                const CVariant* pVar = a_rSyncPointList[k].m_baseVariantsExcluded[m];
                
                if(eSNP == pVar->GetVariantType())
                {
                    if(isOverlap(pVar->GetStart(), pVar->GetEnd(), tmpVarList[0]->GetStart(), tmpVarList[0]->GetEnd()))
                    {
                        if(pVar->m_genotype[0] != 0 && pVar->m_genotype[1] != 0)
                        {
                            bIsCompliant = false;
                            break;
                        }
                    }
                }
                
                else if(eINDEL == pVar->GetVariantType())
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
                for(int aa = 0; aa < tmpVarList.size(); aa++)
                    a_pCompliantList.push_back(tmpVarList[aa]);
            }
            
            else
            {
                for(int aa = 0; aa < tmpVarList.size(); aa++)
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
        CheckFor0Path(a_nChrId, true,  a_rOvarList, fatherViolants, fatherCompliants);
        CheckFor0Path(a_nChrId, false, a_rOvarList, motherViolants, motherCompliants);
    }
    else
    {
        CheckFor0PathFor00(a_nChrId, true,  a_rOvarList, fatherViolants, fatherCompliants);
        CheckFor0PathFor00(a_nChrId, false, a_rOvarList, motherViolants, motherCompliants);
    }
    
    std::cout << std::endl;
    std::cout << "00 Check Father Compliant Found:" << fatherCompliants.size() << std::endl;
    std::cout << "00 Check Mother Compliant Found:" << motherCompliants.size() << std::endl;
    std::cout << std::endl;
    
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
    for(int k = 0; k < chromosomeListToProcess.size(); k++)
    {
        chromosomeLists[threadPoolIt].push_back(chromosomeListToProcess[k]);
        threadPoolIt = (threadPoolIt+1) % exactThreadCount;
    }
    
    //Assign divided task to the threads
    for(int k = 0; k < exactThreadCount; k++)
    {
        m_pThreadPool[k] = std::thread(&CMendelianAnalyzer::ProcessChromosome, this, chromosomeLists[k]);
    }
    
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
        
    }

}


void CMendelianAnalyzer::MergeFunc(int a_nChromosomeId)
{
    std::vector<const CVariant*> MendelianViolationVars;
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
                        MendelianViolationVars.push_back(&varMC->GetVariant());
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
    
    
    //Process remaining vars in FatherChild
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
    
    //Process remaining vars in MotherChild
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
    
    CheckFor0Path(a_nChromosomeId, true, check0atFatherSide, violationVarsFrom0CheckFather, compliantVarsFrom0CheckFather);
    CheckFor0Path(a_nChromosomeId, false, check0atMotherSide, violationVarsFrom0CheckMother, compliantVarsFrom0CheckMother);
    
    std::vector<const CVariant*> childVariants = m_provider.GetVariantList(eCHILD, a_nChromosomeId);
    
    std::vector<int>childProcessedArray(childVariants.size());
    for(int elem : childProcessedArray)
        elem = 0;
    
    //Mark mendelian compliant vars as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < MendelianCompliantVars.size(); k++)
    {
        if(childVariants[k]->m_nId == MendelianCompliantVars[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
        
    }
    
    //Mark mendelian violation vars as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < MendelianViolationVars.size(); k++)
    {
        if(childVariants[k]->m_nId == MendelianViolationVars[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    //Mark father child only vars as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < fatherChildOnly.size(); k++)
    {
        if(childVariants[k]->m_nId == fatherChildOnly[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    //Mark mother child only vars as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < motherChildOnly.size(); k++)
    {
        if(childVariants[k]->m_nId == motherChildOnly[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }

    //Mark check 0 condition at father side as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < check0atFatherSide.size(); k++)
    {
        if(childVariants[k]->m_nId == check0atFatherSide[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }

    //Mark check 0 condition at mother side as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < check0atMotherSide.size(); k++)
    {
        if(childVariants[k]->m_nId == check0atMotherSide[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    //Mark child 00 vars as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < check00Child.size(); k++)
    {
        if(childVariants[k]->m_nId == check00Child[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }

    //Mark child 00 vars GT matched as processed
    for(int k = 0, m = 0; k < childVariants.size() && m < check00ChildGTMatched.size(); k++)
    {
        if(childVariants[k]->m_nId == check00ChildGTMatched[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }

    
    
    int noMatch = 0;
    for(int childItr = 0; childItr < childProcessedArray.size(); childItr++)
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
                noMatch++;
            }
        }
    }
    
    std::vector<const CVariant*> compliantVarsFrom00Check;
    std::vector<const CVariant*> violationVarsFrom00Check;

    std::vector<const CVariant*> compliantVarsFrom00CheckGT;
    std::vector<const CVariant*> violationVarsFrom00CheckGT;

    CheckFor00Child(a_nChromosomeId, check00Child, violationVarsFrom00Check, compliantVarsFrom00Check, false);
    CheckFor00Child(a_nChromosomeId, check00ChildGTMatched, violationVarsFrom00CheckGT, compliantVarsFrom00CheckGT, true);
    
    std::vector<const CVariant*> compliants;
    std::vector<const CVariant*> violations;
    
    compliants.insert(std::end(compliants), std::begin(MendelianCompliantVars), std::end(MendelianCompliantVars));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckFather), std::end(compliantVarsFrom0CheckFather));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom0CheckMother), std::end(compliantVarsFrom0CheckMother));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom00Check), std::end(compliantVarsFrom00Check));
    compliants.insert(std::end(compliants), std::begin(compliantVarsFrom00CheckGT), std::end(compliantVarsFrom00CheckGT));
    std::sort(compliants.begin(), compliants.end(), variantCompare);
    
    violations.insert(std::end(violations), std::begin(MendelianViolationVars), std::end(MendelianViolationVars));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckFather), std::end(violationVarsFrom0CheckFather));
    violations.insert(std::end(violations), std::begin(violationVarsFrom0CheckMother), std::end(violationVarsFrom0CheckMother));
    violations.insert(std::end(violations), std::begin(fatherChildOnly), std::end(fatherChildOnly));
    violations.insert(std::end(violations), std::begin(motherChildOnly), std::end(motherChildOnly));
    violations.insert(std::end(violations), std::begin(violationVarsFrom00Check), std::end(violationVarsFrom00Check));
    violations.insert(std::end(violations), std::begin(childUniqueList), std::end(childUniqueList));
    violations.insert(std::end(violations), std::begin(violationVarsFrom00CheckGT), std::end(violationVarsFrom00CheckGT));
    std::sort(violations.begin(), violations.end(), variantCompare);
    
    std::cout << "===================== STATISTICS " << a_nChromosomeId + 1 << "===================" << std::endl;
    std::cout << std::endl;
    std::cout << "Violation Vars Same Allele Match:" << MendelianViolationVars.size() << std::endl;
    std::cout << "Compliant Variants:" << MendelianCompliantVars.size() << std::endl;
    std::cout << std::endl;
    std::cout << "Violation Vars 0 Check Mother:" << violationVarsFrom0CheckMother.size() << std::endl;
    std::cout << "Compliant Vars 0 Check Mother:" << compliantVarsFrom0CheckMother.size() << std::endl;
    std::cout << "0 Check Total: " << check0atMotherSide.size() << " Sum of two:" << violationVarsFrom0CheckMother.size() + compliantVarsFrom0CheckMother.size() << std::endl;
    std::cout << std::endl;
    std::cout << "Violation Vars 0 Check Father:" << violationVarsFrom0CheckFather.size() << std::endl;
    std::cout << "Compliant Vars 0 Check Father:" << compliantVarsFrom0CheckFather.size() << std::endl;
    std::cout << "0 Check Total: " << check0atFatherSide.size() << " Sum of two:" << violationVarsFrom0CheckFather.size() + compliantVarsFrom0CheckFather.size() << std::endl;
    std::cout << std::endl;
    std::cout << "Violation Vars Father Only:"  <<  fatherChildOnly.size() << std::endl;
    std::cout << "Violation Vars Mother Only:"  <<  motherChildOnly.size() << std::endl;
    std::cout << "Violation Vars Child Unique:" <<  childUniqueList.size() << std::endl;
    std::cout << std::endl;
    std::cout << "Violation Vars 00 Check:" << violationVarsFrom00Check.size() << std::endl;
    std::cout << "Compliant Vars 00 Check :" << compliantVarsFrom00Check.size() << std::endl;
    std::cout << "00 Check Total:" << check00Child.size() << " Sum of two:" << violationVarsFrom00Check.size() + compliantVarsFrom00Check.size() << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "Total Compliants:" << compliants.size() << std::endl;
    std::cout << "Total Violations:" << violations.size() << std::endl;
    std::cout << "Child Var Size:" << childVariants.size() << std::endl;
    
    std::ofstream compliantsAll;
    std::string commonPath = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/chr" + std::to_string(a_nChromosomeId + 1);
    
    compliantsAll.open(commonPath + "_CompliantsALL.txt");
    
    for(const CVariant* k : compliants)
        compliantsAll << k->m_nOriginalPos << std::endl;
    compliantsAll.close();

    std::ofstream violationsAll;
    violationsAll.open(commonPath + "_ViolationsALL.txt");
    
    for(const CVariant* k : violations)
        violationsAll << k->m_nOriginalPos << std::endl;
    violationsAll.close();
    
/*
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
        outputFileF.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21Father0sCompliant.txt");
        for(const CVariant* k : compliantVarsFrom0CheckFather)
            outputFileF << k->m_nOriginalPos << std::endl;
        outputFileF.close();

        std::ofstream outputFileFvio;
        outputFileFvio.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21Father0sViolation.txt");
        for(const CVariant* k : violationVarsFrom0CheckFather)
            outputFileFvio << k->m_nOriginalPos << std::endl;
        outputFileFvio.close();
        

        std::ofstream outputFileM;
        outputFileM.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21Mother0sCompliant.txt");
        for(const CVariant* k : compliantVarsFrom0CheckMother)
            outputFileM << k->m_nOriginalPos << std::endl;
        outputFileM.close();
        
        
        std::ofstream outputFileMvio;
        outputFileMvio.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Chr21Mother0sViolation.txt");
        for(const CVariant* k : violationVarsFrom0CheckMother)
            outputFileMvio << k->m_nOriginalPos << std::endl;
        outputFileMvio.close();

        std::ofstream outputFileCC1;
        outputFileCC1.open(commonPath + "_00Compliants.txt");
        for(const CVariant* k : compliantVarsFrom00Check)
            outputFileCC1 << k->m_nOriginalPos << std::endl;
        outputFileCC1.close();
    
        std::ofstream outputFileCC11;
        outputFileCC11.open(commonPath + "_00CompliantsGT.txt");
        for(const CVariant* k : compliantVarsFrom00CheckGT)
            outputFileCC11 << k->m_nOriginalPos << std::endl;
        outputFileCC11.close();

        std::ofstream outputFileCC2;
        outputFileCC2.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/Child00CheckViolations.txt");
        for(const CVariant* k : violationVarsFrom00Check)
            outputFileCC2 << k->m_nOriginalPos << std::endl;
        outputFileCC2.close();

 
        std::ofstream outputFileChildUnique;
        outputFileChildUnique.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/childUniqueWithout00.txt");
        for(int k : childUniqueVarIndexList)
            outputFileChildUnique << childVariants[k]->m_nOriginalPos << std::endl;
        outputFileChildUnique.close();
*/
    
}






