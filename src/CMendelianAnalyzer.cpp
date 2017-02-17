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

#include <fstream>


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
    
    std::vector<const COrientedVariant*> motherChildHeterozygous;
    std::vector<const COrientedVariant*> motherChildHomozygous;
    std::vector<const COrientedVariant*> fatherChildHeterozygous;
    std::vector<const COrientedVariant*> fatherChildHomozygous;

    std::vector<int> childUniqueVarIndexList;

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
            if(varMC->GetVariant().m_bIsHeterozygous)
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
                fatherChildHeterozygous.push_back(varFC);
            else
                fatherChildHomozygous.push_back(varFC);
            
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
                motherChildHeterozygous.push_back(varMC);
            else
                motherChildHomozygous.push_back(varMC);
            
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
            fatherChildHeterozygous.push_back(varFC);
        else
            fatherChildHomozygous.push_back(varFC);
        
        if(!FatherChildVariants.hasNext())
            break;
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
            motherChildHeterozygous.push_back(varMC);
        else
            motherChildHomozygous.push_back(varMC);
        
        if(!MotherChildVariants.hasNext())
            break;
    }
    
    
    std::vector<const CVariant*> childVariants = m_provider.GetVariantList(eCHILD, a_nChromosomeId);
    
    std::vector<int>childProcessedArray(childVariants.size());
    for(int elem : childProcessedArray)
        elem = 0;
    
    for(int k = 0, m = 0; k < childVariants.size() && m < MendelianCompliantVars.size(); k++)
    {
        if(childVariants[k]->m_nId == MendelianCompliantVars[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
        
    }
    
    for(int k = 0, m = 0; k < childVariants.size() && m < MendelianViolationVars.size(); k++)
    {
        if(childVariants[k]->m_nId == MendelianViolationVars[m]->m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    for(int k = 0, m = 0; k < childVariants.size() && m < motherChildHeterozygous.size(); k++)
    {
        if(childVariants[k]->m_nId == motherChildHeterozygous[m]->GetVariant().m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }

    for(int k = 0, m = 0; k < childVariants.size() && m < motherChildHomozygous.size(); k++)
    {
        if(childVariants[k]->m_nId == motherChildHomozygous[m]->GetVariant().m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    for(int k = 0, m = 0; k < childVariants.size() && m < fatherChildHeterozygous.size(); k++)
    {
        if(childVariants[k]->m_nId == fatherChildHeterozygous[m]->GetVariant().m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    for(int k = 0, m = 0; k < childVariants.size() && m < fatherChildHomozygous.size(); k++)
    {
        if(childVariants[k]->m_nId == fatherChildHomozygous[m]->GetVariant().m_nId)
        {
            childProcessedArray[k]++;
            m++;
        }
    }
    
    int noMatch = 0;
    for(int elem : childProcessedArray)
        if(elem == 0)
            noMatch++;
    
   /*
    if(a_nChromosomeId == 0)
    {
    std::ofstream outputFile;
    outputFile.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/ORIGINALChr1MendelCompliants.txt");
    
    for(const CVariant* k : MendelianCompliantVars)
        outputFile << k->m_nOriginalPos << std::endl;
    
    outputFile.close();
    }
    */
    
    std::cout << "==="  << a_nChromosomeId + 1 << "===" << std::endl;
    std::cout << "MendelianCompliant:" << MendelianCompliantVars.size() << std::endl;
    std::cout << "MendelianViolation:" << MendelianViolationVars.size() << std::endl;
    std::cout << "Heterozygous Mother Side Only:" << motherChildHeterozygous.size() << std::endl;
    std::cout << "Homozygous Mother Side Only:" << motherChildHomozygous.size() << std::endl;
    std::cout << "Heterozygous Father Side Only:" << fatherChildHeterozygous.size() << std::endl;
    std::cout << "Homozygous Father Side Only:" << fatherChildHomozygous.size() << std::endl;
    std::cout << "ChildUnique:" << noMatch << std::endl;
    std::cout << "MendelianViolation Total:" << MendelianViolationVars.size() + motherChildHomozygous.size() + motherChildHeterozygous.size() + fatherChildHomozygous.size() + fatherChildHeterozygous.size()
    + noMatch << std::endl;
    std::cout << "TOTAL:" << MendelianViolationVars.size() + MendelianCompliantVars.size() + motherChildHomozygous.size() + motherChildHeterozygous.size() + fatherChildHomozygous.size() + fatherChildHeterozygous.size()
    + noMatch << std::endl;
    std::cout << "Child Var Original:" << childVariants.size() << std::endl;
    
}






