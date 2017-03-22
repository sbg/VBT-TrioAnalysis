//
//  test.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 1/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include <stdio.h>
#include <fstream>
#include <vector>


#include "CVcfReader.h"
#include <iostream>
#include "CVcfAnalyzer.h"
#include "COrientedVariant.h"

void testGa4ghOutput()
{
    CVcfReader vcfeval;
    CVcfReader sbg;
    
    CVariant variantVcfEval;
    CVariant variantSbg;
    
    SConfig config;
    
    int id = 0;
    
    int unmatchCnt = 0;
    int matchCnt = 0;
    
    vcfeval.Open("/Users/c1ms21p6h3qk/Desktop/RTG/rtg-core/Inputs/rtgoutput/output.vcf.gz");
    sbg.Open("/Users/c1ms21p6h3qk/Desktop/VcfCompOutput/ga4ghOutput.vcf");
    
    vcfeval.SelectSample("TRUTH");
    sbg.SelectSample("TRUTH");
    
    
    bool valExistsSBG = sbg.GetNextRecord(&variantSbg, id, config);
    bool valExistsVCFEVAL = vcfeval.GetNextRecord(&variantVcfEval, id, config);
    
    while(valExistsSBG && valExistsVCFEVAL)
    {
        int vcfEvalStart = variantVcfEval.m_nOriginalPos;
        int sbgStart = variantSbg.m_nOriginalPos;
        
        if(vcfEvalStart == sbgStart &&
           0 == variantVcfEval.m_allelesStr.compare(variantSbg.m_allelesStr))
        {
            id++;
            matchCnt++;
            valExistsSBG = sbg.GetNextRecord(&variantSbg, id, config);
            valExistsVCFEVAL = vcfeval.GetNextRecord(&variantVcfEval, id, config);
            continue;
        }
        
        else if (vcfEvalStart > sbgStart)
        {
            unmatchCnt++;
            valExistsSBG = sbg.GetNextRecord(&variantSbg, id, config);
            continue;
        }
        
        else if (sbgStart > vcfEvalStart)
        {
            unmatchCnt++;
            valExistsVCFEVAL = vcfeval.GetNextRecord(&variantVcfEval, id, config);
            continue;
        }
        
        else
        {
            std::cout << id << "\t\t" << variantSbg.m_nStartPos << "\t\t" << variantSbg.m_allelesStr << "\t\t\t" << variantVcfEval.m_allelesStr << std::endl;
            unmatchCnt++;
            id++;
            valExistsSBG = sbg.GetNextRecord(&variantSbg, id, config);
            valExistsVCFEVAL = vcfeval.GetNextRecord(&variantVcfEval, id, config);
            
            continue;
        }
    }
    
    std::cout << "Match:" << matchCnt << std::endl;
    std::cout << "Unmatch:" << unmatchCnt << std::endl;

}

/*
void testRefOverlapOutput(int argc, char** argv)
{
    CVcfAnalyzer analyzer1;
    analyzer1.Run(argc, argv);
    
    strcpy(argv[4], std::string("/Users/c1ms21p6h3qk/Desktop/VcfCompOpt2").c_str());
    strcpy(argv[argc-1], std::string("-ref-overlap").c_str());
    
    CVcfAnalyzer analyzer2;
    analyzer2.Run(argc, argv);
    
    std::vector<const COrientedVariant*> set1IncludedBase = analyzer1.m_aBestPaths[20].m_baseSemiPath.GetIncludedVariants();
    std::vector<const COrientedVariant*> set2IncludedBase = analyzer2.m_aBestPaths[20].m_baseSemiPath.GetIncludedVariants();
    
    std::vector<const COrientedVariant*> set1IncludedCalled = analyzer1.m_aBestPaths[20].m_calledSemiPath.GetIncludedVariants();
    std::vector<const COrientedVariant*> set2IncludedCalled = analyzer2.m_aBestPaths[20].m_calledSemiPath.GetIncludedVariants();

    std::vector<const COrientedVariant*>::iterator set1IncludedBaseIterator = set1IncludedBase.begin();
    std::vector<const COrientedVariant*>::iterator set2IncludedBaseIterator = set2IncludedBase.begin();
    std::vector<const COrientedVariant*>::iterator set1IncludedCalledIterator = set1IncludedCalled.begin();
    std::vector<const COrientedVariant*>::iterator set2IncludedCalledIterator = set2IncludedCalled.begin();
    
    
    int matchCount = 0;
    int unmatchCount = 0;
    
//    std::vector<const CVariant&> unmatchedVariants1;
//    std::vector<const CVariant&> unmatchedVariants2;
    
    //TEST INCLUDED BASE
    while(true)
    {
        if(set1IncludedBaseIterator == set1IncludedBase.end() || set2IncludedBaseIterator == set2IncludedBase.end())
            break;
        
        if((*set1IncludedBaseIterator)->GetVariant().GetOriginalPos() == (*set2IncludedBaseIterator)->GetVariant().GetOriginalPos()
           &&
           (*set1IncludedBaseIterator)->GetVariant().m_allelesStr == (*set2IncludedBaseIterator)->GetVariant().m_allelesStr)
        {
            matchCount++;
            set1IncludedBaseIterator++;
            set2IncludedBaseIterator++;
        }
        else if((*set1IncludedBaseIterator)->GetVariant().GetOriginalPos() > (*set2IncludedBaseIterator)->GetVariant().GetOriginalPos())
        {
            unmatchCount++;
            std::cout << "SET2:" << (*set2IncludedBaseIterator)->GetVariant().ToString() << std::endl;
            std::cout << "===" << std::endl;

            set2IncludedBaseIterator++;
        }
        else if((*set1IncludedBaseIterator)->GetVariant().GetOriginalPos() < (*set2IncludedBaseIterator)->GetVariant().GetOriginalPos())
        {
            unmatchCount++;
            std::cout << "SET1:" << (*set1IncludedBaseIterator)->GetVariant().ToString() << std::endl;
            std::cout << "===" << std::endl;

            set1IncludedBaseIterator++;
        }
        else
        {
            unmatchCount++;
            set1IncludedBaseIterator++;
            set2IncludedBaseIterator++;
        }
        
    }
    
    std::cout << "Match Count:" << matchCount << std::endl;
    std::cout << "Unmatch Count:" << unmatchCount << std::endl;
    
    matchCount = 0;
    unmatchCount = 0;
    
    //TEST INCLUDED CALLED
    while(true)
    {
        if(set1IncludedCalledIterator == set1IncludedCalled.end() || set2IncludedCalledIterator == set2IncludedCalled.end())
            break;
        
        if((*set1IncludedCalledIterator)->GetVariant().GetOriginalPos() == (*set2IncludedCalledIterator)->GetVariant().GetOriginalPos()
           &&
           (*set1IncludedCalledIterator)->GetVariant().m_allelesStr == (*set2IncludedCalledIterator)->GetVariant().m_allelesStr)
        {
            matchCount++;
            set1IncludedCalledIterator++;
            set2IncludedCalledIterator++;
        }
        else if((*set1IncludedCalledIterator)->GetVariant().GetOriginalPos() > (*set2IncludedCalledIterator)->GetVariant().GetOriginalPos())
        {
            unmatchCount++;
            std::cout << "SET2:" << (*set2IncludedCalledIterator)->GetVariant().ToString() << std::endl;
            std::cout << "===" << std::endl;

            set2IncludedCalledIterator++;
        }
        else if((*set1IncludedCalledIterator)->GetVariant().GetOriginalPos() < (*set2IncludedCalledIterator)->GetVariant().GetOriginalPos())
        {
            unmatchCount++;
            std::cout << "SET1:" << (*set1IncludedCalledIterator)->GetVariant().ToString() << std::endl;
            std::cout << "===" << std::endl;

            set1IncludedCalledIterator++;
        }
        else
        {
            std::cout << "SET1:" << (*set1IncludedCalledIterator)->GetVariant().ToString() << std::endl;
            std::cout << "SET2:" << (*set2IncludedCalledIterator)->GetVariant().ToString() << std::endl;
            std::cout << "===" << std::endl;
            
            unmatchCount++;
            set1IncludedCalledIterator++;
            set2IncludedCalledIterator++;
        }
        
    }
    
    std::cout << "Match Count:" << matchCount << std::endl;
    std::cout << "Unmatch Count:" << unmatchCount << std::endl;
    
}
*/


void compareTrimmedVariables()
{
    std::vector<std::string> sbgLines;
    std::vector<std::string> vcfevalLines;
    
    std::ifstream sbgFILE("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::ifstream vcfevalFILE("/Users/c1ms21p6h3qk/Desktop/outECLIPSE.txt");
    
    std::string line;
    
    while(std::getline(sbgFILE, line))
        sbgLines.push_back(line);
    
    while(std::getline(vcfevalFILE, line))
        vcfevalLines.push_back(line);

//    std::sort(sbgLines.begin(), sbgLines.end());
//    std::sort(vcfevalLines.begin(), vcfevalLines.end());
    
    int unmatchCnt = 0;
    int sbgExtra = 0;
    int vcfevalExtra = 0;
    std::vector<std::string> unmatchedIds;
    std::vector<std::string> sbgExtraIds;
    std::vector<std::string> vcfevalExtraIds;
    
    int k=0;
    int p=0;
    
    while(k != sbgLines.size() && p!= vcfevalLines.size())
    {
        if(sbgLines[k].compare(vcfevalLines[p]) == 0)
        {
            k++;
            p++;
        }
        else if(sbgLines[k].compare(vcfevalLines[p]) < 0)
        {
            unmatchCnt++;
            sbgExtra++;
            unmatchedIds.push_back(sbgLines[k] + "    SBG");
            sbgExtraIds.push_back(sbgLines[k]);
            k++;
        }
        
        else
        {
            unmatchCnt++;
            unmatchedIds.push_back(vcfevalLines[p] + "    VCFEVAL");
            vcfevalExtra++;
            vcfevalExtraIds.push_back(vcfevalLines[p]);
            p++;
        }
    }
    
    std::cout << "UnmatchCnt: " << unmatchCnt << std::endl;
    
    sbgFILE.close();
    vcfevalFILE.close();
}


void testTripletReader2()
{
    const int TEST_CHR = 1;
    
    CVcfReader father;
    CVcfReader mother;
    CVcfReader child;
    
    
    father.Open("/Users/c1ms21p6h3qk/Desktop/MendelianInput/Father_HG00096.vcf");
    mother.Open("/Users/c1ms21p6h3qk/Desktop/MendelianInput/Mother_HG00171.vcf");
    child.Open("/Users/c1ms21p6h3qk/Desktop/MendelianInput/Child_HG00096xHG00171.vcf");

    std::vector<CVariant> motherVariants;
    std::vector<CVariant> fatherVariants;
    std::vector<CVariant> childVariants;

    
    CVariant variant;
    
    int id = 0;
    SConfig config;
    
    while(father.GetNextRecord(&variant, id, config))
    {
        if(variant.m_nChrId == TEST_CHR)
        {
            fatherVariants.push_back(variant);
            id++;
        }
    }

    id = 0;
    while(mother.GetNextRecord(&variant, id, config))
    {
        if(variant.m_nChrId == TEST_CHR)
        {
            motherVariants.push_back(variant);
            id++;
        }
    }
    
    id = 0;
    while(child.GetNextRecord(&variant, id, config))
    {
        if(variant.m_nChrId == TEST_CHR)
        {
            childVariants.push_back(variant);
            id++;
        }
    }
    
    std::cout << "Mother variants Size:" << motherVariants.size() << std::endl;
    std::cout << "Father variants Size:" << fatherVariants.size() << std::endl;
    std::cout << "Child variants Size:" << childVariants.size() << std::endl;
    
    
    std::vector<CVariant*> motherChildCommon;
    std::vector<EVariantMatch> motherChildMatchTypes;
    std::vector<int> motherChildMatchAllele;
    std::vector<CVariant*> motherUnique;
    std::vector<CVariant*> childUniqueMC;

    std::vector<CVariant*> fatherChildCommon;
    std::vector<EVariantMatch> fatherChildMatchTypes;
    std::vector<int> fatherChildMatchAllele;
    std::vector<CVariant*> fatherUnique;
    std::vector<CVariant*> childUniqueFC;

    //MOTHER-CHILD COMPARISON
    int k=0, m=0;
    for(; k < motherVariants.size() && m < childVariants.size();)
    {
        if(motherVariants[k].m_nOriginalPos == 5021073)
        {
            int aas = 0;
            aas++;
        }
        
        if(motherVariants[k].m_nOriginalPos == childVariants[m].m_nOriginalPos)
        {
            int t1,t2,t3,t4;
            
            t1 = motherVariants[k].m_alleles[0].m_sequence == childVariants[m].m_alleles[0].m_sequence ? 1 : 0;
            t2 = motherVariants[k].m_alleles[1].m_sequence == childVariants[m].m_alleles[1].m_sequence ? 1 : 0;
            
            t3 = motherVariants[k].m_alleles[0].m_sequence == childVariants[m].m_alleles[1].m_sequence ? 1 : 0;
            t4 = motherVariants[k].m_alleles[1].m_sequence == childVariants[m].m_alleles[0].m_sequence ? 1 : 0;
            
            
            if(t1+t2 == 2 || t3+t4 == 2)
            {
                motherChildMatchTypes.push_back(eGENOTYPE_MATCH);
                motherChildMatchAllele.push_back(0);
            }
            else if(t1 + t2 == 1 || t3+t4 == 1)
            {
                motherChildMatchTypes.push_back(eALLELE_MATCH);
                motherChildMatchAllele.push_back(t1+t4 > 0 ? 0 : 1);
            }
            else
            {
                m++;
                k++;
                continue;
            }
        
            motherChildCommon.push_back(&motherVariants[k]);
            m++;
            k++;
        
        }
        else if(motherVariants[k].m_nOriginalPos > childVariants[m].m_nOriginalPos)
        {
            childUniqueMC.push_back(&childVariants[m]);
            m++;
        }
        else
        {
            motherUnique.push_back(&motherVariants[k]);
            k++;
        }
    }
    
    for(int x = k; x < motherVariants.size(); x++)
        motherUnique.push_back(&motherVariants[x]);
    for(int x = m; x < childVariants.size(); x++)
        childUniqueMC.push_back(&childVariants[x]);
    
        

    
    //FATHER-CHILD COMPARISON
    for(k=0, m=0; k < fatherVariants.size() && m < childVariants.size();)
    {
        if(fatherVariants[k].m_nOriginalPos == 5021073)
        {
            int aas = 0;
            aas++;
        }

        
        if(fatherVariants[k].m_nOriginalPos == childVariants[m].m_nOriginalPos)
        {
            int t1,t2,t3,t4;
            
            t1 = fatherVariants[k].m_alleles[0].m_sequence == childVariants[m].m_alleles[0].m_sequence ? 1 : 0;
            t2 = fatherVariants[k].m_alleles[1].m_sequence == childVariants[m].m_alleles[1].m_sequence ? 1 : 0;

            t3 = fatherVariants[k].m_alleles[0].m_sequence == childVariants[m].m_alleles[1].m_sequence ? 1 : 0;
            t4 = fatherVariants[k].m_alleles[1].m_sequence == childVariants[m].m_alleles[0].m_sequence ? 1 : 0;
 
            
            if(t1+t2 == 2 || t3+t4 == 2)
            {
                fatherChildMatchTypes.push_back(eGENOTYPE_MATCH);
                fatherChildMatchAllele.push_back(0);
            }
            else if(t1 + t2 == 1 || t3+t4 == 1)
            {
                fatherChildMatchTypes.push_back(eALLELE_MATCH);
                fatherChildMatchAllele.push_back(t1+t4 > 0 ? 0 : 1);
            }
            else
            {
                m++;
                k++;
                continue;
            }
            
            fatherChildCommon.push_back(&fatherVariants[k]);
            m++;
            k++;
        }
        
        
        else if(fatherVariants[k].m_nOriginalPos > childVariants[m].m_nOriginalPos)
        {
            childUniqueFC.push_back(&childVariants[m]);
            m++;
        }
        else
        {
            fatherUnique.push_back(&fatherVariants[k]);
            k++;
        }
    }

    for(int x = k; x < fatherVariants.size(); x++)
        fatherUnique.push_back(&fatherVariants[x]);
    for(int x = m; x < childVariants.size(); x++)
        childUniqueFC.push_back(&childVariants[x]);
    
    
    
    std::vector<int> mendelianViolations;
    std::vector<int> mendelianCompliants;
    
    
    for(k = 0, m = 0; k < motherChildCommon.size() && m < fatherChildCommon.size();)
    {
        if(motherChildCommon[k]->m_nOriginalPos == fatherChildCommon[m]->m_nOriginalPos)
        {
            if(motherChildCommon[k]->m_nOriginalPos == 5021073)
            {
                int aas = 0;
                aas++;
            }

            if(motherChildMatchAllele[k] != fatherChildMatchAllele[m])
            {
                mendelianCompliants.push_back(motherChildCommon[k]->m_nOriginalPos);
                k++;
                m++;
            }
            
            else if(motherChildMatchTypes[k] == eGENOTYPE_MATCH || fatherChildMatchTypes[m] == eGENOTYPE_MATCH)
            {
                mendelianCompliants.push_back(motherChildCommon[k]->m_nOriginalPos);
                k++;
                m++;
            }
            
            else
            {
                mendelianViolations.push_back(motherChildCommon[k]->m_nOriginalPos);
                k++;
                m++;
            }
        }
        else if(motherChildCommon[k]->m_nOriginalPos > fatherChildCommon[m]->m_nOriginalPos)
        {
            mendelianViolations.push_back(fatherChildCommon[m]->m_nOriginalPos);
            m++;
        }
        else
        {
            mendelianViolations.push_back(motherChildCommon[k]->m_nOriginalPos);
            k++;
        }
    }
    
    for(int x = k; x < motherChildCommon.size(); x++)
        mendelianViolations.push_back(motherChildCommon[x]->m_nOriginalPos);
    for(int x = m; x < fatherChildCommon.size(); x++)
        mendelianViolations.push_back(fatherChildCommon[x]->m_nOriginalPos);
    
    for(k = 0, m=0; k < childUniqueMC.size() && m <  childUniqueFC.size();)
    {
        if(childUniqueMC[k]->m_nOriginalPos == childUniqueFC[m]->m_nOriginalPos)
        {
            mendelianViolations.push_back(childUniqueMC[k]->m_nOriginalPos);
            k++;
            m++;
        }
        else if(childUniqueMC[k]->m_nOriginalPos > childUniqueFC[m]->m_nOriginalPos)
        {
            m++;
        }
        else
        {
            k++;
        }
    }
    
    std::cout << "Mother-Child Common:" << motherChildCommon.size() << std::endl;
    std::cout << "Mother Unique : " << motherUnique.size() << std::endl;
    std::cout << "Child Unique MC: "  << childUniqueMC.size() << std:: endl;
    
    std::cout << "Father-Child Common:" << fatherChildCommon.size() << std::endl;
    std::cout << "Father Unique : " << fatherUnique.size() << std::endl;
    std::cout << "Child Unique FC: "  << childUniqueFC.size() << std:: endl;
    
    std::cout << "Mendelian Compliant Vars:" << mendelianCompliants.size() << std::endl;
    std::cout << "Mendelian Violation Vars:" << mendelianViolations.size() << std::endl;
    
    std::ofstream outputFile;
    outputFile.open("/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TestChr1MendelCompliants.txt");
    
    for(int k : mendelianCompliants)
        outputFile << k << std::endl;
    
    outputFile.close();
    
    
    

}

struct SCompareResult
{
    int common;
    int unMatchOrig;
};

SCompareResult compare2List(const std::string& a_pathTest, const std::string& a_pathQuery)
{
    std::ifstream testfile;
    std::ifstream origfile;
    
    testfile.open(a_pathTest.c_str());
    origfile.open(a_pathQuery.c_str());
        
    std::vector<int> testIndexes;
    std::vector<int> origIndexes;
    
    std::string line;
    
    while(std::getline(testfile, line))
        testIndexes.push_back(atoi(line.c_str()));
    
    while(std::getline(origfile, line))
        origIndexes.push_back(atoi(line.c_str()));

    
    std::vector<int> commonVariants;
    std::vector<int> unmatchTest;
    std::vector<int> unmatchOrig;
    
    std::cout << "Test File Size: " << testIndexes.size() << std::endl;
    std::cout << "Orig File Size: " << origIndexes.size() << std::endl;
    
    int k,m;
    
    for(k = 0, m=0; k < origIndexes.size() && m <  testIndexes.size();)
    {
        if(origIndexes[k] == testIndexes[m])
        {
            commonVariants.push_back(origIndexes[k]);
            k++;
            m++;
        }
        else if(origIndexes[k] > testIndexes[m])
        {
            unmatchTest.push_back(testIndexes[m]);
            m++;
        }
        else
        {
            unmatchOrig.push_back(origIndexes[k]);
            k++;
        }
    }
    
    for(int x = k; x < origIndexes.size(); x++)
        unmatchOrig.push_back(origIndexes[x]);
    for(int x = m; x < testIndexes.size(); x++)
        unmatchTest.push_back(testIndexes[x]);
    
    std::cout << "Common Cnt:"   << commonVariants.size() << std::endl;
    std::cout << "Unmatch Test:" << unmatchTest.size()  << std::endl;
    std::cout << "Unmatch Orig:" << unmatchOrig.size()  << std::endl;
    
    std::cout << std::endl << std::endl;
    
    SCompareResult res;
    res.common = (int)commonVariants.size();
    res.unMatchOrig = (int)unmatchOrig.size();
    
    //for(int k : commonVariants)
    //    std::cout << k+1 << std::endl;
    
    std::cout << "Unmatch Test Variants:" << std::endl;
    for(int k : unmatchTest)
         std::cout << k+1 << std::endl;
    
    //std::cout << std::endl << std::endl;
    
    std::cout << "Unmatch Query Variants" << std::endl;
    for(int k : unmatchOrig)
        std::cout << k << std::endl;

    
    return res;
}


void UnitTestTrioComparison(int a_nChrNumber, bool a_bIsFilterOverlap, bool a_bIsFilter00)
{
    std::string vcfFileName;
    
    std::string queryCompliantSetFileName;
    std::string truthCompliantSetFileName;
    
    std::string queryViolationSetFileName;
    std::string truthViolationSetFileName;
    
    
    if(a_bIsFilter00 && a_bIsFilterOverlap)
    {
        vcfFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered_NoChild00.vcf";
        truthCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETFilteredOverlap00/Chr" + std::to_string(a_nChrNumber) + "_CompliantsTRUTH.txt";
        truthViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETFilteredOverlap00/Chr" + std::to_string(a_nChrNumber) + "_ViolationsTRUTH.txt";
        
        queryCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/chr" + std::to_string(a_nChrNumber) + "_CompliantsALL.txt";
        queryViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/chr" + std::to_string(a_nChrNumber) + "_ViolationsALL.txt";
    }
    
    else if(a_bIsFilterOverlap)
    {
        vcfFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered.vcf";
        
        truthCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETFilteredOverlap/Chr" + std::to_string(a_nChrNumber) + "_CompliantsTRUTH.txt";
        truthViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETFilteredOverlap/Chr" + std::to_string(a_nChrNumber) + "_ViolationsTRUTH.txt";
        
        queryCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/NoParent00Test/chr" + std::to_string(a_nChrNumber) + "_CompliantsALL.txt";
        queryViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/NoParent00Test/chr" + std::to_string(a_nChrNumber) + "_ViolationsALL.txt";
        
    }
    
    else
    {
        vcfFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat.vcf";
        
        truthCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETUnfiltered/Chr" + std::to_string(a_nChrNumber) + "_CompliantsTRUTH.txt";
        truthViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETUnfiltered/Chr" + std::to_string(a_nChrNumber) + "_ViolationsTRUTH.txt";
        
        queryCompliantSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/chr" + std::to_string(a_nChrNumber) + "_CompliantsALL.txt";
        queryViolationSetFileName = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/chr" + std::to_string(a_nChrNumber) + "_ViolationsALL.txt";
        
    }

    std::cout << "##################       CHR " << a_nChrNumber << "        ##################" << std::endl;
    SCompareResult resCompliant;
    SCompareResult resViolation;
    
    std::cout << "===================COMPLIANT COMPARISON=======================" << std::endl << std::endl;
    //COMPLIANTS COMPARISON
    resCompliant = compare2List(truthCompliantSetFileName, queryCompliantSetFileName);
    std::cout << "===================VIOLATION COMPARISON=======================" << std::endl << std::endl;
    //VIOLATIONS COMPARISON
    resViolation = compare2List(truthViolationSetFileName, queryViolationSetFileName);
    
    std::cout << "TOTAL COMMON: " << resCompliant.common + resViolation.common << std::endl;
    std::cout << "TOTAL UNMATCH:" << resViolation.unMatchOrig + resCompliant.unMatchOrig << std::endl;
    
}



void GenerateTruthSetsMaria(std::string a_rFilename, bool a_bIsFilterOverlap, bool a_bIsFilter00)
{
    CVcfReader triplet;
    
    std::string modeSubPath = "";
    if(a_bIsFilterOverlap)
        modeSubPath = modeSubPath + "Overlap";
    if(a_bIsFilter00)
        modeSubPath = modeSubPath + "00";
    
    triplet.Open(a_rFilename.c_str());
    
    int CHR_ITERATOR = 23;
    
    std::string commonPath = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/TruthSETFiltered" + modeSubPath + "/Chr";
    std::ofstream outputFileCompliants;
    std::ofstream outputFileViolations;

    std::vector<CVariant> motherVariants;
    std::vector<CVariant> fatherVariants;
    std::vector<CVariant> childVariants;
    
    
    CVariant variants[3];
    bool isAdd;
    
    bool endLoop = false;
    
    while(triplet.GetNextRecordMultiSample(variants))
    {
        
        if(variants[0].m_nChrId == -1 || CHR_ITERATOR > 23)
        {
            endLoop = true;
        }
            
        if(variants[0].m_nChrId != CHR_ITERATOR)
        {
            std::vector<int> mendelianViolations;
            std::vector<int> mendelianCompliants;
            
            for(int k = 0; k < motherVariants.size(); k++)
            {
                bool c1 = motherVariants[k].m_genotype[0] == childVariants[k].m_genotype[0] || motherVariants[k].m_genotype[1] == childVariants[k].m_genotype[0];
                bool c2 = motherVariants[k].m_genotype[0] == childVariants[k].m_genotype[1] || motherVariants[k].m_genotype[1] == childVariants[k].m_genotype[1];
                bool c3 = fatherVariants[k].m_genotype[0] == childVariants[k].m_genotype[0] || fatherVariants[k].m_genotype[1] == childVariants[k].m_genotype[0];
                bool c4 = fatherVariants[k].m_genotype[0] == childVariants[k].m_genotype[1] || fatherVariants[k].m_genotype[1] == childVariants[k].m_genotype[1];
                
                if((c1 == true && c4 == true) || (c2 == true && c3 == true))
                {
                    mendelianCompliants.push_back(childVariants[k].m_alleles[0].m_nStartPos);
                }
                else
                    mendelianViolations.push_back(childVariants[k].m_alleles[0].m_nStartPos);
            }
            
            outputFileCompliants.open(commonPath + std::to_string(CHR_ITERATOR) + "_CompliantsTRUTH.txt");
            for(int k : mendelianCompliants)
                outputFileCompliants << k << std::endl;
            outputFileCompliants.close();
            
            outputFileViolations.open(commonPath + std::to_string(CHR_ITERATOR) + "_ViolationsTRUTH.txt");
            for(int k : mendelianViolations)
                outputFileViolations << k << std::endl;
            outputFileViolations.close();
            
            motherVariants.clear();
            fatherVariants.clear();
            childVariants.clear();
            mendelianCompliants.clear();
            mendelianViolations.clear();
            CHR_ITERATOR++;
            
            if(endLoop == true)
                break;
        }
        
        isAdd = true;
        for(int k = 0; k < variants[0].m_nAlleleCount; k++)
        {
            if(variants[0].m_alleles[k].m_sequence == "*")
            {
                variants[0].m_genotype[0] = 0;
                variants[0].m_genotype[1] = 0;
            }
            
        }

        for(int k = 0; k < variants[1].m_nAlleleCount; k++)
        {
            if(variants[1].m_alleles[k].m_sequence == "*")
            {
                variants[1].m_genotype[0] = 0;
                variants[1].m_genotype[1] = 0;
            }
        }

        bool isSkip = false;
        
        for(int k = 0; k < variants[2].m_nAlleleCount; k++)
        {
            if(variants[2].m_alleles[k].m_sequence == "*")
                isSkip = true;
        }

        if(!isSkip)
        {
            motherVariants.push_back(variants[0]);
            fatherVariants.push_back(variants[1]);
            childVariants.push_back(variants[2]);
        }

    }
}











