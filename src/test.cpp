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







