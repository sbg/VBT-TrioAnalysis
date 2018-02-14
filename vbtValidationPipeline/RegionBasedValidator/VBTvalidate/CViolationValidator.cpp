/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  Created by Berke Cagkan Toptas
 *
 */

#include "CViolationValidator.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "CLBLComparisonTool.h"

using namespace vbtvalidator;

bool IsInside(const SInterval& a_rInterval, int a_nViolationPos)
{
    return a_rInterval.m_nStart <= a_nViolationPos && a_rInterval.m_nEnd > a_nViolationPos;
}

bool IsOverlap(const SInterval& lhs, const SInterval& rhs)
{
    if(lhs.m_nStart == rhs.m_nStart)
        return true;
    if (lhs.m_nStart == lhs.m_nEnd)
        return (rhs.m_nEnd > lhs.m_nStart && rhs.m_nStart <= lhs.m_nStart);
    else if(rhs.m_nStart == rhs.m_nEnd)
        return (lhs.m_nEnd > rhs.m_nStart && lhs.m_nStart <= rhs.m_nStart);
    else
        return std::min(rhs.m_nEnd, lhs.m_nEnd) - std::max(lhs.m_nStart, rhs.m_nStart) > 0;
}

int CViolationValidator::Run(int argc, const char* argv[])
{
    std::string fastaFile;
    std::string inputTrio;
    std::string originalTrio;
    std::string motherChildIntervalFile;
    std::string fatherChildIntervalFile;
    std::string pedigreeFile;
    bool bCountPerRegionVariant = false;
    bool bIsTrimBeginningFirst = true;
    
    for(int k = 2; k < argc;)
    {
        if(strcmp(argv[k], "-input-vbt-trio") == 0)
        {
            inputTrio = std::string(argv[k+1]);
            k+=2;
        }
        
        else if(strcmp(argv[k], "-input-original-trio") == 0)
        {
            originalTrio = std::string(argv[k+1]);
            k+=2;
        }
        
        else if(strcmp(argv[k], "-reference") == 0)
        {
            fastaFile = std::string(argv[k+1]);
            k+=2;
        }
        
        else if(strcmp(argv[k], "-mother-child-interval") == 0)
        {
            motherChildIntervalFile = std::string(argv[k+1]);
            k+=2;
        }
        
        else if(strcmp(argv[k], "-father-child-interval") == 0)
        {
            fatherChildIntervalFile = std::string(argv[k+1]);
            k+=2;
        }
        
        else if(strcmp(argv[k], "--count-per-region-variant") == 0)
        {
            bCountPerRegionVariant = true;
            k++;
        }
        
        else if(strcmp(argv[k], "--trim-endings-first") == 0)
        {
            bIsTrimBeginningFirst = false;
            k++;
        }
        
        else if(strcmp(argv[k], "-ped-file") == 0)
        {
            pedigreeFile = std::string(argv[k+1]);
            k+=2;
        }
        
        else
        {
            std::cerr << "Wrong Parameter Detected! Terminating Program" << std::endl;
            return -1;
        }
    }
    
    //Read and merge mother-child and father-child syncpoints
    SetIntervalFiles(motherChildIntervalFile, fatherChildIntervalFile);
    
    if(bCountPerRegionVariant == true)
    {
        std::cerr << "[STDERR] Counting Total Variant Counts Per Region" << std::endl;
        SetInputFiles(inputTrio, originalTrio, fastaFile, pedigreeFile, bIsTrimBeginningFirst);
        
        CalculatePerRegionVariantCounts(inputTrio);
    }
    
    else
    {
        std::cerr << "[STDERR] Validating VBT and LBL Tool Decisions" << std::endl;
        
        //Set name of annotated trio vcf
        std::string annotatedLBLTrioPath = originalTrio.substr(0, originalTrio.length()-4) + "_LBLannotated.vcf";
        
        //**** OPEN HERE TO GENERATE NAIVE DECISION ANNOTATED VCF FILES ****
        //Generating Naive decision (LBL) annotated vcf file
        /*
            std::cerr << "Generating Naive decision (LBL) annotated vcf file" << std::endl;
            CLBLComparisonTool LBLtool;
            LBLtool.SetInputTrio(originalTrio, pedigreeFile, bIsTrimBeginningFirst);
            LBLtool.WriteOutputTrio(annotatedLBLTrioPath);
        */
         
        SetInputFiles(inputTrio, originalTrio, fastaFile, pedigreeFile, bIsTrimBeginningFirst);
        ValidateRegions();
    }
    
    return 0;
}

void CViolationValidator::ValidateRegions()
{
    int vbtViolationMissedRegion = 0;
    int vbtExtraViolationRegion = 0;
    int vbtTrueRegion = 0;

    int lblViolationMissedRegion = 0;
    int lblExtraViolationRegion = 0;
    int lblTrueRegion = 0;

    int totalRegionCount = 0;
    int complexRegionCount = 0;

    std::string contigNameVBT, contigNameLBL;
    
    bool bVBTHasNext = true;
    bool bLBLHasNext = true;
    
    while(bVBTHasNext && bLBLHasNext)
    {
        bVBTHasNext = m_vbtVariantProvider.ReadNextChromosome(contigNameVBT);
        bLBLHasNext = m_lblVariantProvider.ReadNextChromosome(contigNameLBL);
        
        std::cerr << "Processing chromosome " << contigNameVBT << std::endl;
        
        std::vector<SInterval> ViolationIntervals;
        
        //If read chromosomes are different then break
        if(contigNameLBL != contigNameVBT)
            break;
        
        //Read New Chromosome Contig from Fasta
        std::cerr << "Reading contig from FASTA" << std::endl;
        m_contig.Clean();
        m_fastaReader.FetchNewChromosome(contigNameVBT, m_contig);
        
        //Get violation intervals for VBT and LBL
        m_vbtVariantProvider.GetViolationIntervals(m_intervalsListMap[contigNameVBT], ViolationIntervals);
        m_lblVariantProvider.GetViolationIntervals(m_intervalsListMap[contigNameLBL], ViolationIntervals);
        
        //Remove Duplicates
        std::sort(ViolationIntervals.begin(), ViolationIntervals.end(), [](const SInterval& left, const SInterval& right) {return left.m_nStart < right.m_nStart;});
        ViolationIntervals.erase(std::unique(ViolationIntervals.begin(),
                                             ViolationIntervals.end(),
                                             [](const SInterval& left, const SInterval& right) {return left.m_nStart == right.m_nStart && left.m_nEnd == right.m_nEnd;}),
                                 ViolationIntervals.end());
        
        //Update total region count
        totalRegionCount += ViolationIntervals.size();
        
        std::vector<EIntervalDecision> lblIntervalDecisions(ViolationIntervals.size());
        std::vector<EIntervalDecision> vbtIntervalDecisions(ViolationIntervals.size());
        std::vector<int> lblViolationCounts(ViolationIntervals.size());
        std::vector<int> vbtViolationCounts(ViolationIntervals.size());

        
        //STEP 1: for each selected region, run Interval validator (VBT and Other tool separately)
        for(int m = 0; m < ViolationIntervals.size(); m++)
        {
            //std::cerr << "Checking " << m << " of " << ViolationIntervals.size() << " intervals.." << std::endl;
            lblIntervalDecisions[m] = m_lblIntervalValidator.IntervalTestMV(contigNameLBL, ViolationIntervals[m], lblViolationCounts[m]);
            vbtIntervalDecisions[m] = m_vbtIntervalValidator.IntervalTestMV(contigNameVBT, ViolationIntervals[m], vbtViolationCounts[m]);
        }
        
        std::vector<SInterval*> vbtWrongIntervals;
        std::vector<SInterval*> vbtWrongIntervalsFromConsistent;
        std::vector<SInterval*> lblWrongIntervals;
        
        //STEP 2: for each region, compare result of VBT and Other tool
        for(int m = 0; m < ViolationIntervals.size(); m++)
        {
            //Discard complex regions
            if(vbtIntervalDecisions[m] == EIntervalDecision::eInterval_Complex || lblIntervalDecisions[m] == EIntervalDecision::eInterval_Complex)
            {
                complexRegionCount++;
                continue;
            }
            
            //Both is wrong
            else if (lblIntervalDecisions[m] == EIntervalDecision::eInterval_WrongAssessed && vbtIntervalDecisions[m] == EIntervalDecision::eInterval_WrongAssessed)
            {
                vbtViolationMissedRegion++;
                lblViolationMissedRegion++;
                vbtWrongIntervals.push_back(&ViolationIntervals[m]);
                lblWrongIntervals.push_back(&ViolationIntervals[m]);
            }
            
            // Naive is correct VBT is wrong
            else if (lblIntervalDecisions[m] == EIntervalDecision::eInterval_CorrectlyAssessed && vbtIntervalDecisions[m] == EIntervalDecision::eInterval_WrongAssessed)
            {
                vbtViolationMissedRegion++;
                lblTrueRegion++;
                vbtWrongIntervals.push_back(&ViolationIntervals[m]);
            }
            
            // Naive is wrong VBT is correct
            else if (lblIntervalDecisions[m] == EIntervalDecision::eInterval_WrongAssessed && vbtIntervalDecisions[m] == EIntervalDecision::eInterval_CorrectlyAssessed)
            {
                vbtTrueRegion++;
                lblViolationMissedRegion++;
                lblWrongIntervals.push_back(&ViolationIntervals[m]);
            }
            
            // Both tool give correct
            else if (lblIntervalDecisions[m] == EIntervalDecision::eInterval_CorrectlyAssessed && vbtIntervalDecisions[m] == EIntervalDecision::eInterval_CorrectlyAssessed)
            {
                
                if(lblViolationCounts[m] == vbtViolationCounts[m])
                {
                    vbtTrueRegion++;
                    lblTrueRegion++;
                }
                else if(lblViolationCounts[m] > vbtViolationCounts[m])
                {
                    vbtTrueRegion++;
                    lblExtraViolationRegion++;
                }
                else
                {
                    vbtWrongIntervalsFromConsistent.push_back(&ViolationIntervals[m]);
                    vbtExtraViolationRegion++;
                    lblTrueRegion++;
                }
            }
        }
        
        std::cerr << "VBT Wrong Intervals: (" << vbtWrongIntervals.size() << "/" << ViolationIntervals.size() << ")" << std::endl;
        for(int m = 0; m < vbtWrongIntervals.size(); m++)
            std::cerr << vbtWrongIntervals[m]->m_nStart << " " << vbtWrongIntervals[m]->m_nEnd << std::endl;
        std::cerr << std::endl;
        std::cerr << "VBT Wrong Intervals From consistent: (" << vbtWrongIntervalsFromConsistent.size() << "/" << ViolationIntervals.size() << ")" << std::endl;
        for(int m = 0; m < vbtWrongIntervalsFromConsistent.size(); m++)
            std::cerr << vbtWrongIntervalsFromConsistent[m]->m_nStart << " " << vbtWrongIntervalsFromConsistent[m]->m_nEnd << std::endl;
        std::cerr << std::endl;
        std::cerr << "LBL Wrong Intervals: (" << lblWrongIntervals.size() << "/" << ViolationIntervals.size() << ")" << std::endl;
        for(int m = 0; m < lblWrongIntervals.size(); m++)
            std::cerr << lblWrongIntervals[m]->m_nStart << " " << lblWrongIntervals[m]->m_nEnd << std::endl;
        std::cerr << std::endl;

    }
    
    std::cout << "Total Region Count   : " << totalRegionCount << std::endl;
    std::cout << std::endl;
    std::cout << "VBT Correct Interval                     : " << vbtTrueRegion << std::endl;
    std::cout << "VBT Wrong   Interval (Missing Violation) : " << vbtViolationMissedRegion << std::endl;
    std::cout << "VBT Wrong   Interval (Extra Violation  ) : " << vbtExtraViolationRegion << std::endl;
    std::cout << std::endl;
    std::cout << "LBL Correct Interval                     : " << lblTrueRegion << std::endl;
    std::cout << "LBL Wrong   Interval (Missing Violation) : " << lblViolationMissedRegion << std::endl;
    std::cout << "LBL Wrong   Interval (Extra Violation  ) : " << lblExtraViolationRegion << std::endl;
    std::cout << std::endl;
    std::cout << "VBT Precision : " << (float)vbtTrueRegion / (float)(vbtTrueRegion + vbtExtraViolationRegion) << std::endl;
    std::cout << "VBT Recall    : " << (float)vbtTrueRegion / (float)(vbtTrueRegion + vbtViolationMissedRegion) << std::endl;
    std::cout << std::endl;
    std::cout << "LBL Precision : " << (float)lblTrueRegion / (float)(lblTrueRegion + lblExtraViolationRegion) << std::endl;
    std::cout << "LBL Recall    : " << (float)lblTrueRegion / (float)(lblTrueRegion + lblViolationMissedRegion) << std::endl;
    std::cout << std::endl;
    std::cout << "VBT Accuracy : " << (float)vbtTrueRegion / (float)(totalRegionCount - complexRegionCount) << std::endl;
    std::cout << "LBL Accuracy : " << (float)lblTrueRegion / (float)(totalRegionCount - complexRegionCount) << std::endl;
}

void CViolationValidator::CalculatePerRegionVariantCounts(std::string a_rTrioPath)
{
    CVcfReader TrioReader;
    CVariant trioVars[3];
    TrioReader.Open(a_rTrioPath.c_str());
    SConfig cfg;
    cfg.m_bIsRefOverlap = true;
    
    std::map<int, int> RegionCountPerLength;
    
    std::string currentChromosome = "1";
    int intervalListMapIterator = 0;
    int currentVarCount = 0;
    
    while (TrioReader.GetNextRecordMultiSample(trioVars, bIsTrimBeginningFirst))
    {
        if(currentChromosome != trioVars[0].m_chrName)
        {
            RegionCountPerLength[currentVarCount]++;
            currentVarCount = 0;
            intervalListMapIterator = 0;
            currentChromosome = trioVars[0].m_chrName;
        }
        
        else if(m_intervalsListMap[currentChromosome][intervalListMapIterator].m_nEnd <= trioVars[0].m_nStartPos)
        {
            RegionCountPerLength[currentVarCount]++;
            currentVarCount = 0;
            intervalListMapIterator++;
        }
        
        for(int k = 0; k < 3; k++)
        {
            if(trioVars[k].m_genotype[0] != 0 || trioVars[k].m_genotype[1] != 0)
                currentVarCount++;
        }
    }
    
    TrioReader.Close();
    
    int totalMissedVarCount = 0;
    
    std::cout << "Region counts per length : (starting from length 0):" << std::endl;
    
    for(int k = 0; k <= 100; k++)
    {
        std::cout << RegionCountPerLength[k] << std::endl;
        if(k > 14)
            totalMissedVarCount += (RegionCountPerLength[k] * k);
    }
    
    std::cout << "Total Missed Var: " << totalMissedVarCount << std::endl;
}

void CViolationValidator::SetInputFiles(const std::string& a_rTrioPathVBT,
                                        const std::string& a_rTrioPathOriginal,
                                        const std::string& a_rFastaPath,
                                        const std::string& a_rPedigreeFile,
                                        bool a_bIsBeginTrimmingFirst)
{
    //Initialize fasta reader
    m_fastaReader.OpenFastaFile(a_rFastaPath.c_str());
    
    std::string MotherSample, FatherSample, ChildSample;
    
    //Initialize variant providers
    m_lblVariantProvider.SetTrioPath(a_rTrioPathOriginal, a_rPedigreeFile, a_bIsBeginTrimmingFirst);
    m_vbtVariantProvider.SetTrioPath(a_rTrioPathVBT, "", a_bIsBeginTrimmingFirst);
    m_lblTrioPath = a_rTrioPathVBT;
    m_vbtTrioPath = a_rTrioPathOriginal;
    
    
    
    //Set variant provider to validators
    m_lblIntervalValidator.SetVariantProvider(&m_lblVariantProvider);
    m_vbtIntervalValidator.SetVariantProvider(&m_vbtVariantProvider);
    
    //Set fasta reader to validators
    m_lblIntervalValidator.SetFastaReader(m_fastaReader, &m_contig);
    m_vbtIntervalValidator.SetFastaReader(m_fastaReader, &m_contig);
    
    //Set the initial trimming choice
    bIsTrimBeginningFirst = a_bIsBeginTrimmingFirst;
}


void CViolationValidator::MergeIntervals(const std::string& a_rMotherChildIntervalsPath, const std::string& a_rFatherChildIntervalsPath)
{
    
    std::ifstream motherIntervalFile;
    std::ifstream fatherIntervalFile;
    
    motherIntervalFile.open(a_rMotherChildIntervalsPath.c_str());
    fatherIntervalFile.open(a_rFatherChildIntervalsPath.c_str());
    
    SInterval nextMotherInterval;
    std::string chrMother;
    SInterval nextFatherInterval;
    std::string chrFather;
    
    //initialize father and mother chromosomes and first intervals
    fatherIntervalFile >> chrFather >> nextFatherInterval.m_nStart >> nextFatherInterval.m_nEnd;
    motherIntervalFile >> chrMother >> nextMotherInterval.m_nStart >> nextMotherInterval.m_nEnd;
    
    while(!motherIntervalFile.eof() && !fatherIntervalFile.eof())
    {
        if(chrFather == chrMother)
        {
            if(nextFatherInterval.m_nEnd == nextMotherInterval.m_nEnd)
            {
                //father and mother sync point is equal, push it to the merged interval list
                m_intervalsListMap[chrFather].push_back(nextFatherInterval);
                fatherIntervalFile >> chrFather >> nextFatherInterval.m_nStart >> nextFatherInterval.m_nEnd;
                motherIntervalFile >> chrMother >> nextMotherInterval.m_nStart >> nextMotherInterval.m_nEnd;
            }
            else
            {
                //sync points are differ, iterate over mother and father until a common sync point will be found
                SInterval tmpInterval;
                tmpInterval.m_nStart = nextFatherInterval.m_nStart;
                
                while(true)
                {
                    while(!motherIntervalFile.eof() && chrFather == chrMother && nextMotherInterval.m_nEnd < nextFatherInterval.m_nEnd)
                        motherIntervalFile >> chrMother >> nextMotherInterval.m_nStart >> nextMotherInterval.m_nEnd;
                    
                    while(!fatherIntervalFile.eof() && chrFather == chrMother && nextFatherInterval.m_nEnd < nextMotherInterval.m_nEnd)
                        fatherIntervalFile >> chrFather >> nextFatherInterval.m_nStart >> nextFatherInterval.m_nEnd;
                    
                    //Parsed until end of chromosome and no common sync point found, end searching
                    if(chrFather != chrMother)
                        break;
                    
                    //Find a common sync point, add it to the merged interval list
                    if(nextFatherInterval.m_nEnd == nextMotherInterval.m_nEnd)
                    {
                        tmpInterval.m_nEnd = nextFatherInterval.m_nEnd;
                        m_intervalsListMap[chrFather].push_back(tmpInterval);
                        fatherIntervalFile >> chrFather >> nextFatherInterval.m_nStart >> nextFatherInterval.m_nEnd;
                        motherIntervalFile >> chrMother >> nextMotherInterval.m_nStart >> nextMotherInterval.m_nEnd;
                        break;
                    }
                }
            }
        }
        
        else if(std::stoi(chrFather) > std::stoi(chrMother))
        {
            //There are more sync points in mother side, iterate over them until we advance to the next chromosome
            while(!motherIntervalFile.eof() && chrFather != chrMother)
                motherIntervalFile >> chrMother >> nextMotherInterval.m_nStart >> nextMotherInterval.m_nEnd;
        }
        
        else
        {
            //There are more sync points in father side, iterate over them until we advance to the next chromosome
            while(!fatherIntervalFile.eof() && chrFather != chrMother)
                fatherIntervalFile >> chrFather >> nextFatherInterval.m_nStart >> nextFatherInterval.m_nEnd;
        }
        
    }
    
    motherIntervalFile.close();
    fatherIntervalFile.close();
}


void CViolationValidator::SetIntervalFiles(const std::string& a_rMotherChildIntervalsPath, const std::string& a_rFatherChildIntervalsPath)
{
    std::cerr << "[stderr] Merging Intervals..." << std::endl;
    MergeIntervals(a_rMotherChildIntervalsPath, a_rFatherChildIntervalsPath);
    std::cerr << "[stderr] Intervals are merged" << std::endl;
}






