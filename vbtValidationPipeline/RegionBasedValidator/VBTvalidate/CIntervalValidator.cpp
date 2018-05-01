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

#include "CVariant.h"
#include "CIntervalValidator.h"
#include "CFastaParser.h"
#include "CVariantProvider.h"
#include <vector>

using namespace vbtvalidator;

void CIntervalValidator::SetVariantProvider(CVariantProvider *a_pProvider)
{
    m_pVariantProvider = a_pProvider;
}

void CIntervalValidator::SetFastaReader(CFastaParser& a_rFastaReader, SContig* a_pContig)
{
    m_pFastaReader = &a_rFastaReader;
    m_pContig = a_pContig;
}

void CIntervalValidator::GeneratePhases(std::deque<SPhase>& a_rPhaseVector, std::vector<const CVariant*>& a_rVariants) const
{
    SPhase initialPhase;
    
    //Initialize phase list with null SPhase object
    a_rPhaseVector.push_back(initialPhase);
    
    //Get All possible phasing list for variants
    while(true)
    {
        if(a_rPhaseVector.front().m_phaseVector.size() == a_rVariants.size())
            break;
        
        SPhase curSeq1(a_rPhaseVector.front());
        SPhase curSeq2(a_rPhaseVector.front());
        a_rPhaseVector.pop_front();
        
        curSeq1.m_phaseVector.push_back(1);
        a_rPhaseVector.push_back(curSeq1);
        
        if(a_rVariants[curSeq2.m_phaseVector.size()]->m_bIsHeterozygous)
        {
            curSeq2.m_phaseVector.push_back(0);
            a_rPhaseVector.push_back(curSeq2);
        }
    }
}

EIntervalDecision CIntervalValidator::IntervalTestMV_RegionBased(const std::string& a_rChrName, const SInterval& a_rInterval, int& hasViolation) const
{
    //Containers to store variants in given interval
    std::vector<const CVariant*> childVars;
    std::vector<const CVariant*> motherVars;
    std::vector<const CVariant*> fatherVars;
    
    //Containers to store phasing combinations of variants in given interval
    std::deque<SPhase> childSubPhases;
    std::deque<SPhase> motherSubPhases;
    std::deque<SPhase> fatherSubPhases;
    
    
    //Fill mother, father and child variants reside in given interval (will get consistent variants only)
    m_pVariantProvider->GetVariantsAll(a_rInterval, motherVars, fatherVars, childVars, hasViolation);
    
    
    //If there are more than 15 variants in a region, we mark the region as complex and do not process it
    if(motherVars.size() + fatherVars.size() + childVars.size() > 15)
        return eInterval_Complex;
    
    //If there is no consistent variant in region, we mark it as consistent without processing
    if(motherVars.size() + fatherVars.size() + childVars.size() == 0)
        return eInterval_CorrectlyAssessed;
    
    //Some variants overlaps with the sync point to prevent errors, we added 100 base pair to both sides
    SInterval interval;
    interval.m_nStart = a_rInterval.m_nStart - 100;
    interval.m_nEnd = a_rInterval.m_nEnd + 100;
    
    //Clip the substring from FASTA referance
    std::string intervalString = m_pFastaReader->GetSubSequence(*m_pContig, interval.m_nStart, interval.m_nEnd);
    
    //Get All possible phasing list for Child variants
    GeneratePhases(childSubPhases, childVars);
    GeneratePhases(motherSubPhases, motherVars);
    GeneratePhases(fatherSubPhases, fatherVars);
    
    //List of sequences after variants are applied to the clipped reference in all possible combinations
    std::deque<SSequence> childSequences;
    std::deque<SSequence> motherSequences;
    std::deque<SSequence> fatherSequences;
    
    //Create all replayed sequence given child phase list
    for(int k = 0; k < childSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(childSubPhases[k], interval, intervalString, childVars);
        if(seq.m_haplotypeA != "NONE")
            childSequences.push_back(seq);
    }
    
    //Create all replayed sequence given mother phase list
    for(int k = 0; k < motherSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(motherSubPhases[k], interval, intervalString, motherVars);
        if(seq.m_haplotypeA != "NONE")
            motherSequences.push_back(seq);
    }
    
    //Create all replayed sequence given father phase list
    for(int k = 0; k < fatherSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(fatherSubPhases[k], interval, intervalString, fatherVars);
        if(seq.m_haplotypeA != "NONE")
            fatherSequences.push_back(seq);
    }
    
    //Compare all generated sequences and check if there is a child sequence that matches to both mother and father
    for(int c = 0; c < childSequences.size(); c++)
    {
        for(int m = 0; m < motherSequences.size(); m++)
        {
            for(int f = 0; f < fatherSequences.size(); f++)
            {
                if(childSequences[c].m_haplotypeA == motherSequences[m].m_haplotypeA
                   &&
                   childSequences[c].m_haplotypeB == fatherSequences[f].m_haplotypeA)
                {
                    //A consistent combination found!
                    if(hasViolation > 0)
                        return eInterval_WrongAssessed;
                    else
                        return eInterval_CorrectlyAssessed;
                }
            }
        }
    }
    
    //No consistent combination found!
    if(hasViolation > 0)
        return eInterval_CorrectlyAssessed;
    else
        return eInterval_WrongAssessed;
}

EIntervalDecision CIntervalValidator::IntervalTestMV(const std::string& a_rChrName, const SInterval& a_rInterval, int& a_rViolationCount) const
{
    //Containers to store variants in given interval
    std::vector<const CVariant*> childVars;
    std::vector<const CVariant*> motherVars;
    std::vector<const CVariant*> fatherVars;
    bool bIsViolationOverlapWithConsistent = false;
    
    //Containers to store phasing combinations of variants in given interval
    std::deque<SPhase> childSubPhases;
    std::deque<SPhase> motherSubPhases;
    std::deque<SPhase> fatherSubPhases;
    
    
    //Fill mother, father and child variants reside in given interval (will get consistent variants only)
    m_pVariantProvider->GetVariants(a_rInterval, motherVars, fatherVars, childVars, bIsViolationOverlapWithConsistent, a_rViolationCount);
    
    //If there is a violation that overlaps with a consistent variant, then mark interval decision as wrong
    if(true == bIsViolationOverlapWithConsistent)
       return eInterval_WrongAssessed;
    
    //If there are more than 15 variants in a region, we mark the region as complex and do not process it
    if(motherVars.size() + fatherVars.size() + childVars.size() > 15)
        return eInterval_Complex;
    
    //If there is no consistent variant in region, we mark it as consistent without processing
    if(motherVars.size() + fatherVars.size() + childVars.size() == 0)
        return eInterval_CorrectlyAssessed;
    
    //Some variants overlaps with the sync point to prevent errors, we added 100 base pair to both sides
    SInterval interval;
    interval.m_nStart = a_rInterval.m_nStart - 100;
    interval.m_nEnd = a_rInterval.m_nEnd + 100;
    
    //Clip the substring from FASTA referance
    std::string intervalString = m_pFastaReader->GetSubSequence(*m_pContig, interval.m_nStart, interval.m_nEnd);
    
    //Get All possible phasing list for Child variants
    GeneratePhases(childSubPhases, childVars);
    GeneratePhases(motherSubPhases, motherVars);
    GeneratePhases(fatherSubPhases, fatherVars);
    
    //List of sequences after variants are applied to the clipped reference in all possible combinations
    std::deque<SSequence> childSequences;
    std::deque<SSequence> motherSequences;
    std::deque<SSequence> fatherSequences;
    
    //Create all replayed sequence given child phase list
    for(int k = 0; k < childSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(childSubPhases[k], interval, intervalString, childVars);
        if(seq.m_haplotypeA != "NONE")
            childSequences.push_back(seq);
    }
    
    //Create all replayed sequence given mother phase list
    for(int k = 0; k < motherSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(motherSubPhases[k], interval, intervalString, motherVars);
        if(seq.m_haplotypeA != "NONE")
            motherSequences.push_back(seq);
    }
    
    //Create all replayed sequence given father phase list
    for(int k = 0; k < fatherSubPhases.size(); k++)
    {
        SSequence seq = ApplyPhasing(fatherSubPhases[k], interval, intervalString, fatherVars);
        if(seq.m_haplotypeA != "NONE")
            fatherSequences.push_back(seq);
    }
    
    //Compare all generated sequences and check if there is a child sequence that matches to both mother and father
    for(int c = 0; c < childSequences.size(); c++)
    {
        for(int m = 0; m < motherSequences.size(); m++)
        {
            for(int f = 0; f < fatherSequences.size(); f++)
            {
                if(childSequences[c].m_haplotypeA == motherSequences[m].m_haplotypeA
                   &&
                   childSequences[c].m_haplotypeB == fatherSequences[f].m_haplotypeA)
                {
                    //A consistent combination found! Return consistent
                    return eInterval_CorrectlyAssessed;
                }
            }
        }
    }
    
    //No consistent combination found. Return violation
    return eInterval_WrongAssessed;
}

SSequence CIntervalValidator::ApplyPhasing(const SPhase& a_rPhase, const vbtvalidator::SInterval& a_rInterval, const std::string& a_rReferenceSeq, std::vector<const CVariant*>& a_rVariantList) const
{
    
    SSequence seq;
    
    int shift = a_rInterval.m_nStart;
    seq.m_haplotypeA = "";
    seq.m_haplotypeB = "";
    
    int lastHapAIndex = 0;
    int lastHapBIndex = 0;

    int refLength = static_cast<int>(a_rReferenceSeq.length());
    
    for(int p = 0; p < a_rPhase.m_phaseVector.size(); p++)
    {
        int phase = a_rPhase.m_phaseVector[p];
        int revphase = 1 - a_rPhase.m_phaseVector[p];
        
        //====PLAY HAPLOTYPE A====
        
        //If the allele is not reference, apply to the haplotype
        if(!a_rVariantList[p]->m_alleles[phase].m_bIsIgnored)
        {
            int refBasePairSizeToAdd_A = (a_rVariantList[p]->m_alleles[phase].m_nStartPos - shift) - lastHapAIndex;
            
            //This phasing vector cannot be played
            if(refLength <= lastHapAIndex || refBasePairSizeToAdd_A > refLength || refBasePairSizeToAdd_A < 0)
            {
                seq.m_haplotypeA = "NONE";
                seq.m_haplotypeB = "NONE";
                return seq;
            }
            
            else
            {
                //Add basepairs between variants from reference
                seq.m_haplotypeA += a_rReferenceSeq.substr(lastHapAIndex, refBasePairSizeToAdd_A);
                //Add basepairs of variant
                seq.m_haplotypeA += a_rVariantList[p]->m_alleles[phase].m_sequence;
                
                lastHapAIndex = a_rVariantList[p]->m_alleles[phase].m_nEndPos - shift;
            }
        }
        
        //====PLAY HAPLOTYPE B====
        
        
        //If the allele is not reference, apply to the haplotype
        if(!a_rVariantList[p]->m_alleles[revphase].m_bIsIgnored)
        {
            int refBasePairSizeToAdd_B = (a_rVariantList[p]->m_alleles[revphase].m_nStartPos - shift) - lastHapBIndex;
        
            //This phasing vector cannot be played
            if(refLength <= lastHapBIndex || refBasePairSizeToAdd_B > refLength || refBasePairSizeToAdd_B < 0)
            {
                seq.m_haplotypeA = "NONE";
                seq.m_haplotypeB = "NONE";
                return seq;
            }
            
            else
            {
                //Add basepairs between variants from reference
                seq.m_haplotypeB += a_rReferenceSeq.substr(lastHapBIndex, refBasePairSizeToAdd_B);
                //Add basepairs of variant
                seq.m_haplotypeB += a_rVariantList[p]->m_alleles[revphase].m_sequence;
                
                lastHapBIndex = a_rVariantList[p]->m_alleles[revphase].m_nEndPos - shift;
            }
        }
    }
    
    //Add remaining base pairs from last variant end position to end of interval
    int lengthA = 1 + a_rInterval.m_nEnd - shift - lastHapAIndex;
    if(lastHapAIndex + lengthA <= refLength && lengthA > 0)
        seq.m_haplotypeA += a_rReferenceSeq.substr(lastHapAIndex, lengthA);
    else
    {
        seq.m_haplotypeA = "NONE";
        seq.m_haplotypeB = "NONE";
        return seq;
    }
    
    int lengthB = 1 + a_rInterval.m_nEnd - shift - lastHapBIndex;
    if(lastHapBIndex + lengthB <= refLength && lengthB > 0)
        seq.m_haplotypeB += a_rReferenceSeq.substr(lastHapBIndex, lengthB);
    else
    {
        seq.m_haplotypeA = "NONE";
        seq.m_haplotypeB = "NONE";
        return seq;
    }

    return seq;
}

