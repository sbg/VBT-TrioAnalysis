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

#include "CVariantProvider.h"
#include "CIntervalValidator.h"
#include <algorithm>
#include <iostream>

using namespace vbtvalidator;

bool CVariantProvider::CompareVariants(const CVariant& var1, const CVariant& var2)
{
    if(var1.m_nStartPos != var2.m_nStartPos)
        return var1.m_nStartPos < var2.m_nStartPos;
    else if(var1.m_nEndPos != var2.m_nEndPos)
        return var1.m_nEndPos < var2.m_nEndPos;
    else
        return var1.m_nId < var2.m_nId;
}

//Checks if the given two range is overlapping
bool isOverlap(int left1, int right1, int left2, int right2)
{
    if(left1 == left2)
        return true;
    if(right1 == left1)
        return (right2 > left1 && left2 <= left1);
    else if(right2 == left2)
        return (right1 > left2 && left1 <= left2);
    else
        return std::min(right1, right2) - std::max(left1, left2) > 0;
}

void CVariantProvider::SetTrioPath(const std::string& a_rTrioPath,
                                   const std::string& a_rPedigreePath,
                                   bool a_bIsTrimBeginningFirst)
{
    
    //Open Pedigree file
    m_pedParser.ParsePedigree(a_rPedigreePath);
    
    CVcfReader temp;
    temp.Open(a_rTrioPath.c_str());
    std::vector<std::string> sampleNames;
    temp.GetSampleNames(sampleNames);
    
    std::vector<std::string> mfcSamplesOrdered;
    
    if(a_rPedigreePath == "")
    {
        mfcSamplesOrdered.push_back("MOTHER");
        mfcSamplesOrdered.push_back("FATHER");
        mfcSamplesOrdered.push_back("CHILD");
    }
    else
    {
        mfcSamplesOrdered = m_pedParser.GetIdsMFC(sampleNames[0], sampleNames[1], sampleNames[2]);
    }
    
    //Open Vcf for child read
    m_vcfReaderChild.Open(a_rTrioPath.c_str());
    m_vcfReaderChild.SelectSample(mfcSamplesOrdered[2]);
    
    //Open Vcf for father read
    m_vcfReaderFather.Open(a_rTrioPath.c_str());
    m_vcfReaderFather.SelectSample(mfcSamplesOrdered[1]);
    
    //Open Vcf for mother read
    m_vcfReaderMother.Open(a_rTrioPath.c_str());
    m_vcfReaderMother.SelectSample(mfcSamplesOrdered[0]);
    
    //Set Trimming decision
    m_bIsTrimBeginningFirst = a_bIsTrimBeginningFirst;
}

bool CVariantProvider::ReadNextChromosomeForSample(std::string& a_rChrName, EVcfName a_uFrom)
{
    //Clear iterators
    m_nMotherItr = 0;
    m_nFatherItr = 0;
    m_nChildItr = 0;
    
    CVcfReader* pVcfReader = 0;
    std::vector<CVariant>* pVariantList = 0;
    
    switch (a_uFrom)
    {
        case eCHILD:
            pVcfReader = &m_vcfReaderChild;
            pVariantList = &m_aChildVariants;
            break;
            
        case eFATHER:
            pVcfReader = &m_vcfReaderFather;
            pVariantList = &m_aFatherVariants;
            break;
            
        case eMOTHER:
            pVcfReader = &m_vcfReaderMother;
            pVariantList = &m_aMotherVariants;
            break;
            
        default:
            break;
    }
    
    //Clear old data
    pVariantList->clear();
    
    bool hasNext = false;
    
    SConfig config;
    config.m_bIsRefOverlap = true;
    config.m_bIsFilterEnabled = false;
    config.m_bTrimBeginningFirst = m_bIsTrimBeginningFirst;
    
    bcf1_t* recordPointer = pVcfReader->GetRecordPointer();
    bcf_hdr_t* headerPointer = pVcfReader->GetHeaderPointer();
    
    CVariant variant;
    int id = 0;
    std::string preChrId = "";
    
    std::vector<CVariant> multiTrimmableVarList;
    
    while(pVcfReader->GetNextRecord(&variant, id, config))
    {
        variant.m_mendelianDecision = GetMendelianDecision(recordPointer, headerPointer);
        
        if(preChrId == "")
        {
            preChrId = variant.m_chrName;
            a_rChrName = variant.m_chrName;
        }
        
        if(preChrId != variant.m_chrName)
        {
            hasNext = true;
            break;
        }
        
        if(!variant.m_bIsNoCall && IsHomRef(variant))
            continue;
        
        else if(IsStructuralVariant(variant, 1000))
            continue;
        
        else if(true == variant.m_bHaveMultipleTrimOption)
        {
            multiTrimmableVarList.push_back(variant);
            id++;
        }
        
        else
        {
            pVariantList->push_back(variant);
            id++;
        }
    }
    
    preChrId = "";
    id = 0;
    
    //Find optimal trimmings for variants where multiple trimming option is available
    FindOptimalTrimmings(multiTrimmableVarList, a_uFrom);
    AppendTrimmedVariants(multiTrimmableVarList, a_uFrom);

    //Sort variants according to the start and original position (and IDs)
    std::sort(pVariantList->begin(), pVariantList->end(), CVariantProvider::CompareVariants);

    return hasNext;
}

bool CVariantProvider::ReadNextChromosome(std::string &a_rChrName)
{
    bool bHasNext = true;

    //Read Next chromosome from father
    bHasNext = ReadNextChromosomeForSample(a_rChrName, eFATHER);
    //Read Next chromosome from mother
    bHasNext = ReadNextChromosomeForSample(a_rChrName, eMOTHER);
    //Read Next chromosome from child
    bHasNext = ReadNextChromosomeForSample(a_rChrName, eCHILD);
    
    return bHasNext;
}


void CVariantProvider::FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom)
{
    if(a_rVariantList.size() == 0)
        return;
    
    std::vector<CVariant>* allVarList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            allVarList = &m_aChildVariants;
            break;
        case eFATHER:
            allVarList = &m_aFatherVariants;
            break;
        case eMOTHER:
            allVarList = &m_aMotherVariants;
            break;
        default:
            allVarList = 0;
            break;
    }
    
    unsigned int varItr[2];
    varItr[0] = 0;
    varItr[1] = 0;
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
    {
        for(int i = 0; i < 2; i++)
        {
            if(a_rVariantList[k].m_alleles[i].m_bIsIgnored || a_rVariantList[k].m_alleles[i].m_bIsTrimmed)
                continue;
            
            //The maximum possible trimming nucleotide count from start and end of each allele
            unsigned int canTrimStart, canTrimEnd;
            a_rVariantList[k].GetMaxTrimStartEnd(i, canTrimStart, canTrimEnd);
            
            std::vector<CVariant> tmpoverlapVariants;
            
            while(varItr[i] < (*allVarList).size() && a_rVariantList[k].m_alleles[i].m_nStartPos > (*allVarList)[varItr[i]].m_nEndPos)
                varItr[i]++;
            
            if(varItr[i] == (*allVarList).size())
                continue;
            
            unsigned int secondItr = varItr[i];
            
            while(secondItr < (*allVarList).size() && a_rVariantList[k].m_alleles[i].m_nEndPos >= (*allVarList)[secondItr].m_nStartPos)
            {
                if(isOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
                             a_rVariantList[k].m_alleles[i].m_nEndPos,
                             (*allVarList)[secondItr].m_nStartPos,
                             (*allVarList)[secondItr].m_nEndPos))
                    tmpoverlapVariants.push_back((*allVarList)[secondItr]);
                secondItr++;
            }
            
            //Trim variants as standard if there is no overlap
            if(tmpoverlapVariants.size() == 0)
                a_rVariantList[k].TrimVariant(i, m_bIsTrimBeginningFirst);
            
            else
            {
                for(unsigned int ovarItr = 0; ovarItr < tmpoverlapVariants.size(); ovarItr++)
                {
                    //Check each allele of overlapping variant
                    for(int tmpItr = 0; tmpItr < 2; tmpItr++)
                    {
                        //If the allele does not overlap, continue
                        if(!isOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
                                      a_rVariantList[k].m_alleles[i].m_nEndPos,
                                      tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos,
                                      tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos))
                            continue;
                        
                        unsigned int overlapStart = std::max(a_rVariantList[k].m_alleles[i].m_nStartPos, tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos);
                        unsigned int overlapEnd = std::min(a_rVariantList[k].m_alleles[i].m_nEndPos, tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos);
                        
                        //Trim from beginning
                        if(a_rVariantList[k].m_alleles[i].m_nStartPos == (int)overlapStart)
                        {
                            if(a_rVariantList[k].m_alleles[i].m_nStartPos + (int)canTrimStart >= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos)
                            {
                                int toClip = overlapEnd - a_rVariantList[k].m_alleles[i].m_nStartPos;
                                if(toClip > 0 && toClip <= (int)canTrimStart)
                                {
                                    a_rVariantList[k].TrimVariant(i, toClip, 0);
                                    canTrimStart -= toClip;
                                    continue;
                                }
                            }
                        }
                        
                        //Trim from beginning or end --
                        else if(a_rVariantList[k].m_alleles[i].m_nEndPos > (int) overlapEnd)
                        {
                            //Try to trim from end
                            if(a_rVariantList[k].m_alleles[i].m_nEndPos - (int)canTrimEnd <= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos)
                            {
                                int toClip = a_rVariantList[k].m_alleles[i].m_nEndPos - overlapStart;
                                if(toClip > 0 && toClip <= (int)canTrimEnd)
                                {
                                    a_rVariantList[k].TrimVariant(i, 0, toClip);
                                    canTrimEnd -= toClip;
                                    continue;
                                }
                            }
                            
                            //Try to trim from front
                            if(a_rVariantList[k].m_alleles[i].m_nStartPos + (int)canTrimStart >= tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nEndPos)
                            {
                                int toClip = overlapEnd - a_rVariantList[k].m_alleles[i].m_nStartPos;
                                if(toClip > 0 && toClip <= (int)canTrimStart)
                                {
                                    a_rVariantList[k].TrimVariant(i, toClip, 0);
                                    canTrimStart -= toClip;
                                    continue;
                                }
                            }
                            
                        }
                        
                        //Trim from end
                        else
                        {
                            if(a_rVariantList[k].m_alleles[i].m_nEndPos - (int)canTrimEnd < tmpoverlapVariants[ovarItr].m_alleles[tmpItr].m_nStartPos)
                            {
                                int toClip = a_rVariantList[k].m_alleles[i].m_nEndPos - overlapStart;
                                if(toClip > 0 && toClip <= (int)canTrimEnd)
                                {
                                    a_rVariantList[k].TrimVariant(i, 0, toClip);
                                    canTrimEnd -= toClip;
                                    continue;
                                }
                            }
                        }
                    }
                }
                
                if(!a_rVariantList[k].m_alleles[i].m_bIsTrimmed)
                    a_rVariantList[k].TrimVariant(i, m_bIsTrimBeginningFirst);
                else
                    a_rVariantList[k].TrimVariant(i, canTrimStart, canTrimEnd);
            }
        }
        
        
        //Update Variant Range
        int minStart = INT_MAX;
        int maxEnd = -1;
        
        for(int i =0; i < a_rVariantList[k].m_nAlleleCount; i++)
        {
            if(!a_rVariantList[k].m_alleles[i].m_bIsIgnored)
            {
                maxEnd = std::max(maxEnd, static_cast<int>(a_rVariantList[k].m_alleles[i].m_nEndPos));
                minStart = std::min(minStart, static_cast<int>(a_rVariantList[k].m_alleles[i].m_nStartPos));
            }
        }
        
        a_rVariantList[k].m_nEndPos = maxEnd == -1 ? a_rVariantList[k].m_nOriginalPos : maxEnd;
        a_rVariantList[k].m_nStartPos = minStart == INT_MAX ? a_rVariantList[k].m_nOriginalPos : minStart;
        
    }
    
}

void CVariantProvider::AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom)
{
    std::vector<CVariant>* variantList;
    
    switch (a_uFrom)
    {
        case eCHILD:
            variantList = &m_aChildVariants;
            break;
        case eFATHER:
            variantList = &m_aFatherVariants;
            break;
        case eMOTHER:
            variantList = &m_aMotherVariants;
            break;
        default:
            variantList = 0;
            break;
    }
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
        (*variantList).push_back(a_rVariantList[k]);
}

bool CVariantProvider::IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const
{
    std::size_t found;
    
    for(int k = 0; k < a_rVariant.m_nAlleleCount; k++)
    {
        std::string allele = a_rVariant.m_alleles[k].m_sequence;
        
        if((int)allele.size() > a_nMaxLength)
            return true;
        found = allele.find('[');
        if(found != std::string::npos)
            return true;
        found = allele.find('<');
        if(found != std::string::npos)
            return true;
        found = allele.find('*');
        if(found != std::string::npos)
            return true;
        found = allele.find('.');
        if(found != std::string::npos)
            return true;
    }
    return false;
}

bool CVariantProvider::IsHomRef(const CVariant& a_rVariant) const
{
    bool res = true;
    
    for(int k = 0; k < a_rVariant.m_nZygotCount; k++)
    {
        if(a_rVariant.m_genotype[k] != 0)
        {
            res = false;
            break;
        }
    }
    return res;
}


void CVariantProvider::GetVariants(const SInterval& a_rInterval,
                                   std::vector<const CVariant*>& a_rMotherVars,
                                   std::vector<const CVariant*>& a_rFatherVars,
                                   std::vector<const CVariant*>& a_rChildVars,
                                   bool& a_rIsViolationOverlapConsistent,
                                   int& a_rTotalViolationCount)
{
    std::vector<EVcfName> violationSamples;
    std::vector<const CVariant*> violationVars;
    a_rIsViolationOverlapConsistent = false;
    
    //Skip irrelevant variants
    while(m_nMotherItr < m_aMotherVariants.size() && m_aMotherVariants[m_nMotherItr].m_nEndPos < a_rInterval.m_nStart)
        m_nMotherItr++;
    while(m_nFatherItr < m_aFatherVariants.size() && m_aFatherVariants[m_nFatherItr].m_nEndPos < a_rInterval.m_nStart)
        m_nFatherItr++;
    while(m_nChildItr < m_aChildVariants.size() && m_aChildVariants[m_nChildItr].m_nEndPos < a_rInterval.m_nStart)
        m_nChildItr++;

    //Temporary iterators to not lose starting indexes
    int tmpFatherItr = m_nFatherItr;
    int tmpMotherItr = m_nMotherItr;
    int tmpChildItr = m_nChildItr;
    
    
    //Push all consistent variants in given interval
    while(tmpMotherItr < m_aMotherVariants.size() && m_aMotherVariants[tmpMotherItr].m_nStartPos < a_rInterval.m_nEnd)
    {
        if(m_aMotherVariants[tmpMotherItr].m_mendelianDecision == eCONSISTENT
           &&
           true == isOverlap(a_rInterval.m_nStart, a_rInterval.m_nEnd, m_aMotherVariants[tmpMotherItr].m_nStartPos, m_aMotherVariants[tmpMotherItr].m_nEndPos))
        {
            a_rMotherVars.push_back(&m_aMotherVariants[tmpMotherItr]);
        }
        else if(m_aMotherVariants[tmpMotherItr].m_mendelianDecision == eVIOLATION)
        {
            violationVars.push_back(&m_aMotherVariants[tmpMotherItr]);
            violationSamples.push_back(eMOTHER);
        }
        
        tmpMotherItr++;
    }
    
    //Push all consistent variants in given interval
    while(tmpFatherItr < m_aFatherVariants.size() && m_aFatherVariants[tmpFatherItr].m_nStartPos < a_rInterval.m_nEnd)
    {
        if(m_aFatherVariants[tmpFatherItr].m_mendelianDecision == eCONSISTENT
           &&
           true == isOverlap(a_rInterval.m_nStart, a_rInterval.m_nEnd, m_aFatherVariants[tmpFatherItr].m_nStartPos, m_aFatherVariants[tmpFatherItr].m_nEndPos))
        {
            a_rFatherVars.push_back(&m_aFatherVariants[tmpFatherItr]);
        }
        else if(m_aFatherVariants[tmpFatherItr].m_mendelianDecision == eVIOLATION)
        {
            violationVars.push_back(&m_aFatherVariants[tmpFatherItr]);
            violationSamples.push_back(eFATHER);
        }
        
        tmpFatherItr++;
    }

    //Push all consistent variants in given interval
    while(tmpChildItr < m_aChildVariants.size() && m_aChildVariants[tmpChildItr].m_nStartPos < a_rInterval.m_nEnd)
    {
        if(m_aChildVariants[tmpChildItr].m_mendelianDecision == eCONSISTENT
           &&
           true == isOverlap(a_rInterval.m_nStart, a_rInterval.m_nEnd, m_aChildVariants[tmpChildItr].m_nStartPos, m_aChildVariants[tmpChildItr].m_nEndPos))
        {
            a_rChildVars.push_back(&m_aChildVariants[tmpChildItr]);
        }
        
        else if(m_aChildVariants[tmpChildItr].m_mendelianDecision == eVIOLATION)
        {
            violationVars.push_back(&m_aChildVariants[tmpChildItr]);
            violationSamples.push_back(eCHILD);
        }

        tmpChildItr++;
    }
    
    //Check if violations overlap with consistent variants
    for(int k = 0; k < violationVars.size(); k++)
    {
        for(int m = 0; m < a_rMotherVars.size(); m++)
        {
            if(true == isOverlap(a_rMotherVars[m]->m_nStartPos, a_rMotherVars[m]->m_nEndPos, violationVars[k]->m_nStartPos, violationVars[k]->m_nEndPos))
            {
                a_rIsViolationOverlapConsistent = true;
                if(violationSamples[k] == eCHILD)
                    
                break;
            }
        }
        
        if(true == a_rIsViolationOverlapConsistent)
            break;
        
        for(int m = 0; m < a_rFatherVars.size(); m++)
        {
            if(true == isOverlap(a_rFatherVars[m]->m_nStartPos, a_rFatherVars[m]->m_nEndPos, violationVars[k]->m_nStartPos, violationVars[k]->m_nEndPos))
            {
                a_rIsViolationOverlapConsistent = true;
                break;
            }
        }
        
        if(true == a_rIsViolationOverlapConsistent)
            break;

        for(int m = 0; m < a_rChildVars.size(); m++)
        {
            if(true == isOverlap(a_rChildVars[m]->m_nStartPos, a_rChildVars[m]->m_nEndPos, violationVars[k]->m_nStartPos, violationVars[k]->m_nEndPos))
            {
                a_rIsViolationOverlapConsistent = true;
                break;
            }
        }
    }
    
    a_rTotalViolationCount = (int)(violationVars.size());
    
}

EMendelianDecision CVariantProvider::GetMendelianDecision(bcf1_t *a_pRecordPtr, bcf_hdr_t *a_pHeaderPtr)
{
    EMendelianDecision decision;
    
    int* mdarr = NULL;
    int nmdarr = 0;
    bcf_get_info_int32(a_pHeaderPtr, a_pRecordPtr,"MD", &mdarr, &nmdarr);
    
    switch (mdarr[0])
    {
        case 0:
            decision = eUNKNOWN;
            break;
        case 1:
            decision = eCONSISTENT;
            break;
        case 2:
            decision = eVIOLATION;
            break;
        case 3:
            decision = eNOCALL_PARENT;
            break;
        case 4:
            decision = eNOCALL_CHILD;
            break;
        default:
            decision = eUNKNOWN;
            break;
    }
    
    delete[] mdarr;
    return decision;
}

void CVariantProvider::GetViolationIntervals(const std::vector<SInterval>& a_rIntervalList, std::vector<SInterval>& a_rViolationInterval)
{
    int motherItr = 0;
    int fatherItr = 0;
    int childItr = 0;
    
    for(int k = 0; k < a_rIntervalList.size(); k++)
    {
        bool bIsViolationInterval = false;
        
        while(motherItr < m_aMotherVariants.size() && m_aMotherVariants[motherItr].m_nEndPos <= a_rIntervalList[k].m_nEnd)
        {
            if(m_aMotherVariants[motherItr].m_mendelianDecision == eVIOLATION
               &&
               true == isOverlap(a_rIntervalList[k].m_nStart, a_rIntervalList[k].m_nEnd, m_aMotherVariants[motherItr].m_nStartPos, m_aMotherVariants[motherItr].m_nEndPos))
            {
                a_rViolationInterval.push_back(a_rIntervalList[k]);
                bIsViolationInterval = true;
                motherItr++;
                break;
            }
            motherItr++;
        }
        
        if(true == bIsViolationInterval)
            continue;
        
        while(fatherItr < m_aFatherVariants.size() && m_aFatherVariants[fatherItr].m_nEndPos <= a_rIntervalList[k].m_nEnd)
        {
            if(m_aFatherVariants[fatherItr].m_mendelianDecision == eVIOLATION
               &&
               true == isOverlap(a_rIntervalList[k].m_nStart, a_rIntervalList[k].m_nEnd, m_aFatherVariants[fatherItr].m_nStartPos, m_aFatherVariants[fatherItr].m_nEndPos))
            {
                a_rViolationInterval.push_back(a_rIntervalList[k]);
                bIsViolationInterval = true;
                fatherItr++;
                break;
            }
            fatherItr++;
        }

        if(true == bIsViolationInterval)
            continue;
        
        while(childItr < m_aChildVariants.size() && m_aChildVariants[childItr].m_mendelianDecision == eVIOLATION
              &&
              m_aChildVariants[childItr].m_nEndPos <= a_rIntervalList[k].m_nEnd)
        {
            if(true == isOverlap(a_rIntervalList[k].m_nStart, a_rIntervalList[k].m_nEnd, m_aChildVariants[childItr].m_nStartPos, m_aChildVariants[childItr].m_nEndPos))
            {
                a_rViolationInterval.push_back(a_rIntervalList[k]);
                childItr++;
                break;
            }
            childItr++;
        }
    }
}


