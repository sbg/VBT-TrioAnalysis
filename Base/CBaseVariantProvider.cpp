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
 *  CBaseVariantProvider.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 12/7/17.
 *
 */

#include "CBaseVariantProvider.h"
#include "CVariant.h"
#include "COrientedVariant.h"
#include "CUtils.h"

void CBaseVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
    {
        if(pVar->m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pVar->m_variantStatus != eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pVar-> m_variantStatus == eNOT_ASSESSED)
            pVar->m_variantStatus = a_status;
        else
            continue;
    }
}

void CBaseVariantProvider::SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const core::COrientedVariant* pOVar : a_rVariantList)
    {
        if(pOVar->GetVariant().m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pOVar->GetVariant().m_variantStatus != eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pOVar->GetVariant().m_variantStatus == eNOT_ASSESSED)
            pOVar->GetVariant().m_variantStatus = a_status;
        else
            continue;
    }
}

bool CBaseVariantProvider::ReadContig(std::string a_chrId, SContig& a_rContig)
{
    return m_referenceFasta.FetchNewChromosome(a_chrId, a_rContig);
}


void CBaseVariantProvider::FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, std::vector<std::vector<CVariant>>* a_pAllVarList, const SConfig& a_rConfig)
{
    if(a_rVariantList.size() == 0)
        return;
    
    unsigned int varItr[2];
    varItr[0] = 0;
    varItr[1] = 0;
    
    int currentChrId = a_rVariantList[0].m_nChrId;
    
    for(unsigned int k = 0; k < a_rVariantList.size(); k++)
    {
        for(int i = 0; i < a_rVariantList[k].m_nZygotCount; i++)
        {
            if(a_rVariantList[k].m_alleles[i].m_bIsIgnored || a_rVariantList[k].m_alleles[i].m_bIsTrimmed)
                continue;
            
            //The maximum possible trimming nucleotide count from start and end of each allele
            unsigned int canTrimStart, canTrimEnd;
            a_rVariantList[k].GetMaxTrimStartEnd(i, canTrimStart, canTrimEnd);
            
            std::vector<CVariant> tmpoverlapVariants;
            
            if(currentChrId != a_rVariantList[k].m_nChrId)
            {
                varItr[0] = 0;
                varItr[1] = 0;
                currentChrId = a_rVariantList[k].m_nChrId;
            }
            
            while(varItr[i] < (*a_pAllVarList)[currentChrId].size() && a_rVariantList[k].m_alleles[i].m_nStartPos > (*a_pAllVarList)[currentChrId][varItr[i]].m_nEndPos)
                varItr[i]++;
            
            if(varItr[i] == (*a_pAllVarList)[currentChrId].size())
                continue;
            
            unsigned int secondItr = varItr[i];
            
            while(secondItr < (*a_pAllVarList)[currentChrId].size() && a_rVariantList[k].m_alleles[i].m_nEndPos >= (*a_pAllVarList)[currentChrId][secondItr].m_nStartPos)
            {
                if(CUtils::IsOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
                                     a_rVariantList[k].m_alleles[i].m_nEndPos,
                                     (*a_pAllVarList)[currentChrId][secondItr].m_nStartPos,
                                     (*a_pAllVarList)[currentChrId][secondItr].m_nEndPos))
                    tmpoverlapVariants.push_back((*a_pAllVarList)[currentChrId][secondItr]);
                secondItr++;
            }
            
            //Trim variants as standard if there is no overlap
            if(tmpoverlapVariants.size() == 0)
                a_rVariantList[k].TrimVariant(i, a_rConfig.m_bTrimBeginningFirst);
            
            else
            {
                for(unsigned int ovarItr = 0; ovarItr < tmpoverlapVariants.size(); ovarItr++)
                {
                    //Check each allele of overlapping variant
                    for(int tmpItr = 0; tmpItr < tmpoverlapVariants[ovarItr].m_nZygotCount; tmpItr++)
                    {
                        //If the allele does not overlap, continue
                        if(!CUtils::IsOverlap(a_rVariantList[k].m_alleles[i].m_nStartPos,
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
                    a_rVariantList[k].TrimVariant(i, a_rConfig.m_bTrimBeginningFirst);
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
