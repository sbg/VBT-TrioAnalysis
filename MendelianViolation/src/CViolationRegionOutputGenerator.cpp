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
 *  CViolationRegionOutputGenerator.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke.Toptas on 4/18/18.
 *
 */

#include "CViolationRegionOutputGenerator.h"
#include "CMendelianVariantProvider.h"
#include "CMendelianResultLog.h"
#include "SChrIdTriplet.h"
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace mendelian;

CViolationRegionOutputGenerator::CViolationRegionOutputGenerator(const std::vector<core::CPath>& a_aBestPathsFatherChildGT,
                                const std::vector<core::CPath>& a_aBestPathsFatherChildAM,
                                const std::vector<core::CPath>& a_aBestPathsMotherChildGT,
                                const std::vector<core::CPath>& a_aBestPathsMotherChildAM,
                                const CMendelianVariantProvider& a_rProvider,
                                CMendelianResultLog& a_rResultLog) :
m_aBestPathsFatherChildGT(a_aBestPathsFatherChildGT),
m_aBestPathsFatherChildAM(a_aBestPathsFatherChildAM),
m_aBestPathsMotherChildGT(a_aBestPathsMotherChildGT),
m_aBestPathsMotherChildAM(a_aBestPathsMotherChildAM),
m_provider(a_rProvider),
m_resultLog(a_rResultLog)
{
    m_bIsOutputEnabled = false;
}

void CViolationRegionOutputGenerator::OpenBedForWrite(const std::string& a_rFilePath)
{
    m_bIsOutputEnabled = true;
    m_outputBEDfile.open(a_rFilePath);
}



void CViolationRegionOutputGenerator::GenerateSyncPoints(const SChrIdTriplet& a_triplet, std::vector<int>& intersectedSyncPoints)
{
    //Generate the sync point list
    std::vector<int> syncPointsFather(m_aBestPathsFatherChildGT[a_triplet.m_nTripleIndex].m_aSyncPointList);
    std::vector<int> syncPointsFatherAM(m_aBestPathsFatherChildAM[a_triplet.m_nTripleIndex].m_aSyncPointList);
    //Concat AM + GT sync points
    syncPointsFather.insert(syncPointsFather.end(), syncPointsFatherAM.begin(), syncPointsFatherAM.end());
    //Sort sync points
    std::sort(syncPointsFather.begin(), syncPointsFather.end());
    syncPointsFather.erase( unique( syncPointsFather.begin(), syncPointsFather.end() ), syncPointsFather.end());
    
    
    //Generate the sync point list
    std::vector<int> syncPointsMother(m_aBestPathsMotherChildGT[a_triplet.m_nTripleIndex].m_aSyncPointList);
    std::vector<int> syncPointsMotherAM(m_aBestPathsMotherChildAM[a_triplet.m_nTripleIndex].m_aSyncPointList);
    //Concat AM + GT sync points
    syncPointsMother.insert(syncPointsMother.end(), syncPointsMotherAM.begin(), syncPointsMotherAM.end());
    //Sort sync points
    std::sort(syncPointsMother.begin(), syncPointsMother.end());
    syncPointsMother.erase( unique( syncPointsMother.begin(), syncPointsMother.end() ), syncPointsMother.end());
    
    //Take the intersection of mother and father side syncpoints
    std::set_intersection(syncPointsMother.begin(),syncPointsMother.end(), syncPointsFather.begin(), syncPointsFather.end(), back_inserter(intersectedSyncPoints));
    
    //Add the first and last sync point
    intersectedSyncPoints.insert(intersectedSyncPoints.begin(), 0);
    intersectedSyncPoints.push_back(m_provider.GetContig(a_triplet.m_chrName).length);
    
    //Get variant list
    int syncPointItr = 1;
    std::vector<const CVariant*> motherVariants = m_provider.GetVariantList(eMOTHER, a_triplet.m_nMid);
    std::vector<const CVariant*> fatherVariants = m_provider.GetVariantList(eFATHER, a_triplet.m_nFid);
    std::vector<const CVariant*> childVariants = m_provider.GetVariantList(eCHILD, a_triplet.m_nCid);
    
    //Merge consecutive regions if a mother variant overlaps with boundaries of the region
    for(int k = 0; k < (int)motherVariants.size(); k++)
    {
        while(intersectedSyncPoints[syncPointItr] <= motherVariants[k]->m_nStartPos)
            syncPointItr++;
        
        while(motherVariants[k]->m_nEndPos > intersectedSyncPoints[syncPointItr])
            intersectedSyncPoints.erase(intersectedSyncPoints.begin() + syncPointItr);
    }
    
    //Merge consecutive regions if a father variant overlaps with boundaries of the region
    syncPointItr = 1;
    for(int k = 0; k < (int)fatherVariants.size(); k++)
    {
        while(intersectedSyncPoints[syncPointItr] <= fatherVariants[k]->m_nStartPos)
            syncPointItr++;
        
        while(fatherVariants[k]->m_nEndPos > intersectedSyncPoints[syncPointItr])
            intersectedSyncPoints.erase(intersectedSyncPoints.begin() + syncPointItr);
    }
    
    //Merge consecutive regions if a child variant overlaps with boundaries of the region
    syncPointItr = 1;
    for(int k = 0; k < (int)childVariants.size(); k++)
    {
        while(intersectedSyncPoints[syncPointItr] <= childVariants[k]->m_nStartPos)
            syncPointItr++;
        
        while(childVariants[k]->m_nEndPos > intersectedSyncPoints[syncPointItr])
            intersectedSyncPoints.erase(intersectedSyncPoints.begin() + syncPointItr);
    }
    
}


void CViolationRegionOutputGenerator::GenerateViolationRegions(SChrIdTriplet& a_triplet, std::vector<SVcfRecord>& a_rRecordList, std::vector<EMendelianDecision>& a_rDecisionList)
{
    std::vector<int> intersectionlist;
    GenerateSyncPoints(a_triplet, intersectionlist);
 
    int recordItr = 0;
    
    int violationRegionCount = 0;
    int nocallParentCount = 0;
    int nocallChildCount = 0;
    int consistentRegionCount = 0;
    
    for(int k = 1; k < (int)intersectionlist.size(); k++)
    {
        bool bFoundViolationFlag = false;
        bool bFoundNocallParentFlag = false;
        bool bFoundNocallChildFlag = false;
        bool bHasVariant = false;
        
        while(recordItr < (int)a_rRecordList.size() && a_rRecordList[recordItr].left < intersectionlist[k-1])
            recordItr++;
        
        
        while(recordItr < (int)a_rRecordList.size() && a_rRecordList[recordItr].right <= intersectionlist[k])
        {
            bHasVariant = true;
            if(a_rDecisionList[recordItr] == eViolation)
                bFoundViolationFlag = true;
            else if(a_rDecisionList[recordItr] == eNoCallParent)
                bFoundNocallParentFlag = true;
            else if(a_rDecisionList[recordItr] == eNoCallChild)
                bFoundNocallChildFlag = true;
            recordItr++;
        }
        
        if(!bFoundNocallChildFlag && !bFoundNocallParentFlag && bFoundViolationFlag)
        {
            if(true == m_bIsOutputEnabled)
                m_outputBEDfile << a_triplet.m_chrName << "\t" << intersectionlist[k-1] << "\t" << intersectionlist[k] << std::endl;
            violationRegionCount++;
        }
        else if(bFoundNocallChildFlag)
            nocallChildCount++;
        else if(bFoundNocallParentFlag)
            nocallParentCount++;
        else if(bHasVariant)
            consistentRegionCount++;
    }
    
    m_resultLog.LogRegionBasedCounts(consistentRegionCount, violationRegionCount, nocallParentCount, nocallChildCount);
}

void CViolationRegionOutputGenerator::CloseBed()
{
    m_outputBEDfile.close();
}


