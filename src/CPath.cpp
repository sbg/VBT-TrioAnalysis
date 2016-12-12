//  Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
//  Copyright (c) 2014. Real Time Genomics Limited.
//
//  CPath.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CPath.h"
#include <iostream>
#include "CSyncPoint.h"
#include "CVariantProvider.h"

CPath::CPath()
{}

CPath::CPath(const char* a_aRefSequence, int a_nRefSize)
: m_baseSemiPath(a_aRefSequence, a_nRefSize, eBASE),
  m_calledSemiPath(a_aRefSequence, a_nRefSize, eCALLED),
  m_nCSinceSync(0),
  m_nBSinceSync(0)
{
    m_nPathId = -1;
}

CPath::CPath(const CPath& a_rObj)
: m_baseSemiPath(a_rObj.m_baseSemiPath),
  m_calledSemiPath(a_rObj.m_calledSemiPath)
{
    m_aSyncPointList = std::vector<int>(a_rObj.m_aSyncPointList);
    m_nCSinceSync = a_rObj.m_nCSinceSync;
    m_nBSinceSync = a_rObj.m_nBSinceSync;
    
    m_nPathId = a_rObj.m_nPathId;
}

CPath::CPath(const CPath& a_rObj, int  a_nSyncPointToPush)
: m_baseSemiPath(a_rObj.m_baseSemiPath),
  m_calledSemiPath(a_rObj.m_calledSemiPath)
{
    m_aSyncPointList = std::vector<int>(a_rObj.m_aSyncPointList);
    m_aSyncPointList.push_back(a_nSyncPointToPush);
    m_nCSinceSync = a_rObj.m_nCSinceSync;
    m_nBSinceSync = a_rObj.m_nBSinceSync;
    
    m_nPathId = a_rObj.m_nPathId;
}

bool CPath::IsEqual(const CPath& a_rObj) const
{
    return (CompareTo(a_rObj) == 0);
}

int CPath::CompareTo(const CPath& a_rObj) const
{
    int res = m_calledSemiPath.CompareTo(a_rObj.m_calledSemiPath);
    
    if(res != 0)
        return res;
    else
        return m_baseSemiPath.CompareTo(a_rObj.m_baseSemiPath);
}

CPath& CPath::Exclude(EVcfName a_nVCF, const CVariant& a_rVariant, int a_nVariantIndex)
{      
    switch(a_nVCF)
    {
        case eBASE:
            m_baseSemiPath.ExcludeVariant(a_rVariant, a_nVariantIndex);
            break;
        case eCALLED:
            m_calledSemiPath.ExcludeVariant(a_rVariant, a_nVariantIndex);
        break;
    }

    return *this;
}

CPath& CPath::Include(EVcfName a_nVCF, COrientedVariant& a_rVariant, int a_nVariantIndex)
{
    switch(a_nVCF)
    {
        case eBASE:
            m_baseSemiPath.IncludeVariant(a_rVariant, a_nVariantIndex);
            m_nBSinceSync++;
            break;
        case eCALLED:
            m_calledSemiPath.IncludeVariant(a_rVariant, a_nVariantIndex);
            m_nCSinceSync++;
            break;
    }

    return *this;
}

int CPath::AddVariant(CPathContainer* a_pPathList, EVcfName a_nVcfName, CVariantProvider* a_pVariantProvider, int a_nVariantIndex, int a_nChromosomeId)
{
    int pathCount = 0;
    
    const CVariant* pNextVariant = a_pVariantProvider->GetVariant(a_nVcfName, a_nChromosomeId, a_nVariantIndex);
    COrientedVariant* Ovar1 = a_pVariantProvider->GetOrientedVariant(a_nVcfName, a_nChromosomeId, a_nVariantIndex, true);
    COrientedVariant* Ovar2 = a_pVariantProvider->GetOrientedVariant(a_nVcfName, a_nChromosomeId, a_nVariantIndex, false);
    
    if (true == InSync())
    {
        m_nCSinceSync = 0;
        m_nBSinceSync = 0;

        // Create a path extension that excludes this variant
        a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
        a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
        //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +150;
        pathCount++;
    
        // Create a path extension that includes this variant in the possible phases
        if (!pNextVariant->IsHeterozygous())
        {
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
            const CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar1(a_rVariant, true);
            //Make sure variant is not overlap with the previous one
            if(p->IsNew(*Ovar1))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +300;
                pathCount++;
            }
        }
        else
        {
            //Include with ordered genotype
            a_pPathList[pathCount].m_pPath =  std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
            CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar1(a_rVariant, true);
            //Make sure variant is not overlap with the previous one
            if(p->IsNew(*Ovar1))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +300;
                pathCount++;
            }
            
            //Include with unordered genotype
            a_pPathList[pathCount].m_pPath =  std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
            p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar2(a_rVariant, false);
            //Make sure variant is not overlap with the previous one
            if(p->IsNew(*Ovar2))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar2, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +450;
                pathCount++;
            }
        }
        
    }
    
    else
    {
        // Create a path extension that excludes this variant
        a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
        a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
        //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +150;
        pathCount++;
        
        // Create a path extension that includes this variant in the possible phases
        if (!pNextVariant->IsHeterozygous())
        {
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
            CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar1(a_rVariant, true);
            //Make sure variant is not overlap with the previous one
            if(p->IsNew(*Ovar1))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId +300;
                pathCount++;
            }
        }
        else
        {
            //Include with ordered genotype
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
            CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar1(a_rVariant, true);
            if(p->IsNew(*Ovar1))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId + 300;
                pathCount++;
            }
            
            //Include with unordered genotype
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
            p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
            //COrientedVariant Ovar2(a_rVariant, false);
            if(p->IsNew(*Ovar2))
            {
                a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar2, a_nVariantIndex);
                //a_pPathList[pathCount].m_pPath->m_nPathId = this->m_nPathId + 450;
                pathCount++;
            }
        }
    }
    
    
    return pathCount;
}

bool CPath::InSync() const
{
    if (m_calledSemiPath.CompareHaplotypePositions() != 0)
      return false;
    else if (m_baseSemiPath.CompareHaplotypePositions() != 0)
      return false;
    else if (m_calledSemiPath.GetPosition() != m_baseSemiPath.GetPosition()) 
      return false;
    else if (m_calledSemiPath.GetPosition() < m_calledSemiPath.GetVariantEndPosition())
      return false;
    else if (m_baseSemiPath.GetPosition() < m_baseSemiPath.GetVariantEndPosition())
      return false;
    else if (!m_calledSemiPath.IsOnTemplate())
      return false;
    else if (!m_baseSemiPath.IsOnTemplate())
      return false;
    else
      return true;
}


bool CPath::HasFinished() const
{
    return m_baseSemiPath.HasFinished() && m_calledSemiPath.HasFinished();
}

bool CPath::HasNoOperation() const
{
    return (m_nBSinceSync == 0 && m_nCSinceSync > 0) || (m_nCSinceSync == 0 && m_nBSinceSync > 0);
}


void CPath::Step()
{
    if(m_calledSemiPath.CompareHaplotypePositions() > 0)
    {
        //make haplotype B catch up haplotype A
        m_calledSemiPath.StepHaplotypeB();
        m_baseSemiPath.StepHaplotypeB();
    }
    else if(m_calledSemiPath.CompareHaplotypePositions() < 0)
    {
        //make haplotype A catch up haplotype B
        m_calledSemiPath.StepHaplotypeA();
        m_baseSemiPath.StepHaplotypeA();
    }
    else
    {
        //Step both
        m_calledSemiPath.StepHaplotypeA();
        m_calledSemiPath.StepHaplotypeB();
        m_baseSemiPath.StepHaplotypeA();
        m_baseSemiPath.StepHaplotypeB();
    }

}

void CPath::MoveForward(int a_nPosition)
{
    m_calledSemiPath.MoveForward(a_nPosition);
    m_baseSemiPath.MoveForward(a_nPosition);
}

bool CPath::Matches()
{
    return m_calledSemiPath.Matches(m_baseSemiPath);
}


std::vector<CSyncPoint> GetSyncPointsList(const std::vector<int>& syncpoints,const std::vector<COrientedVariant*>& baseLine, const std::vector<COrientedVariant*>& called)
{
    std::vector<CSyncPoint> list;
    int basePos = 0;
    int callPos = 0;
    
    for (int loc : syncpoints)
    {
        int baseLineCount = 0;
        int calledCount = 0;
        while (basePos < baseLine.size() && baseLine[basePos]->GetVariant().GetStart() <= loc)
        {
            baseLineCount++;
            basePos++;
        }
        
        while (callPos < called.size() && called[callPos]->GetVariant().GetStart() <= loc)
        {
            calledCount++;
            callPos++;
        }
        
        list.push_back(CSyncPoint(loc, calledCount, baseLineCount));
    }
    return list;
}


void CPath::CalculateWeights()
{
    assert(m_aSyncPointList.size() >= 1);
    std::vector<CSyncPoint> syncpoints = GetSyncPointsList(m_aSyncPointList,m_baseSemiPath.GetIncludedVariants(), m_calledSemiPath.GetIncludedVariants());
    assert(syncpoints.size() == m_aSyncPointList.size());
    
    std::vector<CSyncPoint>::iterator it = syncpoints.begin();
    int syncStart = 0;
    CSyncPoint syncpoint = *it;

    for (const COrientedVariant* v : m_calledSemiPath.GetIncludedVariants())
    {
        while (syncpoint.GetPosition() < v->GetVariant().GetStart() && it != syncpoints.end())
        {
            // Jump to sync point entry containing this variant
            syncStart = syncpoint.GetPosition();
            syncpoint = *it;
            it++;
        }
        if (syncpoint.m_nBaselineTPCount == 0)
        {
            //A rare condition where variants cancels out each other
            v->SetWeight(0);
        }
        else
        {
            // Inside count could still be 0 if the only TP are outside evaluation regions
            v->SetWeight( static_cast<double>(syncpoint.m_nBaselineTPCount) / static_cast<double>(syncpoint.m_nCalledTPCount));
        }
    }
}


void CPath::Print() const
{
    std::cout << "-----" << std::endl;
    std::cout << "Sync:" << m_nBSinceSync << " " << m_nCSinceSync << std::endl;
    std::cout << "Sync Points: ";
    for (int i = (int)m_aSyncPointList.size()-1; i>=0; i--)
        std::cout << m_aSyncPointList[i] << " ";
    std::cout<< std::endl;
    std::cout << "--Base Semipath--" << std::endl;
    m_baseSemiPath.Print();
    std::cout << "--Called Semipath--" << std::endl;
    m_calledSemiPath.Print();    
}
