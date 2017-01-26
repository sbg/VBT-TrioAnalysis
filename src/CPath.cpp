//  Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
//  Copyright (c) 2014. Real Time Genomics Limited.
//
//  CPath.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//

#include "CPath.h"
#include <iostream>
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

CPath& CPath::Include(EVcfName a_nVCF, const COrientedVariant& a_rVariant, int a_nVariantIndex)
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

int CPath::AddVariant(CPathContainer* a_pPathList,
                      EVcfName a_nVcfName,
                      const std::vector<const CVariant*>& a_pVariantList,
                      const std::vector<const COrientedVariant*>& a_pOVariantList,
                      int a_nVariantIndex,
                      bool a_bIsGenotypeMatch)
{
    int pathCount = 0;

    if(a_bIsGenotypeMatch)
    {
        const CVariant* pNextVariant = a_pVariantList[a_nVariantIndex];
        const COrientedVariant* Ovar1 = a_pOVariantList[2* a_nVariantIndex];
        const COrientedVariant* Ovar2 = a_pOVariantList[2* a_nVariantIndex + 1];
        
        if (true == InSync())
        {
            m_nCSinceSync = 0;
            m_nBSinceSync = 0;
            
            // Create a path extension that excludes this variant
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
            a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
            pathCount++;
            
            // Create a path extension that includes this variant in the possible phases
            if (!pNextVariant->IsHeterozygous())
            {
                a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
                const CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                //Make sure variant is not overlap with the previous one
                if(p->IsNew(*Ovar1))
                {
                    a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                    pathCount++;
                }
            }
            else
            {
                //Include with ordered genotype
                a_pPathList[pathCount].m_pPath =  std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
                CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                //Make sure variant is not overlap with the previous one
                if(p->IsNew(*Ovar1))
                {
                    a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
                    pathCount++;
                }
                //Include with unordered genotype
                a_pPathList[pathCount].m_pPath =  std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
                p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                //Make sure variant is not overlap with the previous one
                if(p->IsNew(*Ovar2))
                {
                    a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar2, a_nVariantIndex);
                    pathCount++;
                }
            }
            
        }
        
        else
        {
            // Create a path extension that excludes this variant
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
            a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
            pathCount++;
            
            // Create a path extension that includes this variant in the possible phases
            if (!pNextVariant->IsHeterozygous())
            {
                a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
                CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                //Make sure variant is not overlap with the previous one
                if(p->IsNew(*Ovar1))
                {
                    a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovar1, a_nVariantIndex);
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
    }
    
    else
    {
        const CVariant* pNextVariant = a_pVariantList[a_nVariantIndex];
        const COrientedVariant* Ovars[] = {a_pOVariantList[2* a_nVariantIndex], a_pOVariantList[2* a_nVariantIndex + 1]};
        
        if(true == InSync())
        {
            m_nCSinceSync = 0;
            m_nBSinceSync = 0;
            
            // Create a path extension that excludes this variant
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
            a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
            pathCount++;
            
            for(int k = 0; k < pNextVariant->m_nZygotCount; k++)
            {
                if(pNextVariant->m_genotype[k] != 0)
                {
                    a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this, m_calledSemiPath.GetPosition()));
                    const CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                    //Make sure variant is not overlap with the previous one
                    if(p->IsNew(*Ovars[k]))
                    {
                        a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovars[k], a_nVariantIndex);
                        pathCount++;
                    }
                }
            }
            
        }
        
        else
        {
            // Create a path extension that excludes this variant
            a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
            a_pPathList[pathCount].m_pPath->Exclude(a_nVcfName, *pNextVariant, a_nVariantIndex);
            pathCount++;

            for(int k = 0; k < pNextVariant->m_nZygotCount; k++)
            {
                if(pNextVariant->m_genotype[k] != 0)
                {
                    a_pPathList[pathCount].m_pPath = std::shared_ptr<CPath>(new CPath(*this));
                    CSemiPath* p = a_nVcfName == eBASE ? &a_pPathList[pathCount].m_pPath->m_baseSemiPath : &a_pPathList[pathCount].m_pPath->m_calledSemiPath;
                    //Make sure variant is not overlap with the previous one
                    if(p->IsNew(*Ovars[k]))
                    {
                        a_pPathList[pathCount].m_pPath->Include(a_nVcfName, *Ovars[k], a_nVariantIndex);
                        pathCount++;
                    }
                }
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

void CPath::ClearIncludedVariants()
{
    m_calledSemiPath.ClearIncludedVariants();
    m_baseSemiPath.ClearIncludedVariants();
}

void CPath::AddIncludedVariants(std::vector<const COrientedVariant*>& a_rIncludedVarListCalled, std::vector<const COrientedVariant*>& a_rIncludedVarListBase)
{
    m_calledSemiPath.AddIncludedVariants(a_rIncludedVarListCalled);
    m_baseSemiPath.AddIncludedVariants(a_rIncludedVarListBase);
}

void CPath::ClearExcludedVariants()
{
    m_calledSemiPath.ClearExcludedVariants();
    m_baseSemiPath.ClearExcludedVariants();
}

void CPath::AddExcludedVariants(std::vector<int>& a_rExcludedVarListCalled, std::vector<int>& a_rExcludedVarListBase)
{
    m_baseSemiPath.AddExcludedVariants(a_rExcludedVarListBase);
    m_calledSemiPath.AddExcludedVariants(a_rExcludedVarListCalled);
}

void CPath::ClearSyncPointList()
{
    m_aSyncPointList.clear();
}

void CPath::AddSyncPointList(std::vector<int>& a_rSyncPointArray)
{
    m_aSyncPointList = a_rSyncPointArray;
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
