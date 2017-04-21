//  Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
//  Copyright (c) 2016 Seven Bridges Genomics
//
//  CPathReplay.cpp
//  VCFComparison
//
//  Created by Berke Toptas
//

#include "CPathReplay.h"
#include <iostream>
#include <vector>
#include <cassert>
#include "CPath.h"


CPathReplay::CPathReplay(std::vector<const CVariant*>& a_aVarListBase,
            std::vector<const CVariant*>& a_aVarListCalled,
            std::vector<const COrientedVariant*>& a_aOVarListBase,
            std::vector<const COrientedVariant*>& a_aOvarlistCalled)
: m_aVariantListBase(a_aVarListBase),
  m_aVariantListCalled(a_aVarListCalled),
  m_aOrientedVariantListBase(a_aOVarListBase),
 m_aOrientedVariantListCalled(a_aOvarlistCalled)
{
    m_nMaxPathSize = DEFAULT_MAX_PATH_SIZE;
    m_nMaxIterationCount = DEFAULT_MAX_ITERATION_SIZE;
}

void CPathReplay::SetMaxPathAndIteration(int a_nMaxPathSize, int a_nMaxIterationCount)
{
    m_nMaxIterationCount = a_nMaxIterationCount;
    m_nMaxPathSize = a_nMaxPathSize;
}


CPath CPathReplay::FindBestPath(SContig a_contig, bool a_bIsGenotypeMatch)
{
    CPathContainer initialPath(a_contig.m_pRefSeq, a_contig.m_nRefLength);
    m_pathList.Add(initialPath);
    CPathContainer best(initialPath);
    CPathContainer lastSyncPath;
    int maxPaths = 0;
    int currentIterations = 0;
    int currentMaxIterations = 0;
    int currentMax = 0;
    m_nCurrentPosition = 0;
    int lastSyncPos = 0;
    int complexRegionCount = 0;
    int totalSkippedVariantCount = 0;
    
    CPathContainer processedPath;
    
    while(!m_pathList.Empty())
    {
        currentMax = std::max(currentMax, m_pathList.Size());
        currentMaxIterations = std::max(currentMaxIterations, currentIterations++);
        m_pathList.GetLeastAdvanced(processedPath);
        
        //std::cout << "====" << TestID++ << "====" << "PathID:" << processedPath.m_pPath->m_nPathId << std::endl;
        //std::cout << "Size:" << m_pathList.Size() + 1 << " Range:" << lastSyncPos + 1 << "-" << m_nCurrentPosition + 1 << " LocalIterations:" << currentIterations << std::endl;
       
        if(m_pathList.Size() == 0)
        {
            const std::vector<const COrientedVariant*>* calledIncluded;
            const std::vector<int>* calledExcluded;
            calledIncluded = &(processedPath.m_pPath->m_calledSemiPath.GetIncludedVariants());
            calledExcluded = &(processedPath.m_pPath->m_calledSemiPath.GetExcluded());
            
            const std::vector<const COrientedVariant*>* baseIncluded;
            const std::vector<int>* baseExcluded;
            baseIncluded = &(processedPath.m_pPath->m_baseSemiPath.GetIncludedVariants());
            baseExcluded = &(processedPath.m_pPath->m_baseSemiPath.GetExcluded());

            if(calledIncluded->size() > 0 || baseIncluded->size() > 0)
            {
                for(int k = 0; k < (int)calledIncluded->size(); k++)
                    m_IncludedVariantsCalledBest.push_back((*calledIncluded)[k]);
                for(int k = 0; k < (int)baseIncluded->size(); k++)
                    m_IncludedVariantsBaselineBest.push_back((*baseIncluded)[k]);
                
                for(int k = 0; k < (int)calledExcluded->size(); k++)
                    m_ExcludedVariantsCalledBest.push_back((*calledExcluded)[k]);
                for(int k = 0; k < (int)baseExcluded->size(); k++)
                    m_ExcludedVariantsBaselineBest.push_back((*baseExcluded)[k]);
                
                for(int k = 0; k < (int)processedPath.m_pPath->m_aSyncPointList.size(); k++)
                    m_SyncPointsBest.push_back(processedPath.m_pPath->m_aSyncPointList[k]);
                
                processedPath.m_pPath->ClearSyncPointList();
                processedPath.m_pPath->ClearIncludedVariants();
                processedPath.m_pPath->ClearExcludedVariants();
            }

            //std::cout << TestID << "Included Size:" << processedPath.m_pPath->m_baseSemiPath.GetIncludedVariants().size() << " Excluded Size:" << processedPath.m_pPath->m_baseSemiPath.GetExcluded().size() << std::endl;
            
            int currentSyncPos = processedPath.m_pPath->m_calledSemiPath.GetPosition();
            if(currentMax > maxPaths)
            {
                maxPaths = currentMax;
                //maxPathRegion =  "chr21:" + std::to_string(lastSyncPos + 1) + "-" + std::to_string(currentSyncPos + 1);
                //std::cout << "Maximum path complexity now " << maxPaths << ", at " << maxPathRegion << " with "  << currentIterations << " iterations" << std::endl;
            }
            currentMax = 0;
            currentMaxIterations = std::max(currentIterations, currentMaxIterations);
            currentIterations = 0;
            lastSyncPos = currentSyncPos;
            lastSyncPath = processedPath;
        }
        else if(m_pathList.Size() >  m_nMaxPathSize || currentIterations > m_nMaxIterationCount)
        {
            complexRegionCount++;
            std::cerr << "Evaluation is too complex!";
            std::cerr << " There are " << m_pathList.Size() << " unresolved paths, " << currentIterations << " iterations at reference region ";
            std::cerr << a_contig.m_chromosome << ":" << (lastSyncPos + 1) << "-" << (m_nCurrentPosition + 2) << std::endl;

            //Drop all paths currently in play
            m_pathList.Clear();
            currentIterations = 0;
            // Create new head containing path up until last sync point
            processedPath = lastSyncPath;
            //Ignore variants until Current Position
            totalSkippedVariantCount += SkipVariantsTo(*processedPath.m_pPath, a_contig, m_nCurrentPosition+1);
        }

        if(processedPath.m_pPath->HasFinished())
        {
            //std::cout << "processed path is finished" << std::endl;
            //Path is done. Update the Best Path if it is better
            CPathContainer processedCopy(*processedPath.m_pPath, processedPath.m_pPath->m_calledSemiPath.GetPosition());
            best = FindBetter(best, processedCopy) ? best : processedCopy;
            continue;
        }

        if(EnqueueVariant(*processedPath.m_pPath, eCALLED, a_contig.m_nChrId, a_bIsGenotypeMatch))
        {
            //std::cout << "Called semipath enqueued" << std::endl;
            continue;
        }
        
        if(EnqueueVariant(*processedPath.m_pPath, eBASE, a_contig.m_nChrId, a_bIsGenotypeMatch))
        {
            //std::cout << "Base semipath enqueued" << std::endl;
            continue;
        }

        processedPath.m_pPath->Step();
        
        //processedPath.m_pPath->Print();

        
        if(processedPath.m_pPath->InSync())
        {
            SkipToNextVariant(*processedPath.m_pPath, a_contig);
            //std::cout << "In Sync/ skip to next variant" << std::endl;
        }
        else
        {
            //std::cout << "Not in sync" << std::endl;
        }

        if(processedPath.m_pPath->Matches())
        {
            //std::cout << "Head matches, keeping" << std::endl;
            AddIfBetter(processedPath);
            //m_pathList.Print();
        }
        
        else
        {
           //std::cout << "Head mismatch, discard" << std::endl;
        }
        
    }
    
    
    const std::vector<const COrientedVariant*>* calledIncluded;
    const std::vector<int>* calledExcluded;
    calledIncluded = &(best.m_pPath->m_calledSemiPath.GetIncludedVariants());
    calledExcluded = &(best.m_pPath->m_calledSemiPath.GetExcluded());
    
    const std::vector<const COrientedVariant*>* baseIncluded;
    const std::vector<int>* baseExcluded;
    baseIncluded = &(best.m_pPath->m_baseSemiPath.GetIncludedVariants());
    baseExcluded = &(best.m_pPath->m_baseSemiPath.GetExcluded());

    for(int k = 0; k < (int)calledIncluded->size(); k++)
        m_IncludedVariantsCalledBest.push_back((*calledIncluded)[k]);
    for(int k = 0; k < (int)baseIncluded->size(); k++)
        m_IncludedVariantsBaselineBest.push_back((*baseIncluded)[k]);
    for(int k = 0; k < (int)calledExcluded->size(); k++)
        m_ExcludedVariantsCalledBest.push_back((*calledExcluded)[k]);
    for(int k = 0; k < (int)baseExcluded->size(); k++)
        m_ExcludedVariantsBaselineBest.push_back((*baseExcluded)[k]);
    for(int k = 0; k < (int)best.m_pPath->m_aSyncPointList.size(); k++)
        m_SyncPointsBest.push_back(best.m_pPath->m_aSyncPointList[k]);
    
    best.m_pPath->ClearSyncPointList();
    best.m_pPath->AddSyncPointList(m_SyncPointsBest);
    best.m_pPath->ClearIncludedVariants();
    best.m_pPath->AddIncludedVariants(m_IncludedVariantsCalledBest, m_IncludedVariantsBaselineBest);
    best.m_pPath->ClearExcludedVariants();
    best.m_pPath->AddExcludedVariants(m_ExcludedVariantsCalledBest, m_ExcludedVariantsBaselineBest);
    
    std::cerr << "FINISHED " << a_contig.m_nChrId + 1 << ": Complex Region: " << complexRegionCount;
    std::cerr << " Skipped Variant Count :" << totalSkippedVariantCount;
    std::cerr << " Maximum path complexity is " << maxPaths << ", with "  << currentMaxIterations << " iterations " << std::endl;
    return *best.m_pPath;
}

void CPathReplay::AddIfBetter(const CPathContainer& a_path)
{
    std::set<CPathContainer>::iterator it;

    int cnt = m_pathList.Size();
    
    it = m_pathList.Find(a_path);
    if(it != m_pathList.End())
    {
        CPathContainer other = m_pathList.floor(a_path);
        CPathContainer best =  FindBetter(a_path, other) ? a_path : other;
    
        if(true == best.m_pPath->IsEqual(*a_path.m_pPath))
        {
            m_pathList.Erase(other);
            m_pathList.Add(best);
            assert(cnt == m_pathList.Size());
        }
    }
    
    else
    {
        int cnt = m_pathList.Size();
        m_pathList.Add(a_path);
        assert(cnt < m_pathList.Size());
    }
}

bool CPathReplay::FindBetter2(const CPathContainer& lhs, const CPathContainer& rhs)
{
    // See if we have obvious no-ops we would rather drop
    bool lhsSync = lhs.m_pPath->InSync() || lhs.m_pPath->HasFinished();
    bool rhsSync = rhs.m_pPath->InSync() || rhs.m_pPath->HasFinished();
    
    if(lhsSync && rhsSync)
    {
        if(lhs.m_pPath->HasNoOperation())
        {
            return false;
        }
        else if(rhs.m_pPath->HasNoOperation())
        {
            return true;
        }
    }

    const int lhsVariantCountCalled = static_cast<int>(lhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size());
    const int rhsVariantCountCalled = static_cast<int>(rhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size());
    
    if(lhsVariantCountCalled != rhsVariantCountCalled)
        return lhsVariantCountCalled > rhsVariantCountCalled;
    
    const int lhsVariantCountBaseline = static_cast<int>(lhs.m_pPath->m_baseSemiPath.GetIncludedVariants().size());
    const int rhsVariantCountBaseline = static_cast<int>(rhs.m_pPath->m_baseSemiPath.GetIncludedVariants().size());
    
    if(lhsVariantCountBaseline != rhsVariantCountBaseline)
        return lhsVariantCountBaseline > rhsVariantCountBaseline;

    // Prefer solutions that minimize discrepencies between baseline and call counts since last sync point
    const int lhsDelta = abs(lhs.m_pPath->m_nBSinceSync - lhs.m_pPath->m_nCSinceSync);
    const int rhsDelta = abs(rhs.m_pPath->m_nBSinceSync - rhs.m_pPath->m_nCSinceSync);
    if(lhsDelta != rhsDelta)
        return lhsDelta < rhsDelta ? true : false;
    
    // Prefer solutions that sync more regularly (more likely to be "simpler")
    const int syncDelta = (lhs.m_pPath->m_aSyncPointList.size() == 0 ? 0 : lhs.m_pPath->m_aSyncPointList.back()) - (rhs.m_pPath->m_aSyncPointList.size() == 0 ? 0 : rhs.m_pPath->m_aSyncPointList.back());
    if(syncDelta != 0)
        return syncDelta > 0 ? true : false;
    
    
    const std::vector<const COrientedVariant*> lhsIncludedVariants = (lhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? lhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
                                                                                                                                         lhs.m_pPath->m_calledSemiPath.GetIncludedVariants();
    const std::vector<const COrientedVariant*> rhsIncludedVariants = (rhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? rhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
                                                                                                                                         rhs.m_pPath->m_calledSemiPath.GetIncludedVariants();

    // At this point break ties arbitrarily based on allele ordering
    return (lhsIncludedVariants.back()->GetAlleleIndex() < rhsIncludedVariants.back()->GetAlleleIndex()) ? true : false;

}





bool CPathReplay::FindBetter(const CPathContainer& lhs, const CPathContainer& rhs)
{
    // See if we have obvious no-ops we would rather drop
    bool lhsSync = lhs.m_pPath->InSync() || lhs.m_pPath->HasFinished();
    bool rhsSync = rhs.m_pPath->InSync() || rhs.m_pPath->HasFinished();
    
    if(lhsSync && rhsSync)
    {
        if(lhs.m_pPath->HasNoOperation())
        {
            return false;
        }
        else if(rhs.m_pPath->HasNoOperation())
        {
            return true;
        }
    }
    
    // Prefer paths that maximise total number of included variants (baseline + called)

    const std::vector<const COrientedVariant*> lhsIncludedVariants = (lhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? lhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
                                                                                                                                   lhs.m_pPath->m_calledSemiPath.GetIncludedVariants();
    const std::vector<const COrientedVariant*> rhsIncludedVariants = (rhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? rhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
                                                                                                                                   rhs.m_pPath->m_calledSemiPath.GetIncludedVariants();
    
    const int lhsVariantCount = static_cast<int>(lhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() + lhs.m_pPath->m_baseSemiPath.GetIncludedVariants().size());
    const int rhsVariantCount = static_cast<int>(rhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() + rhs.m_pPath->m_baseSemiPath.GetIncludedVariants().size());
    
    if(lhsVariantCount == rhsVariantCount)
    {
        //Tie break equivalently scoring paths for greater aesthetics
        if(lhsIncludedVariants.size() != 0 && rhsIncludedVariants.size() != 0)
        {
            
            // Prefer solutions that minimize discrepencies between baseline and call counts since last sync point
            const int lhsDelta = abs(lhs.m_pPath->m_nBSinceSync - lhs.m_pPath->m_nCSinceSync);
            const int rhsDelta = abs(rhs.m_pPath->m_nBSinceSync - rhs.m_pPath->m_nCSinceSync);
            if(lhsDelta != rhsDelta)
                return lhsDelta < rhsDelta ? true : false;

            // Prefer solutions that sync more regularly (more likely to be "simpler")
            const int syncDelta = (lhs.m_pPath->m_aSyncPointList.size() == 0 ? 0 : lhs.m_pPath->m_aSyncPointList.back()) - (rhs.m_pPath->m_aSyncPointList.size() == 0 ? 0 : rhs.m_pPath->m_aSyncPointList.back());
            if(syncDelta != 0)
                return syncDelta > 0 ? true : false;
            
            // At this point break ties arbitrarily based on allele ordering
            return (lhsIncludedVariants.back()->GetAlleleIndex() < rhsIncludedVariants.back()->GetAlleleIndex()) ? true : false;
        }
    }
    
    return lhsVariantCount > rhsVariantCount ? true : false;
}


bool CPathReplay::EnqueueVariant(CPath& a_rPathToPlay, EVcfName a_uVcfSide, int a_nChromosomeId, bool a_bIsGenotypeMatch)
{
    const CSemiPath* pSemiPath = 0;

    switch(a_uVcfSide)
    {
        case eBASE:
            pSemiPath = &a_rPathToPlay.m_baseSemiPath;
            break;
        case eCALLED:
            pSemiPath = &a_rPathToPlay.m_calledSemiPath; 
            break;
    }

    int nVariantId = GetNextVariant(*pSemiPath);
    
    if(nVariantId != -1)
    {
        const CVariant* pNext = a_uVcfSide == eBASE ? m_aVariantListBase[nVariantId] : m_aVariantListCalled[nVariantId];
        //std::cout << "Add alternatives to " << ((a_uVcfSide == eBASE) ? "BASE " : "CALLED ") << pNext->ToString() << std::endl;
        
        m_nCurrentPosition = std::max(m_nCurrentPosition, pNext->GetStart());
        CPathContainer paths[3];
        int pathCount = a_rPathToPlay.AddVariant(paths,
                                                  a_uVcfSide,
                                                 (a_uVcfSide == eBASE ? m_aVariantListBase : m_aVariantListCalled),
                                                 (a_uVcfSide == eBASE ? m_aOrientedVariantListBase : m_aOrientedVariantListCalled),
                                                 nVariantId,
                                                 a_bIsGenotypeMatch);
        
        for(int k=0; k < pathCount; k++)
        {
            AddIfBetter(paths[k]);
        }
        return true;
    }

    return false;
}

void CPathReplay::SkipToNextVariant(CPath& a_rProcessedPath, const SContig& a_rContig)
{
    int aNext = FutureVariantPosition(a_rProcessedPath.m_calledSemiPath, eCALLED, a_rContig);
    int bNext = FutureVariantPosition(a_rProcessedPath.m_baseSemiPath, eBASE, a_rContig);
    int lastTemplatePos = a_rContig.m_nRefLength - 1;
    // -1 because we want to be before the position
    int nextPos = std::min(std::min(aNext,bNext), lastTemplatePos) -1;

    //std::cout << "Next Position is:" << nextPos << std::endl;

    assert (a_rProcessedPath.m_calledSemiPath.GetPosition() == a_rProcessedPath.m_baseSemiPath.GetPosition());
    
    if(nextPos > a_rProcessedPath.m_calledSemiPath.GetPosition())
       a_rProcessedPath.MoveForward(nextPos);

}

int CPathReplay::FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, const SContig& a_rContig) const
{
    int nextIdx = a_rSemiPath.GetVariantIndex() + 1;
    int varListSize = a_uVcfName == eBASE ? static_cast<int>(m_aVariantListBase.size()) : static_cast<int>(m_aVariantListCalled.size());
    
    if (nextIdx >= varListSize)
    {
        return a_rContig.m_nRefLength -1;
    } 

    else 
    {
        int nextVarStart = a_uVcfName == eBASE ? m_aVariantListBase[nextIdx]->GetStart() : m_aVariantListCalled[nextIdx]->GetStart();
        return nextVarStart;
    }
}

int CPathReplay::GetNextVariant(const CSemiPath& a_rSemiPath) const
{
    int nextId = a_rSemiPath.GetVariantIndex() + 1;
    int varListSize = a_rSemiPath.GetVcfName() == eBASE ? static_cast<int>(m_aVariantListBase.size()) : static_cast<int>(m_aVariantListCalled.size());
    
    if(nextId >= varListSize)
        return -1;

    const CVariant* nextVar = a_rSemiPath.GetVcfName() == eBASE ? m_aVariantListBase[nextId] : m_aVariantListCalled[nextId];
    
    if(nextVar->GetStart() <= (a_rSemiPath.GetPosition() + 1))
        return nextId;
    
    if(a_rSemiPath.WantsFutureVariantBases() && (nextVar->GetStart() <= a_rSemiPath.GetVariantEndPosition()))
        return nextId;

    return -1;

}

int CPathReplay::SkipVariantsTo(CPath& a_rPath, const SContig& a_rContig, int a_nMaxPos)
{
    //BASE SEMIPATH
    int varIndex = a_rPath.m_baseSemiPath.GetVariantIndex();
    int baseSkippedCount = 0;
    int calledSkippedCount = 0;
    
    while(varIndex < (int)m_aVariantListBase.size() && (varIndex == -1  || m_aVariantListBase[varIndex]->GetStart() < a_nMaxPos))
    {
        varIndex++;
        baseSkippedCount++;
    }
    varIndex--;
    a_rPath.m_baseSemiPath.SetVariantIndex(varIndex);
    a_rPath.m_baseSemiPath.MoveForward(std::min(a_nMaxPos, a_rContig.m_nRefLength -1));
    
    //std::cerr << "Baseline skipped Variant Count:" << baseSkippedCount << std::endl;
    
    //CALLED SEMIPATH
    varIndex = a_rPath.m_calledSemiPath.GetVariantIndex();
    
    while(varIndex < (int)m_aVariantListCalled.size() && (varIndex == -1  || m_aVariantListCalled[varIndex]->GetStart() < a_nMaxPos))
    {
        varIndex++;
        calledSkippedCount++;
    }
    varIndex--;
    a_rPath.m_calledSemiPath.SetVariantIndex(varIndex);
    a_rPath.m_calledSemiPath.MoveForward(std::min(a_nMaxPos, a_rContig.m_nRefLength-1));
    
    //std::cerr << "Called skipped Variant Count:" <<  calledSkippedCount << std::endl;

    
    return calledSkippedCount + baseSkippedCount;
}

void CPathReplay::Clear()
{
    m_pathList.Clear();
    m_nCurrentPosition = 0;
    m_SyncPointsBest.clear();
    m_ExcludedVariantsCalledBest.clear();
    m_IncludedVariantsCalledBest.clear();
    m_ExcludedVariantsBaselineBest.clear();
    m_IncludedVariantsBaselineBest.clear();
}






