//  Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
//  Copyright (c) 2014. Real Time Genomics Limited.
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

void CPathReplay::InitializeReaders(const SConfig& a_rConfig)
{
    //Initialize refseq provider
    m_refFASTA.Open(a_rConfig.m_pFastaFileName);
    m_refFASTA.ReadContig();
    
    //Initialize variant provider
    m_variantProvider.SetFastaReader(m_refFASTA);
    m_variantProvider.InitializeReaders(a_rConfig);

}

CPath CPathReplay::FindBestPath(int a_nChrId)
{
    CPathContainer initialPath(m_refFASTA.GetRefSeq(),m_refFASTA.GetRefSeqSize());
    m_pathList.Add(initialPath);
    CPathContainer best(initialPath);
    CPathContainer lastSyncPath;
    int maxPaths = 0;
    int currentIterations = 0;
    int currentMaxIterations = 0;
    int currentMax = 0;
    m_nCurrentPosition = 0;
    int lastSyncPos = 0;
    //int TestID = 0;
    //std::string maxPathRegion;
    
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
            int currentSyncPos = processedPath.m_pPath->m_calledSemiPath.GetPosition();
            if(currentMax > maxPaths)
            {
                maxPaths = currentMax;
                //maxPathRegion =  "chr21:" + std::to_string(lastSyncPos + 1) + "-" + std::to_string(currentSyncPos + 1);
                //std::cout << "Maximum path complexity now " << maxPaths << ", at " << maxPathRegion << " with "  << currentIterations << " iterations" << std::endl;
            }
            currentMax = 0;
            currentIterations = 0;
            lastSyncPos = currentSyncPos;
            lastSyncPath = processedPath;
        }
        else if(m_pathList.Size() >  MAX_COMPLEXITY || currentIterations > MAX_ITERATION)
        {
            std::cout << "Evaluation is too complex!!" << std::endl;
            //Drop all paths currently in play
            m_pathList.Clear();
            currentIterations = 0;
            // Create new head containing path up until last sync point
            processedPath = lastSyncPath;
            //Ignore variants until Current Position
            SkipVariantsTo(*processedPath.m_pPath, a_nChrId, m_nCurrentPosition+1);
        }

        if(processedPath.m_pPath->HasFinished())
        {
            //std::cout << "processed path is finished" << std::endl;
            //Path is done. Update the Best Path if it is better
            CPathContainer processedCopy(*processedPath.m_pPath, processedPath.m_pPath->m_calledSemiPath.GetPosition());
            best = FindBetter(best, processedCopy) ? best : processedCopy;
            continue;
        }

        if(EnqueueVariant(*processedPath.m_pPath, eCALLED, a_nChrId))
        {
            //std::cout << "Called semipath enqueued" << std::endl;
            continue;
        }
        
        if(EnqueueVariant(*processedPath.m_pPath, eBASE, a_nChrId))
        {
            //std::cout << "Base semipath enqueued" << std::endl;
            continue;
        }

        processedPath.m_pPath->Step();
        
        //processedPath.m_pPath->Print();

        
        if(processedPath.m_pPath->InSync())
        {
            SkipToNextVariant(*processedPath.m_pPath, a_nChrId);
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
    
     //std::cout << "Best Path Found" << std::endl;
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

    const std::vector<COrientedVariant*> lhsIncludedVariants = (lhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? lhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
                                                                                                                                   lhs.m_pPath->m_calledSemiPath.GetIncludedVariants();
    const std::vector<COrientedVariant*> rhsIncludedVariants = (rhs.m_pPath->m_calledSemiPath.GetIncludedVariants().size() == 0) ? rhs.m_pPath->m_baseSemiPath.GetIncludedVariants() :
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

bool CPathReplay::EnqueueVariant(CPath& a_rPathToPlay, EVcfName a_uVcfSide, int a_nChromosomeId)
{
    const CSemiPath* pSemiPath;

    switch(a_uVcfSide)
    {
        case eBASE:
            pSemiPath = &a_rPathToPlay.m_baseSemiPath;
            break;
        case eCALLED:
            pSemiPath = &a_rPathToPlay.m_calledSemiPath; 
            break;
    }

    int nVariantId = GetNextVariant(*pSemiPath, a_nChromosomeId);
    const CVariant* pNext = m_variantProvider.GetVariant(a_uVcfSide, a_nChromosomeId, nVariantId);
    
    if(nVariantId != -1)
    {
        //std::cout << "Add alternatives to " << ((a_uVcfSide == eBASE) ? "BASE " : "CALLED ") << pNext->ToString() << std::endl;
        
        m_nCurrentPosition = std::max(m_nCurrentPosition, pNext->GetStart());
        CPathContainer paths[3];
        int pathCount = a_rPathToPlay.AddVariant(paths, a_uVcfSide, &m_variantProvider, nVariantId, a_nChromosomeId);
        
        for(int k=0; k < pathCount; k++)
        {
            AddIfBetter(paths[k]);
        }
        return true;
    }

    return false;
}

void CPathReplay::SkipToNextVariant(CPath& a_rProcessedPath, int a_nChromosomeId)
{
    int aNext = FutureVariantPosition(a_rProcessedPath.m_calledSemiPath, eCALLED, a_nChromosomeId);
    int bNext = FutureVariantPosition(a_rProcessedPath.m_baseSemiPath, eBASE, a_nChromosomeId);
    int lastTemplatePos = m_refFASTA.GetRefSeqSize() -1;
    // -1 because we want to be before the position
    int nextPos = std::min(std::min(aNext,bNext), lastTemplatePos) -1;

    //std::cout << "Next Position is:" << nextPos << std::endl;

    assert (a_rProcessedPath.m_calledSemiPath.GetPosition() == a_rProcessedPath.m_baseSemiPath.GetPosition());
    
    if(nextPos > a_rProcessedPath.m_calledSemiPath.GetPosition())
       a_rProcessedPath.MoveForward(nextPos);

}

int CPathReplay::FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, int a_nChromosomeId) const
{
    int nextIdx = a_rSemiPath.GetVariantIndex() + 1;
    
    if (nextIdx >= m_variantProvider.GetVariantListSize(a_uVcfName, a_nChromosomeId)) 
    {
      return m_refFASTA.GetRefSeqSize() -1;
    } 

    else 
    {
      return m_variantProvider.GetVariant(a_uVcfName, a_nChromosomeId, nextIdx)->GetStart();
    }
}

int CPathReplay::GetNextVariant(const CSemiPath& a_rSemiPath, int a_nChromosomeId) const
{
    int nextId = a_rSemiPath.GetVariantIndex() + 1;
    
    if(nextId >= m_variantProvider.GetVariantListSize(a_rSemiPath.GetVcfName(), a_nChromosomeId))
        return -1;

    const CVariant* nextVar = m_variantProvider.GetVariant(a_rSemiPath.GetVcfName(), a_nChromosomeId, nextId);
    
    if(nextVar->GetStart() <= (a_rSemiPath.GetPosition() + 1))
        return nextId;
    
    if(a_rSemiPath.WantsFutureVariantBases() && (nextVar->GetStart() <= a_rSemiPath.GetVariantEndPosition()))
        return nextId;

    return -1;

}

void CPathReplay::SkipVariantsTo(CPath& a_rPath, int a_nChromosomeId, int a_nMaxPos)
{
    //BASE SEMIPATH
    int varIndex = a_rPath.m_baseSemiPath.GetVariantIndex();
    
    while(varIndex < m_variantProvider.GetVariantListSize(eBASE, a_nChromosomeId) && (varIndex == -1  || m_variantProvider.GetVariant(eBASE, a_nChromosomeId, varIndex)->GetStart() < a_nMaxPos))
    {
        varIndex++;
    }
    varIndex--;
    a_rPath.m_baseSemiPath.SetVariantIndex(varIndex);
    a_rPath.m_baseSemiPath.MoveForward(std::min(a_nMaxPos, m_refFASTA.GetRefSeqSize()-1));
    
    
    //CALLED SEMIPATH
    varIndex = a_rPath.m_calledSemiPath.GetVariantIndex();
    
    while(varIndex < m_variantProvider.GetVariantListSize(eCALLED, a_nChromosomeId) && (varIndex == -1  || m_variantProvider.GetVariant(eCALLED, a_nChromosomeId, varIndex)->GetStart() < a_nMaxPos))
    {
        varIndex++;
    }
    varIndex--;
    a_rPath.m_calledSemiPath.SetVariantIndex(varIndex);
    a_rPath.m_calledSemiPath.MoveForward(std::min(a_nMaxPos, m_refFASTA.GetRefSeqSize()-1));
    
}







