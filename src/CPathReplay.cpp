#include "CPathReplay.h"
#include <iostream>
#include <vector>

inline int max(int a, int b)
{
    return a > b ? a : b;
}

inline int min(int a, int b)
{
    return a < b ? a : b;
}

void CPathReplay::InitializeReaders(const char* a_pBaseVcfFile,
                       const char* a_pCalledVcfFile,
                       const char* a_pRefFastaFile)
{
    //Initialize variant provider
    m_variantProvider.InitializeReaders(a_pBaseVcfFile, a_pCalledVcfFile);

    //Initialize refseq provider
    m_refFASTA.Open(a_pRefFastaFile);
    m_refFASTA.ReadContig();
}

CPath CPathReplay::FindBestPath(int a_nChrId)
{
    CPath initialPath(m_refFASTA.GetRefSeq(),m_refFASTA.GetRefSeqSize());
    m_pathList.Add(initialPath);
    CPath* best = 0;
    CPath* lastSyncPath = 0;
    int maxPaths = 0;
    int currentIterations = 0;
    int currentMaxIterations = 0;
    int currentMax = 0;
    m_nCurrentPosition = 0;
    int lastSyncPos = 0;
    std::string maxPathRegion;

    while(!m_pathList.Empty())
    {
        std::cout << "Iteration No:" << currentIterations << std::endl;
        currentMax = max(currentMax, m_pathList.Size());
        currentMaxIterations = max(currentMaxIterations, currentIterations++);
        CPath processedPath = m_pathList.GetLeastAdvanced();
        if(m_pathList.Size() == 0)
        {
            std::cout << "Pathlist size is 0" << std::endl;
            int currentSyncPos = processedPath.m_calledSemiPath.GetPosition();
            if(currentMax > maxPaths)
            {
                maxPaths = currentMax;
                maxPathRegion = ":" + std::to_string(lastSyncPos + 1) + "-" + std::to_string(currentSyncPos + 1);
                std::cout << "Maximum path complexity now " << maxPaths << ", at " << maxPathRegion << " with "  << currentIterations << " iterations" << std::endl;
            }
            currentMax = 0;
            currentIterations = 0;
            lastSyncPos = currentSyncPos;
            lastSyncPath = &processedPath;
        }
        else if(m_pathList.Size() >  MAX_COMPLEXITY || currentIterations > MAX_ITERATION)
        {
            std::cout << "Evaluation is too complex!!" << std::endl;
            //Drop all paths currently in play
            m_pathList.Clear();
            currentIterations = 0;
            // Create new head containing path up until last sync point
            processedPath = *lastSyncPath;
        }

        if(processedPath.HasFinished())
        {
            std::cout << "processed path is finished" << std::endl;
            //Path is done. Update the Best Path if it is better
            std::list<int> syncPoints(processedPath.m_aSyncPointList);
            syncPoints.push_front(processedPath.m_calledSemiPath.GetPosition());
            CPath processedCopy(processedPath, syncPoints);
            best = FindBetter(*best, processedCopy) ? best : &processedCopy; 
            continue;
        }

        if(EnqueueVariant(processedPath, eCALLED, a_nChrId))
        {
            std::cout << "Called semipath enqueued" << std::endl;
            continue;
        }
        
        if(EnqueueVariant(processedPath, eBASE, a_nChrId))
        {
            std::cout << "Base semipath enqueued" << std::endl;
            continue;
        }

        processedPath.Step();
        
        if(processedPath.InSync())
        {
            SkipToNextVariant(processedPath, a_nChrId);
            std::cout << "In Sync/ skip to next variant" << std::endl;
        }
        else
        {
            std::cout << "Not in sync" << std::endl;
        }

        if(processedPath.Matches())
        {
            std::cout << "Head matches, keeping" << std::endl;
            AddIfBetter(processedPath);
        }
        else
        {
            std::cout << "Head mismatch, discard" << std::endl;
        }


    }
    
    std::cout << "Best Path Found" << std::endl;
     return *best;
}

void CPathReplay::AddIfBetter(std::vector<CPath> a_pathsToAdd)
{
    for(int k=0; k< a_pathsToAdd.size(); k++)
    {
        if(m_pathList.Contains(a_pathsToAdd[k]))
        {
            CPath other = m_pathList.floor(a_pathsToAdd[k]);           
            CPath best =  FindBetter(a_pathsToAdd[k], other) ? a_pathsToAdd[k] : other;

            if(true == best.IsEqual(a_pathsToAdd[k]))
            {
                m_pathList.Erase(other);
                m_pathList.Add(best);
            }
        }
        else
        {
            m_pathList.Add(a_pathsToAdd[k]);
        }
    }
}

void CPathReplay::AddIfBetter(const CPath& a_path)
{
    std::vector<CPath> pathstoAdd;
    pathstoAdd.push_back(a_path);
    AddIfBetter(pathstoAdd);
}


bool CPathReplay::FindBetter(const CPath& lhs, const CPath& rhs)
{
    //TODO: Implement
    return false;
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
        m_nCurrentPosition = max(m_nCurrentPosition, pNext->GetStart());
        std::vector<CPath> paths = a_rPathToPlay.AddVariant(a_uVcfSide, *pNext, nVariantId);
        for(int k = 0; k < paths.size(); k++)
        {
            std::cout<< paths[k].m_nBSinceSync << " " << paths[k].m_nCSinceSync << " " << std::endl;
            std::cout<< paths[k].m_baseSemiPath.GetPosition() << " " << paths[k].m_baseSemiPath.GetVariantEndPosition() << " " << paths[k].m_baseSemiPath.GetVariantIndex() << std::endl;
            std::cout<< paths[k].m_calledSemiPath.GetPosition() << " " << paths[k].m_calledSemiPath.GetVariantEndPosition() << " " << paths[k].m_calledSemiPath.GetVariantIndex() << std::endl;
        }
        AddIfBetter(paths);
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
    int nextPos = min(min(aNext,bNext), lastTemplatePos) -1;

    if(nextPos > a_rProcessedPath.m_calledSemiPath.GetPosition())
       a_rProcessedPath.MoveForward(nextPos);

}

int CPathReplay::FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, int a_nChromosomeId)
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


int CPathReplay::GetNextVariant(const CSemiPath& a_rSemiPath, int a_nChromosomeId)
{
    int nextId = a_rSemiPath.GetVariantIndex() + 1;
    
    if(nextId >= m_variantProvider.GetVariantListSize(a_rSemiPath.GetVcfName(), a_nChromosomeId))
        return -1;

    const CVariant* nextVar = m_variantProvider.GetVariant(a_rSemiPath.GetVcfName(), a_nChromosomeId, nextId);
    
    if(nextVar->GetStart() <= a_rSemiPath.GetPosition() + 1)
        return nextId;
    
    if(a_rSemiPath.WantsFutureVariantBases() && nextVar->GetStart() <= a_rSemiPath.GetVariantEndPosition())
        return nextId;

    return -1;

}
