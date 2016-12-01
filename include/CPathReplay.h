/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_PATH_REPLAY_H_
#define _C_PATH_REPLAY_H_

#include "CVariantProvider.h"
#include "CThreadSafePathList.h"
#include "CPath.h"
#include "CVariant.h"
#include "EReplayChoice.h"
#include "CFastaReader.h"
#include "SConfig.h"

const int MAX_ITERATION = 10000000;
const int MAX_COMPLEXITY = 50000;

class CPathReplay
{
    public:

        // Open the base and called VCF files and Reference Fasta file
        void InitializeReaders(const SConfig& a_rConfig);

        // Finds the best path by generating all possible paths for the given chromosome
        CPath FindBestPath(int a_nChrId);

    private:

        // Add the paths to the sorted path list if there is no better path
        void AddIfBetter(std::vector<CPath> a_pathList);
        void AddIfBetter(const CPath& a_pathList);

        //Move the path to just before the next variant for either side
        void SkipToNextVariant(CPath& a_rProcessedPath, int a_nChromosomeId);
    
        //Gets the next upstream variant position
        int FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, int a_nChromosomeId) const;
    
        //Compare the two paths and find the one that maximize TP count
        bool FindBetter(const CPath& lhs, const CPath& rhs);
    
        //Process next variant for the input path
        bool EnqueueVariant(CPath& a_rPathToPlay, EVcfName a_uVcfSide, int a_nChromosomeId);
    
        //Gets the index of the next variant if it should be enqueued to the supplied HalfPath at the current position,
        //or -1 if there is none to be enqueued at the current position
        int GetNextVariant(const CSemiPath& a_rSemiPath, int a_nChromosomeId) const;
    
        //Move the path to the specified position, ignoring any intervening variants.
        void SkipVariantsTo(CPath& a_rPath, int a_nChromosomeId, int a_nMaxPos);
    
        //Path list to store generated paths
        CThreadSafePathList m_pathList;
    
        //FASTA reader
        CFastaReader m_refFASTA;
    
        //Variant Provider
        CVariantProvider m_variantProvider;

        int m_nCurrentPosition;

};




#endif // _C_PATH_REPLAY_H_
