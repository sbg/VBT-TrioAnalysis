/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
 * Copyright (c) 2014. Real Time Genomics Limited.
 
 *
 * Moved to C++ by Berke Cagkan Toptas
 */

#ifndef _C_PATH_REPLAY_H_
#define _C_PATH_REPLAY_H_

#include "CVariantProvider.h"
#include "CThreadSafePathList.h"
#include "CVariant.h"

const int MAX_ITERATION = 10000000;
const int MAX_COMPLEXITY = 50000;

class CPath;
class CPathContainer;

class CPathReplay
{
    public:

        //Sets variant provider object instance
        void SetVariantProvider(const CVariantProvider& a_rVariantProvider);
    
        // Finds the best path by generating all possible paths for the given chromosome
        CPath FindBestPath(SContig a_contig);

    private:

        // Add the paths to the sorted path list if there is no better path
        void AddIfBetter(const CPathContainer& a_path);

        //Move the path to just before the next variant for either side
        void SkipToNextVariant(CPath& a_rProcessedPath, const SContig& a_rContig);
    
        //Gets the next upstream variant position
        int FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, const SContig& a_rContig) const;
    
        //Compare the two paths and find the one that maximize TP count
        bool FindBetter(const CPathContainer& lhs, const CPathContainer& rhs);
    
        //Process next variant for the input path
        bool EnqueueVariant(CPath& a_rPathToPlay, EVcfName a_uVcfSide, int a_nChromosomeId);
    
        //Gets the index of the next variant if it should be enqueued to the supplied HalfPath at the current position,
        //or -1 if there is none to be enqueued at the current position
        int GetNextVariant(const CSemiPath& a_rSemiPath, int a_nChromosomeId) const;
    
        //Move the path to the specified position, ignoring any intervening variants.
        void SkipVariantsTo(CPath& a_rPath, const SContig& a_rContig, int a_nMaxPos);
    
        //Path list to store generated paths
        CThreadSafePathList m_pathList;
        
        //Access to Variant Provider
        const CVariantProvider* m_variantProvider;

        int m_nCurrentPosition;
    
        std::vector<const COrientedVariant*> m_IncludedVariantsBaselineBest;
        std::vector<int> m_ExcludedVariantsBaselineBest;
    
        std::vector<const COrientedVariant*> m_IncludedVariantsCalledBest;
        std::vector<int> m_ExcludedVariantsCalledBest;
    
        std::vector<int> m_SyncPointsBest;
};




#endif // _C_PATH_REPLAY_H_
