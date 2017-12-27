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
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/PathFinder.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 * Licensed under the Simplified BSD License: https://github.com/RealTimeGenomics/rtg-tools/blob/master/LICENSE.txt
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 *
 */

#ifndef _C_PATH_REPLAY_H_
#define _C_PATH_REPLAY_H_

#include "CVariantProvider.h"
#include "CPathSet.h"
#include "CVariant.h"

namespace core
{

class CPath;
class CPathContainer;

/**
 * @brief Core Variant Comparison class that compares given 2 variant sets
 *
 * CPathReplay is the core variant comparison class. It takes 4 variant sets (Variant and Oriented variant list for base/call vcf file)
 */
class CPathReplay
{
    public:

        ///Sets variant for comparison
        CPathReplay(std::vector<const CVariant*>& a_aVarListBase,
                    std::vector<const CVariant*>& a_aVarListCalled,
                    std::vector<const COrientedVariant*>& a_aOVarListBase,
                    std::vector<const COrientedVariant*>& a_aOvarlistCalled);
    
    
        ///Sets maximum pathsize and maximum path iteration count
        void SetMaxPathAndIteration(int a_nMaxPathSize, int a_nMaxIterationCount);
    
        ///Clears variants belong to best path
        void Clear();
    
        /**
         * @brief Finds the best path by generating all possible paths for the given chromosome
         *
         * Main comparison function
         *
         * @param a_contig Chromosome to be processed
         * @param a_bIsGenotypeMatch comparison mode (true is genotype matching - ga4gh method3, and false is allele matching - ga4gh method2)
         */
        CPath FindBestPath(SContig a_contig, bool a_bIsGenotypeMatch);

    private:

        ///Add the paths to the sorted path list if there is no better path
        void AddIfBetter(const CPathContainer& a_path);

        ///Move the path to just before the next variant for either side
        void SkipToNextVariant(CPath& a_rProcessedPath, const SContig& a_rContig);
    
        ///Gets the next upstream variant position
        int FutureVariantPosition(const CSemiPath& a_rSemiPath, EVcfName a_uVcfName, const SContig& a_rContig) const;
    
        ///Compare the two paths and find the one that maximize TP count
        bool FindBetter(const CPathContainer& lhs, const CPathContainer& rhs);
    
        ///Process next variant for the input path
        bool EnqueueVariant(CPath& a_rPathToPlay, EVcfName a_uVcfSide, bool a_bIsGenotypeMatch);
    
        /**
         *Gets the index of the next variant if it should be enqueued to the supplied HalfPath at the current position,
         *or -1 if there is none to be enqueued at the current position
         */
        int GetNextVariant(const CSemiPath& a_rSemiPath) const;
    
        ///Move the path to the specified position, ignoring any intervening variants. Returns the skipped variant count
        int SkipVariantsTo(CPath& a_rPath, const SContig& a_rContig, int a_nMaxPos);
    
        ///Path list to store generated paths
        CPathSet m_pathList;
        

        int m_nCurrentPosition;
    
        std::vector<const COrientedVariant*> m_IncludedVariantsBaselineBest;
        std::vector<int> m_ExcludedVariantsBaselineBest;
        std::vector<const COrientedVariant*> m_IncludedVariantsCalledBest;
        std::vector<int> m_ExcludedVariantsCalledBest;
        std::vector<int> m_SyncPointsBest;
    
        //Between each sync points, there is only one path left in our search tree during the process. For this reason, once
        //we have 1 path left in the search tree, we can copy the content of it another list and clear the path data. By doing
        //this we can keep the path size small.

        std::vector<const CVariant*>& m_aVariantListBase;
        std::vector<const CVariant*>& m_aVariantListCalled;
        std::vector<const COrientedVariant*>& m_aOrientedVariantListBase;
        std::vector<const COrientedVariant*>& m_aOrientedVariantListCalled;
    
        ///Cutoff path size to fit in memory
        int m_nMaxPathSize;
        ///Cutoff iteration count without enqueing any variant to the pathlist
        int m_nMaxIterationCount;
};

}


#endif // _C_PATH_REPLAY_H_
