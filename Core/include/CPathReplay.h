/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/PathFinder.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 * Copyright (c) 2016 Seven Bridges Genomics
 *               2017 SBGD Inc.
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
