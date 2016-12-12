/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
 * Copyright (c) 2014. Real Time Genomics Limited.
 
 *
 * Moved to C++ by Berke Cagkan Toptas
 */

#ifndef _C_PATH_H_
#define _C_PATH_H_

#include "SResult.h"
#include "CSemiPath.h"
#include "EReplayChoice.h"
#include "EVcfName.h"

class CPathContainer;
class CVariantProvider;

class CPath
{
    public:
    
    CPath();
    CPath(const char* a_aRefSequence, int a_nRefSize);
    CPath(const CPath& a_rObj, int a_nSyncPointToPush);
    
    //Copy constructor
    CPath(const CPath& a_rObj);
    
    // Check if the two semipaths are synchronized
    bool InSync() const;
    
    //Check if the given path is equal to this
    bool IsEqual(const CPath& a_rObj) const;
    
    //Compare two the given path with this
    int CompareTo(const CPath& a_rObj) const;
    
    //Exclude the variant to the given side
    CPath& Exclude(EVcfName a_nVCF, const CVariant& a_rVariant, int a_nVariantIndex);
    
    //Include variant to the given side
    CPath& Include(EVcfName a_nVCF, const COrientedVariant& a_rVariant, int a_nVariantIndex);
    
    //Add variant to the given side of path and return the path count
    int AddVariant(CPathContainer* a_pPathList, EVcfName a_nVcfName,const CVariantProvider* a_pVariantProvider, int a_nVariantIndex, int a_nChromosomeId);
    
    bool operator<(const CPath& a_rObj) const
    {
        return CompareTo(a_rObj) < 0;
    }
    
    //
    void Step();
    
    //Force move haplotypes to the given position
    void MoveForward(int a_nPosition);

    //[FOR TEST] Print the values of path
    void Print() const;

    //Check if baseline semipath matches with the called semipath
    bool Matches();

    //return If the path has no called or based variant after last sync
    bool HasNoOperation() const;
    
    //Check if the path has finished
    bool HasFinished() const;

    //Find a weighting for all the TP calls in a path. this is done by sync points, within each SyncPoint
    //this will assure that the total number of TP we output will always reflect number of TP in baseline file
    //if there are any call TP without corresponding baseline TP, these are simply assigned a weight of 0.
    //(this can happen when two calls cancel each other out when replayed, although the default path finding now avoids this)
    void CalculateWeights();
    
    //Semi path object for base 
    CSemiPath m_baseSemiPath;
    
    //Semi path object for called
    CSemiPath m_calledSemiPath;
    
    //Position index list of the syncronisation points
    std::vector<int> m_aSyncPointList;
    
    //Added variant count to called since last sync 
    int m_nCSinceSync;
    
    //Added variant count to called since last sync
    int m_nBSinceSync;
    
    //TEST Purpose
    int m_nPathId;
    
};

class CPathContainer
{
public:
    
    CPathContainer()
    {
        m_pPath = 0;
    }
    
    CPathContainer(const CPathContainer& a_rObj)
    {
        m_pPath = a_rObj.m_pPath;
    }
    
    CPathContainer(const char* a_aRefSequence, int a_nRefSize)
    {
        m_pPath = std::shared_ptr<CPath>(new CPath(a_aRefSequence,a_nRefSize));
    }
    
    CPathContainer(const CPath& a_rPath, int a_nSyncPointToPush)
    {
        m_pPath = std::shared_ptr<CPath>(new CPath(a_rPath, a_nSyncPointToPush));
    }
    
    bool operator<(const CPathContainer& a_rObj) const
    {
        return m_pPath->CompareTo(*a_rObj.m_pPath) < 0;
    }
    
    //Pointer to path object
    std::shared_ptr<CPath> m_pPath;
    
};




#endif // _C_PATH_H_













