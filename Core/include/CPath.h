/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval
 * Copyright (c) 2016 Seven Bridges Genomics

 
 *
 * Moved to C++ by Berke Cagkan Toptas
 */

#ifndef _C_PATH_H_
#define _C_PATH_H_

#include "CSemiPath.h"
#include "EVcfName.h"
#include <memory>

namespace core
{

class CPathContainer;
class CSyncPoint;
class CMendelianVariantProvider;

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
    int AddVariant(CPathContainer* a_pPathList,
                   EVcfName a_nVcfName,
                   const std::vector<const CVariant*>& a_pVariantList,
                   const std::vector<const COrientedVariant*>& a_pOVariantList,
                   int a_nVariantIndex,
                   bool a_bIsGenotypeMatch);
    
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
    
    //Delete all Included variants
    void ClearIncludedVariants();
    //Add variants to the included variant list
    void AddIncludedVariants(std::vector<const COrientedVariant*>& a_rIncludedVarListCalled,
                             std::vector<const COrientedVariant*>& a_rIncludedVarListBase);
    
    //Delete all excluded variant indexes
    void ClearExcludedVariants();
    //Add variant indexes to the excluded variant list
    void AddExcludedVariants(std::vector<int>& a_rIncludedVarListCalled, std::vector<int>& a_rIncludedVarListBase);
    
    //Delete all sync point list
    void ClearSyncPointList();
    //Add sync points to the sync point list
    void AddSyncPointList(std::vector<int>& a_rSyncPointArray);
    
    //Sorts included variants (baseline and called) according to variant ids
    void SortIncludedVariants();
    
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

}


#endif // _C_PATH_H_













