/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_PATH_H_
#define _C_PATH_H_

#include "SResult.h"
#include "CSemiPath.h"
#include "EReplayChoice.h"
#include "EVcfName.h"
#include <list>

class CPath
{
    public:

    CPath(const char* a_aRefSequence, int a_nRefSize);
    CPath(const CPath& a_rObj, std::list<int>& a_syncPoints);
    
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
    
    //Add variant to the given side of path
    std::vector<CPath> AddVariant(EVcfName a_nVcfName, const CVariant& a_rVariant, int a_nVariantIndex);
    
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

    //Semi path object for base 
    CSemiPath m_baseSemiPath;
    
    //Semi path object for called
    CSemiPath m_calledSemiPath;
    
    //Position index list of the syncronisation points
    std::list<int> m_aSyncPointList;
    
    //Added variant count to called since last sync 
    int m_nCSinceSync;
    //Added variant count to called since last sync
    int m_nBSinceSync;

};


#endif // _C_PATH_H_
