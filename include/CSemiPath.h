/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_SEMI_PATH_H_
#define _C_SEMI_PATH_H_

#include "CHaplotypeSequence.h"
#include "EVcfName.h"

class CSemiPath
{
    public:

    CSemiPath(const char* a_aRefSequence, int a_nRefSize, EVcfName a_uVcfName);
    CSemiPath(const CSemiPath& a_rObj);

    //Include the given variant to this semipath
    void IncludeVariant(const COrientedVariant& a_rVariant, int a_nVariantIndex);
    
    //Exclude the given variant to this semipath
    void ExcludeVariant(const CVariant& a_rVariant, int a_nVariantIndex);

    //Compare and returns the template position diffence of hapA and hapB
    int CompareHaplotypePositions() const;

    //Return vcf name of the semi path
    EVcfName GetVcfName() const;
    
    //Gets the end position the semipath (max of hapA and hapB)
    int GetPosition() const;
    
    //Return the end position of last variant added
    int GetVariantEndPosition() const;
    
    //Return the end position of last included variant
    int GetIncludedVariantEndPosition() const;
    
    //Return the index of last variant added
    int GetVariantIndex() const;
    
    //Check whether this half path is fully on the template (i.e. no haplotypes are within a variant)
    bool IsOnTemplate() const;
    
    //Compare this half path with the given half path
    int CompareTo(const CSemiPath& a_rObj) const;
    
    //Check whether this half path is equal to the given half path
    bool IsEqual(const CSemiPath& a_rObj) const;

    //Checks if both haplotype A and B is finished
    bool HasFinished() const;

    //Force semipath to move the given position
    void MoveForward(int a_nPosition);
    
    //Chech whether this half path matches with the given half path
    bool Matches(const CSemiPath& a_rOther);

    //Test whether a deficit of variant bases are upstream in the queue in order to perform a step.
    //return false indicates that no variants need to be immediately enqueued
    bool WantsFutureVariantBases() const;
  
    bool FinishedHaplotypeA() const;
    bool FinishedHaplotypeB() const;

    //Return the next base on the B haplotype
    char NextHaplotypeBBase() const;
    
    //Return the next base on the A haplotype
    char NextHaplotypeABase() const;

    void StepHaplotypeA();
    void StepHaplotypeB();

    private:

    // name of the semipath
    EVcfName m_uVcfName;
    // Index of last variant added
    int m_nVariantIndex;
    // End of last variant added
    int m_nVariantEndPosition;
    // Last variant included
    int m_nIncludedVariantEndPosition;

    std::vector<COrientedVariant> m_aIncludedVariants;
    std::vector<CVariant> m_aExcludedVariants; 

    CHaplotypeSequence m_haplotypeA;
    CHaplotypeSequence m_haplotypeB;

    bool m_bFinishedHapA;
    bool m_bFinishedHapB;

};


#endif // _C_SEMI_PATH_H_
