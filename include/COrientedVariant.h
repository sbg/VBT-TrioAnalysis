/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef _C_ORIENTED_VARIANT_H_
#define _C_ORIENTED_VARIANT_H_

#include "CVariant.h"


class COrientedVariant
{
    public: 
    COrientedVariant();
    
    //Initialize Oriented with a variant and orientation selection
    COrientedVariant(const CVariant& a_rObj, bool a_bIsOrderOfGenotype);
    
    //Copy constructor
    COrientedVariant(const COrientedVariant& a_rObj);    
    
    //Get the allele alt string
    const SAllele& GetAllele() const;
    
    //Compare variants according to start/end position
    int CompareTo(const COrientedVariant& a_rObj) const;
    
    //Gets the start position of the variant
    int GetStartPos() const;
    
    //Gets the index of the allele
    int GetAlleleIndex() const;
    
    //Gets the end position of the allele
    int GetEndPos() const;
    
    //Get the Allele in the reverse side of chosen orientation
    COrientedVariant Other() const;
    
    //Return if the variant is null
    bool IsNull() const;
    
    //Return true if hap A is used and return false if hap B is used
    bool IsOrderOfGenotype() const;
    
    //Return the variant
    const CVariant& GetVariant() const;
    
    //Set variant to null
    void SetToNull();
    
    //Set the weight(proportion of variant inside 2 consecutive syncpoints) of the variant
    void SetWeight(double a_fWeight) const;
    
    //[For Test purpose] print the oriented variant
    void Print() const;

    private:
    //Index of the selected allele of this variant
    int m_nAlleleIndex;
    //Index of the other allele of this variant
    int m_nOtherAlleleIndex;
    //If the selected allele is first number or not (eg.  a/b   a-> true b-> false)
    bool m_bIsOrderOfGenotype;
    //Weight of the variant
    mutable double m_fWeight;
    bool m_bIsNull;
    
    const CVariant* m_variant;
    
};

#endif //_C_ORIENTED_VARIANT_H_
