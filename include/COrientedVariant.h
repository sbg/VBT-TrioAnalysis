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
    std::string GetAlleleString() const;
    
    //Compare variants according to start/end position
    int CompareTo(const COrientedVariant& a_rObj) const;
    
    //Gets the start position of the allele
    int GetAlleleStartPos() const;
    
    //Gets the index of the allele
    int GetAlleleIndex() const;
    
    //Gets the end position of the allele
    int GetAlleleEndPos() const;    
    
    //Get the Allele in the reverse side of chosen orientation
    COrientedVariant Other() const;
    
    //Return if the variant is null
    bool IsNull() const;
    
    //Return the variant
    const CVariant& GetVariant() const;
    
    //Set variant to null
    void SetToNull();

    private:
    //Index of the selected allele of this variant
    int m_nAlleleIndex;
    //Index of the other allele of this variant
    int m_nOtherAlleleIndex;
    //If the selected allele is first number or not (eg.  a/b   a-> true b-> false)
    bool m_bIsOrderOfGenotype;
    CVariant m_variant;

};

#endif //_C_ORIENTED_VARIANT_H_
