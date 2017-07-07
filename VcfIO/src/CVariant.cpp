//
//  CVariant.cpp
//  VCFComparison
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#include "CVariant.h"
#include <iostream>
#include <sstream>
#include "Constants.h"


CVariant::CVariant(): m_nVcfId(-1),
            m_nChrId(-1),
            m_chrName(),
            m_bIsPhased(false),
            m_nStartPos(-1)

{
    m_nId = -1;
    m_nAlleleCount = 0;
    m_bIsFirstNucleotideTrimmed = false;
    
    m_allelesStr = "";
    m_nZygotCount = 0;
    m_variantStatus = eNOT_ASSESSED;
    m_genotype[0] = -1;
    m_genotype[1] = -1;
    m_bIsNoCall = false;
}

CVariant::CVariant(const CVariant& a_rObj)
{
    m_nVcfId = a_rObj.m_nVcfId;
    m_nChrId = a_rObj.m_nChrId;
    m_chrName = a_rObj.m_chrName;
    m_bIsPhased = a_rObj.m_bIsPhased;
    m_nAlleleCount = a_rObj.m_nAlleleCount;
    m_alleles[0].m_nEndPos = a_rObj.m_alleles[0].m_nEndPos;
    m_alleles[0].m_nStartPos = a_rObj.m_alleles[0].m_nStartPos;
    m_alleles[0].m_sequence = a_rObj.m_alleles[0].m_sequence;
    m_alleles[0].m_bIsIgnored = a_rObj.m_alleles[0].m_bIsIgnored;
    m_alleles[1].m_nEndPos = a_rObj.m_alleles[1].m_nEndPos;
    m_alleles[1].m_nStartPos = a_rObj.m_alleles[1].m_nStartPos;
    m_alleles[1].m_sequence = a_rObj.m_alleles[1].m_sequence;
    m_alleles[1].m_bIsIgnored = a_rObj.m_alleles[1].m_bIsIgnored;
    
    m_nStartPos = a_rObj.m_nStartPos;
    m_nEndPos = a_rObj.m_nEndPos;
    m_nId = a_rObj.m_nId;
    m_refSequence = a_rObj.m_refSequence;
    m_bIsHeterozygous = a_rObj.m_bIsHeterozygous;
    m_bIsFirstNucleotideTrimmed = a_rObj.m_bIsFirstNucleotideTrimmed;
    
    m_filterString = a_rObj.m_filterString;
    m_allelesStr = a_rObj.m_allelesStr;
    m_nZygotCount = a_rObj.m_nZygotCount;
    m_genotype[0] = a_rObj.m_genotype[0];
    m_genotype[1] = a_rObj.m_genotype[1];
    m_variantStatus = a_rObj.m_variantStatus;
    m_nOriginalPos = a_rObj.m_nOriginalPos;
    m_bIsNoCall = a_rObj.m_bIsNoCall;
    
}

int CVariant::CompareTo(const CVariant& a_rObj) const
{
    if(GetStart() < a_rObj.GetStart())
        return -1;
    else if(GetStart() > a_rObj.GetStart())
        return 1;
    else
        return GetEnd() - a_rObj.GetEnd();
}

bool CVariant::Clear()
{
    m_nVcfId = -1;
    m_nChrId = -1;
    m_nStartPos = -1;
    m_chrName.clear();
    m_nAlleleCount = 0;
    m_nZygotCount = 0;
    m_allelesStr.clear();
    m_filterString.clear();
    m_bIsFirstNucleotideTrimmed = false;
    m_variantStatus = eNOT_ASSESSED;
    m_nOriginalPos = -1;
    m_alleles[0].m_bIsIgnored = false;
    m_alleles[1].m_bIsIgnored = false;
    m_genotype[0] = -1;
    m_genotype[1] = -1;
    m_bIsNoCall = false;
    return true;
}

bool CVariant::IsHeterozygous() const
{
    return m_bIsHeterozygous;
}

void CVariant::Print() const
{
    if(m_nVcfId == 0)
        std::cout << "Belongs to  : Baseline" << std::endl;
    else
        std::cout << "Belongs to  : Called" << std::endl;
    std::cout <<     "Ref : " << m_refSequence << std::endl;
    for(int k = 0; k < m_nAlleleCount; k++)
    {
        std::cout << "Alt" << k << ": " << m_alleles[k].m_sequence << std::endl;
    }
}

int CVariant::GetId() const
{
    return m_nId;
}

int CVariant::GetStart() const
{
    return m_nStartPos;
}

int CVariant::GetEnd() const
{
    return m_nEndPos;
}

int CVariant:: GetOriginalPos() const
{
    return m_nOriginalPos;
}

bool CVariant::IsPhased() const
{
    return m_bIsPhased;
}

std::string CVariant::GetRefSeq() const
{
    return m_refSequence;
}

SAllele CVariant::GetAllele(int a_nAlleleId) const
{
    return m_alleles[a_nAlleleId];
}

std::string CVariant::GetOriginalAlleleStr(int a_nAlleleIndex) const
{
    if(m_bIsNoCall == true)
        return "";
    
    std::vector<std::string> alleleList;
    std::stringstream ss(m_allelesStr);
    std::string alel;
    while (std::getline(ss, alel, ','))
    {
        alleleList.push_back(alel);
    }
    
    assert(a_nAlleleIndex < (int)alleleList.size());
    
    return alleleList[m_genotype[a_nAlleleIndex]];
}

bool CVariant::IsNull() const
{
    if(m_nVcfId == -1)
        return true;
    else
        return false;
}

std::string CVariant::ToString() const
{
    std::string toRet = "";
    
    toRet = m_chrName + ":" + std::to_string(GetStart() + 1) + "-" + std::to_string(GetEnd() + 1) + " (";
    
    for(int k=0; k < m_nAlleleCount; k++)
    {
        if(k > 0)
            toRet = toRet + ":";
        
        toRet = toRet + (m_alleles[k].m_bIsIgnored ? "*" : m_alleles[k].m_sequence);
    }
    toRet = toRet + ")";
    
    return toRet;
}

bool CVariant::IsFilterPASS() const
{
    return m_bIsFilterPASS;
}

EVariantType CVariant::GetVariantType() const
{
    //SNP CASE
    if(m_refSequence.length() == 1 && m_alleles[0].m_sequence.length() == 1 && m_alleles[1].m_sequence.length() == 1)
        return eSNP;
    
    //SV CASE
    std::size_t found;
    for(int k = 0; k < m_nAlleleCount; k++)
    {
        std::string allele = m_alleles[k].m_sequence;

        found = allele.find('[');
        if(found != std::string::npos)
            return eSV;
        found = allele.find('<');
        if(found != std::string::npos)
            return eSV;
        found = allele.find('*');
        if(found != std::string::npos)
            return eSV;
        found = allele.find('.');
        if(found != std::string::npos)
            return eSV;
    }
    
    //INDEL OTHERWISE
    return eINDEL;
}


EVariantCategory CVariant::GetVariantCategory() const
{
    //Split alleles string to vector of allele strings
    std::vector<std::string> alleles;
    std::stringstream ss(m_allelesStr);
    
    std::string al;
    while(std::getline(ss, al, ','))
    {
        alleles.push_back(al);
    }
    
    //SNP CASE
    if(alleles[0].length() == 1 && alleles[m_genotype[0]].length() == 1 && alleles[m_genotype[1]].length() == 1)
        return EVariantCategory::eSNP;

    //INDEL- INSERTION CASE
    else if(alleles[0].length() == 1 && (alleles[m_genotype[0]].length() > 1 || alleles[m_genotype[1]].length() > 1))
    {
        int varLength = (int)std::max(alleles[m_genotype[0]].length(), alleles[m_genotype[1]].length());
        
        if(varLength < 5)
            return EVariantCategory::eINDEL_INSERT_SMALL;
        else if(varLength < 15)
            return EVariantCategory::eINDEL_INSERT_MEDIUM;
        else
            return EVariantCategory::eINDEL_INSERT_LARGE;
    }

    //INDEL - DELETION CASE
    else if(alleles[0].length() > 1 && (alleles[m_genotype[0]].length() == 1 ||  alleles[m_genotype[1]].length() == 1))
    {
        int varLength = (int)alleles[0].length();
    
        if(varLength < 5)
            return EVariantCategory::eINDEL_DELETE_SMALL;
        else if(varLength < 15)
            return EVariantCategory::eINDEL_DELETE_MEDIUM;
        else
            return EVariantCategory::eINDEL_DELETE_LARGE;
    }
 
    //COMPLEX
    else
    {
        int varLength = (int)std::max(alleles[0].length(), std::max(alleles[m_genotype[0]].length(), alleles[m_genotype[1]].length()));

        if(varLength < SMALL_VARIANT_SIZE)
            return EVariantCategory::eINDEL_COMPLEX_SMALL;
        else if(varLength < MEDIUM_VARIANT_SIZE)
            return EVariantCategory::eINDEL_COMPLEX_MEDIUM;
        else
            return EVariantCategory::eINDEL_COMPLEX_LARGE;
    }

}












