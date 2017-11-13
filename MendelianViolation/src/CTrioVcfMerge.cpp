//
//  CTrioVcfMerge.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 11/8/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CTrioVcfMerge.h"
#include "CVariant.h"
#include <algorithm>


using namespace mendelian;

///Constructor. Set the variant sets
CTrioVcfMerge::CTrioVcfMerge(std::vector<const CVariant*>& a_rMotherVariants, std::vector<const CVariant*>& a_rFatherVariants, std::vector<const CVariant*>& a_rChildVariants) :
m_motherVariants(a_rMotherVariants),
m_fatherVariants(a_rFatherVariants),
m_childVariants(a_rChildVariants)
{
    motherItr = 0;
    fatherItr = 0;
    childItr = 0;
}

///Return the next record variants. If a sample is 0/0 null pointer will be returned for that sample. false will be return when iterator comes to the end
bool CTrioVcfMerge::GetNext(const CVariant** a_rMotherVariant, const CVariant** a_rFatherVariant, const CVariant** a_rChildVariant)
{
    bool bIsSuccess = true;
    
    //Refill our variant buffer
    if(motherVariantBuffer.size() == 0 && fatherVariantBuffer.size() == 0 && childVariantBuffer.size() == 0)
         bIsSuccess = RefillBuffers();
    
        
    //We have multiple variants starting at same position. Best combination will be selected as a new record
    if(motherVariantBuffer.size() > 1 || fatherVariantBuffer.size() > 1 || childVariantBuffer.size() > 1)
    {
        int motherStart = motherVariantBuffer.size() == 0 ? INT_MAX : motherVariantBuffer[0]->m_nStartPos;
        int fatherStart = fatherVariantBuffer.size() == 0 ? INT_MAX : fatherVariantBuffer[0]->m_nStartPos;
        int childStart = childVariantBuffer.size() == 0 ? INT_MAX : childVariantBuffer[0]->m_nStartPos;
        
        int minimumStart = std::min({fatherStart, motherStart, childStart});
        
        if(motherStart == minimumStart)
        {
            *a_rMotherVariant = motherVariantBuffer[0];
            motherVariantBuffer.erase(motherVariantBuffer.begin());
        }
        else
            *a_rMotherVariant = 0;
            
        if(fatherStart == minimumStart)
        {
            *a_rFatherVariant = fatherVariantBuffer[0];
            fatherVariantBuffer.erase(fatherVariantBuffer.begin());
        }
        else
            *a_rFatherVariant = 0;

        
        if(childStart == minimumStart)
        {
            *a_rChildVariant = childVariantBuffer[0];
            childVariantBuffer.erase(childVariantBuffer.begin());
        }
        else
            *a_rChildVariant = 0;

    }
    
    //We have at most 1 variant per sample, send them directly to the output pointers and clear buffers
    else
    {
        (*a_rMotherVariant) = motherVariantBuffer.size() == 0 ? 0 : motherVariantBuffer[0];
        (*a_rFatherVariant) = fatherVariantBuffer.size() == 0 ? 0 : fatherVariantBuffer[0];
        (*a_rChildVariant)  = childVariantBuffer.size() == 0 ? 0 : childVariantBuffer[0];
        
        //Clear buffers
        motherVariantBuffer.clear();
        fatherVariantBuffer.clear();
        childVariantBuffer.clear();
    }
    
    return bIsSuccess;
}

bool CTrioVcfMerge::RefillBuffers()
{
    //Has variant left to fill buffers
    if(childItr == (int)m_childVariants.size() && motherItr == (int)m_motherVariants.size() && fatherItr == (int)m_fatherVariants.size())
        return false;
    
    smallestVariantPosition = std::min({childItr  == (int)m_childVariants.size()  ?  INT_MAX : m_childVariants[childItr]->m_nOriginalPos,
                                        motherItr == (int)m_motherVariants.size() ?  INT_MAX : m_motherVariants[motherItr]->m_nOriginalPos,
                                        fatherItr == (int)m_fatherVariants.size() ?  INT_MAX :  m_fatherVariants[fatherItr]->m_nOriginalPos});
    
    while(childItr < (int)m_childVariants.size() && smallestVariantPosition == m_childVariants[childItr]->m_nOriginalPos)
        childVariantBuffer.push_back(m_childVariants[childItr++]);
    
    while(motherItr < (int)m_motherVariants.size() && smallestVariantPosition == m_motherVariants[motherItr]->m_nOriginalPos)
        motherVariantBuffer.push_back(m_motherVariants[motherItr++]);
    
    while(fatherItr < (int)m_fatherVariants.size() && smallestVariantPosition == m_fatherVariants[fatherItr]->m_nOriginalPos)
        fatherVariantBuffer.push_back(m_fatherVariants[fatherItr++]);
    
    return true;
}
