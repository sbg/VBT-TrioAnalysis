//
//  CBaseVariantProvider.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/7/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#include "CBaseVariantProvider.h"
#include "CVariant.h"
#include "COrientedVariant.h"

void CBaseVariantProvider::SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const CVariant* pVar : a_rVariantList)
    {
        if(pVar->m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pVar->m_variantStatus != eGENOTYPE_MATCH)
            pVar->m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pVar-> m_variantStatus == eNOT_ASSESSED)
            pVar->m_variantStatus = a_status;
        else
            continue;
    }
}

void CBaseVariantProvider::SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const
{
    for(const core::COrientedVariant* pOVar : a_rVariantList)
    {
        if(pOVar->GetVariant().m_variantStatus == eCOMPLEX_SKIPPED)
            continue;
        
        if(a_status == eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eALLELE_MATCH && pOVar->GetVariant().m_variantStatus != eGENOTYPE_MATCH)
            pOVar->GetVariant().m_variantStatus = a_status;
        else if(a_status == eNO_MATCH && pOVar->GetVariant().m_variantStatus == eNOT_ASSESSED)
            pOVar->GetVariant().m_variantStatus = a_status;
        else
            continue;
    }
}

bool CBaseVariantProvider::ReadContig(std::string a_chrId, SContig& a_rContig)
{
    return m_referenceFasta.FetchNewChromosome(a_chrId, a_rContig);
}
