/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  CBaseVariantProvider.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Toptas on 12/7/17.
 *
 */

#ifndef _C_BASE_VARIANT_PROVIDER_H_
#define _C_BASE_VARIANT_PROVIDER_H_

#include "EVariantMatch.h"
#include "CVcfReader.h"
#include "CFastaParser.h"
#include <vector>

class CVariant;

namespace core
{
    class COrientedVariant;
}

class CBaseVariantProvider
{
    
public:
    
    ///Set the status of each variant in the given list
    void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
    
    ///Set the status of each variant in the given list
    void SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
            
    ///Read contig given by the chromosome id
    bool ReadContig(std::string a_chrId, SContig& a_rContig);
    
protected:

    ///Find the optimal Trimming for variant list that have more than 1 trimming options. (See Readme under 'core' folder)
    void FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, std::vector<std::vector<CVariant>>* a_pAllVarList, const SConfig& a_rConfig);
    
    //REFERENCE FASTA
    CFastaParser m_referenceFasta;


};



#endif /* _C_BASE_VARIANT_PROVIDER_H_ */
