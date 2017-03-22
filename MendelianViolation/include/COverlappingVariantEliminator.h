//
//  COverlapingVariantEliminator.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/2/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "htslib/vcf.h"
#include "CVcfReader.h"
#include <string>


class COverlappingVariantEliminator
{
    
public:
    
    void FilterOverlaps(const std::string& a_rFileName, bool a_bIsFilterOverlap, bool a_bIsFilter00, int a_00filterId);

private:
    
    htsFile*    m_pHtsFileFiltered;
    bcf_hdr_t * m_pHeaderFiltered;
    bcf1_t *    m_pRecordFiltered;
    
    //Vcf reader instance
    CVcfReader m_vcfReader;
    
    //Output vcf file name
    std::string m_fileName;
};
