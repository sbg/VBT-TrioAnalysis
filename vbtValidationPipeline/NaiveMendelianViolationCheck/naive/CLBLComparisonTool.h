//
//  CLBLComparisonTool.h
//  gatk_sbgd_trio_test
//
//  Created by Berke.Toptas on 11/17/17.
//  Copyright © 2017 SBGD. All rights reserved.
//


#ifndef CLBLComparisonTool_h
#define CLBLComparisonTool_h

#include "CVcfReader.h"
#include "CVcfWriter.h"

class CLBLComparisonTool
{
    
public:
    
    //Set Input trio path
    void SetInputTrio(const std::string& a_rInputTrioPath, const std::string& a_rPedigreeFile);
    
    //Write annotated trio to the given path
    void WriteOutputTrio(const std::string& a_rOutputTrioPath);
    
private:
    
    //Class to read trio
    CVcfReader m_vcfReader;
    
    //Mother, father and child indexes in vcf (0-based Sample Order in VCF )
    int m_motherIndex, m_fatherIndex, m_childIndex;
    
    //Choice of allele trimming order
    bool m_bIsTrimBeginningFirst;
    
    //Class to output annotated trio
    CVcfWriter m_vcfWriter;
    
    std::string pedigreeFile;

};



#endif /* CLBLComparisonTool_h */
