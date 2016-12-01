//  Created by Berke.Toptas on 11/11/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

struct SConfig
{
    //Any configuration that will be used send via this config object
    
    
    //Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize;
    
    //Base Vcf file
    const char* m_pBaseVcfFileName;
    
    //Called Vcf file
    const char* m_pCalledVcfFileName;
    
    //Fasta file
    const char* m_pFastaFileName;
    
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName;
    
    //Threshold value of the variants
    float m_fQualityThreshold;
};


#endif // _S_CONFIG_H
