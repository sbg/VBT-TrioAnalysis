//  Created by Berke.Toptas on 11/11/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

//Any configuration that will be used send via this config object
struct SConfig
{
    //Indicates if enough parameters are received from command line
    bool m_bIsOK;
    
    //Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize = 1000;
    
    //Base Vcf file
    const char* m_pBaseVcfFileName;
    
    //Called Vcf file
    const char* m_pCalledVcfFileName;
    
    //Fasta file
    const char* m_pFastaFileName;
    
    //If filtering enabled
    bool m_bIsFilterEnabled = true;
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName = "PASS";
   
    //If Quality threshold is enabled
    bool m_bIsQualityThresholdEnabled;
    //Threshold value of the variants
    float m_fQualityThreshold;
    
};


#endif // _S_CONFIG_H
