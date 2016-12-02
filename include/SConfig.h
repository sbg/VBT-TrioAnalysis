//  Created by Berke.Toptas on 11/11/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

//Any configuration that will be used send via this config object
struct SConfig
{
    SConfig(int argc, char** argv)
    {
        const char* PARAM_BASE = "-b";
        const char* PARAM_CALLED = "-c";
        const char* PARAM_FILTER = "-f";

        int it = 1;
        
        while(it <= argc)
        {
            if(0 == strcmp(argv[it], PARAM_BASE))
                m_pBaseVcfFileName = argv[it+1];
            else if(0 == strcmp(argv[it], PARAM_CALLED))
                m_pCalledVcfFileName = argv[it+1];
            else if(0 == strcmp(argv[it], PARAM_FILTER))
                m_pFilterName = argv[it+1];
        
            it += 2;
        }
    }
    
    SConfig()
    {}
    
    
    //Indicates if enough parameters are received from command line
    bool m_bIsOK;
    
    //Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize;
    
    //Base Vcf file
    const char* m_pBaseVcfFileName;
    
    //Called Vcf file
    const char* m_pCalledVcfFileName;
    
    //Fasta file
    const char* m_pFastaFileName;
    
    //If filtering enabled
    bool m_bIsFilterEnabled;
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName;
   
    //If Quality threshold is enabled
    bool m_bIsQualityThresholdEnabled;
    //Threshold value of the variants
    float m_fQualityThreshold;
    
};


#endif // _S_CONFIG_H
