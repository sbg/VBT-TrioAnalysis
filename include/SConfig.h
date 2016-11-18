//  Created by Berke.Toptas on 11/11/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

struct SConfig
{
    //Any configuration that will be used send via this config object
    
    
    //Set true to read only the PASS variants to the variant list
    bool m_bIsFilterPASS;
    
    //Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize;
    
    //Base Vcf file
    const char* m_pBaseVcfFileName;
    
    //Called Vcf file
    const char* m_pCalledVcfFileName;
    
    //Fasta file
    const char* m_pFastaFileName;
    
};


#endif // _S_CONFIG_H
