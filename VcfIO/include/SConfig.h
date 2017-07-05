//  Created by Berke.Toptas on 11/11/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

#include "Constants.h"

//Any configuration that will be used send via this config object
struct SConfig
{
    //Indicates if enough parameters are received from command line
    bool m_bIsOK;
    
    //Base Vcf file
    const char* m_pBaseVcfFileName;
    
    //Called Vcf file
    const char* m_pCalledVcfFileName;
    
    //Fasta file
    const char* m_pFastaFileName;
    
    //PED file
    const char* m_pPedigreeFileName;
    bool m_bInitializeFromPED = false;
    
    //Output Directory
    const char* m_pOutputDirectory;
    
    //Output Mode
    const char* m_pOutputMode = "SPLIT";
    
    //Comparison engine mode
    bool m_bIsGenotypeMatch = true;
    
    //If filtering enabled
    bool m_bIsFilterEnabled = true;
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName = "PASS";
   
    //VBT mendelian mode output names
    const char* m_output_vcf_name = "VBT_mendelian_output.vcf";
    const char* m_output_log_name = "VBT_mendelian_output_detailedLog.txt";
    
    //If Quality threshold is enabled
    bool m_bIsQualityThresholdEnabled = false;
    //Threshold value of the variants
    float m_fQualityThreshold;
    
    //If Base sample name choosing is enabled
    bool m_bBaseSampleEnabled = false;
    //Name of the base sample to be selected
    const char* m_pBaseSample;
    
    //If Called sample name choosing is enabled
    bool m_bCalledSampleEnabled = false;
    //Name of the called sample to be selected
    const char* m_pCalledSample;
    
    //If set true, check only autosomes
    bool m_bAutosomeOnly = false;
    
    //Process only SNPs when true
    bool m_bSNPOnly = false;
    
    //Process only INDELs when true
    bool m_bINDELOnly = false;
    
    //Enable REF Overlap mode
    bool m_bIsRefOverlap = false;
    
    //Number of thread to use during execution
    int m_nThreadCount = DEFAULT_THREAD_COUNT;
    
    //Maximum number of path that variant comparison core can store at a time [History Table for Dynamic Programming]
    int m_nMaxPathSize = DEFAULT_MAX_PATH_SIZE;

    //Maximum number of iteration to resolve a variant (max iteration count for a variant to give TP/FP/FN decision)
    int m_nMaxIterationCount = DEFAULT_MAX_ITERATION_SIZE;
    
    //Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize = DEFAULT_MAX_BP_LENGTH;
    
};


#endif // _S_CONFIG_H
