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
 *  Created by Berke Cagkan Toptas
 *
 */

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
    
    //Output Directory
    const char* m_pOutputDirectory;
    
    //If filtering enabled
    bool m_bIsFilterEnabled = true;
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName = "PASS";
   
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
    
    //Process only SNPs when true
    bool m_bSNPOnly = false;
    
    //Process only INDELs when true
    bool m_bINDELOnly = false;
    
    //Enable REF Overlap mode
    bool m_bIsRefOverlap = false;
    
    ///When a variant allele has multiple trimming option, TRUE: clips from beginning first FALSE: clips from ending first
    bool m_bTrimBeginningFirst = true;
    
    //Run in Platform mode [1 Thread per chromosome]
    bool m_bIsPlatformMode = false;
    
};


#endif // _S_CONFIG_H
