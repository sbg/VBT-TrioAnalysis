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
 *  SConfig.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 11/11/16.
 *
 */

#ifndef _S_CONFIG_H_
#define _S_CONFIG_H_

#include "Constants.h"

/**
 * @brief Any configuration that will be used will be distributed over classes via this object
 *
 */
struct SConfig
{
    ///Indicates if enough parameters are received from command line
    bool m_bIsOK;
    
    ///Base Vcf file
    const char* m_pBaseVcfFileName;
    
    ///Called Vcf file
    const char* m_pCalledVcfFileName;
    
    ///Fasta file
    const char* m_pFastaFileName;
    
    ///PED file
    const char* m_pPedigreeFileName;
    bool m_bInitializeFromPED = false;
    
    ///BED file
    const char* m_pBedFileName;
    bool m_bInitializeFromBed = false;
    
    ///Output Directory
    const char* m_pOutputDirectory;
    
    ///Output Mode
    const char* m_pOutputMode = "SPLIT";
    
    ///Comparison engine mode
    bool m_bIsGenotypeMatch = true;
    
    ///If filtering enabled
    bool m_bIsFilterEnabled = true;
    //Filter name (Filter variants according to filter column)
    const char* m_pFilterName = "PASS";
    
    ///VBT mendelian mode PREFIX
    const char* m_output_prefix = "out";
    
    ///If Quality threshold is enabled
    bool m_bIsQualityThresholdEnabled = false;
    //Threshold value of the variants
    float m_fQualityThreshold;
    
    ///If Base sample name choosing is enabled
    bool m_bBaseSampleEnabled = false;
    ///Name of the base sample to be selected
    const char* m_pBaseSample;
    
    ///If Called sample name choosing is enabled
    bool m_bCalledSampleEnabled = false;
    ///Name of the called sample to be selected
    const char* m_pCalledSample;
    
    ///If set true, check only autosomes
    bool m_bAutosomeOnly = false;
    
    ///Process only SNPs when true
    bool m_bSNPOnly = false;
    
    ///Process only INDELs when true
    bool m_bINDELOnly = false;
    
    ///Enable REF Overlap mode
    bool m_bIsRefOverlap = true;
    
    ///When a variant allele has multiple trimming option, TRUE: clips from beginning first FALSE: clips from ending first
    bool m_bTrimBeginningFirst = true;
    
    ///Enable generating syncpoint files which is the intermediate output of core module
    bool m_bGenerateSyncPoints = false;
    
    //Enable generating violation regions as BED file
    bool m_bGenerateViolationRegions = false;
    
    ///Enable reading whole info format into a structure while parsing VCF file
    bool m_bIsReadINFO = false;
    std::string m_infotags;
    
    ///Number of thread to use during execution
    int m_nThreadCount = DEFAULT_THREAD_COUNT;
    
    ///Maximum number of path that variant comparison core can store at a time [History Table for Dynamic Programming]
    int m_nMaxPathSize = DEFAULT_MAX_PATH_SIZE;

    ///Maximum number of iteration to resolve a variant (max iteration count for a variant to give TP/FP/FN decision)
    int m_nMaxIterationCount = DEFAULT_MAX_ITERATION_SIZE;
    
    ///Maximum size of the variant that will be processed by VCF comparison algorithm (Use it to eliminate SVs)
    int m_nMaxVariantSize = DEFAULT_MAX_BP_LENGTH;
    
};


#endif // _S_CONFIG_H
