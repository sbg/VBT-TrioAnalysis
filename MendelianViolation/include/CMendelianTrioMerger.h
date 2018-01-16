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
 *  CMendelianTrioMerger.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 3/9/17.
 *
 */

#ifndef _C_MENDELIAN_TRIO_MERGER_H_
#define _C_MENDELIAN_TRIO_MERGER_H_

#include "CVcfReader.h"
#include "CVcfWriter.h"
#include "CVariant.h"
#include "EMendelianDecision.h"
#include "EMendelianVcfName.h"
#include <vector>
#include "Constants.h"
#include "ENoCallMode.h"
#include "CMendelianResultLog.h"
#include "SChrIdTriplet.h"

namespace mendelian
{

/**
 * @brief Generates output trio vcf annotated with mendelian decisions using vcfs of family members
 *
 */
class CMendelianTrioMerger
{
    
public:
    
    ///Generates the trio vcf by merging parent-child variants into single 3-sample vcf file
    void GenerateTrioVcf(std::vector<SChrIdTriplet>& a_rCommonChromosomes);
    
    ///Set contigs to write output header. Also initializes size of decision and variant arrays according to number of common chromosomes
    void SetContigList(const std::vector<SVcfContig>& a_rContigs,
                       int a_nCommonContigCount,
                       int a_nChildVarListSize,
                       int a_nFatherVarListSize,
                       int a_nMotherVarListSize);
    
    ///Set variants and their decisions for marking MendelianDecision info columns
    void SetDecisionsAndVariants(SChrIdTriplet& a_rTriplet, EMendelianVcfName a_vcfName,
                                 const std::vector<EMendelianDecision>& a_rDecisionList,
                                 const std::vector<const CVariant*>& a_rVarList);

    ///Set No Call mode to decide if no calls will be printed as ./. or 0/0s
    void SetNoCallMode(ENoCallMode a_mode);
    
    ///Set the full path of output trio vcf
    void SetTrioPath(const std::string& a_nTrioPath);
    
    ///Set the access of result log from mendelian vcf analyzer to for detailed logs
    void SetResultLogPointer(CMendelianResultLog* a_pResultLog);
    
    ///Set the filename of child to access INFOs. If this function is called, then all info columns will be copied to the generated trio records
    void SetInfoReadParameters(const std::string& a_rChildInputPath,
                               const std::string& a_rFatherInputPath,
                               const std::string& a_rMotherInputPath);
    
private:

    ///Merge 3 variant set into one trio.vcf that mendelian decisions are marked
    void AddRecords(SChrIdTriplet& a_rTriplet);
    
    ///Write VCF records to the VCF file and fill the log tables
    void WriteRecords(const std::vector<SVcfRecord>& recordList,
                            const std::vector<EVariantCategory>& recordCategoryList,
                            const std::vector<EMendelianDecision>& recordDecisionList);

    ///Fill the header part of the trio vcf
    void FillHeader();
    
    ///Fill the INFO tags of the heaader that are taken from input vcf files
    void FillInfoHeaderLines();
    
    //Merge given three variant to a single record
    void DoMerge(const CVariant* a_pVarMother,
                 const CVariant* a_pVarFather,
                 const CVariant* a_pVarChild,
                 EMendelianDecision a_initDecision,
                 std::vector<SVcfRecord>& a_rRecordList);
   
    //Add unique alleles of given variant to the allele list
    void AddAllele(const CVariant* a_pVariant, int a_nMaxRefSequenceLength, std::vector<std::string>& alleles);
    
    ///Fill a_rSample with the information of variant
    void AddSample(const CVariant* a_pVariant, std::vector<std::string>& a_rAlleles, SPerSampleData& a_rSample);
    
    ///Decide the mendelian type of the given three variant in a row
    EMendelianDecision GetMendelianDecision(const CVariant* a_pVarMother,
                                            const CVariant* a_pVarFather,
                                            const CVariant* a_pVarChild,
                                            SChrIdTriplet& a_rTriplet);
    
    ///Register a line of merged vcf to the detailed report table [updates m_logEntry]
    void RegisterMergedLine(EMendelianDecision a_decision, EVariantCategory a_category);
    
    ///Register the genotype of merged vcf to genotype table [updates m_logGenotypes]
    void RegisterGenotype(const CVariant* a_pMother,
                          const CVariant* a_pFather,
                          const CVariant* a_pChild,
                          EMendelianDecision a_initDecision,
                          EVariantCategory a_category);
    
    ///Register the genotype of merged vcf record to genotype table [updates m_logGenotypes]
    void RegisterGenotype(const SVcfRecord& a_rRecord, EVariantCategory a_uCategory, EMendelianDecision a_uDecision);
    
    ///For the given recordList unify the overlapping variant decisions
    void ProcessRefOverlappedRegions(std::vector<SVcfRecord>&  a_rRecordList, std::vector<EMendelianDecision>& a_rRecordDecisionList);
    
    ///if the record at recordItr and temporaryItr are overlapping change the decision of record at temopraryItr
    bool CheckForOverlap(std::vector<SVcfRecord>&  a_rRecordList, std::vector<EMendelianDecision>& a_rRecordDecisionList, int recordItr, int temporaryItr, EMendelianDecision curVariantDecision);
    
    ///Vcf writer instance
    CVcfWriter m_vcfWriter;
    
    std::vector<std::vector<EMendelianDecision>> m_aChildDecisions;
    std::vector<std::vector<EMendelianDecision>> m_aFatherDecisions;
    std::vector<std::vector<EMendelianDecision>> m_aMotherDecisions;

    std::vector<std::vector<const CVariant*>> m_aChildVariants;
    std::vector<std::vector<const CVariant*>> m_aFatherVariants;
    std::vector<std::vector<const CVariant*>> m_aMotherVariants;
    
    ///No call mode selection of mendelian violation check (implicit, explicit or none)
    ENoCallMode m_noCallMode;
    
    ///Contig list which will be used to fill header part of output vcf
    std::vector<SVcfContig> m_contigs;
    
    ///Output trio full path
    std::string m_trioPath;
    
    //Input VCF file paths (To access INFO column)
    std::string m_childPath;
    std::string m_motherPath;
    std::string m_fatherPath;
    
    //If annotations will be written to output merged vcf
    bool m_bIsAnnotationsON = false;
    
    ///Objects to pass the result Log
    SMendelianDetailedLogEntry m_logEntry;
    SMendelianDetailedLogGenotypes m_logGenotypes;
    
    CMendelianResultLog* m_pResultLog;
};

}


#endif // _C_MENDELIAN_TRIO_MERGER_H_
