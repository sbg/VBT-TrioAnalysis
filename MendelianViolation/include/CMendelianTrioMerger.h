//
//  CMendelianTrioMerger.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

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
    
private:

    ///Merge 3 variant set into one trio.vcf that mendelian decisions are marked
    void AddRecords(SChrIdTriplet& a_rTriplet);
    
    
    ///Fill the header part of the trio vcf
    void FillHeader();
    
    //Merge given three variant to a single record
    void DoMerge(const CVariant* a_pVarMother,
                 const CVariant* a_pVarFather,
                 const CVariant* a_pVarChild,
                 EMendelianDecision a_initDecision,
                 std::vector<SVcfRecord>& a_rRecordList);
    
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
    
    ///Objects to pass the result Log
    SMendelianDetailedLogEntry m_logEntry;
    SMendelianDetailedLogGenotypes m_logGenotypes;
    
    CMendelianResultLog* m_pResultLog;
};

}


#endif // _C_MENDELIAN_TRIO_MERGER_H_
