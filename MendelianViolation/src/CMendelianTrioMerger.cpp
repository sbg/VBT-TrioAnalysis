//
//  CMendelianTrioMerger.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianTrioMerger.h"
#include "CTrioVcfMerge.h"
#include "Constants.h"
#include "CUtils.h"
#include <algorithm>
#include <iostream>
#include <sstream>


using namespace mendelian;

//Checks if the given two range is overlapping
bool IsOverlap(SVcfRecord& rec1, SVcfRecord& rec2)
{
    return CUtils::IsOverlap(rec1.left, rec1.right, rec2.left, rec2.right);
}

void CMendelianTrioMerger::SetTrioPath(const std::string& a_nTrioPath)
{
    m_trioPath = std::string(a_nTrioPath);
}

void CMendelianTrioMerger::SetDecisionsAndVariants(SChrIdTriplet& a_rTriplet, EMendelianVcfName a_vcfName, const std::vector<EMendelianDecision>& a_rDecisionList, const std::vector<const CVariant*>& a_rVarList)
{
    std::vector<EMendelianDecision>* pDecisionList;
    std::vector<const CVariant*>* pVarList;

    assert(a_rDecisionList.size() == a_rVarList.size());
    
    switch (a_vcfName)
    {
        case eCHILD:
            pDecisionList = &m_aChildDecisions[a_rTriplet.m_nTripleIndex];
            pVarList = &m_aChildVariants[a_rTriplet.m_nCid];
            break;
            
        case eFATHER:
            pDecisionList = &m_aFatherDecisions[a_rTriplet.m_nTripleIndex];
            pVarList = &m_aFatherVariants[a_rTriplet.m_nFid];
            break;
            
        case eMOTHER:
            pDecisionList = &m_aMotherDecisions[a_rTriplet.m_nTripleIndex];
            pVarList = &m_aMotherVariants[a_rTriplet.m_nMid];
            break;
            
        default:
            break;
    }

    for(unsigned int k = 0; k < a_rVarList.size(); k++)
    {
        pDecisionList->push_back(a_rDecisionList[k]);
        pVarList->push_back(a_rVarList[k]);
    }
}

void CMendelianTrioMerger::SetNoCallMode(ENoCallMode a_mode)
{
    m_noCallMode = a_mode;
}

void CMendelianTrioMerger::SetResultLogPointer(CMendelianResultLog* a_pResultLog)
{
    m_pResultLog = a_pResultLog;
}

void CMendelianTrioMerger::SetContigList(const std::vector<SVcfContig>& a_rCommonContigs, int a_nCommonContigCount, int a_nChildVarListSize, int a_nFatherVarListSize, int a_nMotherVarListSize)
{
    m_contigs = a_rCommonContigs;
    
    //Allocate decisions vectors
    m_aChildDecisions  = std::vector<std::vector<EMendelianDecision>>(a_nCommonContigCount);
    m_aFatherDecisions = std::vector<std::vector<EMendelianDecision>>(a_nCommonContigCount);
    m_aMotherDecisions = std::vector<std::vector<EMendelianDecision>>(a_nCommonContigCount);
    
    //Allocate variant vectors
    m_aChildVariants  = std::vector<std::vector<const CVariant*>>(a_nChildVarListSize);
    m_aFatherVariants = std::vector<std::vector<const CVariant*>>(a_nFatherVarListSize);
    m_aMotherVariants = std::vector<std::vector<const CVariant*>>(a_nMotherVarListSize);
}

void CMendelianTrioMerger::FillHeader()
{
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= VBT Mendelian Violation Tool " + VBT_VERSION );
    
    //ADD MENDELIAN VIOLATION INFO TYPE
    m_vcfWriter.AddHeaderLine("##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Mendelian Violation Decision. (0)-complex, (1)-compliant, (2)-violation (3)-NoCall Parent (4)-NoCall Child \">");
    
    //ADD GT COLUMN
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    
    //ADD CONTIG IDs
    for(unsigned int k = 0; k < m_contigs.size(); k++)
        m_vcfWriter.AddHeaderLine("##contig=<ID=" + m_contigs[k].name + ",length=" + std::to_string(m_contigs[k].length) + ">");
    
    //ADD REQUIRED SAMPLES
    m_vcfWriter.AddSampleName("MOTHER");
    m_vcfWriter.AddSampleName("FATHER");
    m_vcfWriter.AddSampleName("CHILD");
    
    //CLOSE HEADER
    m_vcfWriter.WriteHeaderToVcf();
    
}

void CMendelianTrioMerger::GenerateTrioVcf(std::vector<SChrIdTriplet>& a_rCommonChromosomes)
{
    //Open Vcf file to write
    m_vcfWriter.CreateVcf(m_trioPath.c_str());
    
    //Fill the header
    FillHeader();
    
    //Initialize log objects
    m_logEntry.clear();
    m_logGenotypes.clear();
    
    for(unsigned int k = 0; k < a_rCommonChromosomes.size(); k++)
    {
        std::cerr << "[stderr] Writing chromosome " << a_rCommonChromosomes[k].m_chrName << std::endl;
        AddRecords(a_rCommonChromosomes[k]);
    }
    
    //Send the logs to the log class
    m_pResultLog->LogDetailedReport(m_logEntry);
    m_pResultLog->LogGenotypeMatrix(m_logGenotypes);
    
    m_vcfWriter.CloseVcf();
}


void CMendelianTrioMerger::AddRecords(SChrIdTriplet &a_rTriplet)
{
    std::vector<SVcfRecord> recordList;
    std::vector<EVariantCategory> recordCategoryList;
    std::vector<EMendelianDecision> recordDecisionList;
    
    CTrioVcfMerge trioMerger(m_aMotherVariants[a_rTriplet.m_nMid], m_aFatherVariants[a_rTriplet.m_nFid], m_aChildVariants[a_rTriplet.m_nCid]);
    
    const CVariant* motherVariant;
    const CVariant* fatherVariant;
    const CVariant* childVariant;
    
    while(true == trioMerger.GetNext(&motherVariant, &fatherVariant, &childVariant))
    {
        //Determine the decision of the record
        EMendelianDecision decision = GetMendelianDecision(motherVariant, fatherVariant, childVariant, a_rTriplet);
        
        //Ignore complex-skipped variants
        if(decision == eSkipped)
            continue;
        
        //Push record decision
        recordDecisionList.push_back(decision);
        
        //Push record category -> TODO: We should do it for final record?
        EVariantCategory category = motherVariant == 0 ? ((fatherVariant == 0) ? childVariant->GetVariantCategory() : fatherVariant->GetVariantCategory()) : motherVariant->GetVariantCategory();
        recordCategoryList.push_back(category);
   
        //Merge variant and push it to the recordList
        DoMerge(motherVariant, fatherVariant, childVariant, decision, recordList);
    }
    
    std::cerr << "Processing Overlapping Regions..." << std::endl;
    
    //Update records for overlapping regions
    ProcessRefOverlappedRegions(recordList, recordDecisionList);
    
    std::cerr << "Writing variants to output vcf..." << std::endl;
    
    //Write the final updated variants to output vcf and logs to the report table
    for(int k = 0; k < static_cast<int>(recordList.size()); k++)
    {
        m_vcfWriter.AddMendelianRecord(recordList[k]);
        RegisterMergedLine(recordDecisionList[k], recordCategoryList[k]);
        RegisterGenotype(recordList[k], recordCategoryList[k], recordDecisionList[k]);
    }
    
}

void CMendelianTrioMerger::AddSample(const CVariant* a_pVariant, std::vector<std::string>& a_rAlleles, SPerSampleData& a_rSampleData)
{
    SPerSampleData dataFather;
    if(a_pVariant != NULL)
    {
        a_rSampleData.m_nHaplotypeCount = a_pVariant->m_nZygotCount;
        a_rSampleData.m_bIsPhased = a_pVariant->m_bIsPhased; // Future Work: Phasings of variants we found can be written to output
        a_rSampleData.m_bIsNoCallVariant = m_noCallMode == eNone ? false : a_pVariant->m_bIsNoCall;
        int additionalBasePairCount = static_cast<int>(a_rAlleles[0].length() - a_pVariant->m_refSequence.length());
        
        for(unsigned int k = 0; k < (unsigned int)a_pVariant->m_nZygotCount; k++)
        {
            std::string currentAllele = a_pVariant->GetOriginalAlleleStr(k) + a_rAlleles[0].substr((int)a_rAlleles[0].length() - additionalBasePairCount);
            
            for(unsigned int m = 0; m < a_rAlleles.size(); m++)
            {
                if(currentAllele == a_rAlleles[m])
                {
                    a_rSampleData.m_aGenotype[k] = (int)m;
                    break;
                }
            }
        }
    }
    else if(m_noCallMode == eImplicitNoCall)
        a_rSampleData.m_bIsNoCallVariant = true;
    else
    {
        a_rSampleData.m_aGenotype[0] = 0;
        a_rSampleData.m_aGenotype[1] = 0;
        a_rSampleData.m_bIsNoCallVariant = false;
    }
}

void CMendelianTrioMerger::DoMerge(const CVariant* a_pVarMother,
                                   const CVariant* a_pVarFather,
                                   const CVariant* a_pVarChild,
                                   EMendelianDecision a_decision,
                                   std::vector<SVcfRecord>& a_rRecordList)
{
    
    SVcfRecord vcfrecord;

    vcfrecord.m_nPosition = a_pVarMother == NULL ? (a_pVarFather != NULL ? a_pVarFather->m_nOriginalPos : a_pVarChild->m_nOriginalPos) : a_pVarMother->m_nOriginalPos;
    vcfrecord.m_chrName = a_pVarMother == NULL ? (a_pVarFather != NULL ? a_pVarFather->m_chrName : a_pVarChild->m_chrName) : a_pVarMother->m_chrName;
    vcfrecord.m_mendelianDecision = std::to_string(static_cast<int>(a_decision));
    
    //Left most position of vcf record
    vcfrecord.left = std::min({a_pVarMother != 0 ? a_pVarMother->m_nStartPos : INT_MAX,
                               a_pVarFather != 0 ? a_pVarFather->m_nStartPos : INT_MAX,
                               a_pVarChild  != 0 ? a_pVarChild->m_nStartPos : INT_MAX});
    
    //Right most position of vcf record
    vcfrecord.right = std::max({a_pVarMother != 0 ? a_pVarMother->m_nEndPos : INT_MIN,
                                a_pVarFather != 0 ? a_pVarFather->m_nEndPos : INT_MIN,
                                a_pVarChild  != 0 ? a_pVarChild->m_nEndPos : INT_MIN});
    
    //Fill the alleles part according to trio
    std::vector<std::string> alleles;
    
    //Detect the longest reference sequence
    int maxRefSequenceLength = std::max({a_pVarMother != 0 ? (int)a_pVarMother->m_refSequence.length() : INT_MIN,
                                         a_pVarFather != 0 ? (int)a_pVarFather->m_refSequence.length() : INT_MIN,
                                         a_pVarChild  != 0 ? (int)a_pVarChild->m_refSequence.length() : INT_MIN});
    
    //Push the longest reference sequence as our reference
    if(a_pVarChild != NULL && (int)a_pVarChild->m_refSequence.length() == maxRefSequenceLength)
        alleles.push_back(a_pVarChild->m_refSequence);
    else if (a_pVarMother != NULL && (int)a_pVarMother->m_refSequence.length() == maxRefSequenceLength)
        alleles.push_back(a_pVarMother->m_refSequence);
    else
        alleles.push_back(a_pVarFather->m_refSequence);
    
    //Add child allele
    if(a_pVarChild != NULL)
    {
        //If child reference is shorter than the the records reference, we complete all of its alleles with the missing Base pairs
        int additionalBasePairCount = maxRefSequenceLength - (int)a_pVarChild->m_refSequence.length();
        
        for(unsigned int k= 0; k < (unsigned int)a_pVarChild->m_nAlleleCount; k++)
        {
            std::string childAllele = a_pVarChild->GetOriginalAlleleStr(k) + alleles[0].substr(maxRefSequenceLength-additionalBasePairCount);
            if(childAllele != "" && std::find(alleles.begin(), alleles.end(), childAllele) == alleles.end())
                alleles.push_back(childAllele);
        }
    }
    
    //Add different mother alleles
    if(a_pVarMother != NULL)
    {
        //If mother reference is shorter than the the records reference, we complete all of its alleles with the missing Base pairs
        int additionalBasePairCount = maxRefSequenceLength - (int)a_pVarMother->m_refSequence.length();
        
        for(unsigned int k= 0; k < (unsigned int)a_pVarMother->m_nAlleleCount; k++)
        {
            std::string motherAllele = a_pVarMother->GetOriginalAlleleStr(k) + alleles[0].substr(maxRefSequenceLength-additionalBasePairCount);
            if(motherAllele != "" && std::find(alleles.begin(), alleles.end(), motherAllele) == alleles.end())
                alleles.push_back(motherAllele);
        }
    }
    
    //Add different father alleles
    if(a_pVarFather != NULL)
    {
        //If father reference is shorter than the the records reference, we complete all of its alleles with the missing Base pairs
        int additionalBasePairCount = maxRefSequenceLength - (int)a_pVarFather->m_refSequence.length();
        
        for(unsigned int k= 0; k < (unsigned int)a_pVarFather->m_nAlleleCount; k++)
        {
            std::string fatherAllele = a_pVarFather->GetOriginalAlleleStr(k) + alleles[0].substr(maxRefSequenceLength-additionalBasePairCount);
            if(fatherAllele != "" && std::find(alleles.begin(), alleles.end(), fatherAllele) == alleles.end())
                alleles.push_back(fatherAllele);
        }
    }
    
    //Write alleles to the vcf record
    std::string alleleString = "";
    for(unsigned int k=0; k < alleles.size(); k++)
    {
        if(k != 0)
            alleleString = alleleString + ",";
        alleleString = alleleString + alleles[k];
    }
    vcfrecord.m_alleles = alleleString;
    
    //Add Samples
    SPerSampleData dataMother, dataFather, dataChild;
    AddSample(a_pVarMother, alleles, dataMother);
    AddSample(a_pVarFather, alleles, dataFather);
    AddSample(a_pVarChild, alleles, dataChild);
    
    //Push trip genotypes to the vcfrecord
    vcfrecord.m_aSampleData.push_back(dataMother);
    vcfrecord.m_aSampleData.push_back(dataFather);
    vcfrecord.m_aSampleData.push_back(dataChild);
    
    //Add record to the trio
    a_rRecordList.push_back(vcfrecord);
}

void CMendelianTrioMerger::RegisterMergedLine(EMendelianDecision a_decision, EVariantCategory a_category)
{
    //Ignore the complex types
    if(a_decision == EMendelianDecision::eUnknown)
        return;
    
    switch (a_category)
    {
        case EVariantCategory::eSNP:
            m_logEntry.m_nSNP[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_INSERT_SMALL:
            m_logEntry.m_nInsertSmall[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_INSERT_MEDIUM:
            m_logEntry.m_nInsertMedium[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_INSERT_LARGE:
            m_logEntry.m_nInsertLarge[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_DELETE_SMALL:
            m_logEntry.m_nDeleteSmall[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_DELETE_MEDIUM:
            m_logEntry.m_nDeleteMedium[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_DELETE_LARGE:
            m_logEntry.m_nDeleteLarge[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_COMPLEX_SMALL:
            m_logEntry.m_nComplexSmall[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_COMPLEX_MEDIUM:
            m_logEntry.m_nComplexMedium[a_decision-1]++;
            break;
        case EVariantCategory::eINDEL_COMPLEX_LARGE:
            m_logEntry.m_nComplexLarge[a_decision-1]++;
            break;
            
        default:
            break;
    }
    
}

EMendelianDecision CMendelianTrioMerger::GetMendelianDecision(const CVariant* a_pVarMother, const CVariant* a_pVarFather, const CVariant* a_pVarChild, SChrIdTriplet& a_rTriplet)
{
    EMendelianDecision decision = eUnknown;
    
    if(a_pVarMother != 0 && a_pVarMother->m_variantStatus == eCOMPLEX_SKIPPED)
        decision = eSkipped;
    else if(a_pVarFather != 0 && a_pVarFather->m_variantStatus == eCOMPLEX_SKIPPED)
        decision = eSkipped;
    else if(a_pVarChild != 0 && a_pVarChild->m_variantStatus == eCOMPLEX_SKIPPED)
        decision = eSkipped;
    
    //Check for nocall - Implicit Mode
    else if(m_noCallMode == ENoCallMode::eImplicitNoCall)
    {
        if(a_pVarChild == 0 || a_pVarChild->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallChild;
        else if(a_pVarMother == 0 || a_pVarMother->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else if(a_pVarFather == 0 || a_pVarFather->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
    }
    
    //Check for nocall - Eplicit Mode
    else if(m_noCallMode == ENoCallMode::eExplicitNoCall)
    {
        if(a_pVarChild !=0 && a_pVarChild->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallChild;
        else if(a_pVarMother != 0 && a_pVarMother->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else if(a_pVarFather != 0 && a_pVarFather->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
    }
    
    //If variant is not no-call
    if(decision == eUnknown)
    {
        EMendelianDecision childDecision = a_pVarChild   != 0 ? m_aChildDecisions[a_rTriplet.m_nTripleIndex][a_pVarChild->m_nId]   : EMendelianDecision::eUnknown;
        EMendelianDecision fatherDecision = a_pVarFather != 0 ? m_aFatherDecisions[a_rTriplet.m_nTripleIndex][a_pVarFather->m_nId] : EMendelianDecision::eUnknown;
        EMendelianDecision motherDecision = a_pVarMother != 0 ? m_aMotherDecisions[a_rTriplet.m_nTripleIndex][a_pVarMother->m_nId] : EMendelianDecision::eUnknown;
        
        if(childDecision == eViolation || motherDecision == eViolation || fatherDecision == eViolation)
            decision = eViolation;
        else
            decision = eCompliant;
    }
    
    return decision;
}

void CMendelianTrioMerger::RegisterGenotype(const SVcfRecord& a_rRecord, EVariantCategory a_uCategory, EMendelianDecision a_uDecision)
{
 
    if(a_uCategory == EVariantCategory::eNONE)
        return;
    
    bool bIsMotherMultiAllelic = !(a_rRecord.m_aSampleData[0].m_aGenotype[0] < 2 && a_rRecord.m_aSampleData[0].m_aGenotype[1] < 2);
    bool bIsFatherMultiAllelic = !(a_rRecord.m_aSampleData[1].m_aGenotype[0] < 2 && a_rRecord.m_aSampleData[1].m_aGenotype[1] < 2);
    bool bIsChildMultiAllelic =  !(a_rRecord.m_aSampleData[2].m_aGenotype[0] < 2 && a_rRecord.m_aSampleData[2].m_aGenotype[1] < 2);
    
    int motherGT = a_rRecord.m_aSampleData[0].m_aGenotype[0] + a_rRecord.m_aSampleData[0].m_aGenotype[1];
    int fatherGT = a_rRecord.m_aSampleData[1].m_aGenotype[0] + a_rRecord.m_aSampleData[1].m_aGenotype[1];
    int childGT  = a_rRecord.m_aSampleData[2].m_aGenotype[0] + a_rRecord.m_aSampleData[2].m_aGenotype[1];
    
    int columnNo = motherGT * 9 + fatherGT * 3 + childGT;
    
    if(a_uDecision == EMendelianDecision::eCompliant)
    {
        if(bIsMotherMultiAllelic || bIsFatherMultiAllelic || bIsChildMultiAllelic)
            m_logGenotypes.genotypesCompliant[static_cast<int>(a_uCategory)][27]++;
        else
            m_logGenotypes.genotypesCompliant[static_cast<int>(a_uCategory)][columnNo]++;
        
    }
    
    if(a_uDecision == EMendelianDecision::eViolation)
    {
        if(bIsMotherMultiAllelic || bIsFatherMultiAllelic || bIsChildMultiAllelic)
            m_logGenotypes.genotypesViolation[static_cast<int>(a_uCategory)][27]++;
        else
            m_logGenotypes.genotypesViolation[static_cast<int>(a_uCategory)][columnNo]++;
    }
    
}

void CMendelianTrioMerger::ProcessRefOverlappedRegions(std::vector<SVcfRecord>&  a_rRecordList, std::vector<EMendelianDecision>& a_rRecordDecisionList)
{

    unsigned int recordItr = 0;
    
    while(recordItr < a_rRecordList.size())
    {
        EMendelianDecision curVariantDecision = static_cast<EMendelianDecision>(std::stoi(a_rRecordList[recordItr].m_mendelianDecision));
        
        //Skip consistent variants
        if(curVariantDecision != eCompliant && curVariantDecision != eUnknown)
        {
            //Go Backward and use a second iterator
            int temporaryItr = (int)recordItr-1;
            while(temporaryItr >= 0)
            {
                if(IsOverlap(a_rRecordList[recordItr], a_rRecordList[temporaryItr]))
                {
                    if(a_rRecordDecisionList[temporaryItr] == eCompliant)
                    {
                        a_rRecordDecisionList[temporaryItr] = curVariantDecision;
                        a_rRecordList[temporaryItr].m_mendelianDecision = a_rRecordList[recordItr].m_mendelianDecision;
                    }
                    
                    else if(a_rRecordDecisionList[temporaryItr] == eViolation && curVariantDecision != eViolation)
                    {
                        a_rRecordDecisionList[temporaryItr] = curVariantDecision;
                        a_rRecordList[temporaryItr].m_mendelianDecision = a_rRecordList[recordItr].m_mendelianDecision;
                    }
                    
                    temporaryItr--;
                }
                
                else
                    break;
            }
            
            //Go Forward with second iterator
            temporaryItr = (int)recordItr + 1;
            while(recordItr < a_rRecordList.size())
            {
                if(IsOverlap(a_rRecordList[recordItr], a_rRecordList[temporaryItr]))
                {
                    if(a_rRecordDecisionList[temporaryItr] == eCompliant)
                    {
                        a_rRecordDecisionList[temporaryItr] = curVariantDecision;
                        a_rRecordList[temporaryItr].m_mendelianDecision = a_rRecordList[recordItr].m_mendelianDecision;
                    }
                    
                    else if(a_rRecordDecisionList[temporaryItr] == eViolation && curVariantDecision != eViolation)
                    {
                        a_rRecordDecisionList[temporaryItr] = curVariantDecision;
                        a_rRecordList[temporaryItr].m_mendelianDecision = a_rRecordList[recordItr].m_mendelianDecision;
                    }
                    
                    temporaryItr++;
                }
                else
                    break;
            }
        }
        
        recordItr++;
    }
}
