//
//  CMendelianTrioMerger.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianTrioMerger.h"
#include <algorithm>
#include <iostream>


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

    for(int k = 0; k < (int)a_rVarList.size(); k++)
    {
        if(a_rVarList[k]->m_variantStatus != eNOT_ASSESSED)
        {
            pDecisionList->push_back(a_rDecisionList[k]);
            pVarList->push_back(a_rVarList[k]);
        }
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
    m_vcfWriter.AddHeaderLine("##source= SBG Mendelian Comparison Tool Ver. 1.0 (Beta), 2017");
    
    //ADD MENDELIAN VIOLATION INFO TYPE
    m_vcfWriter.AddHeaderLine("##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Mendelian Violation Decision. (0)-complex, (1)-compliant, (2)-violation (3)-NoCall Parent (4)-NoCall Child \">");
    //m_vcfWriter.AddHeaderLine("##INFO=<ID=POS,Number=1,Type=Integer,Description=\"Complex variant decision position. Shows the position of child variant that is connected\">");
    
    //ADD GT COLUMN
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    
    //ADD CONTIG IDs
    for(int k = 0; k < (int)m_contigs.size(); k++)
        m_vcfWriter.AddHeaderLine("##contig=<ID=" + m_contigs[k].name + ",length=" + std::to_string(m_contigs[k].length) + ">");
    
    //Instead of hardcoding we are reading contig info from input vcfs
    /* [DEPRECATED]
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr1,length=249250621>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr2,length=243199373>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr3,length=198022430>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr4,length=191154276>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr5,length=180915260>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr6,length=171115067>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr7,length=159138663>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr8,length=146364022>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr9,length=141213431>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr10,length=135534747>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr11,length=135006516>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr12,length=133851895>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr13,length=115169878>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr14,length=107349540>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr15,length=102531392>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr16,length=90354753>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr17,length=81195210>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr18,length=78077248>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr19,length=59128983>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr20,length=63025520>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr21,length=48129895>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chr22,length=51304566>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrM,length=16571>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrX,length=155270560>");
    m_vcfWriter.AddHeaderLine("##contig=<ID=chrY,length=59373566>");
    */
     
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
    
    for(int k = 0; k < (int)a_rCommonChromosomes.size(); k++)
    {
        std::cout << "Writing chromosome " << a_rCommonChromosomes[k].m_chrName << std::endl;
        AddRecords(a_rCommonChromosomes[k]);
    }
    
    //Send the logs to the log class
    m_pResultLog->LogDetailedReport(m_logEntry);
    m_pResultLog->LogGenotypeMatrix(m_logGenotypes);
    
    m_vcfWriter.CloseVcf();
}



void CMendelianTrioMerger::AddRecords(SChrIdTriplet& a_rTriplet)
{
    
    int childItr = 0, motherItr = 0, fatherItr = 0;
    
    while(childItr < (int)m_aChildVariants[a_rTriplet.m_nCid].size() || motherItr < (int)m_aMotherVariants[a_rTriplet.m_nMid].size() || fatherItr < (int)m_aFatherVariants[a_rTriplet.m_nFid].size())
    {
        bool mergeMotherChild;
        bool mergeFatherChild;
        
        if(childItr == (int)m_aChildVariants[a_rTriplet.m_nCid].size() || motherItr == (int)m_aMotherVariants[a_rTriplet.m_nMid].size())
            mergeMotherChild = false;
        else
            mergeMotherChild = IsMerge(m_aChildVariants[a_rTriplet.m_nCid][childItr], m_aMotherVariants[a_rTriplet.m_nMid][motherItr]);
        if(childItr == (int)m_aChildVariants[a_rTriplet.m_nCid].size() || fatherItr == (int)m_aFatherVariants[a_rTriplet.m_nFid].size())
            mergeFatherChild = false;
        else
            mergeFatherChild = IsMerge(m_aChildVariants[a_rTriplet.m_nCid][childItr], m_aFatherVariants[a_rTriplet.m_nFid][fatherItr]);
        
        if(mergeFatherChild && mergeMotherChild)
        {
            EMendelianDecision decision = GetMendelianDecision(m_aMotherVariants[a_rTriplet.m_nMid][motherItr],
                                                               m_aFatherVariants[a_rTriplet.m_nFid][fatherItr],
                                                               m_aChildVariants[a_rTriplet.m_nCid][childItr],
                                                               m_aChildDecisions[a_rTriplet.m_nTripleIndex][childItr]);
            
            //Register the line
            EVariantCategory category = RegisterMergedLine(m_aChildVariants[a_rTriplet.m_nCid][childItr], decision);
            
            //Register the genotypes at the 27x27 genotype matrix
            RegisterGenotype(m_aMotherVariants[a_rTriplet.m_nMid][motherItr],
                             m_aFatherVariants[a_rTriplet.m_nFid][fatherItr],
                             m_aChildVariants[a_rTriplet.m_nCid][childItr],
                             decision,
                             category);
            
            DoTripleMerge(a_rTriplet, childItr, fatherItr, motherItr, decision);
            continue;
        }
        
        else if(mergeMotherChild)
        {
            //CHECK IF FATHER IS SMALLER
            if(fatherItr != (int)m_aFatherVariants[a_rTriplet.m_nFid].size() && m_aFatherVariants[a_rTriplet.m_nFid][fatherItr]->m_nOriginalPos < m_aChildVariants[a_rTriplet.m_nCid][childItr]->m_nOriginalPos)
            {
                EMendelianDecision decision = GetMendelianDecision(0, m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], 0, m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr]);
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], decision);

                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(0, m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], 0, decision, category);
                
                DoSingleVar(a_rTriplet, fatherItr, eFATHER, decision);
                
            }
            //DO MOTHER CHILD MERGE
            else
            {
                EMendelianDecision decision = GetMendelianDecision(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], 0, m_aChildVariants[a_rTriplet.m_nCid][childItr], m_aChildDecisions[a_rTriplet.m_nTripleIndex][childItr]);

                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aChildVariants[a_rTriplet.m_nCid][childItr], decision);

                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(m_aMotherVariants[a_rTriplet.m_nMid][motherItr],
                                 0,
                                 m_aChildVariants[a_rTriplet.m_nCid][childItr],
                                 decision,
                                 category);
                
                DoDoubleMerge(a_rTriplet, motherItr, childItr, eMOTHER, eCHILD, decision);
                
            }
            continue;
        }
        
        else if(mergeFatherChild)
        {
            //CHECK IF MOTHER IS SMALLER
            if(motherItr != (int)m_aMotherVariants[a_rTriplet.m_nMid].size() && m_aMotherVariants[a_rTriplet.m_nMid][motherItr]->m_nOriginalPos < m_aChildVariants[a_rTriplet.m_nCid][childItr]->m_nOriginalPos)
            {
                EMendelianDecision decision = GetMendelianDecision(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], 0, 0, m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr]);
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], decision);

                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], 0, 0, decision, category);
                
                
                DoSingleVar(a_rTriplet, motherItr, eMOTHER, decision);
                
            }
            //DO FATHER CHILD MERGE
            else
            {
                EMendelianDecision decision = GetMendelianDecision(0, m_aFatherVariants[a_rTriplet.m_nFid][fatherItr] , m_aChildVariants[a_rTriplet.m_nCid][childItr], m_aChildDecisions[a_rTriplet.m_nTripleIndex][childItr]);

                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aChildVariants[a_rTriplet.m_nCid][childItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(0, m_aFatherVariants[a_rTriplet.m_nFid][fatherItr] , m_aChildVariants[a_rTriplet.m_nCid][childItr], decision, category);
                
    
                DoDoubleMerge(a_rTriplet, fatherItr, childItr, eFATHER, eCHILD, decision);
                
            }
            
            continue;
        }
        
        // No Child variant is merged check for father-mother merge
        bool mergeMotherFather;
        if(motherItr == (int)m_aMotherVariants[a_rTriplet.m_nMid].size() || fatherItr == (int)m_aFatherVariants[a_rTriplet.m_nFid].size())
            mergeMotherFather = false;
        else
            mergeMotherFather = IsMerge(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], m_aFatherVariants[a_rTriplet.m_nFid][fatherItr]);
        
        if(mergeMotherFather)
        {
            //CHECK IF CHILD IS SMALLER
            if(childItr != (int)m_aChildVariants[a_rTriplet.m_nCid].size() && m_aChildVariants[a_rTriplet.m_nCid][childItr]->m_nOriginalPos < m_aMotherVariants[a_rTriplet.m_nMid][motherItr]->m_nOriginalPos)
            {
                EMendelianDecision decision = GetMendelianDecision(0, 0, m_aChildVariants[a_rTriplet.m_nCid][childItr], m_aChildDecisions[a_rTriplet.m_nTripleIndex][childItr]);
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aChildVariants[a_rTriplet.m_nCid][childItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(0, 0, m_aChildVariants[a_rTriplet.m_nCid][childItr], decision, category);
                
                
                DoSingleVar(a_rTriplet, childItr, eCHILD, decision);
            }
            
            //DO MOTHER FATHER MERGE
            else
            {
                EMendelianDecision decision;
                
                if(m_noCallMode == ENoCallMode::eImplicitNoCall)
                    decision = eNoCallChild;
                else if(m_noCallMode == ENoCallMode::eExplicitNoCall && (m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr] == eNoCallParent || m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr] == eNoCallParent))
                    decision = eNoCallParent;
                else if(m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr] == eViolation || m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr] == eViolation)
                    decision = eViolation;
                else if(m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr] == eCompliant || m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr] == eCompliant)
                    decision = eCompliant;
                else if(m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr] == eNoCallParent || m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr] == eNoCallParent)
                    decision = eNoCallParent;
                else
                    decision = eUnknown;
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], 0, decision, category);
                
                DoDoubleMerge(a_rTriplet, motherItr, fatherItr, eMOTHER, eFATHER, decision);
                
            }
            continue;
        }
        
        //There is no merge between three variant print the smallest one
        else
        {
            int motherPos = motherItr != (int)m_aMotherVariants[a_rTriplet.m_nMid].size() ? m_aMotherVariants[a_rTriplet.m_nMid][motherItr]->m_nOriginalPos : INT_MAX;
            int fatherPos = fatherItr != (int)m_aFatherVariants[a_rTriplet.m_nFid].size() ? m_aFatherVariants[a_rTriplet.m_nFid][fatherItr]->m_nOriginalPos : INT_MAX;
            int childPos  = childItr  != (int)m_aChildVariants[a_rTriplet.m_nCid].size()  ? m_aChildVariants[a_rTriplet.m_nCid][childItr]->m_nOriginalPos   : INT_MAX;
        
            if(motherPos <= fatherPos && motherPos <= childPos)
            {
                EMendelianDecision decision = m_noCallMode == ENoCallMode::eImplicitNoCall ? EMendelianDecision::eNoCallChild : m_aMotherDecisions[a_rTriplet.m_nTripleIndex][motherItr];
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(m_aMotherVariants[a_rTriplet.m_nMid][motherItr], 0, 0, decision, category);
                
                DoSingleVar(a_rTriplet, motherItr, eMOTHER, decision);
            }
            else if(fatherPos <= motherPos && fatherPos <= childPos)
            {
                EMendelianDecision decision = m_noCallMode == ENoCallMode::eImplicitNoCall ? EMendelianDecision::eNoCallChild : m_aFatherDecisions[a_rTriplet.m_nTripleIndex][fatherItr];
                
                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(0, m_aFatherVariants[a_rTriplet.m_nFid][fatherItr], 0, decision, category);
                
                DoSingleVar(a_rTriplet, fatherItr, eFATHER, decision);

            }
            else
            {
                EMendelianDecision decision = m_noCallMode == ENoCallMode::eImplicitNoCall ? EMendelianDecision::eNoCallParent : m_aChildDecisions[a_rTriplet.m_nTripleIndex][childItr];

                //Register the line
                EVariantCategory category = RegisterMergedLine(m_aChildVariants[a_rTriplet.m_nCid][childItr], decision);
                
                //Register the genotypes at the 27x27 genotype matrix
                RegisterGenotype(0, 0, m_aChildVariants[a_rTriplet.m_nCid][childItr], decision, category);

                DoSingleVar(a_rTriplet, childItr, eCHILD, decision);
            }
        }
    
    }
}


bool CMendelianTrioMerger::IsMerge(const CVariant* a_pVar1, const CVariant* a_pVar2)
{
    if(a_pVar1->m_nOriginalPos == a_pVar2->m_nOriginalPos && a_pVar1->m_refSequence == a_pVar2->m_refSequence)
        return true;
    else
        return false;
}

void CMendelianTrioMerger::DoTripleMerge(SChrIdTriplet& a_rTriplet, int& a_nChildItr, int& a_nFatherItr, int& a_nMotherItr, EMendelianDecision a_decision)
{
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_nOriginalPos;
    vcfrecord.m_chrName = m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_chrName;
    vcfrecord.m_mendelianDecision = std::to_string(a_decision);
    
    std::vector<std::string> alleles;
    
    //Push reference sequence
    alleles.push_back(m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_refSequence);
    
    //Add non ref child alleles
    for(int k= 0; k < m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_nAlleleCount; k++)
    {
        if(m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_genotype[k] > 0)
            alleles.push_back(m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->GetOriginalAlleleStr(k));
    }
    
    //Add different mother alleles
    for(int k= 0; k < m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_nAlleleCount; k++)
    {
        if(std::find(alleles.begin(), alleles.end(), m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->GetOriginalAlleleStr(k)) == alleles.end())
            alleles.push_back(m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->GetOriginalAlleleStr(k));
    }
    
    //Add different father alleles
    for(int k= 0; k < m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->m_nAlleleCount; k++)
    {
        if(std::find(alleles.begin(), alleles.end(), m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->GetOriginalAlleleStr(k)) == alleles.end())
            alleles.push_back(m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->GetOriginalAlleleStr(k));
    }
    
    std::string alleleString = "";
    for(int k=0; k < (int)alleles.size(); k++)
    {
        if(k != 0)
            alleleString = alleleString + ",";
        alleleString = alleleString + alleles[k];
    }
    vcfrecord.m_alleles = alleleString;
    
    //ADD MOTHER
    SPerSampleData dataMother;
    dataMother.m_nHaplotypeCount = m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_nZygotCount;
    dataMother.m_bIsPhased = m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_bIsPhased; // TODO: This should be altered
    dataMother.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_bIsNoCall;
    for(int k = 0; k < m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_nZygotCount; k++)
    {
        for(int m = 0; m < (int)alleles.size(); m++)
        {
            if(m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->GetOriginalAlleleStr(k) == alleles[m])
            {
                dataMother.m_aGenotype[k] = m;
                break;
            }
        }
    }
    
    //ADD FATHER
    SPerSampleData dataFather;
    dataFather.m_nHaplotypeCount = m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->m_nZygotCount;
    dataFather.m_bIsPhased = m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->m_bIsPhased; // TODO: This should be altered
    dataFather.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->m_bIsNoCall;
    for(int k = 0; k < m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->m_nZygotCount; k++)
    {
        for(int m = 0; m < (int)alleles.size(); m++)
        {
            if(m_aFatherVariants[a_rTriplet.m_nFid][a_nFatherItr]->GetOriginalAlleleStr(k) == alleles[m])
            {
                dataFather.m_aGenotype[k] = m;
                break;
            }
        }
    }
    
    //ADD CHILD
    SPerSampleData dataChild;
    dataChild.m_nHaplotypeCount = m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_nZygotCount;
    dataChild.m_bIsPhased = m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_bIsPhased; // TODO: This should be altered
    dataChild.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->m_bIsNoCall;
    for(int k = 0; k < m_aMotherVariants[a_rTriplet.m_nMid][a_nMotherItr]->m_nZygotCount; k++)
    {
        for(int m = 0; m < (int)alleles.size(); m++)
        {
            if(m_aChildVariants[a_rTriplet.m_nCid][a_nChildItr]->GetOriginalAlleleStr(k) == alleles[m])
            {
                dataChild.m_aGenotype[k] = m;
                break;
            }
        }
    }
    
    //Push trip genotypes to the vcfrecord
    vcfrecord.m_aSampleData.push_back(dataMother);
    vcfrecord.m_aSampleData.push_back(dataFather);
    vcfrecord.m_aSampleData.push_back(dataChild);
    
    //Add record to the trio
    m_vcfWriter.AddMendelianRecord(vcfrecord);
    
    a_nChildItr++;
    a_nFatherItr++;
    a_nMotherItr++;
}


void CMendelianTrioMerger::DoDoubleMerge(SChrIdTriplet& a_rTriplet, int& a_nItr1, int& a_nItr2, EMendelianVcfName a_name1, EMendelianVcfName a_name2, EMendelianDecision a_decision)
{
    const CVariant* pVarMother = (a_name1 == eMOTHER) ? m_aMotherVariants[a_rTriplet.m_nMid][a_nItr1] : ((a_name2 == eMOTHER) ? m_aMotherVariants[a_rTriplet.m_nMid][a_nItr2] : NULL);
    const CVariant* pVarFather = (a_name1 == eFATHER) ? m_aFatherVariants[a_rTriplet.m_nFid][a_nItr1] : ((a_name2 == eFATHER) ? m_aFatherVariants[a_rTriplet.m_nFid][a_nItr2] : NULL);
    const CVariant* pVarChild  = (a_name1 == eCHILD)  ? m_aChildVariants[a_rTriplet.m_nCid][a_nItr1]  : ((a_name2 == eCHILD)  ? m_aChildVariants[a_rTriplet.m_nCid][a_nItr2]  : NULL);
    
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = (pVarChild != NULL) ? pVarChild->m_nOriginalPos : pVarMother->m_nOriginalPos;
    vcfrecord.m_chrName = (pVarChild != NULL) ? pVarChild->m_chrName : pVarMother->m_chrName;
    vcfrecord.m_mendelianDecision = std::to_string(a_decision);
    
    //Fill the alleles part according to trio
    std::vector<std::string> alleles;
    
    //Push reference sequence
    if(pVarChild != NULL)
        alleles.push_back(pVarChild->m_refSequence);
    else
        alleles.push_back(pVarMother->m_refSequence);
    
    
    //Add child allele
    if(pVarChild != NULL)
    {
        for(int k= 0; k < pVarChild->m_nAlleleCount; k++)
        {
            if(std::find(alleles.begin(), alleles.end(), pVarChild->GetOriginalAlleleStr(k)) == alleles.end())
                alleles.push_back(pVarChild->GetOriginalAlleleStr(k));
        }
    }
    
    //Add different mother alleles
    if(pVarMother != NULL)
    {
        for(int k= 0; k < pVarMother->m_nAlleleCount; k++)
        {
            if(std::find(alleles.begin(), alleles.end(), pVarMother->GetOriginalAlleleStr(k)) == alleles.end())
                alleles.push_back(pVarMother->GetOriginalAlleleStr(k));
        }
    }
    
    //Add different father alleles
    if(pVarFather != NULL)
    {
        for(int k= 0; k < pVarFather->m_nAlleleCount; k++)
        {
            if(std::find(alleles.begin(), alleles.end(), pVarFather->GetOriginalAlleleStr(k)) == alleles.end())
                alleles.push_back(pVarFather->GetOriginalAlleleStr(k));
        }
    }
    
    std::string alleleString = "";
    for(int k=0; k < (int)alleles.size(); k++)
    {
        if(k != 0)
            alleleString = alleleString + ",";
        alleleString = alleleString + alleles[k];
    }
    vcfrecord.m_alleles = alleleString;

    //ADD MOTHER
    SPerSampleData dataMother;
    if(pVarMother != NULL)
    {
        dataMother.m_nHaplotypeCount = pVarMother->m_nZygotCount;
        dataMother.m_bIsPhased = pVarMother->m_bIsPhased; // TODO: This should be altered
        dataMother.m_bIsNoCallVariant = m_noCallMode == eNone ? false : pVarMother->m_bIsNoCall;
        for(int k = 0; k < pVarMother->m_nZygotCount; k++)
        {
            for(int m = 0; m < (int)alleles.size(); m++)
            {
                if(pVarMother->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataMother.m_aGenotype[k] = m;
                    break;
                }
            }
        }
    }
    else if(m_noCallMode == eImplicitNoCall)
        dataMother.m_bIsNoCallVariant = true;
    else
    {
        dataMother.m_aGenotype[0] = 0;
        dataMother.m_aGenotype[1] = 0;
        dataMother.m_bIsNoCallVariant = false;
    }
    
    //ADD FATHER
    SPerSampleData dataFather;
    if(pVarFather != NULL)
    {
        dataFather.m_nHaplotypeCount = pVarFather->m_nZygotCount;
        dataFather.m_bIsPhased = pVarFather->m_bIsPhased; // TODO: This should be altered
        dataFather.m_bIsNoCallVariant = m_noCallMode == eNone ? false : pVarFather->m_bIsNoCall;
        for(int k = 0; k < pVarFather->m_nZygotCount; k++)
        {
            for(int m = 0; m < (int)alleles.size(); m++)
            {
                if(pVarFather->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataFather.m_aGenotype[k] = m;
                    break;
                }
            }
        }
    }
    else if(m_noCallMode == eImplicitNoCall)
        dataFather.m_bIsNoCallVariant = true;
    else
    {
        dataFather.m_aGenotype[0] = 0;
        dataFather.m_aGenotype[1] = 0;
        dataFather.m_bIsNoCallVariant = false;
    }
    
    //ADD CHILD
    SPerSampleData dataChild;
    if(pVarChild != NULL)
    {
        dataChild.m_nHaplotypeCount = pVarChild->m_nZygotCount;
        dataChild.m_bIsPhased = pVarChild->m_bIsPhased; // TODO: This should be altered
        dataChild.m_bIsNoCallVariant = m_noCallMode == eNone ? false : pVarChild->m_bIsNoCall;
        for(int k = 0; k < pVarChild->m_nZygotCount; k++)
        {
            for(int m = 0; m < (int)alleles.size(); m++)
            {
                if(pVarChild->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataChild.m_aGenotype[k] = m;
                    break;
                }
            }
        }
    }
    else if(m_noCallMode == eImplicitNoCall)
        dataChild.m_bIsNoCallVariant = true;
    else
    {
        dataChild.m_aGenotype[0] = 0;
        dataChild.m_aGenotype[1] = 0;
        dataChild.m_bIsNoCallVariant = false;
    }
    
    //Push trip genotypes to the vcfrecord
    vcfrecord.m_aSampleData.push_back(dataMother);
    vcfrecord.m_aSampleData.push_back(dataFather);
    vcfrecord.m_aSampleData.push_back(dataChild);
    
    //Add record to the trio
    m_vcfWriter.AddMendelianRecord(vcfrecord);
    
    //Increment counters
    a_nItr1++;
    a_nItr2++;
    
}


void CMendelianTrioMerger::DoSingleVar(SChrIdTriplet& a_rTriplet, int& a_nItr, EMendelianVcfName a_name, EMendelianDecision a_decision)
{
    const CVariant* pVariant = (a_name == eMOTHER) ? m_aMotherVariants[a_rTriplet.m_nMid][a_nItr] : ((a_name == eFATHER) ? m_aFatherVariants[a_rTriplet.m_nFid][a_nItr] : m_aChildVariants[a_rTriplet.m_nCid][a_nItr]);
    
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = pVariant->m_nOriginalPos;
    vcfrecord.m_chrName = pVariant->m_chrName;
    vcfrecord.m_mendelianDecision = std::to_string(a_decision);
    vcfrecord.m_alleles = pVariant->m_allelesStr;
    
    SPerSampleData dataMother;
    SPerSampleData dataFather;
    SPerSampleData dataChild;
    
    if(a_name == eMOTHER)
    {
        dataMother.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataMother.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        dataMother.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : pVariant->m_bIsNoCall;
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataMother.m_aGenotype[k] = pVariant->m_genotype[k];
    }
    else if(m_noCallMode == ENoCallMode::eNone || m_noCallMode == ENoCallMode::eExplicitNoCall)
    {
        dataMother.m_aGenotype[0] = 0;
        dataMother.m_aGenotype[1] = 0;
        dataMother.m_bIsNoCallVariant = false;
    }
    else
        dataMother.m_bIsNoCallVariant = true;
    
    if(a_name == eFATHER)
    {
        dataFather.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataFather.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        dataFather.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : pVariant->m_bIsNoCall;
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataFather.m_aGenotype[k] = pVariant->m_genotype[k];
    
    }
    else if(m_noCallMode == ENoCallMode::eNone || m_noCallMode == ENoCallMode::eExplicitNoCall)
    {
        dataFather.m_aGenotype[0] = 0;
        dataFather.m_aGenotype[1] = 0;
        dataFather.m_bIsNoCallVariant = false;
    }
    else
        dataFather.m_bIsNoCallVariant = true;

    if(a_name == eCHILD)
    {
        dataChild.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataChild.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        dataChild.m_bIsNoCallVariant = m_noCallMode == ENoCallMode::eNone ? false : pVariant->m_bIsNoCall;
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataChild.m_aGenotype[k] = pVariant->m_genotype[k];
    }
    else if(m_noCallMode == ENoCallMode::eNone || m_noCallMode == ENoCallMode::eExplicitNoCall)
    {
        dataChild.m_aGenotype[0] = 0;
        dataChild.m_aGenotype[1] = 0;
        dataChild.m_bIsNoCallVariant = false;
    }
    else
        dataChild.m_bIsNoCallVariant = true;
    
    //Push trip genotypes to the vcfrecord
    vcfrecord.m_aSampleData.push_back(dataMother);
    vcfrecord.m_aSampleData.push_back(dataFather);
    vcfrecord.m_aSampleData.push_back(dataChild);
    
    //Add record to the trio
    m_vcfWriter.AddMendelianRecord(vcfrecord);
    
    //Increment the counter
    a_nItr++;
    
}

EVariantCategory CMendelianTrioMerger::RegisterMergedLine(const CVariant* a_pVariant, EMendelianDecision a_decision)
{
    //Ignore the complex types
    if(a_decision == EMendelianDecision::eUnknown)
        return EVariantCategory::eNONE;
    
    EVariantCategory category = a_pVariant->GetVariantCategory();
    
    switch (category)
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
    
    return category;
}

EMendelianDecision CMendelianTrioMerger::GetMendelianDecision(const CVariant* a_pVarMother, const CVariant* a_pVarFather, const CVariant* a_pVarChild, EMendelianDecision a_initDecision)
{
    EMendelianDecision decision;

    if(m_noCallMode == ENoCallMode::eImplicitNoCall)
    {
        if(a_pVarChild == 0 || a_pVarChild->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallChild;
        else if(a_pVarMother == 0 || a_pVarMother->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else if(a_pVarFather == 0 || a_pVarFather->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else
            decision = a_initDecision;
    }
    
    else if(m_noCallMode == ENoCallMode::eExplicitNoCall)
    {
        if(a_pVarChild !=0 && a_pVarChild->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallChild;
        else if(a_pVarMother != 0 && a_pVarMother->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else if(a_pVarFather != 0 && a_pVarFather->m_bIsNoCall)
            decision = EMendelianDecision::eNoCallParent;
        else
            decision = a_initDecision;
    }
    
    else // m_noCallMode == eNone
    {
        decision = a_initDecision;
    }
    
    return decision;
}


void CMendelianTrioMerger::RegisterGenotype(const CVariant* a_pMother,
                                            const CVariant* a_pFather,
                                            const CVariant* a_pChild,
                                            EMendelianDecision a_decision,
                                            EVariantCategory a_category)
{
    if(a_category == EVariantCategory::eNONE)
        return;
    
    bool bIsMotherMultiAllelic = a_pMother == 0 ? false : !(a_pMother->m_genotype[0] < 2 && a_pMother->m_genotype[1] < 2);
    bool bIsFatherMultiAllelic = a_pFather == 0 ? false : !(a_pFather->m_genotype[0] < 2 && a_pFather->m_genotype[1] < 2);
    bool bIsChildMultiAllelic = a_pChild == 0 ? false   : !(a_pChild->m_genotype[0]  < 2 && a_pChild->m_genotype[1]  < 2);
    
    int motherGT = a_pMother == 0 ? 0 : a_pMother->m_genotype[0] + a_pMother->m_genotype[1];
    int fatherGT = a_pFather == 0 ? 0 : a_pFather->m_genotype[0] + a_pFather->m_genotype[1];
    int childGT = a_pChild == 0 ? 0   : a_pChild->m_genotype[0]  + a_pChild->m_genotype[1];
    
    int columnNo = motherGT * 9 + fatherGT * 3 + childGT;

    if(a_decision == EMendelianDecision::eCompliant)
    {
        if(bIsMotherMultiAllelic || bIsFatherMultiAllelic || bIsChildMultiAllelic)
            m_logGenotypes.genotypesCompliant[static_cast<int>(a_category)][27]++;
        else
            m_logGenotypes.genotypesCompliant[static_cast<int>(a_category)][columnNo]++;
        
    }
            
    if(a_decision == EMendelianDecision::eViolation)
    {
        if(bIsMotherMultiAllelic || bIsFatherMultiAllelic || bIsChildMultiAllelic)
            m_logGenotypes.genotypesViolation[static_cast<int>(a_category)][27]++;
        else
            m_logGenotypes.genotypesViolation[static_cast<int>(a_category)][columnNo]++;
    }
    
}










