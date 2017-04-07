//
//  CMendelianTrioMerger.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#include "CMendelianTrioMerger.h"
#include <algorithm>


void CMendelianTrioMerger::SetTrioPath(const std::string& a_nTrioPath)
{
    m_trioPath = std::string(a_nTrioPath);
}

void CMendelianTrioMerger::SetVariants(int a_nChromosomeId, EMendelianVcfName a_vcfName, const std::vector<const CVariant*>& a_rVarList)
{
    std::vector<const CVariant*>* pVarList;
    
    switch (a_vcfName)
    {
        case eCHILD:
            pVarList = &m_aChildVariants[a_nChromosomeId];
            break;
        case eFATHER:
            pVarList = &m_aFatherVariants[a_nChromosomeId];
            break;
        case eMOTHER:
            pVarList = &m_aMotherVariants[a_nChromosomeId];
            break;
        default:
            break;
    }
    
    for(const CVariant* pVar : a_rVarList)
        pVarList->push_back(pVar);
}

void CMendelianTrioMerger::SetDecisions(int a_nChromosomeId, EMendelianVcfName a_vcfName, const std::vector<EMendelianDecision>& a_rDecisionList)
{
    std::vector<EMendelianDecision>* pDecisionList;
    
    switch (a_vcfName)
    {
        case eCHILD:
            pDecisionList = &m_aChildDecisions[a_nChromosomeId];
            break;
        case eFATHER:
            pDecisionList = &m_aFatherDecisions[a_nChromosomeId];
            break;
        case eMOTHER:
            pDecisionList = &m_aMotherDecisions[a_nChromosomeId];
            break;
        default:
            break;
    }
    
    for(EMendelianDecision decision : a_rDecisionList)
        pDecisionList->push_back(decision);
}


void CMendelianTrioMerger::FillHeader()
{
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= SBG Mendelian Comparison Tool Ver. 1.0 (Beta), 2017");
    
    //ADD MENDELIAN VIOLATION INFO TYPE
    m_vcfWriter.AddHeaderLine("##INFO=<ID=MD,Number=1,Type=String,Description=\"Mendelian Violation Decision.\">");
    
    //ADD GT COLUMN
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    
    //ADD CONTIG IDs
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
    
    //ADD REQUIRED SAMPLES
    m_vcfWriter.AddSampleName("MOTHER");
    m_vcfWriter.AddSampleName("FATHER");
    m_vcfWriter.AddSampleName("CHILD");
    
    //CLOSE HEADER
    m_vcfWriter.WriteHeaderToVcf();
    
}

void CMendelianTrioMerger::GenerateTrioVcf()
{
    m_vcfWriter.CreateVcf(m_trioPath.c_str());
    
    //Fill the header
    FillHeader();
    
    for(int k = 0; k < CHROMOSOME_COUNT; k++)
    {
        if(m_aChildVariants[k].size() != 0)
            AddRecords(k);
    }

    m_vcfWriter.CloseVcf();
}



void CMendelianTrioMerger::AddRecords(int a_nChromosomeId)
{
    
    int childItr = 0, motherItr = 0, fatherItr = 0;
    
    while(childItr < m_aChildVariants[a_nChromosomeId].size() || motherItr < m_aMotherVariants[a_nChromosomeId].size() || fatherItr < m_aFatherVariants[a_nChromosomeId].size())
    {
        bool mergeMotherChild;
        bool mergeFatherChild;
        
        if(childItr == m_aChildVariants[a_nChromosomeId].size() || motherItr == m_aMotherVariants[a_nChromosomeId].size())
            mergeMotherChild = false;
        else
            mergeMotherChild = IsMerge(m_aChildVariants[a_nChromosomeId][childItr], m_aMotherVariants[a_nChromosomeId][motherItr]);
        if(childItr == m_aChildVariants[a_nChromosomeId].size() || fatherItr == m_aFatherVariants[a_nChromosomeId].size())
            mergeFatherChild = false;
        else
            mergeFatherChild = IsMerge(m_aChildVariants[a_nChromosomeId][childItr], m_aFatherVariants[a_nChromosomeId][fatherItr]);
        
        if(mergeFatherChild && mergeMotherChild)
        {
            DoTripleMerge(a_nChromosomeId, childItr, fatherItr, motherItr);
            continue;
        }
        
        else if(mergeMotherChild)
        {
            //CHECK IF FATHER IS SMALLER
            if(fatherItr != m_aFatherVariants[a_nChromosomeId].size() && m_aFatherVariants[a_nChromosomeId][fatherItr]->m_nOriginalPos < m_aChildVariants[a_nChromosomeId][childItr]->m_nOriginalPos)
                DoSingleVar(a_nChromosomeId, fatherItr, eFATHER, m_aFatherDecisions[a_nChromosomeId][fatherItr]);
            //DO MOTHER CHILD MERGE
            else
                DoDoubleMerge(a_nChromosomeId, motherItr, childItr, eMOTHER, eCHILD, m_aChildDecisions[a_nChromosomeId][childItr]);
            
            continue;
        }
        
        else if(mergeFatherChild)
        {
            //CHECK IF MOTHER IS SMALLER
            if(motherItr != m_aMotherVariants[a_nChromosomeId].size() && m_aMotherVariants[a_nChromosomeId][motherItr]->m_nOriginalPos < m_aChildVariants[a_nChromosomeId][childItr]->m_nOriginalPos)
                DoSingleVar(a_nChromosomeId, motherItr, eMOTHER, m_aMotherDecisions[a_nChromosomeId][motherItr]);
            //DO MOTHER CHILD MERGE
            else
                DoDoubleMerge(a_nChromosomeId, fatherItr, childItr, eFATHER, eCHILD, m_aChildDecisions[a_nChromosomeId][childItr]);
            
            continue;
        }
        
        // No Child variant is merged check for father-mother merge
        bool mergeMotherFather;
        if(motherItr == m_aMotherVariants[a_nChromosomeId].size() || fatherItr == m_aFatherVariants[a_nChromosomeId].size())
            mergeMotherFather = false;
        else
            mergeMotherFather = IsMerge(m_aMotherVariants[a_nChromosomeId][motherItr], m_aFatherVariants[a_nChromosomeId][fatherItr]);
        
        if(mergeMotherFather)
        {
            //CHECK IF CHILD IS SMALLER
            if(childItr != m_aChildVariants[a_nChromosomeId].size() && m_aChildVariants[a_nChromosomeId][childItr]->m_nOriginalPos < m_aMotherVariants[a_nChromosomeId][motherItr]->m_nOriginalPos)
                DoSingleVar(a_nChromosomeId, childItr, eCHILD, m_aChildDecisions[a_nChromosomeId][childItr]);
            //DO MOTHER CHILD MERGE
            else
            {
                EMendelianDecision dec;
                if(m_aMotherDecisions[a_nChromosomeId][motherItr] == eViolation || m_aFatherDecisions[a_nChromosomeId][fatherItr] == eViolation)
                    dec = eViolation;
                else if(m_aMotherDecisions[a_nChromosomeId][motherItr] == eCompliant || m_aFatherDecisions[a_nChromosomeId][fatherItr] == eCompliant)
                    dec = eCompliant;
                else
                    dec = eUnknown;
                
                DoDoubleMerge(a_nChromosomeId, motherItr, fatherItr, eMOTHER, eFATHER, dec);
            }
            continue;
        }
        
        //There is no merge between three variant print the smallest one
        else
        {
            int motherPos = motherItr != m_aMotherVariants[a_nChromosomeId].size() ? m_aMotherVariants[a_nChromosomeId][motherItr]->m_nOriginalPos : INT_MAX;
            int fatherPos = fatherItr != m_aFatherVariants[a_nChromosomeId].size() ? m_aFatherVariants[a_nChromosomeId][fatherItr]->m_nOriginalPos : INT_MAX;
            int childPos  = childItr  != m_aChildVariants[a_nChromosomeId].size()  ? m_aChildVariants[a_nChromosomeId][childItr]->m_nOriginalPos   : INT_MAX;
        
            if(motherPos <= fatherPos && motherPos <= childPos)
                DoSingleVar(a_nChromosomeId, motherItr, eMOTHER, m_aMotherDecisions[a_nChromosomeId][motherItr]);
            else if(fatherPos <= motherPos && fatherPos <= childPos)
                DoSingleVar(a_nChromosomeId, fatherItr, eFATHER, m_aFatherDecisions[a_nChromosomeId][fatherItr]);
            else
                DoSingleVar(a_nChromosomeId, childItr, eCHILD, m_aChildDecisions[a_nChromosomeId][childItr]);
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

void CMendelianTrioMerger::DoTripleMerge(int a_nChromosomeId, int& a_nChildItr, int& a_nFatherItr, int& a_nMotherItr)
{
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_nOriginalPos;
    vcfrecord.m_nChrId = m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_nChrId;
    vcfrecord.m_mendelianDecision = m_aChildDecisions[a_nChromosomeId][a_nChildItr] == eCompliant ? "compliant" : (m_aChildDecisions[a_nChromosomeId][a_nChildItr] == eViolation ? "violation" : "complex");
    
    std::vector<std::string> alleles;
    
    //Push reference sequence
    alleles.push_back(m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_refSequence);
    
    //Add non ref child alleles
    for(int k= 0; k < m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_nAlleleCount; k++)
    {
        if(m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_genotype[k] != 0)
            alleles.push_back(m_aChildVariants[a_nChromosomeId][a_nChildItr]->GetOriginalAlleleStr(k));
    }
    
    //Add different mother alleles
    for(int k= 0; k < m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->m_nAlleleCount; k++)
    {
        if(std::find(alleles.begin(), alleles.end(), m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->GetOriginalAlleleStr(k)) == alleles.end())
            alleles.push_back(m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->GetOriginalAlleleStr(k));
    }
    
    //Add different father alleles
    for(int k= 0; k < m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->m_nAlleleCount; k++)
    {
        if(std::find(alleles.begin(), alleles.end(), m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->GetOriginalAlleleStr(k)) == alleles.end())
            alleles.push_back(m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->GetOriginalAlleleStr(k));
    }
    
    std::string alleleString = "";
    for(int k=0; k < alleles.size(); k++)
    {
        if(k != 0)
            alleleString = alleleString + ",";
        alleleString = alleleString + alleles[k];
    }
    vcfrecord.m_alleles = alleleString;
    
    //ADD MOTHER
    SPerSampleData dataMother;
    dataMother.m_nHaplotypeCount = m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->m_nZygotCount;
    dataMother.m_bIsPhased = m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->m_bIsPhased; // TODO: This should be altered
    for(int k = 0; k < m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->m_nZygotCount; k++)
    {
        for(int m = 0; m < alleles.size(); m++)
        {
            if(m_aMotherVariants[a_nChromosomeId][a_nMotherItr]->GetOriginalAlleleStr(k) == alleles[m])
            {
                dataMother.m_aGenotype[k] = m;
                break;
            }
        }
    }
    
    //ADD FATHER
    SPerSampleData dataFather;
    dataFather.m_nHaplotypeCount = m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->m_nZygotCount;
    dataFather.m_bIsPhased = m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->m_bIsPhased; // TODO: This should be altered
    for(int k = 0; k < m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->m_nZygotCount; k++)
    {
        for(int m = 0; m < alleles.size(); m++)
        {
            if(m_aFatherVariants[a_nChromosomeId][a_nFatherItr]->GetOriginalAlleleStr(k) == alleles[m])
            {
                dataFather.m_aGenotype[k] = m;
                break;
            }
        }
    }
    
    //ADD CHILD
    SPerSampleData dataChild;
    dataChild.m_nHaplotypeCount = m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_nZygotCount;
    dataChild.m_bIsPhased = m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_bIsPhased; // TODO: This should be altered
    for(int k = 0; k < m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_nZygotCount; k++)
        dataChild.m_aGenotype[k] = m_aChildVariants[a_nChromosomeId][a_nChildItr]->m_genotype[k];
    
    
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


void CMendelianTrioMerger::DoDoubleMerge(int a_nChromosomeId, int& a_nItr1, int& a_nItr2, EMendelianVcfName a_name1, EMendelianVcfName a_name2, EMendelianDecision a_decision)
{
    const CVariant* pVarMother = (a_name1 == eMOTHER) ? m_aMotherVariants[a_nChromosomeId][a_nItr1] : ((a_name2 == eMOTHER) ? m_aMotherVariants[a_nChromosomeId][a_nItr2] : NULL);
    const CVariant* pVarFather = (a_name1 == eFATHER) ? m_aFatherVariants[a_nChromosomeId][a_nItr1] : ((a_name2 == eFATHER) ? m_aFatherVariants[a_nChromosomeId][a_nItr2] : NULL);
    const CVariant* pVarChild  = (a_name1 == eCHILD)  ? m_aChildVariants[a_nChromosomeId][a_nItr1]  : ((a_name2 == eCHILD)  ? m_aChildVariants[a_nChromosomeId][a_nItr2]  : NULL);
    
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = (pVarChild != NULL) ? pVarChild->m_nOriginalPos : pVarMother->m_nOriginalPos;
    vcfrecord.m_nChrId = (pVarChild != NULL) ? pVarChild->m_nChrId : pVarMother->m_nChrId;
    vcfrecord.m_mendelianDecision = a_decision == eCompliant ? "compliant" : (a_decision == eViolation ? "violation" : "complex");
    
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
    for(int k=0; k < alleles.size(); k++)
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
        for(int k = 0; k < pVarMother->m_nZygotCount; k++)
        {
            for(int m = 0; m < alleles.size(); m++)
            {
                if(pVarMother->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataMother.m_aGenotype[k] = m;
                    break;
                }
            }
        }
    }
    
    //ADD FATHER
    SPerSampleData dataFather;
    if(pVarFather != NULL)
    {
        dataFather.m_nHaplotypeCount = pVarFather->m_nZygotCount;
        dataFather.m_bIsPhased = pVarFather->m_bIsPhased; // TODO: This should be altered
        for(int k = 0; k < pVarFather->m_nZygotCount; k++)
        {
            for(int m = 0; m < alleles.size(); m++)
            {
                if(pVarFather->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataFather.m_aGenotype[k] = m;
                    break;
                }
            }
        }
    }
    
    //ADD CHILD
    SPerSampleData dataChild;
    if(pVarChild != NULL)
    {
        dataChild.m_nHaplotypeCount = pVarChild->m_nZygotCount;
        dataChild.m_bIsPhased = pVarChild->m_bIsPhased; // TODO: This should be altered
        for(int k = 0; k < pVarChild->m_nZygotCount; k++)
        {
            for(int m = 0; m < alleles.size(); m++)
            {
                if(pVarChild->GetOriginalAlleleStr(k) == alleles[m])
                {
                    dataChild.m_aGenotype[k] = m;
                    break;
                }
            }
        }
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


void CMendelianTrioMerger::DoSingleVar(int a_nChromosomeId, int& a_nItr, EMendelianVcfName a_name, EMendelianDecision a_decision)
{
    const CVariant* pVariant = (a_name == eMOTHER) ? m_aMotherVariants[a_nChromosomeId][a_nItr] : ((a_name == eFATHER) ? m_aFatherVariants[a_nChromosomeId][a_nItr] : m_aChildVariants[a_nChromosomeId][a_nItr]);
    
    SVcfRecord vcfrecord;
    
    vcfrecord.m_nPosition = pVariant->m_nOriginalPos;
    vcfrecord.m_nChrId = pVariant->m_nChrId;
    vcfrecord.m_mendelianDecision = a_decision == eCompliant ? "compliant" : (a_decision == eViolation ? "violation" : "complex");
    vcfrecord.m_alleles = pVariant->m_allelesStr;
    
    SPerSampleData dataMother;
    SPerSampleData dataFather;
    SPerSampleData dataChild;
    
    if(a_name == eMOTHER)
    {
        dataMother.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataMother.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataMother.m_aGenotype[k] = pVariant->m_genotype[k];
    }
    
    else if(a_name == eFATHER)
    {
        dataFather.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataFather.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataFather.m_aGenotype[k] = pVariant->m_genotype[k];
    
    }
    
    else
    {
        dataChild.m_nHaplotypeCount = pVariant->m_nZygotCount;
        dataChild.m_bIsPhased = pVariant->m_bIsPhased; // TODO: This should be altered
        for(int k = 0; k < pVariant->m_nZygotCount; k++)
            dataChild.m_aGenotype[k] = pVariant->m_genotype[k];
    }
    
    //Push trip genotypes to the vcfrecord
    vcfrecord.m_aSampleData.push_back(dataMother);
    vcfrecord.m_aSampleData.push_back(dataFather);
    vcfrecord.m_aSampleData.push_back(dataChild);
    
    //Add record to the trio
    m_vcfWriter.AddMendelianRecord(vcfrecord);
    
    //Increment the counter
    a_nItr++;
}























