//
//  CLBLComparisonTool.cpp
//  gatk_sbgd_trio_test
//
//  Created by Berke.Toptas on 11/17/17.
//  Copyright © 2017 SBGD. All rights reserved.
//

#include "CLBLComparisonTool.h"
#include <iostream>
#include "CSimplePEDParser.h"



void CLBLComparisonTool::SetInputTrio(const std::string& a_rInputTrioPath, const std::string& a_rPedigreeFile)
{
    //Parse pedigree
    CSimplePEDParser pedParser;
    pedParser.ParsePedigree(a_rPedigreeFile);
    
    //Read vcf from given input path
    m_vcfReader.Open(a_rInputTrioPath.c_str());
    std::vector<std::string> sampleNames;
    
    //Get sample names ordered as mother father and child
    m_vcfReader.GetSampleNames(sampleNames);
    std::vector<std::string> orderedSamples = pedParser.GetIdsMFC(sampleNames[0], sampleNames[1], sampleNames[2]);
    
    if(orderedSamples[0] == sampleNames[0])
        m_motherIndex = 0;
    if(orderedSamples[0] == sampleNames[1])
        m_motherIndex = 1;
    if(orderedSamples[0] == sampleNames[2])
        m_motherIndex = 2;
    
    if(orderedSamples[1] == sampleNames[0])
        m_fatherIndex = 0;
    if(orderedSamples[1] == sampleNames[1])
        m_fatherIndex = 1;
    if(orderedSamples[1] == sampleNames[2])
        m_fatherIndex = 2;
    
    if(orderedSamples[2] == sampleNames[0])
        m_childIndex = 0;
    if(orderedSamples[2] == sampleNames[1])
        m_childIndex = 1;
    if(orderedSamples[2] == sampleNames[2])
        m_childIndex = 2;
    
    m_bIsTrimBeginningFirst = false;
}

void CLBLComparisonTool::WriteOutputTrio(const std::string& a_rOutputTrioPath)
{
    //Create VCF at the given output path
    m_vcfWriter.CreateVcf(a_rOutputTrioPath.c_str());
    
    //INIT VCF HEADER
    m_vcfWriter.InitHeader();
    m_vcfWriter.AddHeaderLine("##source= VBT Validation - Naive Mendelian Decision Annotator");
    
    //ADD GT AND MD COLUMN
    m_vcfWriter.AddHeaderLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    m_vcfWriter.AddHeaderLine("##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Mendelian Violation Decision. (0)-complex, (1)-compliant, (2)-violation (3)-NoCall Parent (4)-NoCall Child \">");
    
    //GET THE CONTIGS FROM VCF FILE
    std::vector<SVcfContig> contigs = m_vcfReader.GetContigs();
    
    //ADD CONTIG IDs
    for(int k = 0; k < (int)contigs.size(); k++)
        m_vcfWriter.AddHeaderLine("##contig=<ID=" + contigs[k].name + ",length=" + std::to_string(contigs[k].length) + ">");
    
    //ADD SAMPLE IDs
    std::vector<std::string> sampleNames;
    m_vcfReader.GetSampleNames(sampleNames);
    for(int k = 0; k < sampleNames.size(); k++)
        m_vcfWriter.AddSampleName(sampleNames[k]);
    
    //WRITE HEADER TO VCF
    m_vcfWriter.WriteHeaderToVcf();

    CVariant truthVars[3];
        
    while(m_vcfReader.GetNextRecordMultiSample(truthVars, m_bIsTrimBeginningFirst))
    {
        //Skip complex variants
        if(truthVars[m_motherIndex].m_alleles[0].m_sequence == "*" || truthVars[m_motherIndex].m_alleles[1].m_sequence == "*")
            continue;
        
        if(truthVars[m_fatherIndex].m_alleles[0].m_sequence == "*" || truthVars[m_fatherIndex].m_alleles[1].m_sequence == "*")
            continue;
        
        if(truthVars[m_childIndex].m_alleles[0].m_sequence == "*" || truthVars[m_childIndex].m_alleles[1].m_sequence == "*")
            continue;
        
        bool c1 = truthVars[m_motherIndex].m_genotype[0] == truthVars[m_childIndex].m_genotype[0] || truthVars[m_motherIndex].m_genotype[1] == truthVars[m_childIndex].m_genotype[0];
        bool c2 = truthVars[m_motherIndex].m_genotype[0] == truthVars[m_childIndex].m_genotype[1] || truthVars[m_motherIndex].m_genotype[1] == truthVars[m_childIndex].m_genotype[1];
        bool c3 = truthVars[m_fatherIndex].m_genotype[0] == truthVars[m_childIndex].m_genotype[0] || truthVars[m_fatherIndex].m_genotype[1] == truthVars[m_childIndex].m_genotype[0];
        bool c4 = truthVars[m_fatherIndex].m_genotype[0] == truthVars[m_childIndex].m_genotype[1] || truthVars[m_fatherIndex].m_genotype[1] == truthVars[m_childIndex].m_genotype[1];
        
        std::string decision = "0";
        
        //Set MD as nocall
        if(truthVars[0].m_bIsNoCall || truthVars[1].m_bIsNoCall || truthVars[2].m_bIsNoCall)
            decision = "4";
        
        //Set MD as consistent
        else if((c1 == true && c4 == true) || (c2 == true && c3 == true))
            decision = "1";
        
        //Set MD as violation
        else
            decision = "2";
        
        
        SVcfRecord record;
        SPerSampleData data[3];
        
        
        
        record.m_chrName = truthVars[0].m_chrName;
        record.m_alleles = truthVars[0].m_allelesStr;
        record.m_nPosition = truthVars[0].m_nOriginalPos;
        record.m_mendelianDecision = decision;
        
        for(int k = 0; k < 3; k++)
        {
            data[k].m_nHaplotypeCount = truthVars[k].m_nZygotCount;
            data[k].m_bIsPhased = false;
            for(int m = 0; m < truthVars[k].m_nZygotCount; m++)
                data[k].m_aGenotype[m] = truthVars[k].m_genotype[m];
            data[k].m_bIsNoCallVariant = truthVars[k].m_bIsNoCall;
            
            record.m_aSampleData.push_back(data[k]);
        }
        
        m_vcfWriter.AddMendelianRecord(record);
        
    }
    
    m_vcfWriter.CloseVcf();
}

