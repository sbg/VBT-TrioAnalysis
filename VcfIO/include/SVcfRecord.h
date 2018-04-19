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
 *  SVcfRecord.h
 *  VCFComparison
 *
 * Created by Berke Cagkan Toptas on 4/19/18.
 *
 */

#ifndef _S_VCF_RECORD_H_
#define _S_VCF_RECORD_H_

#include "SInfo.h"
#include <vector>
#include <string>

/**
 * @brief Stores FORMAT column of vcf record
 *
 */
struct SPerSampleData
{
    SPerSampleData()
    {
        m_decisionBD = bcf_str_missing;
        m_matchTypeBK = bcf_str_missing;
        m_nHaplotypeCount = 2;
        m_aGenotype[0] = -1;
        m_aGenotype[1] = -1;
    }
    
    ///FORMAT :  Decison of variant (TP/FP/FN/N)
    std::string m_decisionBD;
    ///FORMAT : Match type of variant (gt/allele match)
    std::string m_matchTypeBK;
    ///Haplotype count of the variant
    int m_nHaplotypeCount;
    ///Genotype of the variant
    int m_aGenotype[2];
    ///Is the genotype is phased
    bool m_bIsPhased = false;
    ///Is the genotype Nocall(./.)
    bool m_bIsNoCallVariant = false;
};

/**
 * @brief Container to store necessary information of vcf record for CVcfWriter
 *
 * SVcfRecord stores mandatory information of vcf record + some info columns to output results
 */
struct SVcfRecord
{
    ///Boundaries of Vcf record for overlapping variant calculation
    int left,right;
    ///Position of the variant (0 based)
    int m_nPosition;
    ///Quality of the variant
    float m_fQuality = 0.0f;
    ///Mendelian Decision INFO (Used for mendelian comparison feature)
    std::string m_mendelianDecision = "";
    ///Filter string of the variant (eg. "PASS")
    std::vector<std::string> m_aFilterString;
    ///Alleles string separated by comma of the variant (eg. m_alleles = "AT,G")
    std::string m_alleles;
    ///Chromosome Name of the variant
    std::string m_chrName;
    ///Sample Data (Data to store for each sample)
    std::vector<SPerSampleData> m_aSampleData;
    ///All info data that is stored for VCF Record (here is the priority from which sample the info columns are taken(for 3 single-sample vcf input): child > father > mother)
    SInfo* m_pInfo = NULL;
};


#endif /* _S_VCF_RECORD_H_ */
