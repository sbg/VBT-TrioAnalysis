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

#ifndef CVariantProvider_h
#define CVariantProvider_h

#include "CVariant.h"
#include "CFastaParser.h"
#include "CVcfReader.h"
#include "CSimplePEDParser.h"


//Forward decleration
namespace vbtvalidator
{
struct SInterval;
}


enum EVcfName
{
    eMOTHER,
    eFATHER,
    eCHILD
};

class CVariantProvider
{
    
public:
    
    //Set the trio vcf path to parse
    void SetTrioPath(const std::string& a_rTrioPath,
                     const std::string& a_rPedigreePath,
                     bool a_bIsTrimBeginningFirst);
    
    //Reads next chromosome from given trio
    bool ReadNextChromosome(std::string& a_rChrName);
    
    //Get the variants according to given interval
    void GetVariants(const vbtvalidator::SInterval& a_rInterval,
                std::vector<const CVariant*>& a_rMotherVars,
                std::vector<const CVariant*>& a_rFatherVars,
                std::vector<const CVariant*>& a_rChildVars,
                bool& a_rIsViolationOverlapConsistent,
                int& a_rTotalViolationCount);

    //Returns all the intervals to a_rViolationInterval that at least 1 violation variant exist
    void GetViolationIntervals(const std::vector<vbtvalidator::SInterval>& a_rIntervalList, std::vector<vbtvalidator::SInterval>& a_rViolationInterval);
    
private:
    
    //Reads next chromosome variants from given vcf reader
    bool ReadNextChromosomeForSample(std::string& a_rChrName, EVcfName a_uFrom);
    
    //Merge trimmed variants with the original variant list
    void AppendTrimmedVariants(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom);

    //Find the optimal trimmings for given variant list if there are more than one trimming option available
    void FindOptimalTrimmings(std::vector<CVariant>& a_rVariantList, EVcfName a_uFrom);
    
    //Checks if the given variant is an SV
    bool IsStructuralVariant(const CVariant& a_rVariant, int a_nMaxLength) const;
    
    //Check if the given variant's genotype is homozygous reference
    bool IsHomRef(const CVariant& a_rVariant) const;
    
    //Compare 2 variant for sorting (based on start, end posititions and variant id)
    static bool CompareVariants(const CVariant& var1, const CVariant& var2);
    
    //Get the mendelian decision (Reads MD annotation from input vcf)
    EMendelianDecision GetMendelianDecision(bcf1_t* a_pRecordPtr, bcf_hdr_t* a_pHeaderPtr);
    
    //Variant List to Process [we will process 1 chromosome at a time]
    std::vector<CVariant> m_aChildVariants;
    std::vector<CVariant> m_aMotherVariants;
    std::vector<CVariant> m_aFatherVariants;

    //Vcf readers for mother father and child
    CVcfReader m_vcfReaderChild;
    CVcfReader m_vcfReaderMother;
    CVcfReader m_vcfReaderFather;
    
    //Pedigree Parser
    CSimplePEDParser m_pedParser;
    
    bool m_bIsTrimBeginningFirst;
    
    //Path of input trio VCF
    std::string m_trioPath;
    
    //Internal iterators for GetVariants Function
    int m_nMotherItr;
    int m_nFatherItr;
    int m_nChildItr;
    
};

#endif /* CVariantProvider_h */
