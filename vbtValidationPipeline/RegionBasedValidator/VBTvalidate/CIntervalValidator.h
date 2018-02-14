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

#ifndef CIntervalValidator_h
#define CIntervalValidator_h

#include <deque>
#include <string>
#include <vector>

//forward decleration
class CVariant;
class CVariantProvider;
class CFastaParser;
struct SContig;

namespace vbtvalidator
{

struct SInterval
{
    int m_nStart;
    int m_nEnd;
};

struct SSequence
{
    std::string m_haplotypeA;
    std::string m_haplotypeB;
};

struct SPhase
{
    std::vector<unsigned short> m_phaseVector;
};

enum EIntervalDecision
{
    eInterval_WrongAssessed,
    eInterval_CorrectlyAssessed,
    eInterval_Complex
};

class CIntervalValidator
{
    
public:
    //For given interval test whether Mendelian decisions are correct or wrong
    EIntervalDecision IntervalTestMV(const std::string& a_rChrName, const SInterval& a_rInterval, int& a_rViolationCount) const;
    
    //Get access to the variant provider
    void SetVariantProvider(CVariantProvider* a_pProvider);

    //Get access to the FASTA reader
    void SetFastaReader(CFastaParser& a_rFastaReader, SContig* a_pContig);
    
    //Get access to the current FASTA contig
    void SetContig(SContig& a_rContig);
    
private:

    //Generate phases of given set of variants
    void GeneratePhases(std::deque<SPhase>& a_rPhaseVector, std::vector<const CVariant*>& a_rVariants) const;

    //Apply variant set with given phasing vector to the clipped reference sequence
    SSequence ApplyPhasing(const SPhase& a_rPhase, const SInterval& a_rInterval, const std::string& a_rReferenceSeq, std::vector<const CVariant*>& a_rVariantList) const;
    
    //Access to variant provider
    CVariantProvider* m_pVariantProvider;
    
    //Access to fasta reader
    CFastaParser* m_pFastaReader;
    SContig* m_pContig;
    
};

};

#endif /* CIntervalValidator_h */

