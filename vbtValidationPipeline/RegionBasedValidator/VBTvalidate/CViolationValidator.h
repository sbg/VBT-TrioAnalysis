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


#ifndef CViolationValidator_h
#define CViolationValidator_h

#include "CIntervalValidator.h"
#include "CVariantProvider.h"
#include "CFastaParser.h"
#include <unordered_map>

namespace vbtvalidator
{
    
class CViolationValidator
{
    
public:
    
    int Run(int argc, const char* argv[]);
    
private:

    //Do the main process
    void ValidateRegionsVariantBased();
    
    void ValidateRegionsRegionBased();
    
    void GenerateViolationIntervals(const std::string& a_rContigName, std::vector<SInterval>& ViolationIntervals);
    
    //Takes vcf names as input and outputs total region counts per length
    void CalculatePerRegionVariantCounts(std::string a_rTrioPath);
    
    //Initialize input files of VBT and other tool to be compared
    void SetInputFiles(const std::string& a_rTrioPathVBT,
                       const std::string& a_rTrioPathOriginal,
                       const std::string& a_rFastaPath,
                       const std::string& a_rPedigreeFile,
                       bool a_bIsBeginTrimmingFirst);
    
    //Read the interval files and merge them
    void SetIntervalFiles(const std::string& a_rMotherChildIntervalsPath, const std::string& a_rFatherChildIntervalsPath);
    void MergeIntervals(const std::string& contig);
    
    void ReadIntervals(const std::string& a_rMotherChildIntervalsPath, const std::string& a_rFatherChildIntervalsPath);
    
    //Merged Interval List Map
    std::unordered_map<std::string, std::vector<SInterval>> m_intervalsListMap;
    
    std::unordered_map<std::string, std::vector<int>> m_motherChildSyncPoints;
    std::unordered_map<std::string, std::vector<int>> m_fatherChildSyncPoints;
    
    
    //Classes for VBT
    CVariantProvider m_vbtVariantProvider;
    CIntervalValidator m_vbtIntervalValidator;
    std::string m_vbtTrioPath;
    bool bIsTrimBeginningFirst;
    
    //Classes for LBL
    CVariantProvider m_lblVariantProvider;
    CIntervalValidator m_lblIntervalValidator;
    std::string m_lblTrioPath;
    
    //Fasta Reader
    CFastaParser m_fastaReader;
    SContig m_contig;    
};
    
};


#endif /* CViolationValidator_h */

