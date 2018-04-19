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
 *  CViolationRegionOutputGenerator.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke.Toptas on 4/18/18.
 *
 */

#ifndef CViolationRegionOutputGenerator_h
#define CViolationRegionOutputGenerator_h

#include "CPath.h"
#include <fstream>
#include "EMendelianDecision.h"
#include "SVcfRecord.h"

namespace mendelian
{
    class CMendelianVariantProvider;
    class CMendelianResultLog;
    struct SChrIdTriplet;
  
class CViolationRegionOutputGenerator
{
    
public:
    
    ///Default constructor
    CViolationRegionOutputGenerator(const std::vector<core::CPath>& a_aBestPathsFatherChildGT,
                                    const std::vector<core::CPath>& a_aBestPathsFatherChildAM,
                                    const std::vector<core::CPath>& a_aBestPathsMotherChildGT,
                                    const std::vector<core::CPath>& a_aBestPathsMotherChildAM,
                                    const CMendelianVariantProvider& a_rProvider,
                                    CMendelianResultLog& a_rResultLog);
    
    ///Open the output BED file stream
    void OpenBedForWrite(const std::string& a_rFilePath);
    
    //Close the BED file stream
    void CloseBed();
    
    //Generate violation regions at the given output directory. Using Resultlog, also generate region based statistics
    void GenerateViolationRegions(SChrIdTriplet& a_triplet, std::vector<SVcfRecord>& a_rRecordList, std::vector<EMendelianDecision>& a_rDecisionList);
    
private:
    
    //Generates the syncpoints for given chromosome and returns the list as second parameter
    void GenerateSyncPoints(const SChrIdTriplet& a_triplet, std::vector<int>& intersectedSyncPoints);
    
    //Best Paths (with Genotype Match) written by each thread for each unique chromosome exists [Between father and child]
    const std::vector<core::CPath>& m_aBestPathsFatherChildGT;
    //Best Paths (with Allele Match) written by each thread for each unique chromosome exists [Between father and child]
    const std::vector<core::CPath>& m_aBestPathsFatherChildAM;
    
    //Best Paths (with Genotype Match) written by each thread for each unique chromosome exists [Between mother and child]
    const std::vector<core::CPath>& m_aBestPathsMotherChildGT;
    //Best Paths (with Allele Match) written by each thread for each unique chromosome exists [Between mother and child]
    const std::vector<core::CPath>& m_aBestPathsMotherChildAM;
    
    ///Variant provider for parent-child comparison
    const CMendelianVariantProvider& m_provider;
    
    ///Result log for mendelian comparison
    CMendelianResultLog& m_resultLog;

    //Output directory that the resulting BED file will be written
    std::ofstream m_outputBEDfile;
    
    //True if violation regions will be printed
    bool m_bIsOutputEnabled;

    
};

}

#endif /* CViolationRegionOutputGenerator_h */
