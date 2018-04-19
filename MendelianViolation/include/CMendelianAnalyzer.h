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
 *  CMendelianAnalyzer.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 1/31/17.
 *
 */

#ifndef C_MENDELIAN_ANALYZER_H_
#define C_MENDELIAN_ANALYZER_H_


#include "SConfig.h"
#include "CMendelianVariantProvider.h"
#include "CMendelianResultLog.h"
#include "CMendelianTrioMerger.h"
#include "SChrIdTriplet.h"
#include "CMendelianDecider.h"
#include "ENoCallMode.h"
#include "CViolationRegionOutputGenerator.h"
#include <thread>
#include <mutex>

namespace mendelian
{

/**
 * @brief Mendelian Violation Detection Tool Capability Class
 *
 * CMendelianAnalyzer object is a tool itself that performs mendelian violation detection operation.
 *
 */
class CMendelianAnalyzer
{
    
public:
    
    ///Default Constructor
    CMendelianAnalyzer();
    
    /**
     * @brief Executes the mendelian violation analyzer
     * @param argc Total argument count entered via command line during execution
     * @param argv List of arguments entered via command line during execution
     */
    int run(int argc, char** argv);
    
private:
    
    ///Read the parameters if the execution mode is mendelian. If all mandatory parameters are set, return true.
    bool ReadParameters(int argc, char** argv);
    
    ///A Function to process mendelian violation pipeline for given chromosome id
    void ProcessChromosome(const std::vector<SChrIdTriplet>& a_rChromosomeIds);
        
    ///Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
    int AssignJobsToThreads(int a_nThreadCount);
    
    ///Prints the help menu at console
    void PrintHelp() const;
    
    ///Config file for father-child comparison which is passed to m_fatherChildProvider
    SConfig m_fatherChildConfig;
    
    ///Config file for mother-child comparison which is passed to m_motherChildProvider
    SConfig m_motherChildConfig;
    
    ///No call mode of the variant
    ENoCallMode m_noCallMode;
    
    ///Variant provider for parent-child comparison
    CMendelianVariantProvider m_provider;
    
    ///Mendelian decider to assess mendelian decision of variants
    CMendelianDecider m_mendelianDecider;

    ///Trio merger for outputing
    CMendelianTrioMerger m_trioWriter;
    
    ///Result log for mendelian comparison
    CMendelianResultLog m_resultLog;
    
    ///Generates BED file for violation regions together with related statistics
    CViolationRegionOutputGenerator m_violationRegionGenerator;
    
    //Best Paths written by each thread for each unique chromosome exists [Between father and child]
    std::vector<core::CPath> m_aBestPathsFatherChildGT;
    std::vector<core::CPath> m_aBestPathsFatherChildAM;
    
    //Best Paths written by each thread for each unique chromosome exists [Between mother and child]
    std::vector<core::CPath> m_aBestPathsMotherChildGT;
    std::vector<core::CPath> m_aBestPathsMotherChildAM;
    
    //To prevent data race in multi-thread mode
    std::mutex mtx;

};

}

#endif //C_MENDELIAN_ANALYZER_H_
