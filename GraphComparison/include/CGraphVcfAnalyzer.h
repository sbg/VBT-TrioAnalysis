//
//  CGraphVcfAnalyzer.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/10/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_GRAPH_VCF_ANALYZER_H_
#define _C_GRAPH_VCF_ANALYZER_H_

#include "CGraphVariantProvider.h"
#include "SChrIdTuple.h"
#include "SGraphVarContainer.h"
#include "CPath.h"
#include "CSyncPoint.h"
#include <unordered_map>
#include <mutex>

namespace graphcomparison
{
    /**
     * @brief Graph VCF Comparison Tool Capability Class
     *
     * CGraphVcfAnalyzer object is a tool itself that performs graph vcf comparison operation.
     *
     */
    class CGraphVcfAnalyzer
    {
        
    public:
        
        ///Constructor
        CGraphVcfAnalyzer();
        
        ///Run the application
        int run(int argc, char** argv);
        
    private:
        
        //Prints the help menu at console
        void PrintHelp() const;
                
        //Read Parameters from command line. If the mandatory arguments are given, return true.
        bool ReadParameters(int argc, char** argv);
        
        //Divide the jobs between different threads homogeneously for given number of thread count. Return the actual thread count
        int AssignJobsToThreads(int a_nThreadCount);

        //Function that process chromosome in bulk for SPLIT mode (process either genotype or allele match)
        void ThreadFunctionSPLIT(std::vector<duocomparison::SChrIdTuple> a_aTuples);

        //Get the SyncPointList between 2 graph
        void GetSyncPointList(const std::vector<const core::COrientedVariant*> a_rBaseIncluded,
                              const std::vector<const core::COrientedVariant*> a_rCalledIncluded,
                              const std::vector<const CVariant*>& a_rBaseExcluded,
                              const std::vector<const CVariant*>& a_rCalledExcluded,
                              const std::vector<int>& a_rSyncPointCoordinates,
                              std::vector<core::CSyncPoint>& a_rSyncPointList);
        
        //Filter excluded variants for next iteration if there is no included branch within the same synchronization points
        void FilterGuaranteedUniqueVariants(std::vector<core::CSyncPoint>& a_rSyncPointList,
                                            std::vector<const CVariant*>& a_rBaseExcluded,
                                            std::vector<const CVariant*>& a_rCalledExcluded);
        
        //Print the logs of graphcomparison to log.txt
        void PrintLogs();
        
        //Variant provider instance
        CGraphVariantProvider m_provider;
        
        //For each chromosome, stores the unique variant indexes to be printed
        std::unordered_map<std::string, SGraphVarContainer> m_uniqueVariantsListPerChromosome;
        
        //Graph Vcf Analyzer input parameters
        std::string m_baseVcfPath;
        std::string m_calledVcfPath;
        std::string m_referencePath;
        std::string m_outputDirectoryPath;
        std::string m_bedFilePath;
        bool m_bIsPassFilterEnabled;
        int m_nThreadCount;
        int m_nMaxBasePairLength;
        int m_nMaxPathSize;
        int m_nMaxIterationCount;
        
        //lock mechanisim to prevent data racing for FASTA contig reader
        std::mutex mtx;
        
    };
    
    
    
    
    
    

}




#endif // _C_GRAPH_VCF_ANALYZER_H_
