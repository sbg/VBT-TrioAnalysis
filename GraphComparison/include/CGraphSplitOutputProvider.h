//
//  CGraphSplitOutputProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/15/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_GRAPH_SPLIT_OUTPUT_PROVIDER_H_
#define _C_GRAPH_SPLIT_OUTPUT_PROVIDER_H_

#include <unordered_map>
#include "SGraphVarContainer.h"
#include "SChrIdTuple.h"


namespace graphcomparison
{
    class CGraphVariantProvider;

    /**
     * @brief Outputs TPBase, TPCalled, FP and FN vcf files as output after comparison operation
     *
     */
    class CGraphSplitOutputProvider
    {
        public:
        
        ///Set access to variant provider
        void SetVcfNames(const std::string& a_rBaseVcfPath, const std::string& a_rCalledVcfPath);
        
        ///Pass graph comparison parameters to output provider since the input vcfs will be traced again
        void SetParameters(const std::string& a_rBedPath, bool a_bPassFilterEnabled, int a_nMaxBasePairLength);
        
        ///Set access to best path list
        void SetExcludeIndexes(std::unordered_map<std::string, SGraphVarContainer>& a_rExcludedIndexContainers);
        
        ///Set the output vcfs path FOLDER
        void SetOutputVcfFolder(const std::string& a_rVcfPath);
        
        ///Generates 2 vcf files splitting included and excluded decisions
        void GenerateVcfs(const std::string& a_rIncludeName, const std::string& a_rExcludeName, bool a_bIsBase);
        
        private:
        
        //Access to best paths
        std::unordered_map<std::string, SGraphVarContainer> m_excludedIndexContainers;
        
        //Graph comparison parameters
        bool m_bIsBEDFileEnabled;
        std::string m_bedFilePath;
        bool m_bPassFilterEnabled;
        int m_nMaxBasePairLength;
        
        //Path of input vcf files
        std::string m_baseVcfPath;
        std::string m_calledVcfPath;
        
        //Path of output folder where we place vcf files
        std::string m_vcfsFolder;
        
    };

}

#endif /* _C_GRAPH_SPLIT_OUTPUT_PROVIDER_H_ */
