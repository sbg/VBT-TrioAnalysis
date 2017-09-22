//
//  SGraphVarContainer.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/15/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef SGraphVarContainer_h
#define SGraphVarContainer_h

#include <string>
#include <vector>

namespace graphcomparison
{
    /**
     * @brief stores the indexes of variants for each small cluster in graph vcf
     */
    struct SGraphVarContainer
    {
        ///Chromosome name
        std::string chrName;
        ///Included variant indexes for baseline graph
        std::vector<int> baseIncludedIndexes;
        ///Included variant indexes for called graph
        std::vector<int> calledIncludedIndexes;
        ///Excluded variant indexes for baseline graph
        std::vector<int> baseExcludedIndexes;
        ///Excluded variant indexes for called graph
        std::vector<int> calledExcludedIndexes;
    };
}


#endif /* SGraphVarContainer_h */
