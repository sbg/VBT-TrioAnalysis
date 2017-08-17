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
    struct SGraphVarContainer
    {
        std::string chrName;
        std::vector<int> baseIncludedIndexes;
        std::vector<int> calledIncludedIndexes;
        std::vector<int> baseExcludedIndexes;
        std::vector<int> calledExcludedIndexes;
    };
}


#endif /* SGraphVarContainer_h */
