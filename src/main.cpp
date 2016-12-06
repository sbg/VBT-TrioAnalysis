#include <stdlib.h>
#include <iostream>

#include "CPathReplay.h"
#include "SConfig.h"
#include <ctime>
#include <fstream>

int main (int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "Argument Count is wrong!! include a vcf and fasta file." << std::endl;
        return 1;       
    }
    
    std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    SConfig configurations;
    configurations.m_pBaseVcfFileName = argv[2];
    configurations.m_pCalledVcfFileName = argv[2];
    configurations.m_pFastaFileName = argv[3];
    configurations.m_nMaxVariantSize = 1000;
    configurations.m_bIsFilterEnabled = false;
    configurations.m_pFilterName = "PASS";
    
    CPathReplay pathReplay;
    pathReplay.InitializeReaders(configurations);

    std::clock_t start;
    double duration;
    
    start = std::clock();
    CPath bestPath = pathReplay.FindBestPath(21);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout << "Called Included Count: " << bestPath.m_calledSemiPath.GetIncluded().size() << std::endl;
    std::cout << "Called Excluded Count: " << bestPath.m_calledSemiPath.GetExcluded().size() << std::endl;
    std::cout << "Baseline Included Count: " << bestPath.m_baseSemiPath.GetIncluded().size() << std::endl;
    std::cout << "Baseline Excluded Count: " << bestPath.m_baseSemiPath.GetExcluded().size() << std::endl;
    
    
    std::cout << "Program Completed :" << duration << " secs" << std::endl;

}
