#include <stdlib.h>
#include <iostream>
#include "CVcfReader.h"
#include "CFastaReader.h"
#include "CPathReplay.h"
#include "SConfig.h"
#include <fstream>

int main (int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "Argument Count is wrong!! include a vcf and fasta file." << std::endl;
        return 1;       
    }
    
    std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    SConfig configurations;
    
    configurations.m_bIsFilterPASS = true;
    configurations.m_pBaseVcfFileName = argv[1];
    configurations.m_pCalledVcfFileName = argv[2];
    configurations.m_pFastaFileName = argv[3];

    CPathReplay pathReplay;
    pathReplay.InitializeReaders(configurations);
    pathReplay.FindBestPath(21);

}
