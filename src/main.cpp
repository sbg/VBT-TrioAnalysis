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
        std::cout << "More argument is expected. Include a vcf and fasta file." << std::endl;
        return 1;       
    }
    
    std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    SConfig configurations;
    configurations.m_pBaseVcfFileName = argv[1];
    configurations.m_pCalledVcfFileName = argv[2];
    configurations.m_pFastaFileName = argv[3];
    configurations.m_nMaxVariantSize = 1000;
    configurations.m_bIsFilterEnabled = false;
    configurations.m_pFilterName = "PASS";

    
    std::clock_t start;
    double duration;
    
    CPathReplay pathReplay;
    pathReplay.InitializeReaders(configurations);

    start = std::clock();
    CPath bestPath = pathReplay.FindBestPath(21);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout << "Called Included Count: " << bestPath.m_calledSemiPath.GetIncluded().size() << std::endl;
    std::cout << "Called Excluded Count: " << bestPath.m_calledSemiPath.GetExcluded().size() << std::endl;
    std::cout << "Baseline Included Count: " << bestPath.m_baseSemiPath.GetIncluded().size() << std::endl;
    std::cout << "Baseline Excluded Count: " << bestPath.m_baseSemiPath.GetExcluded().size() << std::endl;
    
    std::cout << "Program Completed in " << duration << " secs" << std::endl;

}


//#include <iostream>
//#include <stdlib.h>
//#include <iostream>
//#include "SConfig.h"
//#include "CVcfReader.h"
//#include <fstream>
//
//int main() {
//    SConfig configurations;
//    
//    configurations.m_pBaseVcfFileName = "/Users/c1ms21p6h3qk/Desktop/SBGProjectClean/VCFcomparison/Inputs/sample1.vcf.gz";
//    configurations.m_nMaxVariantSize = 1000;
//    configurations.m_bIsFilterEnabled = false;
//    configurations.m_pFilterName = "PASS";
//    
//    CVcfReader m_baseVCF;
//    bool bIsSuccess = m_baseVCF.Open(configurations.m_pBaseVcfFileName);
//    if(!bIsSuccess)
//        std::cout << "Baseline VCF file is unable to open!: " << configurations.m_pBaseVcfFileName << std::endl;
//
//    const int patientCount = m_baseVCF.GetNumberOfSamples();
//    CVariant* variantList = new CVariant[patientCount];
//    std::vector<std::string> patientNameList;
//    
//    m_baseVCF.GetNextRecordMultiSample(variantList, patientNameList);
//    
//    for(int k= 0; k < patientCount; k++)
//    {
//        std::cout << patientNameList[k] << " " <<variantList[k].ToString() << std::endl;
//    }
//    
//    return 0;
//    
//}
