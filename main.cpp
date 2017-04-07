#include <stdlib.h>
#include <iostream>

#include "CPathReplay.h"
#include "SConfig.h"
#include <ctime>
#include <fstream>

#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include <string>

#include "COverlappingVariantEliminator.h"

//void testGa4ghOutput();
//void testRefOverlapOutput(int argc, char** argv);
//void compareTrimmedVariables();
void testTripletReader2();
void compare2List(const std::string& s1, const std::string& s2);
void GenerateTruthSetsMaria(std::string a_rFilename, bool a_bIsFilterOverlap, bool a_bIsFilter00);
void UnitTestTrioComparison(int a_nChrNumber, bool a_bIsFilterOverlap, bool a_bIsFilter00);

int main (int argc, char** argv)
{
    //std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    //std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    //std::ofstream err("/Users/c1ms21p6h3qk/Desktop/errXCOCD.txt");
    //std::cerr.rdbuf(err.rdbuf()); //redirect std::cerr to error.txt!
    
    
    //COverlappingVariantEliminator elim;
    //elim.FilterOverlaps("/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered.vcf", false, false, 0);
    //elim.FilterOverlaps("/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered.vcf", false, true, 1);
    //elim.FilterOverlaps("/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered.vcf", false, true, 2);
    //return 0;
    
    //testTripletReader2();
    
    //std::string fileQuery = "/Users/c1ms21p6h3qk/Desktop/aaa.txt";
    //std::string fileQuery = "/Users/c1ms21p6h3qk/Desktop/MendelianOutput/CHR1/chr1_ViolationsALL.txt";
    //std::string fileTruth = "/Users/c1ms21p6h3qk/Desktop/gatk_corrected_positions.txt";
    //compare2List(fileTruth, fileQuery);
    //return 0;
    
    //GenerateTruthSetsMaria("/Users/c1ms21p6h3qk/Desktop/MendelianInput/ALL_filteredPedGraph_D1_D2_D3_9.31.cat_Filtered.vcf", true, false);
    //return 0;
    
    //UnitTestTrioComparison(1, true, false);
    //return 0;
    
    //for(int m = 1; m < 23; m++)
    //    UnitTestTrioComparison(m, true, false);
    //return 0;
    
    //strcpy(argv[3],  "/Users/c1ms21p6h3qk/Desktop/TestVcfBuilder/father.vcf");
    //strcpy(argv[5],  "/Users/c1ms21p6h3qk/Desktop/TestVcfBuilder/child.vcf");
    //strcpy(argv[7],  "/Users/c1ms21p6h3qk/Desktop/TestVcfBuilder/mother.vcf");
    //strcpy(argv[13], "/Users/c1ms21p6h3qk/Desktop/TestVcfBuilder/reference.fasta");
    
    
    if(strcmp(argv[1], "mendelian") == 0)
    {
       CMendelianAnalyzer mendelianAnalyzer;
       mendelianAnalyzer.run(argc, argv);
    }
    
    else if(strcmp(argv[1], "varcomp") == 0)
    {
       CVcfAnalyzer analyzer;
       analyzer.Run(argc, argv);
    }
    
    else
    {
        std::cout << "Invalid feature name.Please try either following:" << std::endl;
        std::cout << "./vbt varcomp [PARAMETERS]" << std::endl;
        std::cout << "./vbt mendelian [PARAMETERS]" << std::endl;
    }
    
}


