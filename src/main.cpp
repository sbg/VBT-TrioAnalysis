#include <stdlib.h>
#include <iostream>

#include "CPathReplay.h"
#include "SConfig.h"
#include <ctime>
#include <fstream>

#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include <string>

//void testGa4ghOutput();
//void testRefOverlapOutput(int argc, char** argv);
//void compareTrimmedVariables();


int main (int argc, char** argv)
{
    std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    std::ofstream err("/Users/c1ms21p6h3qk/Desktop/errXCOCD.txt");
    std::cerr.rdbuf(err.rdbuf()); //redirect std::cerr to error.txt!
    
    
    if(strcmp(argv[1], "-mendelian") == 0)
    {
       CMendelianAnalyzer mendelianAnalyzer;
       mendelianAnalyzer.run(argc, argv);
    }
    
    else
    {
       CVcfAnalyzer analyzer;
       analyzer.Run(argc, argv);
    }
    
    //=====TESTS=====
    /*
    testGa4ghOutput();
    testRefOverlapOutput(argc, argv);
    compareTrimmedVariables();
    */
}


