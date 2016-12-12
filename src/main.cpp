#include <stdlib.h>
#include <iostream>

#include "CPathReplay.h"
#include "SConfig.h"
#include <ctime>
#include <fstream>

#include "CVcfAnalyzer.h"


int main (int argc, char** argv)
{
    std::ofstream out("/Users/c1ms21p6h3qk/Desktop/outXCODE.txt");
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    CVcfAnalyzer analyzer;
    analyzer.Run(argc, argv);
}


