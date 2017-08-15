#include <iostream>
#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include "CGraphVcfAnalyzer.h"

int main (int argc, char** argv)
{
    int successNo = 0;
    
    std::cout << "Docker image id of this copy is :  vbtapp:v0.6.0" << std::endl;
    std::cerr << "Docker image id of this copy is :  vbtapp:v0.6.0" << std::endl;
    
    std::cout << "==== VARIANT BENCHMARKING TOOL VERSION 1.0 (Beta) ==== " << std::endl;
    std::cout << "Based on paper: http://biorxiv.org/content/biorxiv/early/2015/08/02/023754.full.pdf" << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgdinc.com)" << std::endl;
    std::cout << "Please notify me if program fails or return unexpected results" << std::endl;
    std::cout << "COPYRIGHT (C) 2017 SBGD INC" << std::endl;
    std::cout << std::endl;
    
    if(argc == 1)
    {
        std::cerr << "You have not entered an input.Please try either following:" << std::endl;
        std::cerr << "./vbt varcomp [PARAMETERS]" << std::endl;
        std::cerr << "./vbt mendelian [PARAMETERS]" << std::endl;
        std::cerr << "Please type ./vbt <select_feature> --help for short info about the parameter structure." << std::endl;
        
        return -1;
    }
    
    else if(strcmp(argv[1], "mendelian") == 0)
    {
        mendelian::CMendelianAnalyzer analyzer;
        successNo = analyzer.run(argc, argv);
    }
    
    else if(strcmp(argv[1], "varcomp") == 0)
    {
       duocomparison::CVcfAnalyzer analyzer;
       analyzer.Run(argc, argv);
    }

    else if(strcmp(argv[1], "graphcomp") == 0)
    {
        graphcomparison::CGraphVcfAnalyzer analyzer;
        successNo = analyzer.run(argc, argv);
    }
    
    else
    {
        std::cerr << "Invalid feature name.Please try either following:" << std::endl;
        std::cerr << "./vbt varcomp [PARAMETERS]" << std::endl;
        std::cerr << "./vbt mendelian [PARAMETERS]" << std::endl;
        std::cerr << "Please type ./vbt <select_feature> --help for short info about the parameter structure." << std::endl;
        
        return -1;
    }
    
    return successNo;
}


