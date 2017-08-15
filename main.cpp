#include <iostream>
#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include "CGraphVcfAnalyzer.h"

int main (int argc, char** argv)
{
    int successNo = 0;
    
    std::cout << "Docker image id of this copy is :  vbtapp:v0.6.0" << std::endl;
    std::cerr << "Docker image id of this copy is :  vbtapp:v0.6.0" << std::endl;
    
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


