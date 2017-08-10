#include <iostream>
#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"

int main (int argc, char** argv)
{
    int successNo = 0;
    
    std::cout << "Docker image of this copy is :  vbtapp:v5.40" << std::endl;
    std::cerr << "Docker image of this copy is :  vbtapp:v5.40" << std::endl;
    
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
       CMendelianAnalyzer mendelianAnalyzer;
       successNo = mendelianAnalyzer.run(argc, argv);
    }
    
    else if(strcmp(argv[1], "varcomp") == 0)
    {
       duocomparison::CVcfAnalyzer analyzer;
       analyzer.Run(argc, argv);
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


