//
//  main.cpp
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#include <iostream>
#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include "Constants.h"

int main (int argc, char** argv)
{
    int successNo = 0;
    
    std::cout << "==== VARIANT BENCHMARKING TOOLS " << VBT_VERSION << " ==== " << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgdinc.com)" << std::endl;
    std::cout << "Please send me an email if there is a problem (e.g. bug reporting, run-time errors)" << std::endl;
    std::cout << "COPYRIGHT (C) 2016 Seven Bridges Genomics" << std::endl;
    std::cout << "          (C) 2017 SBGD Inc" << std::endl;
    std::cout << "All rights reserved." << std::endl;
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


