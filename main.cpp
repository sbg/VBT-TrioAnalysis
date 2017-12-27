/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  main.cpp
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas
 *
 */

#include <iostream>
#include "CVcfAnalyzer.h"
#include "CMendelianAnalyzer.h"
#include "Constants.h"

int main (int argc, char** argv)
{
    int successNo = 0;
    std::cout << std::endl;
    std::cout << "==== VARIANT BENCHMARKING TOOLS " << VBT_VERSION << " ==== " << std::endl;
    std::cout << "Please contact us if you would like to report bugs or encounter run-time errors" << std::endl;
    std::cout << "Author: Berke Cagkan Toptas (berke.toptas@sbgdinc.com)" << std::endl;
    std::cout << "Copyright 2017 Seven Bridges Genomics Inc." << std::endl;
    std::cout << std::endl;
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
