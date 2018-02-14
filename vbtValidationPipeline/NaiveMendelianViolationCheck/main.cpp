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
 *  Created by Berke Cagkan Toptas
 *
 */

#include <iostream>
#include "CLBLComparisonTool.h"

int main(int argc, const char * argv[])
{
    std::cerr << "NAIVE TRIO COMPARISON ANNOTATOR" << std::endl;

    if(argc == 1)
    {
        std::cerr << "Parameters are missing!" << std::endl;
        return -1;
    }
    
    CLBLComparisonTool LBLtool;
    
    std::string originalTrio = std::string(argv[1]);
    std::string pedigreeFile = std::string(argv[2]);
    std::string outputAnnotatedTrio = std::string(argv[3]);
    
    //Generating Naive decision (LBL) annotated vcf file
    std::cerr << "Generating Naive decision (LBL) annotated vcf file" << std::endl;
    LBLtool.SetInputTrio(originalTrio, pedigreeFile);
    LBLtool.WriteOutputTrio(outputAnnotatedTrio);
}
