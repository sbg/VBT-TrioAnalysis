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
#include "CViolationValidator.h"

int main(int argc, const char * argv[])
{
    std::cout << "VBT MENDELIAN DECISION EVALUATOR v2.1" << std::endl;
    std::cout << "Seven Bridges, 2017" << std::endl;
    
    if(argc == 1)
    {
        std::cerr << "Parameters are missing" << std::endl;
        return -1;
    }
    
    std::cerr << "VBT Validator (region based)" << std::endl;
    vbtvalidator::CViolationValidator validator;
    validator.Run(argc, argv);
        
    std::cerr << "[stderr] VBTValidator Exited" << std::endl;
    return 0;
    
}
