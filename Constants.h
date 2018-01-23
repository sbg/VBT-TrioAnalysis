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
 *  Constants.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas 2/1/17.
 *
 */

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <string>

//MAX NUMBER OF THREAD COUNT
const int MAX_THREAD_COUNT = 25;

//LEAST NUMBER OF VARIANT REQUIRED TO PROCESS THE CHROMOSOME
const int LEAST_VARIANT_THRESHOLD = 2;

//THREAD COUNT BY DEFAULT IN CASE PROGRAM WONT WORK IN PLATFORM MODE
const int DEFAULT_THREAD_COUNT = 2;

//DEFAULT BASE PAIR LENGTH OF A VARIANT THAT COMPARISON ENGINE PROCESSS
const int DEFAULT_MAX_BP_LENGTH = 1000;

//DEFAULT SIZE OF PATH SET (For Dynamic Programming result saving. This variable should be increased carefully since it is easy to exceed available memory)
const int DEFAULT_MAX_PATH_SIZE = 150000;

//DEFAULT SIZE OF ITERATION COUNT (For Dynamic Programming result saving. This variable should be increased carefully since it is easy to exceed available memory)
const int DEFAULT_MAX_ITERATION_SIZE = 10000000;

//DEFAULT SIZE OF SMALL VARIANTS FOR MENDELIAN VIOLATION DETECTION
const int SMALL_VARIANT_SIZE = 5;

//DEFAULT SIZE OF MEDIUM VARIANTS FOR MENDELIAN VIOLATION DETECTION
const int MEDIUM_VARIANT_SIZE = 15;

//VBT VERSION AND YEAR TO BE COPIED OUTPUT VCFS
const std::string VBT_VERSION = "v1.0 (2018)";

#endif //_CONSTANTS_H_
