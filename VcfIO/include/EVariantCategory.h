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
 *  EVariantCategory.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 3/28/17.
 *
 */

#ifndef _E_VARIANT_TYPE_H_
#define _E_VARIANT_TYPE_H_

/**
 * @brief ENUM to present variant category based on variant size and type(SNP/INSERT/DELETE)
 *
 */
enum class EVariantCategory
{
    eSNP = 0,
    eINDEL_INSERT_SMALL = 1,   //[1, 5]
    eINDEL_INSERT_MEDIUM = 2,  //[6, 14]
    eINDEL_INSERT_LARGE = 3,   //[15, N]
    eINDEL_DELETE_SMALL = 4,   //[1, 5]
    eINDEL_DELETE_MEDIUM = 5,  //[6, 14]
    eINDEL_DELETE_LARGE = 6,   //[15, N]
    eINDEL_COMPLEX_SMALL = 7,  //[1, 5]
    eINDEL_COMPLEX_MEDIUM = 8, //[6, 14]
    eINDEL_COMPLEX_LARGE = 9,   //[15, N]
    eNONE = 10
};

#endif //_E_VARIANT_TYPE_H_
