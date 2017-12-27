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
 *  EMendelianVcfName.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 3/9/17.
 *
 */

#ifndef _E_MENDELIAN_VCF_NAME_H_
#define _E_MENDELIAN_VCF_NAME_H_

namespace mendelian
{

/**
 * @brief ENUM for indicating the vcf file that is being processed
 *
 */
enum EMendelianVcfName
{
    eFATHER = 0,
    eMOTHER = 1,
    eCHILD = 2
};

}


#endif // _E_MENDELIAN_VCF_NAME_H_
