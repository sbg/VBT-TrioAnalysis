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
 *  ENoCallMode.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 3/27/17.
 *
 */

#ifndef _E_NO_CALL_MODE_H_
#define _E_NO_CALL_MODE_H_


namespace mendelian
{

/**
 * @brief ENUM for handling no call variants
 *
 *Explicit no call variant : variants exist in vcf file with ./. genotype
 *Implicit no call variant : sites that no variant exist at vcf file (Hidden 0/0)
 *
 *None: Behaves all implicit and explicit no call variants as HOMREF (0/0) and mark mendelian decision according to that.
 *ImplicitNoCall: Marks both implicit and explicit no call variants' mendelian decision as NoCall.
 *ExplicitNoCall: Marks only explicit no call variants' mendelian decision as NoCall.
 *
 */
enum ENoCallMode
{
    eNone,
    eImplicitNoCall,
    eExplicitNoCall
};
    
}


#endif // _E_NO_CALL_MODE_H_
