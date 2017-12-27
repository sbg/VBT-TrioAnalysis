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
 *  EVariantMatch.h
 *  VariantBenchmarkingTools
 *
 *  Created by Berke Cagkan Toptas on 12/29/16.
 *
 */

#ifndef _E_VARIANT_MATCH_H_
#define _E_VARIANT_MATCH_H_

/**
 * @brief ENUM of Variant Match status after a variant is processed by Core Comparison module
 *
 */
enum EVariantMatch
{
    eGENOTYPE_MATCH,
    eALLELE_MATCH,
    eNO_MATCH,
    eNOT_ASSESSED,
    eCOMPLEX_SKIPPED
};

#endif // _E_VARIANT_MATCH_H_
