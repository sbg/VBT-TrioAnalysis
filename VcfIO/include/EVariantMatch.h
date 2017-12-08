//
//  EVariantMatch.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 12/29/16.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

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
