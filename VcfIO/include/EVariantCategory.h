//
//  EVariantType.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/28/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _E_VARIANT_TYPE_H_
#define _E_VARIANT_TYPE_H_

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
