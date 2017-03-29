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
    eSNP,
    eINDEL_INSERT_SMALL,   //[1, 5]
    eINDEL_INSERT_MEDIUM,  //[6, 14]
    eINDEL_INSERT_LARGE,   //[15, N]
    eINDEL_DELETE_SMALL,   //[1, 5]
    eINDEL_DELETE_MEDIUM,  //[6, 14]
    eINDEL_DELETE_LARGE,   //[15, N]
    eINDEL_COMPLEX_SMALL,  //[1, 5]
    eINDEL_COMPLEX_MEDIUM, //[6, 14]
    eINDEL_COMPLEX_LARGE   //[15, N]
};

#endif //_E_VARIANT_TYPE_H_
