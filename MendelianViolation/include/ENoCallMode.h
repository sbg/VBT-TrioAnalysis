//
//  ENoCallMode.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/27/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _E_NO_CALL_MODE_H_
#define _E_NO_CALL_MODE_H_

/* Enumeration for handling no call variants

 Explicit no call variant : variants exist in vcf file with ./. genotype
 Implicit no call variant : sites that no variant exist at vcf file (Hidden 0/0)
 
 None: Behaves all implicit and explicit no call variants as HOMREF (0/0) and mark mendelian decision according to that.
 ImplicitNoCall: Marks both implicit and explicit no call variants' mendelian decision as NoCall.
 ExplicitNoCall: Marks only explicit no call variants' mendelian decision as NoCall.
 
*/

namespace mendelian
{

enum ENoCallMode
{
    eNone,
    eImplicitNoCall,
    eExplicitNoCall
};
    
}


#endif // _E_NO_CALL_MODE_H_
