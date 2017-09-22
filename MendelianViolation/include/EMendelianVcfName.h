//
//  EMendelianVcfName.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

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
    eFATHER,
    eMOTHER,
    eCHILD
};

}


#endif // _E_MENDELIAN_VCF_NAME_H_
