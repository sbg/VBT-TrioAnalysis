//
//  EMendelianVcfName.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
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
    eFATHER = 0,
    eMOTHER = 1,
    eCHILD = 2
};

}


#endif // _E_MENDELIAN_VCF_NAME_H_
