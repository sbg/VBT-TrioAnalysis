//
//  SMendelianThreadParam.h
//  VCFComparison
//
//  Created by Berke.Toptas on 2/2/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _S_MENDELIAN_THREAD_PARAM_H_
#define _S_MENDELIAN_THREAD_PARAM_H_

struct SMendelianThreadParam
{
    SMendelianThreadParam(int a_nChrId, bool a_bIsFatherChild)
    {
        m_nChromosomeIndex = a_nChrId;
        m_bIsFatherChildComparison = a_bIsFatherChild;
    }
    
    //Chromosome index to process
    int m_nChromosomeIndex;
    //True if process father-child comparison False if process mother-child comparison
    bool m_bIsFatherChildComparison;
    
};

#endif //_S_MENDELIAN_THREAD_PARAM_H_
