//
//  SChrIdTriplet.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 5/1/17.
//  Copyright © 2016 Seven Bridges Genomics.
//            © 2017 SBGD Inc.
//  All rights reserved.
//

#ifndef _S_CHROMOSOME_ID_TRIPLET_H_
#define _S_CHROMOSOME_ID_TRIPLET_H_

namespace mendelian
{

/**
 * @brief Groups indexes of common chromosomes for mother father and child
 *
 */
struct SChrIdTriplet
{
    SChrIdTriplet(int a_mother, int a_father, int a_child, std::string a_chrName, int a_tripleIndex)
    {
        m_nMid = a_mother;
        m_nFid = a_father;
        m_nCid = a_child;
        m_chrName = a_chrName;
        m_nTripleIndex = a_tripleIndex;
    }
    
    SChrIdTriplet()
    {
        m_nMid = -1;
        m_nFid = -1;
        m_nCid = -1;
        m_chrName = "none";
        m_nTripleIndex = -1;
    }
    
    ///We will use this index to store best paths [After we eliminate uncommon chromosomes, indexes will be shifted so we will allocate the size of common chromosomes for best paths]
    int m_nTripleIndex;
    
    ///Index of mother chromosome in Mother vcf
    int m_nMid;
    ///Index of father chromosome in Father vcf
    int m_nFid;
    ///Index of child chromosome in Child vcf
    int m_nCid;
    ///Chromosome Name
    std::string m_chrName;
};

}
    
#endif /* _S_CHROMOSOME_ID_TRIPLET_H_ */
