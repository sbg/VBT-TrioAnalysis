//
//  SChrIdTriplet.h
//  VCFComparison
//
//  Created by Berke.Toptas on 5/1/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

// This triplet stores the indexes of a common chromosome for trio. Eg. chromosome MT indexes can be 19, 7, 25 for mother father and child vcfs.


#ifndef _S_CHROMOSOME_ID_TRIPLET_H_
#define _S_CHROMOSOME_ID_TRIPLET_H_

namespace mendelian
{

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
    
    //We will use this index to store best paths [After we eliminate uncommon chromosomes, indexes will be shifted so we will allocate the size of common chromosomes for best paths]
    int m_nTripleIndex;
    
    //Mother index
    int m_nMid;
    //Father index
    int m_nFid;
    //Child index
    int m_nCid;
    //Chromosome Name
    std::string m_chrName;
};

}
    
#endif /* _S_CHROMOSOME_ID_TRIPLET_H_ */
