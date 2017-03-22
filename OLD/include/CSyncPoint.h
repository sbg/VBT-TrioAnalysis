//
//  CSyncPoint.h
//  VCFComparison
//
//  Created by Berke.Toptas on 2/22/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_SYNC_POINT_H_
#define _C_SYNC_POINT_H_

#include "CVariant.h"


class CSyncPoint
{
    
public:
    
    //Index of the sync point
    int m_nIndex;
    
    //Range of the regeion
    int m_nStartPosition;
    int m_nEndPosition;
    
    //Included Base variants at the range
    std::vector<const COrientedVariant*> m_baseVariantsIncluded;
    //Excluded Base variants at the range
    std::vector<const CVariant*> m_baseVariantsExcluded;

    //Included Called variants at the range
    std::vector<const COrientedVariant*> m_calledVariantsIncluded;
    //Excluded Called variants at the range
    std::vector<const CVariant*> m_calledVariantsExcluded;
    
};




#endif // _C_SYNC_POINT_H_
