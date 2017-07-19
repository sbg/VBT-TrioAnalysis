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
#include "COrientedVariant.h"

class CSyncPoint
{
    
public:
    
    CSyncPoint(const CSyncPoint& rhs)
    {
        m_nIndex = rhs.m_nIndex;
        m_nStartPosition = rhs.m_nStartPosition;
        m_nEndPosition = rhs.m_nEndPosition;
        m_baseVariantsExcluded = std::vector<const CVariant*>(rhs.m_baseVariantsExcluded);
        m_calledVariantsExcluded = std::vector<const CVariant*>(rhs.m_calledVariantsExcluded);
        m_baseVariantsIncluded = std::vector<const COrientedVariant*>(rhs.m_baseVariantsIncluded);
        m_calledVariantsIncluded = std::vector<const COrientedVariant*>(rhs.m_calledVariantsIncluded);
    }
    
    CSyncPoint()
    {
        m_nIndex = -1;
        m_nStartPosition = -1;
        m_nEndPosition = -1;
    }
    
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
