//
//  CSyncPoint.h
//  VCFComparison
//
//  Created by Berke.Toptas on 11/18/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_SYNC_POINT_H_
#define _C_SYNC_POINT_H_

class CSyncPoint
{
    
    public:
    
    CSyncPoint();

    CSyncPoint(int a_nPosition, int a_nCalledCount, int a_nBaselineCount);

    int GetPosition();
    
    int CompareTo(const CSyncPoint& a_rObj);
    

    int m_nPosition;
    int m_nCalledTPCount;
    int m_nBaselineTPCount;    
};


#endif // _C_SYNC_POINT_H_
