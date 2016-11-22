//
//  CSyncPoint.cpp
//  VCFComparison
//
//  Created by Berke.Toptas on 11/18/16.
//  Copyright © 2016 Seven Bridges Genomics. All rights reserved.
//


#include "CSyncPoint.h"


CSyncPoint::CSyncPoint()
{
    m_nPosition = -1;
    m_nCalledTPCount = 0.0;
    m_nBaselineTPCount = 0.0;
}

CSyncPoint::CSyncPoint(int a_nPos, int a_nCalledCount, int a_nBaselineCount)
{
    m_nPosition = a_nPos;
    m_nCalledTPCount = a_nCalledCount;
    m_nBaselineTPCount = a_nBaselineCount;
}

int CSyncPoint::GetPosition()
{
    return m_nPosition;
}

int CSyncPoint::CompareTo(const CSyncPoint& a_rObj)
{
    if(m_nPosition > a_rObj.m_nPosition)
        return 1;
    else if (m_nPosition == a_rObj.m_nPosition)
        return 0;
    else
        return -1;
}
