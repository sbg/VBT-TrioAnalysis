/*
 *
 * Copyright 2017 Seven Bridges Genomics Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/Path.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 * Licensed under the Simplified BSD License: https://github.com/RealTimeGenomics/rtg-tools/blob/master/LICENSE.txt
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 *
 */

#ifndef _C_SYNC_POINT_H_
#define _C_SYNC_POINT_H_

#include "CVariant.h"
#include "COrientedVariant.h"

namespace core
{

/**
 * @brief Locations that separate variants into small clusters where within each cluster the comparison is independent
 *
 */
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
    
    ///Index of the sync point
    int m_nIndex;
    
    ///0-based start position of the region (inclusive)
    int m_nStartPosition;
    ///0-base end position of the region (exclusive)
    int m_nEndPosition;
    
    ///Included Base variants at the range
    std::vector<const COrientedVariant*> m_baseVariantsIncluded;
    ///Excluded Base variants at the range
    std::vector<const CVariant*> m_baseVariantsExcluded;

    ///Included Called variants at the range
    std::vector<const COrientedVariant*> m_calledVariantsIncluded;
    ///Excluded Called variants at the range
    std::vector<const CVariant*> m_calledVariantsExcluded;
    
};

}


#endif // _C_SYNC_POINT_H_
