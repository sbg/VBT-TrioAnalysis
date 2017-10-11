/*
 * Originated from: https://github.com/RealTimeGenomics/rtg-tools/blob/master/src/com/rtg/vcf/eval/Path.java
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 * Ported to C++ by Berke Cagkan Toptas
 * Copyright (c) 2016 Seven Bridges Genomics
 *               2017 SBGD Inc.
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
