//
//  EMendelianDecision.h
//  VCFComparison
//
//  Created by Berke.Toptas on 3/9/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _E_MENDELIAN_DECISION_H_
#define _E_MENDELIAN_DECISION_H_

enum EMendelianDecision
{
    eUnknown = 0,
    eCompliant = 1,
    eViolation = 2,
    eNoCallParent = 3,
    eNoCallChild = 4,
    eSkipped = 5
};

#endif // _E_MENDELIAN_DECISION_H_
