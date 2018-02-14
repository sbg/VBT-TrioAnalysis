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
 *  Created by Berke Cagkan Toptas
 *
 */

#ifndef CLBLComparisonTool_h
#define CLBLComparisonTool_h

#include "CVcfReader.h"
#include "CVcfWriter.h"

class CLBLComparisonTool
{
    
public:
    
    //Set Input trio path
    void SetInputTrio(const std::string& a_rInputTrioPath, const std::string& a_rPedigreeFile, bool a_bIsTrimBeginningFirst);
    
    //Write annotated trio to the given path
    void WriteOutputTrio(const std::string& a_rOutputTrioPath);
    
private:
    
    //Class to read trio
    CVcfReader m_vcfReader;
    
    //Mother, father and child indexes in vcf (0-based Sample Order in VCF )
    int m_motherIndex, m_fatherIndex, m_childIndex;
    
    //Choice of allele trimming order
    bool m_bIsTrimBeginningFirst;
    
    //Class to output annotated trio
    CVcfWriter m_vcfWriter;
    
    std::string pedigreeFile;

};



#endif /* CLBLComparisonTool_h */
