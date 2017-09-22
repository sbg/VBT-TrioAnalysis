//
//  CGraphVcfReader.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/11/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_GRAPH_VCF_READER_H_
#define _C_GRAPH_VCF_READER_H_

#include "CVariant.h"
#include "SConfig.h"
#include "CVcfReader.h"
#include <map>

namespace graphcomparison
{
    /**
     * @brief Parser for graph vcf file
     *
     */
    class CGraphVcfReader
    {
    public:
        ///Constructor
        CGraphVcfReader();
        ///Destructor
        ~CGraphVcfReader();
    
        ///Open a VCF file
        bool Open(const char * a_pFilename);
        
        ///Close a VCF file
        bool Close();
        
        ///Get next record in the file. a_nId sets the id of variant (no need to set)
        bool GetNextRecord(CVariant* a_pVariant, int a_nId = 0, bool a_bIsReadAll = false);

        ///Get contig size
        inline int GetContigSize() const {return (int)(m_contigs.size());};
        
        ///Returns Header and Record Pointer of htslib
        bcf_hdr_t* GetHeaderPointer();
        bcf1_t* GetRecordPointer();
        
        ///Store chromosome name indexes
        std::map<std::string, int> m_chrIndexMap;
        
    private:
        
        // Trimms the alt string that allowing reference overlap
        void TrimRefOverlap(SAllele& a_rAllele);
        
        std::string m_filename;
        bool m_bIsOpen;
        htsFile *   m_pHtsFile;
        bcf_hdr_t * m_pHeader;
        bcf1_t *    m_pRecord;
        
        std::vector<SVcfContig> m_contigs;
    
    };
}



#endif //_C_GRAPH_VCF_READER_H_
