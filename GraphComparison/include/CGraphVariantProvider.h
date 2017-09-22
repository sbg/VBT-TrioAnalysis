//
//  CGraphVariantProvider.h
//  VariantBenchmarkingTools
//
//  Created by Berke.Toptas on 8/11/17.
//  Copyright © 2017 Seven Bridges Genomics. All rights reserved.
//

#ifndef _C_GRAPH_VARIANT_PROVIDER_H_
#define _C_GRAPH_VARIANT_PROVIDER_H_

#include "CVariant.h"
#include "COrientedVariant.h"
#include "CGraphVcfReader.h"
#include "CFastaParser.h"
#include "EVcfName.h"
#include "SChrIdTuple.h"
#include <deque>

namespace graphcomparison
{
    /**
     * @brief Reads and (partially) stores variants of base and called graphs
     *
     * CGraphVariantProvider contains functions to parse vcf files of family members. 1 chromosome at a time is stored
     * by provider and variants are accessible via this class. After comparison variants are cleared and next chromosome is read
     */
    class CGraphVariantProvider
    {
        
    public:
        
        ///Set required parameters that CGraphVariantProvider needs
        void SetParameters(const std::string& a_rBaseVcf,
                           const std::string& a_rCalledVcf,
                           const std::string& a_rReference,
                           const std::string& a_rBedFilePath,
                           bool a_bIsPassFilterEnabled,
                           int a_nMaxBasePairLength);
        
        ///Initialize the VCF readers for base and called vcf file
        bool InitializeReaders();
        
        ///Return the variant list with the given index list
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes);
        
        ///Return the all the variants belongs to given chromosome
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo);
                
        ///Return all the oriented variants belongs to given chromosome
        std::vector<const core::COrientedVariant*> GetOrientedVariantList(std::deque<core::COrientedVariant>& a_rOvarList);

        ///Return the oriented variants belongs to given chromosome and given variant list
        std::vector<const core::COrientedVariant*> GetOrientedVariantList(std::deque<core::COrientedVariant>& a_rOvarList, const std::vector<const CVariant*>& a_rVariants);
        
        ///Return the index tuples of chromosomes which both contained by baseline and called VCF
        std::vector<duocomparison::SChrIdTuple>& GetChromosomeIdTuples();
        
        ///Get the contig from fasta file with given chromosome name
        void GetContig(std::string a_chrName, SContig& a_rCtg);
        
        ///Since CPath excluded index list depends on the input variant list, we return a list to the absolute variant index list using variant ids
        static std::vector<int> GetExcludedIndexes(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_rExcludedIndexList);
        
        ///Return the list of variants that are included/excluded according to the variant status stored in the provider
        std::vector<int> GetVariantIndexesByStatus(EVcfName a_uFrom, int a_nChrNo, EVariantMatch a_nStatus);
        
        ///Set the status of each variant in the given list
        void SetVariantStatus(const std::vector<const CVariant*>& a_rVariantList, EVariantMatch a_status) const;
        ///Set the status of each variant in the given list
        void SetVariantStatus(const std::vector<const core::COrientedVariant*>& a_rVariantList, EVariantMatch a_status) const;
        
        ///Read through the variant lists and generate oriented variant list for call and base
        void FillOrientedVariantList(const duocomparison::SChrIdTuple& a_rTuple, std::deque<core::COrientedVariant>& a_rBaseOrientedVars, std::deque<core::COrientedVariant>& a_rCalledOrientedVars);
        
    private:
        
        //Read through the variant file and fill the variant lists. It assumes that positions are sorted.
        void FillVariantLists();
        //Finds the tuple index list of chromosome which is contained by both baseline and called vcf
        void SetChromosomeIdTuples();
        
        //VCF Readers
        CGraphVcfReader m_baseVCF;
        std::string m_baseVcfPath;
        CGraphVcfReader m_calledVCF;
        std::string m_calledVcfPath;
        
        //BED File (Optional)
        std::string m_bedFilePath;
        bool m_bIsBEDFileEnabled;
        
        //Reference to the fasta reader object
        CFastaParser m_fastaParser;
        std::string m_fastaPath;
        

        bool m_bPassFilterEnabled;
        int m_nMaxBasePairLength;
        
        //List that stores base Variants in order
        std::deque<std::deque<CVariant>> m_aBaseVariantList;
        //List that stores called Variants in order
        std::deque<std::deque<CVariant>> m_aCalledVariantList;
                
        //Chromosome id tuples for each common chromosome
        std::vector<duocomparison::SChrIdTuple> m_aCommonChrTupleList;
        
    };
    
}




#endif /* CGraphVariantProvider_h */
