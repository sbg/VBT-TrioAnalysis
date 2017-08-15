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

namespace graphcomparison
{
    class CGraphVariantProvider
    {
        
    public:
                
        void SetParameters(const std::string& a_rBaseVcf,
                           const std::string& a_rCalledVcf,
                           const std::string& a_rReference,
                           const std::string& a_rBedFilePath,
                           bool a_bIsPassFilterEnabled,
                           int a_nMaxBasePairLength);
        
        //Initialize the VCF readers for base and called vcf file
        bool InitializeReaders();
        
        //Return the variant list with the given index list
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int>& a_VariantIndexes);
        
        //Return the all the variants belongs to given chromosome
        std::vector<const CVariant*> GetVariantList(EVcfName a_uFrom, int a_nChrNo);
                
        //Return all the oriented variants belongs to given chromosome
        std::vector<const core::COrientedVariant*> GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo);

        //Return the oriented variants belongs to given chromosome and given variant indexes
        std::vector<const core::COrientedVariant*> GetOrientedVariantList(EVcfName a_uFrom, int a_nChrNo, const std::vector<int> a_rVariantIndexes);

        
        //Return the index tuples of chromosomes which both contained by baseline and called VCF
        std::vector<duocomparison::SChrIdTuple>& GetChromosomeIdTuples();
        
        //Get the contig from fasta file with given chromosome name
        void GetContig(std::string a_chrName, SContig& a_rCtg);
        
        //Since CPath excluded index list depends on the input variant list, we return a list to the absolute variant index list using variant ids
        static std::vector<int> GetExcludedIndexes(const std::vector<const CVariant*> a_rVariantList, const std::vector<int>& a_rExcludedIndexList);
        
    private:
        
        //Read through the variant file and fill the variant lists. It assumes that positions are sorted.
        void FillVariantLists();
        //Read through the variant lists and generate oriented variant list for call and base
        void FillOrientedVariantLists();
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
        std::vector<std::vector<CVariant>> m_aBaseVariantList;
        //List that stores called Variants in order
        std::vector<std::vector<CVariant>> m_aCalledVariantList;
        
        //List that store the base Oriented variant tuples (In the order of genotype)
        std::vector<std::vector<core::COrientedVariant>> m_aBaseOrientedVariantList;
        //List that store the called Oriented variant tuples (In the order of genotype)
        std::vector<std::vector<core::COrientedVariant>> m_aCalledOrientedVariantList;
        
        
        //Chromosome id tuples for each common chromosome
        std::vector<duocomparison::SChrIdTuple> m_aCommonChrTupleList;
        
    };
    
}




#endif /* CGraphVariantProvider_h */
