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

#ifndef VCF_READER_H_
#define VCF_READER_H_
#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "CVariant.h"
#include "SConfig.h"
#include <map>

struct SVcfContig
{
    std::string name;
    int length;
};

class CVcfReader
{

public:
    CVcfReader();
    CVcfReader(const char * a_pFilename);
    ~CVcfReader();
    
    // Open a VCF file
    bool Open(const char * a_pFilename);
    
    // Close a VCF file
    bool Close();
    
    // Get next record in the file. a_nId sets the id of variant (no need to set)
    bool GetNextRecord(CVariant* a_pVariant, int a_nId, const SConfig& a_rConfig);
    
    // Get next record which has multiple samples in the file(Generates list of N variants and name of each sample)
    bool GetNextRecordMultiSample(CVariant* a_pVariant, bool a_bIsTrimmingFirst = true, bool a_bIsRefOverlap = false);
    
    // Selects the sample name from multi sample VCF file and ignore other samples
    bool SelectSample(const std::string& a_sampleName);
    
    //Fills the a_PatientList with sample names
    void GetSampleNames(std::vector<std::string>& a_PatientList);
    
    // Get the filename
    std::string GetFilename() const {return m_filename;};
    
    // Get the number of samples in vcf
    int GetNumberOfSamples() const;
    
    // Get contig size
    inline int GetContigSize() const {return (int)(contigs_.size());};
    
    //Return the list of chromosome that vcf file contains
    const std::vector<SVcfContig>& GetContigs() const;
    
    // Get a contig name
    inline std::string GetContigName(const unsigned int& id) const {if (id < contigs_.size()) return contigs_[id].name; else return "";};
    
    // Get a contig length
    inline int GetContigLength(const unsigned int& id) const {if (id < contigs_.size()) return contigs_[id].length; else return -1;};
    
    // Get contig id by name. -1 means cannot find the contig
    int GetContigId(std::string a_name) const;
    
    // Prints the vcf file to the screen (Test Purpose)
    void PrintVCF(const SConfig& a_rConfig);

    // Set an ID to the VCF file
    void setID(int a_nVcfId) {m_nVcfId = a_nVcfId;};

    // Returns the Filtering string for the given filter key
    const char* getFilterString(int a_nFilterKey);
    
    // Returns the key id of given Filter string
    int getFilterKey(const char* a_pFilterValue);
    
    //Returns the filter names and descriptions in the vcf file
    void GetFilterInfo(std::vector<std::string>& a_rFilterNames, std::vector<std::string>& a_rFilterDescriptions);
    
    //Returns Header and Record Pointer of htslib
    bcf_hdr_t* GetHeaderPointer();
    bcf1_t* GetRecordPointer();
    
    //Store chromosome name indexes
    std::map<std::string, int> m_chrIndexMap;
    
private:

    /// Trimms the alt string that contains ref allele
    void TrimAllele(SAllele& a_rAllele);
    
    /// Return whether given allele can be trimmed more than one way
    bool HasMultipleTrimming(SAllele& a_rAllele);
    
    ///Check if the first nucleotide for alleles are redundant (for indels)
    bool HasRedundantFirstNucleotide() const;
    
    //Return the chromosome number [0 to 24]
    int GetChromosomeNumber(const std::string& a_chrName) const;
    
    std::string m_filename;
    bool m_bIsOpen;
    htsFile *   m_pHtsFile;
    bcf_hdr_t * m_pHeader;
    bcf1_t *    m_pRecord;  

    std::vector<SVcfContig> contigs_;

    int m_nVcfId;
  
};

#endif //VCF_READER_H_
