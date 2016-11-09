/*
 * Copyright (c) 2016 Seven Bridges Genomics
 * Author: Berke Cagkan Toptas
 */

#ifndef VCF_READER_H_
#define VCF_READER_H_
#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "CVariant.h"

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
    
    // Get next record in the file
    bool GetNextRecord(CVariant* a_pVariant);
    
    // Get the filename
    std::string GetFilename() const {return m_filename;};
    
    // Get the number of samples in vcf
    int GetNumberOfSamples() const;
    
    // Get contig size
    inline int GetContigSize() const {return (int)(contigs_.size());};
    
    // Get a contig name
    inline std::string GetContigName(const unsigned int& id) const {if (id < contigs_.size()) return contigs_[id].name; else return "";};
    
    // Get a contig length
    inline int GetContigLength(const unsigned int& id) const {if (id < contigs_.size()) return contigs_[id].length; else return -1;};
    
    // Get contig id by name. -1 means cannot find the contig
    int GetContigId(const char* name) const;
    
    // Prints the vcf file to the screen (Test Purpose)
    void PrintVCF();

    // Set an ID to the VCF file
    void setID(int a_nVcfId) {m_nVcfId = a_nVcfId;};


private:
    std::string m_filename;
    bool m_bIsOpen;
    htsFile *   m_pHtsFile;
    bcf_hdr_t * m_pHeader;
    bcf1_t *    m_pRecord;  

    std::vector<SVcfContig> contigs_;

    int m_nVcfId;    
  
};

#endif //VCF_READER_H_
