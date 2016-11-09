#include <stdlib.h>
#include <iostream>
#include "CVcfReader.h"
#include "CFastaReader.h"
#include "CPathReplay.h"

int main (int argc, char** argv)
{
    if (argc != 4) 
    {
        std::cout << "Argument Count is wrong!! include a vcf and fasta file." << std::endl;
        return 1;       
    }
      
    CPathReplay pathReplay;
    pathReplay.InitializeReaders(argv[1], argv[2], argv[3]);
    pathReplay.FindBestPath(21);
       
    SResult finalResult;
}
