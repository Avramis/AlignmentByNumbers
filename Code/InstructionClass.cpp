/*
 * InstructionClass.cpp
 *
 *  Created on: 3 Sep 2015
 *      Author: avraam
 */

#include "InstructionClass.h"

InstructionClass::InstructionClass(int check){
    std::cout << "--Alignment by numbers--\n";
    if(check == 0){
        std::cout << "Wrong parameter is used\n";
    }
    std::cout<< "Usage:\n";
    std::cout << "Parameters:\n";
    std::cout << "-fa [./genome.fasta]              -- (Required parameter) fasta file directory and file name\n\n";
    std::cout << "-fq [./reads.fastq]               -- (Required parameter) fastq file directory and file name\n\n";
    std::cout << "-o [./samoutputdirectory]         -- (Optional parameter) output folder directory to save the alignment sam file. If not provided the sam file is save in the same directory as the fastq file.\n\n";
    std::cout << "-kmer [x]                         -- (Optional parameter) kmer size to use for the analysis [default: 100]\n\n";
    std::cout << "-rep [Representation method]      -- (Optional parameter) Provide the representation method to use for converting nucleotide sequences to numerical sequences [default: Voss_indicators]\n";
    std::cout << "Options:\n";
    std::cout << "01) Atomic_Numbers\n";
    std::cout << "02) Complex_Numbers\n";
    std::cout << "03) Dna_Walk\n";
    std::cout << "04) EIIP_numbers\n";
    std::cout << "05) Integer_numbers\n";
    std::cout << "06) Pair_numbers\n";
    std::cout << "07) Real_numbers\n";
    std::cout << "08) Tetrahedron\n";
    std::cout << "09) Voss_indicators        --     default option\n";
    std::cout << "10) Z_curve\n\n";
    std::cout << "-tra [Transformation method]      -- (Optional parameter) Provide the transformation method to use for compressing data to a lower dimensional space [default: DFT]\n";
    std::cout << "Options:\n";
    std::cout << "01) DFT                    --     default option\n";
    std::cout << "02) DWT \n";
    std::cout << "03) PAA\n\n";
    std::cout << "-clvl [y]                         -- (Optional parameter) compression level to use for the analysis [default: 2]\n\n";
    std::cout << "-knn [kn]                         -- (Optional parameter)the number of kn neighbor data to use for the KNN search [default: 100]\n\n";
    std::cout << "-s [true/false]                   -- (Optional parameter) run a sensitive search [default: fasle]\n\n";
    std::cout << "-bs [x]                           -- (Optional parameter) set the size of the tree [default: 100000]\n\n";
    std::cout << "-rd [./report_output_directory]   -- (Optional parameter) output directory to save the timing report file. If not provided the report text file is save in the same directory as the fastq file.\n\n";
    std::cout << "Example:";
    std::cout << "./Alignment_by_numbers -fa ./genome.fasta -fq ./reads.fastq\n";
    std::cout << "./Alignment_by_numbers -fa ./genome.fasta -fq ./reads.fastq -kmer 300 -rep Tetrahedron -tra DWT -clvl 4 -knn 30\n";
    std::cout << "./Alignment_by_numbers -fa ./genome.fasta -fq ./reads.fastq -kmer 300 -rep Tetrahedron -tra DWT -o ./samfiledirectory -clvl 4 -knn 30 -s true -bs 10000 -rd ./reportsfile\n";
    
}
