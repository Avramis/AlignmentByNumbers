//
//  main.cpp
//  Check_Rev_ref
//
//  Created by Avraam_Tapinos on 21/04/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Alignment_Execution.h"
#include "InstructionClass.h"
#include "FastaParser2.h"
#include "Process2dArray.h"

int main(int argc, const char * argv[]) {
    std::string fastadir, fastqdir, repmethsel = "Voss_indicators",
    tranmethsel = "DFT", outputdir, repdir;
    int tranlev = 2, kmer = 100, knn = 100;
    bool sens = false;
    bool fqc = false, fac = false, outputprov = false, repdirpro = false, sr = false;
    std::vector <std::string> fastalist;
    int civ = (argc - 1) % 2;
    
    
    if (civ || argc - 1 == 0) {
        for (int i = 1; i < argc; i++) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionClass instreuctions(1);
                return 1;
            }
        }
        InstructionClass instreuctions(0);
        return 1;
    }
    else {
        for (int i = 1; i < argc; i += 2) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionClass instreuctions(0);
            }
            else if ((std::string) (argv[i]) == "-fa") {
                fastadir = (std::string) (argv[i + 1]);
                if (FILE *file = fopen(fastadir.c_str(), "r")) {
                    fclose(file);
                    fac = true;
                    fastalist.push_back(fastadir);
                }
                
            }
            else if ((std::string) (argv[i]) == "-fq") {
                fastqdir = (std::string) (argv[i + 1]);
                if (FILE *file = fopen(fastqdir.c_str(), "r")) {
                    fclose(file);
                    fqc = true;
                }
                else {
                    fqc = false;
                }
            }
            else if ((std::string) (argv[i]) == "-kmer") {
                kmer = atoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-rep") {
                repmethsel = (std::string) (argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-tra") {
                tranmethsel = (std::string) (argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-clvl") {
                tranlev = atoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-knn") {
                knn = atoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-o") {
                outputdir = (std::string) (argv[i + 1]) + "/";
                outputprov = true;
            }
            else if ((std::string) (argv[i]) == "-s") {
                if( (std::string) argv[i + 1] == "true"){
                    sens = true;
                }
            }
            else if((std::string) (argv[i]) == "-rd"){
                repdir = (std::string) (argv[i + 1]) + "/";
                repdirpro = true;
            }
            else if((std::string) (argv[i]) == "-rs"){
                if(atoi(argv[i+1]) == 1){
                    sr = true;
                }
            }
            else {
                InstructionClass instreuctions(0);
            }
        }
    }
    
    if (outputprov == false) {
        outputdir = fastqdir;
    }
    
    if (repdirpro == false) {
        repdir = fastqdir;
    }
    
    struct stat sb;
    
    if (!(stat(outputdir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
        mkdir(outputdir.c_str(), 0777);
    }
    
    if ((fac && fqc) == false) {
        InstructionClass instreuctions(0);
        return 1;
    }
    
    
    std::cout << "\n";
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    std::cout << "%%%%%%%%%%                                  Alignment by numbers                                  %%%%%%%%%%\n";
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    Alignment_Execution AE(fastalist, fastqdir, outputdir, repmethsel , tranmethsel, tranlev, kmer, knn, sens, repdir, sr);
    //Alignment_Execution(fastadir, fastqdir, outputdir, repmethsel , tranmethsel, tranlev, kmer, knn, sens, repdir,sr);
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    
    return 0;
}
