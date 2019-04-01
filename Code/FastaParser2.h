//
//  FastaParser2.h
//  Al_By_num_multiple_genomes
//
//  Created by Avraam_Tapinos on 27/05/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_By_num_multiple_genomes__FastaParser2__
#define __Al_By_num_multiple_genomes__FastaParser2__

#include <stdio.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <cassert>
#include <time.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <exception>
#include <map>

#include "NGS.h"
#include "Representation.h"
#include "Transformation.h"
#include "Process2dArray.h"
#include "ReverseCompliment.h"
class FastaParser2
{
private:
    //Variables
    std::map<std::string,std::string> Reference;
    Representation rep;
    Transformation tra;
    Process2dArray<double> CA;
    ReverseCompliment RevCom;
    
public:
    //Destructor
    ~FastaParser2();
    
    //Constractor
    FastaParser2(std::vector<std::string> filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * NuSequences, int kmer);
    
    
    //Return the full genome reference nucleotide sequences
    std::map<std::string,std::string> returnRead();

};
#endif /* defined(__Al_By_num_multiple_genomes__FastaParser2__) */
