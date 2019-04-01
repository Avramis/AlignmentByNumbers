//
//  FastaParser.h
//  Alignment_by_number_v4
//
//  Created by Avraam_Tapinos on 12/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Alignment_by_number_v4__FastaParser__
#define __Alignment_by_number_v4__FastaParser__

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

#include "NGS.h"
#include "Representation.h"
#include "Transformation.h"
#include "Process2dArray.h"
#include "ReverseCompliment.h"
//#include "PrintInfo.h"
class FastaParser
{
private:
    //Variables
    std::string Reference;
    Representation rep;
    Transformation tra;
    Process2dArray<double> CA;
    ReverseCompliment RevCom;

public:
    //Destructor
    ~FastaParser();
    
    //Constractor
    FastaParser(std::string filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * NuSequences, int kmer);
    
    
    //Return the full genome reference nucleotide sequence
    std::string returnRead();
};
#endif /* defined(__Alignment_by_number_v4__FastaParser__) */
