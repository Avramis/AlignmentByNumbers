//
//  FastqParser.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__FastqParser__
#define __Al_by_num_v2_array__FastqParser__

#include <stdio.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <cassert>
#include <time.h>
#include <string>
#include <string.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <exception>

#include "Representation.h"
#include "Transformation.h"
#include "NGS.h"
#include "Process2dArray.h"

/*FastqParser class imports reads from a fastq file and converts reads to numerical representations 
and create the data transformations.
 */
class FastqParser{
private:
    //Variables
    Representation rep;
    Transformation tra;
    Process2dArray<double> CA;
public:
    //Destructor
    ~FastqParser();
    
    //Constractor
    FastqParser(std::string filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int refkmerlen);
    
};
#endif /* defined(__Al_by_num_v2_array__FastqParser__) */
