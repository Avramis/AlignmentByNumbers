//
//  PAATran.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__PAATran__
#define __Al_by_num_v2_array__PAATran__

#include <stdio.h>
#include "Process2dArray.h"

/*
 PAATran class generates a PAA transformation for a nucleotide string
 */


class PAATran{
private:
    //Variables the compression level of the transformation
    int lvl;
public:
    //Destructor
    ~PAATran();
    
    //Constractor
    PAATran();
    
    
    //Caluaclte PAA transformation
    void genTran(double** &rep, double** &tran, int X, int Y, int tlvl);
    
    //Return (int) transformation level
    int returnLvL();
};

#endif /* defined(__Al_by_num_v2_array__PAATran__) */
