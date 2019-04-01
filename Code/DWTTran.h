//
//  DWTTran.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__DWTTran__
#define __Al_by_num_v2_array__DWTTran__

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Process2dArray.h"

/*
 DWTTran class generates a Haar DWT transformation for a nucleotide string
 */

class DWTTran{
private:
    //Variable
    int lvl;
public:
    //Destructor
    ~DWTTran();
    
    //Constractor
    DWTTran();
    
    
    //Caluaclte DWT transformation
    void genTran(double **rep, double** &tran, int X, int Y, int tlvl);
    
    //Return (int) transformation level
    int returnLvL();
};
#endif /* defined(__Al_by_num_v2_array__DWTTran__) */
