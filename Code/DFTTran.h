//
//  DFTTran.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__DFTTran__
#define __Al_by_num_v2_array__DFTTran__

#include <stdio.h>
#include <math.h>
#include "Process2dArray.h"

/*
 DFTTran class generates a Haar DFT transformation for a nucleotide string
 */

class DFTTran{
private:
    //Variable
    int lvl,
        isign = -1;

    
    //Calcualte the DFT Decomposition
    void fourier(double** &data, int nn);

    
public:
    //Destructor
    ~DFTTran();
    
    //Constractor
    DFTTran();
    
    
    //Caluaclte DFT transformation
    void genTran(double** &rep, double** &tran, int X, int Y, int tlvl);
    
    //Return (int) transformation level
    int returnLvL();
};
#endif /* defined(__Al_by_num_v2_array__DFTTran__) */
