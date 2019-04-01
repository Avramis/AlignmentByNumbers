//
//  NoTran.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__NoTran__
#define __Al_by_num_v2_array__NoTran__

#include <stdio.h>
#include "Process2dArray.h"

/*
 NoTran class generates assinged a full resolution represenation as transformatiomn for a nucleotide string
 */

class NoTran{
private:
    //Variable
    int lvl;
public:
    //Destructor
    ~NoTran();
    //Constractor
    NoTran();
    
    //Assinge the numerical represenation as a trasnformation
    void genTran(double** &rep, double** &tran, int X, int Y);
    
    //Return the X dimension transformation level
    int returnLvL();
};

#endif /* defined(__Al_by_num_v2_array__NoTran__) */
