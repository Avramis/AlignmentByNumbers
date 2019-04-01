//
//  NoTran.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "NoTran.h"


//Destructor
NoTran::~NoTran(){};
//Constractor
NoTran::NoTran(){};

//Assinge the numerical represenation as a trasnformation
void NoTran:: genTran(double** &rep, double** &tran, int X, int Y){
    Process2dArray<double> CA;
    CA.CreateArray(tran, Y, X);
    for (int i = 0; i < Y; i++)
    {
        for(int j = 0; j < X; j++){
            tran[i][j] = rep[i][j];
        }
    }
    lvl = X;
};

//Return the X dimension transformation level
int NoTran::returnLvL(){
    return lvl;
};
