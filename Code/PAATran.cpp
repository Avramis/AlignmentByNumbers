//
//  PAATran.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "PAATran.h"
//Destructor
PAATran::~PAATran(){};

//Constractor
PAATran::PAATran(){};


//Caluaclte PAA transformation
void PAATran::genTran(double** &rep, double** &tran, int X, int Y, int tlvl){
    Process2dArray<double> CA;
    
    //If the transformation level is smaller than the X dim the PAA trasnforamtion is equal to the representation
    if(X <= tlvl){
        CA.CreateArray(tran, Y, X);
        for (int i = 0; i < Y; i++)
        {
            for(int j = 0; j < X; j++){
                tran[i][j] = rep[i][j];
            }
        }
        lvl = X;
    }
    //If the transformation level is bigger than the X dim
    else{
        CA.CreateArray(tran, Y, tlvl);
        for (int i = 0; i < Y; i++){
            double paawave = 0.0, remain = 0.0;
            int count1 = -1;
            for(int j = 0; j < tlvl; j++){
                count1++;
                int count2 = 1;
                //Calculate the jth paa coefficient
                while((count2 * tlvl)+ (int)remain < X)
                {
                    paawave += (rep[i][count1]* (double)tlvl);
                    count2++;
                    count1++;
                }
                
                remain = ((count2 * (double)tlvl) + remain) - (double)X;
                paawave += (rep[i][count1] * ((double)tlvl - remain));
                tran[i][j] = (paawave/ (double) X);
                paawave = 0.0;
                paawave += (rep[i][count1] * (double)remain);
            }
        }
        
        lvl = tlvl;
    }
    
};

//Return (int) transformation level
int PAATran::returnLvL(){
    return lvl;
};