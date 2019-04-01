//
//  DWTTran.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "DWTTran.h"


//Destructor
DWTTran::~DWTTran(){};

//Constractor
DWTTran::DWTTran(){};


//Return (int) transformation level
int DWTTran::returnLvL(){
    return lvl;
};



/*
//Caluaclte DWT transformation
void DWTTran::genTran(double **rep, double** &tran, int X, int Y, int tlvl){
    Process2dArray<double> CA;
    
    //Calcualte the DWT cycles to perform
    int complev = (int)(ceil(log((double)X)/log(2.0))-ceil(log((double) tlvl)/log(2.0)));
    
    //If complevel equals 0 then DWT trasnforamtion is equal to the representation
    if(complev == 0){
        CA.CreateArray(tran, Y, X);
        for (int i = 0; i < Y; i++)
        {
            for(int j = 0; j < X; j++){
                tran[i][j] = rep[i][j];
            }
        }
        lvl = X;
    }
    else{
        //If the complevel is bigger than the 0 start perfroming the transformation
        
        //the window size used to calculate the wavelet
        int win = (int)(pow(2,ceil(log((double)X)/log(2.0)))/pow(2, (ceil(log((double) tlvl)/log(2.0)))));
        
        //The transformation X dim
        lvl = (int)floor(X/win);
        //pow(2,(ceil(log((double) tlvl)/log(2.0))));
        
        //Create the transformation double array
        CA.CreateArray(tran, Y, lvl);
        
        for(int i = 0; i < Y; i++){
            for(int j = 0; j < lvl; j++){
                double wav = 0;
                
                //Calcaulte the wavelet j
                for (int z = j*win; z < (j+1)*win; z++){
                    wav += rep[i][z];
                }
                
                //Assinge wavelet value to transforation
                tran[i][j] = wav/(double)win;
            }
        }
    }
};
*/




//Caluaclte DWT transformation
//New Proper approach
void DWTTran::genTran(double **rep, double** &tran, int X, int Y, int tlvl){
    Process2dArray<double> CA;
    
    //Calcualte the DWT cycles to perform
    int complev = (int)(ceil(log((double)X)/log(2.0))-ceil(log((double) tlvl)/log(2.0)));
    
    //If complevel equals 0 then DWT trasnforamtion is equal to the representation
    if(complev == 0){
        CA.CreateArray(tran, Y, X);
        for (int i = 0; i < Y; i++)
        {
            for(int j = 0; j < X; j++){
                tran[i][j] = rep[i][j];
            }
        }
        lvl = X;
    }
    else{
        
        int inip = (int)(pow(2,ceil(log((double)X)/log(2.0))));//pow(2, (ceil(log((double) tlvl)/log(2.0)))));
        
        double **temparray;
        CA.CreateArray(temparray, Y, inip);
        
        for (int i = 0; i < Y; i++){
            for (int j = 0; j < X; j++){
                temparray[i][j] = rep[i][j];
            }
            for(int j = X; j < inip; j++){
                temparray[i][j] = 0.0;
            }
        }
        
        int inip2 = (int) ceil(log((double)X)/log(2.0));
        int finp2 = (int) ceil(log((double)tlvl)/log(2.0));
        
        //for(int i = 0; i < (inip2-(finp2+1)); i++){
        for(int i = 0; i < (inip2-(finp2)); i++){
            //for(int j = 0; j <)
            //if(i == 4){
                //do nothing
            //}
            for (int j = 0; j < Y; j++){
                for (int z = 0; z < (int) pow(2,inip2 - (i)); z+=2){
                    //temparray[j][z/2] = (temparray[j][z] + temparray[j][z+1]);
                    temparray[j][z/2] = (temparray[j][z] + temparray[j][z+1])/2;
                    //temparray[j][z/2] = temparray[j][z] - ((temparray[j][z] + temparray[j][z+1])/2);
                }
            }
        }
        
        CA.CreateArray(tran, Y, (int)pow(2,finp2));
        
        for(int i = 0; i < Y; i ++){
            for (int j = 0; j < (int)pow(2,finp2); j++){
                tran[i][j] = temparray[i][j];
            }
        }
        
        CA.DeleteArray(temparray);
        lvl = (int) pow(2,finp2);
    }
};
