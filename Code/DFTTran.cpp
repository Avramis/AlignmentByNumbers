//
//  DFTTran.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "DFTTran.h"


//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Private:                                         //
//////////////////////////////////////////////////////////////////////////////////////////////

//Function to calcualte the Fourier Frequency decomposition
void DFTTran::fourier(double** &data, int nn){
    int i,
    n,
    m,
    j,
    mmax,
    istep;
    double wtemp,
    wr,
    wpr,
    wpi,
    wi,
    theta,
    tempr,
    tempi;
    //double **ffttran;
    
    n = nn<<1;
    j = 1;
    for (i=1; i < n; i+=2) {
        if (j>i) {
            double tempswapj, tempswapi;
            tempswapj = data[0][j-1];
            tempswapi = data[0][i-1];
            data[0][i-1] = tempswapj;
            data[0][j-1] = tempswapi;
            tempswapj = data[0][j];
            tempswapi = data[0][i];
            data[0][i] = tempswapj;
            data[0][j] = tempswapi;
        }
        m = nn;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = (-isign)*(2.0*M_PI/((double)mmax));
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr*data[0][j-1] - wi*data[0][j];
                tempi = wr * data[0][j] + wi*data[0][j-1];
                
                
                data[0][j-1] = data[0][i-1] - tempr;
                data[0][j] = data[0][i] - tempi;
                data[0][i-1] += tempr;
                data[0][i] += tempi;
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }
    
};



//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Public:                                          //
//////////////////////////////////////////////////////////////////////////////////////////////

//Destructor
DFTTran::~DFTTran(){};

//Constractor
DFTTran::DFTTran(){};


//Caluaclte DFT transformation
void DFTTran::genTran(double** &rep, double** &tran, int X, int Y, int tlvl){
    Process2dArray<double> CA;
    
    //If X smaller than trasnforamtion level then tlvl becames equal to X
    if(X <= tlvl){
        tlvl = X;
    }
    
    //generate trasnforamtion double array
    CA.CreateArray(tran, Y, tlvl*2);
    
    //Variable required for DFT decomposition calculcation
    int np2,
    nn;
    nn= X;
    np2 = (int) (ceil(log2((double)nn)));
    if (nn !=  pow(2,np2 )){
        nn = (int)(pow(2,np2 ));
    }
    
    
    for (int i = 0; i < Y; i++){
        double **ffttran;
        
        //generate ffttran double array
        CA.CreateArray(ffttran, 1, nn*2);
        
        //Fill ffttran array
        for (int j=0; j < nn; j++){
            if(j < X){
                ffttran[0][j*2] = rep[i][j];
            }
            else{
                ffttran[0][j*2] = 0;
            }
            ffttran[0][(j*2)+1] = 0;
        }
        
        
        DFTTran::fourier(ffttran,nn);
        
        //Delete the ffttran array
       
        if(tlvl >= nn){
            tlvl = nn;
        }
        
        //Assinge FFT to Transformation
        for (int z = 0; z < tlvl*2; z++){
            tran[i][z] = ffttran[0][z]/ sqrt((double)nn);
        }
         CA.DeleteArray(ffttran);
    }
    
    
    //Assigne tlvl x 2 to lvl
    lvl =tlvl*2;
};

//Return (int) transformation level
int DFTTran::returnLvL(){
    return lvl;
};