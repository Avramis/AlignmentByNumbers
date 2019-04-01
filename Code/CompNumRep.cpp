//
//  CompNumRep.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "CompNumRep.h"

//Destructor
CompNumRep::~CompNumRep(){};

//Cosntractor
CompNumRep::CompNumRep(){};

//Create the double array using the Process2dArray class
//generate the Complex-Number representation
void CompNumRep::setRep(double** &rep, std::string read){
    Process2dArray<double> CA;
    X = (int)read.size();
    CA.CreateArray(rep, Y, X);
    
    for (int i = 0; i < X; i++){
        if(std::toupper(read[i]) == 'T')
        {
            rep[0][i] = 1.0;
            rep[1][i] = -1.0;
        }
        else if(std::toupper(read[i]) == 'C'){
            rep[0][i] = -1.0;
            rep[1][i] = 1.0;
        }
        else if(std::toupper(read[i]) == 'G'){
            rep[0][i] = -1.0;
            rep[1][i] = -1.0;
        }
        else if(std::toupper(read[i]) == 'A'){
            rep[0][i] = 1.0;
            rep[1][i] = 1.0;
        }
        else{
            rep[0][i] = 0.0;
            rep[1][i] = 0.0;
        }
    }
    
};

//Return the X length (|------->)of the Complex-Number representation
int CompNumRep::returnX(){
    return X;
};

//
//                    ( ^ )
//                    ( | )
//Return the Y length ( | ) of the Complex-Number representation

int CompNumRep::returnY(){
    return Y;
};