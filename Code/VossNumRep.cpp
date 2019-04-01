//
//  VossNumRep.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "VossNumRep.h"

//Destructor:
VossNumRep::~VossNumRep(){};

//Constractor
VossNumRep::VossNumRep(){};

//Create the double array using the Process2dArray class
//generate the Voss Indicators representation
void VossNumRep::setRep(double** &rep, std::string read){
    Process2dArray<double> CA;
    X = (int)read.size();
    CA.CreateArray(rep, Y, X);
    
    for (int i = 0; i < X; i++){
        if(std::toupper(read[i]) == 'T')
        {
            rep[0][i] = 0.0;
            rep[1][i] = 0.0;
            rep[2][i] = 0.0;
            rep[3][i] = (1.0);
        }
        else if(std::toupper(read[i]) == 'C'){
            rep[0][i] = (1.0);
            rep[1][i] = 0.0;
            rep[2][i] = 0.0;
            rep[3][i] = 0.0;
        }
        else if(std::toupper(read[i]) == 'G'){
            rep[0][i] = 0.0;
            rep[1][i] = (1.0);
            rep[2][i] = 0.0;
            rep[3][i] = 0.0;
        }
        else if(std::toupper(read[i]) == 'A'){
            rep[0][i] = 0.0;
            rep[1][i] = 0.0;
            rep[2][i] = (1.0);
            rep[3][i] = 0.0;
        }
        else{
            rep[0][i] = 0.0;
            rep[1][i] = 0.0;
            rep[2][i] = 0.0;
            rep[3][i] = 0.0;
        }
    }
    
};


//Return the X length (|------->)of the Voss Indicators representation
int VossNumRep::returnX(){
    return X;
};


//
//                    ( ^ )
//                    ( | )
//Return the Y length ( | ) of the Voss Indicators representation
int VossNumRep::returnY(){
    return Y;
};
