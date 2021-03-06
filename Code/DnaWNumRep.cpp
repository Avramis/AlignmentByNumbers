//
//  DnaWNumRep.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "DnaWNumRep.h"

//Destructor
DnaWNumRep::~DnaWNumRep(){};

//Constractor
DnaWNumRep::DnaWNumRep(){};

//Create the double array using the Process2dArray class
//generate the DNA-Walk representation
void DnaWNumRep::setRep(double** &rep, std::string read){
    Process2dArray<double> CA;
    X = (int)read.size();
    CA.CreateArray(rep, Y, X);
    
    for (int i = 0; i < X; i++){
        if(std::toupper(read[i]) == 'T')
        {
            rep[0][i] = 0.0;
            rep[1][i] = 1.0;
        }
        else if(std::toupper(read[i]) == 'C'){
            rep[0][i] = 1.0;
            rep[1][i] = 0.0;
        }
        else if(std::toupper(read[i]) == 'G'){
            rep[0][i] = -1.0;
            rep[1][i] = 0.0;
        }
        else if(std::toupper(read[i]) == 'A'){
            rep[0][i] = 0.0;
            rep[1][i] = -1.0;
        }
        else{
            rep[0][i] = 0.0;
            rep[1][i] = 0.0;
        }
    }
    
};

//Return the X length (|------->)of the DNA-Walk representation
int DnaWNumRep::returnX(){
    return X;
};


//
//                    ( ^ )
//                    ( | )
//Return the Y length ( | ) of the DNA-Walk representation
int DnaWNumRep::returnY(){
    return Y;
};