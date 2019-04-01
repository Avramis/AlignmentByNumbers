//
//  AtomNumRep.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "AtomNumRep.h"

//Destructor
AtomNumRep::~AtomNumRep(){};

//Constractor
AtomNumRep::AtomNumRep(){};

//Create the double array using the Process2dArray class
//generate the Atomic-Number representation
void AtomNumRep::setRep(double** &rep, std::string read){
    Process2dArray<double> CA;
    X = (int)read.size();
    CA.CreateArray(rep, Y, X);
    
    for (int i = 0; i < X; i++){
        if(std::toupper(read[i]) == 'T')
        {
            rep[0][i] = 66.0;
        }
        else if(std::toupper(read[i]) == 'C'){
            rep[0][i] = 58.0;
        }
        else if(std::toupper(read[i]) == 'G'){
            rep[0][i] = 78.0;
        }
        else if(std::toupper(read[i]) == 'A'){
            rep[0][i] = 70.0;
        }
        else{
            rep[0][i] = 68.0;
        }
    }
    
};

//Return the X length (|------->)of the Atomic-Number representation
int AtomNumRep::returnX(){
    return X;
};

//
//                    ( ^ )
//                    ( | )
//Return the Y length ( | ) of the Atomic-Number representation
int AtomNumRep::returnY(){
    return Y;
};
