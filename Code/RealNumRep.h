//
//  RealNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 RealNumRep class generates a Real Number representation for a nucleotide string
 */

#ifndef __Al_by_num_v2_array__RealNumRep__
#define __Al_by_num_v2_array__RealNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

class RealNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 1,
    X;
public:
    //Destructor:
    ~RealNumRep();
    
    //Constractor
    RealNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the Real Numebr representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Real Number representation
    int returnX();
    
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Real Number representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__RealNumRep__) */
