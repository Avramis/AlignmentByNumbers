//
//  VossNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 TetrNumRep class generates a Voss Indicators representation for a nucleotide string
 */

#ifndef __Al_by_num_v2_array__VossNumRep__
#define __Al_by_num_v2_array__VossNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

class VossNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 4,
    X;
public:
    //Destructor:
    ~VossNumRep();
    
    //Constractor
    VossNumRep();
    
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Voss Indicators representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Voss Indicators representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__VossNumRep__) */
