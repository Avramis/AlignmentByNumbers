//
//  IntNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__IntNumRep__
#define __Al_by_num_v2_array__IntNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

/*
 IntNumRep class generates a Integer Number representation for a nucleotide string
 */

class IntNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 1,
    X;
public:
    //Destructor:
    ~IntNumRep();
    
    //Constractor
    IntNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the Integer Number representation
    void setRep(double** &rep, std::string read);
    
     //Return the X length (|------->)of the Real Number representation
    int returnX();
    
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Integer Number representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__IntNumRep__) */
