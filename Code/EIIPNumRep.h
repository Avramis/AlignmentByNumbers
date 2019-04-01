//
//  EIIPNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__EIIPNumRep__
#define __Al_by_num_v2_array__EIIPNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

/*
 EIIPNumRep class generates a Z-Curve representation for a nucleotide string
 */

class EIIPNumRep{
private:
    //Variable
    int Y = 1,
        X;
public:
    //Destructor
    ~EIIPNumRep();
    
    //Constractor
    EIIPNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the EIIP Numerical representation
    void setRep(double** &rep, std::string read);
    
    
    //Return the X length (|------->)of the EIIP Numerical representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the EIIP Numerical representation
    int returnY();
    
};

#endif /* defined(__Al_by_num_v2_array__EIIPNumRep__) */
