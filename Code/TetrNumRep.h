//
//  TetrNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 TetrNumRep class generates a Tetrahedron representation for a nucleotide string
 */

#ifndef __Al_by_num_v2_array__TetrNumRep__
#define __Al_by_num_v2_array__TetrNumRep__

#include <stdio.h>
#include <string>
#include <math.h>

#include "Process2dArray.h"

class TetrNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 3,
    X;
public:
    //Destructor:
    ~TetrNumRep();
    
    //Constractor
    TetrNumRep();
    
    
    //Create the double array using the Process2dArray class
    //generate the Tetrahedron representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Tetrahedron representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Tetrahedron representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__TetrNumRep__) */
