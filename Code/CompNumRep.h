//
//  CompNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__CompNumRep__
#define __Al_by_num_v2_array__CompNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"
/*
 CompNumRep to generate a Complex-Number represenation of a nucleotide sequence
 */

class CompNumRep{
private:
    //Variable
    int Y = 2,
        X;
public:
    //Destructor
    ~CompNumRep();
    
    //Cosntractor
    CompNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the Complex-Number representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Complex-Number representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Complex-Number representation

    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__CompNumRep__) */
