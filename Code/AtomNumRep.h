//
//  AtomNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__AtomNumRep__
#define __Al_by_num_v2_array__AtomNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

/*
 AtomNumRep class generate an Atomic-Number representation for a nuceotide sequence
 */

class AtomNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 1,
        X;
public:
    //Destructor
    ~AtomNumRep();
    
    //Constractor
    AtomNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the Atomic-Number representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Atomic-Number representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Atomic-Number representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__AtomNumRep__) */
