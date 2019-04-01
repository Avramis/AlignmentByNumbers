//
//  DnaWNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__DnaWNumRep__
#define __Al_by_num_v2_array__DnaWNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

/*
 Geenrate a DNA-Walk representation for a nucleotide string
 */

class DnaWNumRep{
private:
    //Variables
    int Y = 2,
        X;
public:
    //Destructor
    ~DnaWNumRep();
    
    //Constractor
    DnaWNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the DNA-Walk representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the DNA-Walk representation
    int returnX();
    
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the DNA-Walk representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__DnaWNumRep__) */
