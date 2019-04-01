//
//  ZCurNumRep.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__ZCurNumRep__
#define __Al_by_num_v2_array__ZCurNumRep__

#include <stdio.h>
#include <string>
#include "Process2dArray.h"

class ZCurNumRep{
private:
    //Variables: the X and Y dimensions of the double array
    int Y = 3,
    X;
public:
    
    //Destructor:
    ~ZCurNumRep();
    
    //Constractor
    ZCurNumRep();
    
    //Create the double array using the Process2dArray class
    //generate the Z-Curve representation
    void setRep(double** &rep, std::string read);
    
    //Return the X length (|------->)of the Z-Curve representation
    int returnX();
    
    //
    //                    ( ^ )
    //                    ( | )
    //Return the Y length ( | ) of the Z-Curve representation
    int returnY();
    
};
#endif /* defined(__Al_by_num_v2_array__ZCurNumRep__) */
