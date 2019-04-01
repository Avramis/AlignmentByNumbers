//
//  Representation.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 Representation class generates a numerical sequence representation for a nucleotide string
 */

#ifndef __Al_by_num_v2_array__Representation__
#define __Al_by_num_v2_array__Representation__

#include <stdio.h>
#include <iostream>
#include <string>
#include <map>

#include "AtomNumRep.h"
#include "CompNumRep.h"
#include "DnaWNumRep.h"
#include "EIIPNumRep.h"
#include "IntNumRep.h"
#include "PairNumRep.h"
#include "RealNumRep.h"
#include "TetrNumRep.h"
#include "VossNumRep.h"
#include "ZCurNumRep.h"

#include "Process2dArray.h"
#include <iostream>
#include <cmath>
#include <cfloat>
//#include <math.h>  

class Representation{
private:
    //
    //                                                                               ( ^        )
    //                                                                               ( |        )
    // Variables: a a map list with the data representation methods and the XY levels(Y|        ) of the representation
    //                                                                               ( |------> )
    //                                                                               (     X    )
    std::map <std::string, int> RepsList;
    int X,
        Y;
    
    
public:
    //Destructor
    ~Representation();
    
    //Constructor filling the RepresentationsList map
    Representation();
    
    //Generate the appropoirate representation according to the representation method
    void createRepresentation(std::string method, std::string read, double** &rep);
    
    void accumulate(double** &rep, int Y, int X);
    
    void normalised(double** &rep, int Y, int X);

    //                        (  ^ )
    //                        (  | )
    //Return the Y dimensions (Y | ) of the representation
    //                        (  | )
    int returnY();
    
    //Return the X dimensions (---->) of the transformation
    int returnX();};

#endif /* defined(__Al_by_num_v2_array__Representation__) */
