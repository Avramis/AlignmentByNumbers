//
//  Transformation.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 Transformation class generates a compressed transformation for a nucleotide string
 */

#ifndef __Al_by_num_v2_array__Transformation__
#define __Al_by_num_v2_array__Transformation__

#include <stdio.h>
#include <map>
#include "DFTTran.h"
#include "DWTTran.h"
#include "PAATran.h"
#include "NoTran.h"

class Transformation{
private:
    //Variables: a a map list with the data trasforamtion methods and the lvl indicates the transformations X (---->)level
    std::map <std::string, int> TransList;
    int lvl;
    
public:
    //Destructor
    ~Transformation();
    
    //Constractor: Generate the tranformations list
    Transformation();
    
    //Generate the appropoirate transformation according to the transformation method
    void createTransformation(std::string method, double** &rep, double** &tran, int X, int Y, int tlvl);
    
    //Return the X dimensions (---->) of the transformation
    int returnLvL();
};
#endif /* defined(__Al_by_num_v2_array__Transformation__) */
