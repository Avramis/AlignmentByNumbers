//
//  Transformation.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "Transformation.h"

//Destructor
Transformation::~Transformation(){};

//Constractor: Generate the tranformations list
Transformation::Transformation(){
    TransList["DFT"] = 1;
    TransList["DWT"] = 2;
    TransList["PAA"] = 3;
    TransList["None"] = 4;
};

//Generate the appropoirate transformation according to the transformation method
void Transformation::createTransformation(std::string method, double** &rep, double** &tran, int X, int Y, int tlvl){
    switch(TransList.at(method)){
        case 1:{
            DFTTran t;
            t.genTran(rep, tran, X, Y, tlvl);
            lvl = t.returnLvL();
            break;
        }
        case 2:{
            DWTTran t;
            t.genTran(rep, tran, X, Y, tlvl);
            lvl = t.returnLvL();
            break;
        }
        case 3:{
            PAATran t;
            t.genTran(rep, tran, X, Y, tlvl);
            lvl = t.returnLvL();
            break;
        }
        case 4:{
            NoTran t;
            t.genTran(rep, tran, X, Y);
            lvl = t.returnLvL();
            break;
        }
        default:{
            NoTran t;
            t.genTran(rep, tran, X, Y);
            lvl = t.returnLvL();
            break;
        }
    }
};

//Return the X dimensions (---->) of the transformation
int Transformation::returnLvL(){
    return lvl;
}
