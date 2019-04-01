//
//  PrintInfo.h
//  Alignment_by_number_v4
//
//  Created by Avraam_Tapinos on 12/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Alignment_by_number_v4__PrintInfo__
#define __Alignment_by_number_v4__PrintInfo__

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "NGS.h"

class PrintInfo{
private:
public:
    //Print multivaraite array details
    PrintInfo(std::string name, double **doublearray, int X, int Y);
    PrintInfo(std::string name, double **doublearray, int X, int Y,bool x);    
};

#endif /* defined(__Alignment_by_number_v4__PrintInfo__) */
