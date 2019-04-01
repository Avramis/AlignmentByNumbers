//
//  PrintInfo.cpp
//  Alignment_by_number_v4
//
//  Created by Avraam_Tapinos on 12/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "PrintInfo.h"

//Print multivaraite array details
PrintInfo::PrintInfo(std::string name, double **doublearray, int X, int Y){
    std::cout << "\n" << name << ":\n";
    for(int i = 0; i < Y; i++){
        std::cout << "Row " << i << ": ";
        for (int j = 0; j < X; j++){
            std::cout << doublearray[i][j] << " | ";
        }
        std::cout << "\n";
    }
};

PrintInfo::PrintInfo(std::string name, double **doublearray, int X, int Y,bool x){
    if(x == true){
        std::cout << "\n\n\n\nclear\nclc\n";
    }
    std::cout << "\n%" << name <<"\n%Dimensions Y = " << Y <<" Dimensions X = " << X << "\n";
    std::cout <<  name << " = [";
    for(int i = 0; i < Y; i++){
        std::cout << doublearray[i][0];
        //std::cout << "Row " << i << ": ";
        for (int j = 1; j < X; j++){
            std::cout << ", " <<doublearray[i][j];
        }
        std::cout << ";";
    }
    std::cout<< "];\n";
};