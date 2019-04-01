//
//  Representation.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "Representation.h"

//Destructor
Representation::~Representation(){};

//Constructor filling the RepresentationsList map
Representation::Representation(){
    RepsList["Atomic_Numbers"] = 0;
    RepsList["Complex_Numbers"] = 1;
    RepsList["Dna_Walk"] = 2;
    RepsList["EIIP_numbers"] = 3;
    RepsList["Integer_numbers"] = 4;
    RepsList["Pair_numbers"] = 5;
    RepsList["Real_numbers"] = 6;
    RepsList["Tetrahedron"] = 7;
    RepsList["Voss_indicators"] = 8;
    RepsList["Z_curve"] = 9;
};

//Generate the appropoirate representation according to the representation method
void Representation::createRepresentation(std::string method, std::string read, double** &rep){
    switch (RepsList.at(method)){
        case 0:{
            AtomNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 1:{
            CompNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 2:{
            DnaWNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 3:{
            EIIPNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 4:{
            IntNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 5:{
            PairNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 6:{
            RealNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 7:{
            TetrNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 8:{
            VossNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        case 9:{
            ZCurNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
        default:{
            TetrNumRep r;
            r.setRep(rep, read);
            Y = r.returnY();
            X = r.returnX();
            break;
        }
    }
};

void Representation::accumulate(double** &rep, int Y, int X){
    for (int  i = 0; i < Y; i++){
        for (int j = 1; j < X; j++){
            rep[i][j]+=rep[i][j-1];
        }
    }
};



void Representation::normalised(double** &rep, int Y, int X){
    double **meanv, **stde;
    Process2dArray <double>CA;
    CA.CreateArray(meanv, Y, 1);
    CA.CreateArray(stde, Y, 1);
    

    for (int i = 0; i < Y; i++){
        meanv[i][0] = rep[i][0];
        
        for (int j = 1; j < X; j++){
            meanv[i][0]+= rep[i][j];
            //std::cout << rep[i][j] << "\t||\t"<< meanv[i][0] << "\n";
        }
        meanv[i][0] = meanv[i][0]/((double)X);
    }
    
    
    
    
    for(int  i = 0; i < Y; i++){
        double difx = rep[i][0] - meanv[i][0];
        stde[i][0] = (difx * difx);
        for (int j = 1; j < X; j++){
            difx = rep[i][j] - meanv[i][0];
            stde[i][0] += (difx * difx);
        }
        stde[i][0] = sqrt(stde[i][0]/((double)(X)));
    }
    
    
    
    for (int i = 0; i < Y; i++){
        for (int j = 0; j < X; j++){
            rep[i][j] = (rep[i][j] - meanv[i][0])/stde[i][0];
        }
    }
    CA.DeleteArray(meanv);
    CA.DeleteArray(stde);
};
    //std::cout << "%After normilization:\n an = [";
    //for (int i = 0; i < Y; i++){
    //    std::cout << rep[i][0];
    //    for (int j = 1; j < X; j++){
    //        std::cout << ", " << rep[i][j];
    //    }
    //    std::cout << ";";
    //}
    //std::cout << "];\n";
    //std::cout << "\n";
    /*
    double meanv = 0.0, stde = 0.0;
    //Process2dArray <double>CA;
    //CA.CreateArray(meanv, Y, 1);
    //CA.CreateArray(stde, Y, 1);
    
    for (int i = 0; i < Y; i++){
        for (int j = 0; j < X; j++){
            meanv+= rep[i][j];
            //std::cout << rep[i][j] << "\t||\t"<< meanv[i][0] << "\n";
        }
    }
    meanv = meanv/(X*Y);
    
    for (int i = 0; i < Y; i++){
        //stde[i][0] = pow(rep[i][0] - meanv[i][0],2);
        for (int j = 1; j <X; j++){
            stde = stde+ (pow(rep[i][j] - meanv,2));
        }
        
    }
    stde =  stde/(X*Y);
    
    for (int i = 0; i < Y; i++){
        if(std::isnan(meanv) == true){
            meanv = 0;
        }
        
        if(std::isinf(meanv) == true){
            meanv = 0;
        }
        
        
        
        if(std::isnan(stde) == true){
            stde = 1;
        }
        
        if(std::isinf(stde) == true){
            stde = 1;
        }
        
        
        for (int j = 0; j < X; j++){
            //std::cout << "Before norm rep[" << i << "][" << j << "] = " << rep[i][j] <<"\n";
            rep[i][j] = (rep[i][j]-meanv)/stde;
            //std::cout << "After norm rep[" << i << "][" << j << "] = " << rep[i][j] <<"\n";
        }
    }
    
    //CA.DeleteArray(meanv);
    //CA.DeleteArray(stde);
    
    */
    
    
    
    
//};




//                        (  ^ )
//                        (  | )
//Return the Y dimensions (Y | ) of the representation
//                        (  | )
int Representation::returnY(){
    return Y;
}

//Return the X dimensions (---->) of the transformation
int Representation::returnX(){
    return X;
}