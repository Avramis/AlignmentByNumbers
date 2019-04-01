//
//  ReverseCompliment.cpp
//  Alignment _by_number_v3
//
//  Created by Avraam_Tapinos on 10/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "ReverseCompliment.h"

//Compares nucleotide letters and return the reverse compliment
char ReverseCompliment::NucleotideComparison(char x){
    if(toupper(x) == 'A' )
    {
        x='T';
    }
    else if(toupper(x) == 'T' || toupper(x) == 'U' )
    {
        x='A';
    }
    else if(toupper(x) == 'C' )
    {
        x='G';
    }
    else if(toupper(x) == 'G' )
    {
        x='C';
    }
    else if(toupper(x) == 'N' )
    {
        x='N';
    }
    return x;
}

//Destructor
ReverseCompliment::~ReverseCompliment(){};

//Constructors
//Accepts a read and generates the reverse compliment
ReverseCompliment::ReverseCompliment(std::string read){
    setReverese(read);
};

//Empty constractor
ReverseCompliment::ReverseCompliment(){};


//Generate the reverse complimen of the nucleotide read
void ReverseCompliment::setReverese(std::string read){
    rvread = read;
    std::reverse(rvread.begin(),rvread.end());
    for(int i = 0; i < (int) rvread.size(); i++){
        rvread[i]=NucleotideComparison(rvread[i]);
    }
    
};



//Return the reverse compliment nucleotide
std::string ReverseCompliment::getReverse(){
    return rvread;
}