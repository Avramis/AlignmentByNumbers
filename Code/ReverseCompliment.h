//
//  ReverseCompliment.h
//  Alignment _by_number_v3
//
//  Created by Avraam_Tapinos on 10/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

/*
 ReverseCompliment class generates the reverse compliment of a nucleotide sequence
 */

#ifndef __Alignment__by_number_v3__ReverseCompliment__
#define __Alignment__by_number_v3__ReverseCompliment__

#include <stdio.h>
#include <string>
#include <algorithm>

class ReverseCompliment{
private:
    //Variable
    std::string rvread;
    
    //Compares nucleotide letters and return the reverse compliment
    char NucleotideComparison(char x);
    
public:
    //Destructor
    ~ReverseCompliment();
    
    //Constructors:
    //Accepts a read and generates the reverse compliment
    ReverseCompliment(std::string read);
    
    //Empty constractor
    ReverseCompliment();
    
    
    //Generate the reverse complimen of the nucleotide read
    void setReverese(std::string read);
    
    //Return the reverse compliment nucleotide
    std::string getReverse();
    
};
#endif /* defined(__Alignment__by_number_v3__ReverseCompliment__) */
