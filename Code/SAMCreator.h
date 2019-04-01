/*
 * SAMCreator.h
 *
 *  Created on: 3 Sep 2015
 *      Author: avraam
 */

#ifndef SAMCREATOR_H_
#define SAMCREATOR_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include "NGS.h"
#include "ReverseCompliment.h"


//Inline function is used to sort the reads inside the vector in descending order according to their alinged position
inline
bool NGSComp (NGS r1, NGS r2){
    return r1.returnAlPos() < r2.returnAlPos();
    
};

class SamCreator
{
private:
    
public:
    //Destructor
    ~SamCreator();
    
    //Constructor generates a sam file for the aligned reads
    SamCreator(std::string directory,std::string repmethod, std::string transmethod, std::string kmerlen,std::string translvlstr, std::vector<NGS> *reads, std::map<std::string, std::string> *reference,std::string outputdir, std::string sns, std::string strknn);
    
};

#endif /* SAMCREATOR_H_ */
