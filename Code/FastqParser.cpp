//
//  FastqParser.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "FastqParser.h"

//Destructor
FastqParser::~FastqParser(){
};

//Constractor
FastqParser::FastqParser(std::string filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int refkmerlen){
    //Variables
    std::ifstream ngsFile(filepath); //Read file from file_path
    std::string lineContents; // string holding the information from the line i
    std::string read, uid, qualstr;
    int count = 0, posind = 0, count6 = 0,count66=10000,  readlength = 0;
    count66--;
    //clock_t start = clock();
    
    //Start parsing the fastq file until it reach the end of the file
    while(!ngsFile.eof()){
        getline(ngsFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '@')
        {
            //If line starts with @ char and count = 0 then asinge line content to uid
            uid = *new std::string;
            readstream >> uid;
            uid.erase(0,1);
            count++;
        }
        else {
            
            if(count == 1)
            {
                //If line doesn't starts with @ char and count = 1  asinge line content to read
                readlength = *new int;
                read = *new std::string;
                readstream >> read;
                readlength = (int)read.size();
                
                count++;
            }
            else if(count == 2)
            {
                //If line doesn't starts with @ char and count = 2  do nothing
                count++;
            }
            else if(count == 3)
            {
                //If line doesn't starts with @ char and count = 3  asinge line content to quality string
                readstream >> qualstr;
                int reqlen = readlength;
                
                //Check size of read against the kmer size
                //If read size smaller than kmer then refkmer used equals read size
                if(reqlen > 0)
                {
                    if(reqlen >= refkmerlen){
                        reqlen = refkmerlen;
                    }
                }
                
                //Declare read variable
                NGS n(read, uid, (int) read.size(), reqlen, qualstr);
                //Generate forward strand representation
                rep.createRepresentation(repmeth, read.substr(0,reqlen), n.returnFwRep());
                
                
                //rep.accumulate(n.returnFwRep(), rep.returnY(), rep.returnX());
                //rep.normalised(n.returnFwRep(), rep.returnY(), rep.returnX());
                
                
                //Assinge the representation and transforamtion's Y dim.
                n.setY(rep.returnY());
                
                //Assinge the representation's X-dim.
                n.setKx(rep.returnX());
                
                //Generate forward strand transformation
                tra.createTransformation(tranmeth, n.returnFwRep(), n.returnFwTran(), n.returnKx(), n.returnY(), tranlev);
                
                //Assinge the transforamtion's X-dim.
                n.setTx(tra.returnLvL());
                
                //Generate reverse strand represenattion
                //rep.createRepresentation(repmeth, n.returnRvRead().substr(0,n.returnKx()), n.returnRvRep());
                
                //Generate reverse strand transformation
                //tra.createTransformation(tranmeth, n.returnRvRep(), n.returnRvTran(), n.returnKx(), n.returnY(), tranlev);
                
                //Release Forward strand representation's memmory
                CA.DeleteArray(n.returnFwRep());
                
                //Release Reverse strand representation's memmory
                //CA.DeleteArray(n.returnRvRep());
                
                //Assinge n read to the vector

                NuSequences->push_back(n);
                count6++;
                if((count6%1000) == 0){
                    std::cout << "> Fastq read: " << count6 << "\n";
                }
                posind++;
                count = 0;
                
            }
        }
    }
    
    //Close file
    ngsFile.close();
};
