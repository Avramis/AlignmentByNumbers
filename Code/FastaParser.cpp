//
//  FastaParser.cpp
//  Alignment_by_number_v4
//
//  Created by Avraam_Tapinos on 12/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "FastaParser.h"

//Destructor
FastaParser::~FastaParser(){
    //std::cout << "Destruct FastaParser\n";
};

//Constractor
FastaParser::FastaParser(std::string filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * NuSequences, int kmer){
    //Variables
    std::string uid, read, lineentry;
    int count = 0, readlength = 0, Kmer;
    std::ifstream ngsFile(filepath);
    //clock_t start=clock();
    Kmer = kmer;
    
    //Start parsing the fasta  file until it reach the end of the file
    while(!ngsFile.eof())
    {
        getline(ngsFile, lineentry);
        std::stringstream readstream(lineentry);
        if(count == 0 && lineentry[0] == '>')
        {
            //If line starts with > char and count = 0 then asinge line content to uid
            //uid = *new std::string;
            readstream >> uid;
            uid.erase(0,1);
            count++;
            //read = *new std::string;
            readlength = 0;
        }
        else {
            //If line doesn't starts with > char and count = 1  asinge line content to read
            std::string curstr;
            readstream >> curstr;
            curstr.erase(std::remove(curstr.begin(), curstr.end(), '\n'), curstr.end());
            //std::replace (curstr.begin(),curstr.end(), '\n' , '');
            readlength += (int)curstr.size();
            read+=curstr;
            count=0;
            //std::cout << count << "\n";
        }
        
    }
    //close file
    ngsFile.close();
    
    //Assinge full read to Reference variable
    Reference = read;
    
    //Start extracting k-size subsequences from referecne and create Represenations and trasnformations

    //Process2dArray<double> CA;
    
    for (int i = 0; i < (Reference.size()-Kmer)+1; i++){
        //std::cout << "Subsequence Read: " << i << "\n";
        //Representation rep;
        //Transformation tra;
        
        if((i%1000) == 0){
            std::cout << "> Reference subsequence: " << i <<"\n";
        }
        
        //std::cout << "Fasta: " << i << "\n";
        
        //Declare read variable
        
        NGS n(read.substr (i,Kmer), uid, (int)read.size(), Kmer,i);
        n.setDirection(false);
        //Generate forward strand representation
        rep.createRepresentation(repmeth, read.substr (i,Kmer), n.returnFwRep());
        
        
        
        
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
        
        //Release Forward strand representation's memmory
        CA.DeleteArray(n.returnFwRep());
        
        //Assinge n read to the vector
        NuSequences->push_back(n);
        
        //std::string revrt;
        //std::cout << read.substr (i,Kmer) << "\n";
        //revrt = read.substr((read.size()-i)-Kmer, Kmer);
        //std::cout << revrt << "\n";
        RevCom.setReverese(read.substr (i,Kmer));
        std::string revread= RevCom.getReverse();
        NGS nr(revread, uid+"_r", (int)read.size(), Kmer,i);
        nr.setDirection(true);
        //Generate forward strand representation
        rep.createRepresentation(repmeth, nr.returnRead().substr (0,Kmer), nr.returnFwRep());
        
        
        
        //rep.accumulate(nr.returnFwRep(), rep.returnY(), rep.returnX());
        //rep.normalised(nr.returnFwRep(), rep.returnY(), rep.returnX());
        
        
        
        nr.setY(rep.returnY());
        nr.setKx(rep.returnX());
        tra.createTransformation(tranmeth, nr.returnFwRep(), nr.returnFwTran(), nr.returnKx(), nr.returnY(), tranlev);
        nr.setTx(tra.returnLvL());
        CA.DeleteArray(nr.returnFwRep());
        NuSequences->push_back(nr);
    }
};


//Return the full genome reference nucleotide sequence
std::string FastaParser::returnRead(){
    return Reference;
};