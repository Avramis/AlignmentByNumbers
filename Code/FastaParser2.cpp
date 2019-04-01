//
//  FastaParser2.cpp
//  Al_By_num_multiple_genomes
//
//  Created by Avraam_Tapinos on 27/05/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "FastaParser2.h"

//Destructor
FastaParser2::~FastaParser2(){
    //std::cout << "Destruct FastaParser\n";
};

FastaParser2::FastaParser2(std::vector<std::string> filepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * NuSequences, int kmer){
    std::string uid, read, lineentry;
    int count = 0, readlength = 0, Kmer;
    bool record = false;
    std::map<std::string, int> Reflengths;
    //clock_t start=clock();
    Kmer = kmer;
    for (int fi = 0; fi < filepath.size();fi++){
        read = "";
        uid = "";
        readlength = 0;
        std::ifstream ngsFile(filepath.at(fi));
        while(!ngsFile.eof()){
            getline(ngsFile, lineentry);
            std::stringstream readstream(lineentry);
            if(count == 0 && lineentry[0] == '>'){
                
                if(record == true){
                    
                    Reference[uid] = read;
                    Reflengths[uid] = readlength;
                    read = "";
                    uid = "";
                    readlength = 0;
                    record =false;
                }
                
                readstream >> uid;
                uid.erase(0,1);
                count++;
                //read = *new std::string;
                readlength = 0;
                
                
            }
            else{
                record = true;
                std::string curstr;
                readstream >> curstr;
                curstr.erase(std::remove(curstr.begin(), curstr.end(), '\n'), curstr.end());
                //std::replace (curstr.begin(),curstr.end(), '\n' , '');
                readlength += (int)curstr.size();
                read+=curstr;
                count=0;

            }
        }
        
        if(record == true){
            Reference[uid] = read;
            Reflengths[uid] = readlength;
            record =false;
        }
        ngsFile.close();
    }
    
    int coutre = 0;
    for (std::map<std::string,std::string>::iterator it = Reference.begin(); it!=Reference.end(); it++){
        for (int i = 0; i <  (Reflengths[it->first] - kmer)+1; i++){
            coutre++;
            if((coutre%1000) == 0){
                 std::cout << "> Reference subsequence: " << coutre <<"\n";
                //coutre = 0;
                //std::cout << it->first << "\n";
                //std::cout << Reflengths[it->first] << "\n";
            }
            
            //if(i == 9719){
            //    std::cout << "Time to break\n";
            //    std::cout << "Time to break\n";
            //}
            NGS n(it->second.substr (i,Kmer), it->first, Reflengths[it->first], Kmer,i);
            
            
            n.setDirection(false);
            //Generate forward strand representation
            rep.createRepresentation(repmeth, it->second.substr (i,Kmer), n.returnFwRep());
            
            
            
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
            RevCom.setReverese(it->second.substr (i,Kmer));
            std::string revread= RevCom.getReverse();
            NGS nr(revread, it->first, Reflengths[it->first], Kmer,i);
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
    }
    
    std::cout << "Total number of reads in vector: " << NuSequences->size() <<"\n";
};

//Return the full genome reference nucleotide sequence
std::map<std::string,std::string> FastaParser2::returnRead(){
    return Reference;
};