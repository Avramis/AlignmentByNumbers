//
//  PrintInfo.cpp
//  Alignment_by_number_v4
//
//  Created by Avraam_Tapinos on 12/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "NGS.h"

//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Private:                                         //
//////////////////////////////////////////////////////////////////////////////////////////////

//Void Functions:
//Assinge r to variable read
void NGS::setRead(std::string r){
    read = r;
};

//Assinge id to variable Nuid
void NGS::setID (std::string id){
    Nuid = id;
};

//Assinge q to variable qualstr
void NGS::setQual(std::string q){
    qual = q;
};

//Assinge x to variable fX (Full length of the sequence)
void NGS::setFullX(int x){
    fX = x;
};

//Assinge idxp to variable idx (indexing number of the read in the file)
void NGS::setIDX(int idxp){
    idx = idxp;
};

//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Public:                                          //
//////////////////////////////////////////////////////////////////////////////////////////////


//Constractors:
NGS::NGS(){};
NGS::NGS(std::string R, std::string ID, int X, int K, int ind){
    setRead(R);
    setID(ID);
    setFullX(X);
    //setY(Y);
    //setTx(Z);
    setKx(K);
    setIDX(ind);
};

NGS::NGS(std::string R, std::string ID, int X, int K, std::string qualstr){
    setRead(R);
    setID(ID);
    setFullX(X);
    setKx(K);
    setQual(qualstr);
};

//Void Functions:
//Assinge values to variables
void NGS::setRead(std::string R, std::string ID, int X, int K, int ind){
    setRead(R);
    setID(ID);
    setFullX(X);
    //setY(Y);
    //setTx(Z);
    setKx(K);
    setIDX(ind);
};

void NGS::setRead(std::string R, std::string ID, int X, int Y, int Z, int K, std::string qualstr){
    setRead(R);
    setID(ID);
    setFullX(X);
    setY(Y);
    setTx(Z);
    setKx(K);
    setQual(qualstr);
};


//Assinge value y to variable nY (the Y dimensions of the representation and the transformation)
void NGS::setY(int y){
    nY = y;
};

//Assinge value k to variable nX (the selected kmer dimensions of the representation)
void  NGS::setKx(int k){
    nX = k;
};

//Assinge value z to variable tX (the x dimensions of the transformation)
void NGS::setTx(int z){
    tX = z;
};

//Assinge i to dir (The alingment direction of the read. False equals to forward strand || True equals to reverse strand)
void NGS::setDirection(bool i){
    dir = i;
};


void NGS::setAlDirection(bool i){
    alDir.push_back(i);
};



//Assinge i to ap (The alingment position of the read to the reference genome)
void NGS::setAliPos(int i){
    ap.push_back(i);
};

//Assinge c to cigar (The alignemnt cigar string of the alinged read)
void NGS::setCigar(std::string c){
    cigar.push_back(c);
};

//Assinge alingment Ref Name

void NGS::setAlRef(std::string ref){
    alRef.push_back(ref);
};

//Clear alingment vectors
void NGS::clearAlVec(){
    std::vector<int>().swap(ap);
    std::vector<bool>().swap(alDir);
    std::vector<std::string>().swap(cigar);
    std::vector<std::string>().swap(alRef);
};


void NGS::UniAl(){
    if(ap.size()> 1){
        //std::map <std::pair<std::pair<std::string , int>, std::pair<std::string, bool>>, int> uniqueAl;
        std::map<std::pair<std::pair<std::string , bool>, std::pair<std::string, int> >, int> uniqueAl;
        for(int i = 0; i < ap.size(); i++){
            //std::pair<std::string , int> pair1;
            //std::pair<std::string, bool> pair2;
            std::pair<std::pair<std::string , bool>, std::pair<std::string, int> > paircom;
            //std::cout << alRef.at(i) << "\t" << alDir.at(i) << "\t" << cigar.at(i) << "\t" << ap.at(i) << "\n";
            paircom = std::make_pair(std::make_pair(alRef.at(i), alDir.at(i)), std::make_pair(cigar.at(i),ap.at(i)));
            uniqueAl[paircom] = 1;
        }
        
        clearAlVec();
        for (std::map<std::pair<std::pair<std::string , bool>, std::pair<std::string, int> >, int>::iterator it = uniqueAl.begin(); it != uniqueAl.end(); it++){
            //std::cout << it->first.first.first << "\t" << it->first.first.second << "\t" << it->first.second.first << "\t" << it->first.second.second << "\n";
            alRef.push_back(it->first.first.first);
            alDir.push_back(it->first.first.second);
            cigar.push_back(it->first.second.first);
            ap.push_back(it->first.second.second);
            //std::cout << it->first.first.first << "\n";
            //std::cout << it->first.first.second << "\n";
            //std::cout << it->first.second.first << "\n";
            //std::cout << it->first.second.second << "\n";
        }
    }
};


//Return Statements:
//Return the forward strand read
std::string NGS::returnRead(){
    return read;
};

std::string NGS::returnRead() const{
    return read;
};

//Return the reverse strand read
std::string NGS::returnRvRead(){
    ReverseCompliment rvread;
    rvread.setReverese(NGS::returnRead());
    return rvread.getReverse();
};

//Return the read's id
std::string NGS::returnID(){
    return Nuid;
};


std::string NGS::returnID() const{
    return Nuid;
};

//Return the read's quality string
std::string NGS::returnQual(){
    return qual;
};

//Return the read's cigar string
std::vector<std::string> NGS::returnCigar(){
    return cigar;
};


//Return alignment Ref
std::vector<std::string> NGS::returnAlRef(){
    return alRef;
};


//Return the representation and transformation's Y-dim
int NGS::returnY(){
    return nY;
};

//Return the representation and transformation's Y-dim (const)
int NGS::returnY() const{
    return nY;
};

//Return the transformation's X-dim
int NGS::returnTx(){
    return tX;
};

//Return the transformation's X-dim (const)
int NGS::returnTx() const{
    return tX;
};

//Return the representations's X-dim
int NGS::returnKx(){
    return nX;
};

//Return the representations's X-dim (const)
int NGS::returnKx() const{
    return nX;
};

//Return the read's full X-dim
int NGS::returnFullX(){
    return fX;
};

//Return the read's full X-dim (const)
int NGS::returnFullX() const{
    return fX;
};

//Return the read's alignemnt positiion
std::vector<int> NGS::returnAlPos(){
    return ap;
};

//Return the read's indexing position
int NGS::returnIdx(){
    return idx;
};


int NGS::returnIdx() const{
    return idx;
};

//Return the read's forward strand representation
double** &NGS::returnFwRep(){
    return nFwrep;
};

//Return the read's reverese strand representation
double** &NGS::returnRvRep(){
    return nRvrep;
};

//Return the read's forward strand transformation
double** &NGS::returnFwTran(){
    return nFwtran;
};

//Return the read's forward strand transformation (const)
double** NGS::returnFwTran() const{
    return nFwtran;
};

//Return the read's reverse strand transformation
double** &NGS::returnRvTran(){
    return nRvtran;
};

//Return the read's reverse strand transformation (const)
double** NGS::returnRvTran() const{
    return nRvtran;
};

//Return the read's alignemnt direction
bool NGS::returnDir(){
    return dir;
};

bool NGS::returnDir() const{
    return dir;
};

std::vector<bool> NGS::returnAlDir(){
    return alDir;
}


