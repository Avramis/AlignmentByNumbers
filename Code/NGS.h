//
//  NGS.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__NGS__
#define __Al_by_num_v2_array__NGS__

#include <stdio.h>
#include <string>
#include <vector>
#include <map>

#include "ReverseCompliment.h"
#include "PrintInfo.h"

//NGS class holds all the information required for a nucleotide sequecne.
//E.g, the read name, nulceotide sequence, transformations, representations, trasnformation leveles XY
//dimensions for the transformations and representations, alignment positions, alingment directions etc.

class NGS{
//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Private:                                         //
//////////////////////////////////////////////////////////////////////////////////////////////
private:
    //Variables:
    std::string Nuid,
    read,
    qual;
    
    int nX,
    nY,
    tX,
    fX,
    idx;
    
    double **nFwrep,
           **nRvrep,
           **nFwtran,
           **nRvtran;
    
    bool dir;
    std::vector<int> ap;
    std::vector<bool> alDir;
    std::vector<std::string> cigar;
    std::vector<std::string> alRef;
    
    //Void Functions:
    //Assinge r to variable read
    void setRead(std::string r);
    
    //Assinge id to variable Nuid
    void setID (std::string id);
    
    //Assinge q to variable qualstr
    void setQual(std::string q);;
    
    //Assinge x to variable fX (Full length of the sequence)
    void setFullX(int x);
    
    //Assinge idxp to variable idx (indexing number of the read in the file)
    void setIDX(int idxp);

    //void setFwRep(){
    //};

    //void setRvRep(){
    //};

    //void setFwTran(){
    //};

    //void setRvTran(){
    //};

    
    
//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Public:                                          //
//////////////////////////////////////////////////////////////////////////////////////////////
public:
    //Variables:
    
    //Constractors:
    NGS();
    NGS(std::string R, std::string ID, int X, int K, int ind);
    
    NGS(std::string R, std::string ID, int X, int K, std::string qualstr);
    
    //Void Functions:
    //Assinge values to variables
    void setRead(std::string R, std::string ID, int X, int K, int ind);
    
    void setRead(std::string R, std::string ID, int X, int Y, int Z, int K, std::string qualstr);


    //Assinge value y to variable nY (the Y dimensions of the representation and the transformation)
    void setY(int y);
    
    //Assinge value k to variable nX (the selected kmer dimensions of the representation)
    void  setKx(int k);
    
    //Assinge value z to variable tX (the x dimensions of the transformation)
    void setTx(int z);
    
    //Assinge i to dir (The alingment direction of the read. False equals to forward strand || True equals to reverse strand)
    void setDirection(bool i);
    
    //Assinge i to alingment directions (The alingment direction of the read. False equals to forward strand || True equals to reverse strand)
    void setAlDirection(bool i);
    
    //Assinge i to ap (The alingment position of the read to the reference genome)
    void setAliPos(int i);
    
     //Assinge c to cigar (The alignemnt cigar string of the alinged read)
    void setCigar(std::string c);
    
    //Assinge alingment Ref Name
    void setAlRef(std::string ref);
    
    //Clear alingment vectors
    void clearAlVec();
    
    //Check for dupliacate alignments
    void UniAl();
    
    
    //Return Statements:
    //Return the forward str;
    std::string returnRead();
    std::string returnRead() const;
    
    //Return the reverse strand read
    std::string returnRvRead();
    
    //Return the read's id
    std::string returnID();
    
    std::string returnID() const;
    
    //Return the read's quality string
    std::string returnQual();
    
    //Return the read's cigar string
    std::vector<std::string> returnCigar();
    
    //Return alignment Ref
    std::vector<std::string> returnAlRef();

    //Return the representation and transformation's Y-dim
    int returnY();
    
    //Return the representation and transformation's Y-dim (const)
    int returnY() const;

    //Return the transformation's X-dim
    int returnTx();
    
    //Return the transformation's X-dim (const)
    int returnTx() const;
    
    //Return the representations's X-dim
    int returnKx();

    //Return the representations's X-dim (const)
    int returnKx() const;

    //Return the read's full X-dim
    int returnFullX();

    //Return the read's full X-dim (const)
    int returnFullX() const;
    
    //Return the read's alignemnt positiion
    std::vector<int> returnAlPos();
    
    //Return the read's indexing position
    int returnIdx();
    
    int returnIdx() const;
    
    //Return the read's forward strand representation
    double** &returnFwRep();

    //Return the read's reverese strand representation
    double** &returnRvRep();

    //Return the read's forward strand transformation
    double** &returnFwTran();

    //Return the read's forward strand transformation (const)
    double** returnFwTran() const;

    //Return the read's reverse strand transformation
    double** &returnRvTran();

    //Return the read's reverse strand transformation (const)
    double** returnRvTran() const;

    //Return the read's direction
    bool returnDir();
    
    bool returnDir()const;
    
    //Return the reads alignment direction
    std::vector<bool> returnAlDir();

};

#endif /* defined(__Al_by_num_v2_array__NGS__) */
