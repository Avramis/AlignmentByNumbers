//
//  Alignment_Execution.h
//  Al_by_num_v2_array
//
//  Created by Avraam Tapinos on 20/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__Alignment_Execution__
#define __Al_by_num_v2_array__Alignment_Execution__

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>



#include <sstream>
#include <cmath>
#include <numeric>
#include <cassert>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <exception>
#include <map>
#include <cstring>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

//#include "VpTree.h"
//#include "vptree.h"
#include "vptree.hh"

//#include "VpTree1.h"
//#include "VpTree2.h"
//#include "MVpTree.h"
//#include "MVpTree2.h"
#include "NGS.h"
#include "ssw_cpp.h"
#include "SAMCreator.h"
#include "FastaParser.h"
#include "FastaParser2.h"
#include "FastqParser.h"
#include "Representation.h"
#include "Transformation.h"
#include "ReverseCompliment.h"
#include "Process2dArray.h"
#include "ReverseCompliment.h"



/**/
 //Case 1 : ||ΣΣ(n1ij - n2ij)^2||
class ngstree: public VPTree<NGS>
{
public:
    virtual double distance(NGS const &nr1, NGS const &nr2){
        double** n1;
        double** n2;
        n1 = nr1.returnFwTran();
        n2 = nr2.returnFwTran();
        double nudist=0.0;
        int x, y;
        y = (int)nr1.returnY();
        if((int)nr1.returnTx() <= nr2.returnTx()){
            x = nr1.returnTx();
        }
        else{
            x = nr2.returnTx();
        }
        
        for (int i = 0; i < y; i++){
            for(int j = 0; j < x; j++){
                double nucdist =n1[i][j] - n2[i][j];
                nudist+=(nucdist*nucdist);
                //nudist+= pow((n1.returnFwTran()[i][j] - n2.returnFwTran()[i][j]),2);
            }
        }
        return sqrt(nudist);
        //return nudist;
    };
};
/**/


/*
//Case 2 : Σ||Σ(n1ij - n2ij)^2||
class ngstree: public VPTree<NGS>
{
public:
    virtual double distance(NGS const &nr1, NGS const &nr2){
        double** n1;
        double** n2;
        n1 = nr1.returnFwTran();
        n2 = nr2.returnFwTran();
        double nudist=0.0;
        int x, y;
        y = (int)nr1.returnY();
        
        if((int)nr1.returnTx() <= nr2.returnTx()){
            x = nr1.returnTx();
        }
        else{
            x = nr2.returnTx();
        }
        for (int i = 0; i < y; i++){
            double nudist2 = 0;
            for(int j = 0; j < x; j++){
                double nucdist =n1[i][j] - n2[i][j];
                nudist2+=(nucdist*nucdist);
                //nudist+= pow((n1.returnFwTran()[i][j] - n2.returnFwTran()[i][j]),2);
            }
            nudist += sqrt(nudist2);
            //nudist += (nudist2);
        }
        return nudist;
    };
};
*/

/*
 //Case 3 : max||Σ(n1ij - n2ij)^2||
 class ngstree: public VPTree<NGS>
 {
 public:
 virtual double distance(NGS const &nr1, NGS const &nr2){
 double** n1;
 double** n2;
 n1 = nr1.returnFwTran();
 n2 = nr2.returnFwTran();
 double nudist=0.0;
 
 int x, y;
 y = (int)nr1.returnY();
 if((int)nr1.returnTx() <= nr2.returnTx()){
 x = nr1.returnTx();
 }
 else{
 x = nr2.returnTx();
 }
 
 for (int i = 0; i < y; i++){
 double nudist2 = 0;
 for(int j = 0; j < x; j++){
 double nucdist =n1[i][j] - n2[i][j];
 nudist2+=(nucdist*nucdist);
 }
 nudist2 = sqrt(nudist2);
 if(nudist < nudist2){
 nudist = nudist2;
 }
 }
 return nudist;
 };
 
 };
 */




/*
 Alignment_Execution class initiate the alignment of short reads to a reference genome
 */

inline
double Fwdist(NGS const n1, NGS const n2){
    double nudist=0.0;
    int x, y;
    y = (int)n1.returnY();
    if((int)n1.returnTx() <= n2.returnTx()){
        x = n1.returnTx();
    }
    else{
        x = n2.returnTx();
    }
    
    for (int i = 0; i < y; i++){
        for(int j = 0; j < x; j++){
            double nucdist =n1.returnFwTran()[i][j] - n2.returnFwTran()[i][j];
            nudist+=(nucdist*nucdist);
            //nudist+= pow((n1.returnFwTran()[i][j] - n2.returnFwTran()[i][j]),2);
        }
    }
    
    return sqrt(nudist);
    //return nudist;
}


inline
double Fwdist1(double** n1, double** n2, int x, int y){
    double nudist=0.0;
    for (int i = 0; i < y; i++){
        for(int j = 0; j < x; j++){
            double nucdist =n1[i][j] - n2[i][j];
            nudist+=(nucdist*nucdist);
        }
    }
    return sqrt(nudist);
}


/*
inline
double Fwdist1(double** n1, double** n2, int x, int y){
    
    
    double nudist=0.0;
    for (int i = 0; i < y; i++){
        double val = 0;
        for(int j = 0; j < x; j++){
            //val+= ((n1[i][j]-n2[i][j])*(n1[i][j]-n2[i][j]));
            val+= ((n1[i][j]*n2[i][j]));
        }
        nudist+= (val);
    }
    //std::cout << "test before sq: " << nudist <<"\n";
    if(nudist < 0){

        nudist = sqrt(nudist*nudist);
    }
    //
    //test = sqrt(test*test);
    //std::cout <<"after sqrt: "<< nudist << "\n";
    return sqrt(nudist);
}
*/






inline
double Nudist(double** &refrep, double** &qfrep, int y, int x){
    double nudist=0.0;
    for (int i = 0; i < y; i++){
        for(int j = 0; j < x; j++){
            double nucdist = refrep[i][j] - qfrep[i][j];
            nudist+=(nucdist*nucdist);
            //nudist+=pow((refrep[i][j] - qfrep[i][j]),2);
        }
    }
    
    
    return nudist;
}







class Alignment_Execution{
private:
    //Variable
    std::map<std::string, std::string> reference;
    ReverseCompliment RevCo;
    //Call FastaParser to extract information from fats a file
    //void executeFastaParser(std::string fafilepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * RefReads, int kmer);
    void executeFastaParser(std::vector<std::string> fafilepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * RefReads, int kmer);
    
    //Call FastqParser to extract information from fats a file
    void executeFastqParser(std::string fqfilepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int kmer);
    
    public:
    //Destractor
    ~Alignment_Execution();
    
    //Constractor
    //Alignment_Execution(std::string fafilepath, std::string fqfilepath, std::string outdir, std::string repmeth ,std::string tranmeth, int tranlev, int kmer, int knn, bool sens, std::string repdir, bool sr);
    
    
    Alignment_Execution(std::vector<std::string> fafilepath, std::string fqfilepath, std::string outdir, std::string repmeth ,std::string tranmeth, int tranlev, int kmer, int knn, bool sens, std::string repdir, bool sr);
    
    
    
};
#endif /* defined(__Al_by_num_v2_array__Alignment_Execution__) */
