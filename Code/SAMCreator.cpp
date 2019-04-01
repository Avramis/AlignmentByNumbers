/*
 * SAMCreator.cpp
 *
 *  Created on: 3 Sep 2015
 *      Author: avraam
 */

/*
 SAMCreator class generates a sam file for the aligned reads
 */

#include "SAMCreator.h"
//Destructor
SamCreator::~SamCreator(){};

//Constructor generates a sam file for the aligned reads
SamCreator::SamCreator(std::string directory,std::string repmethod, std::string transmethod, std::string kmerlen,std::string translvlstr, std::vector<NGS> *reads, std::map<std::string, std::string> *reference,std::string outputdir, std::string sns, std::string strknn){
    
    std::size_t found = directory.find_last_of(std::string("/"));
    std::string dir = outputdir;
    
    struct stat sb;
    
    if (!(stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
    {
        mkdir(dir.c_str(), 0777);
        //int statval = mkdir(dir.c_str());
        //std::cout << dir <<" directory does not exist.\n";
    }
    
    //Generate the file name
    dir+=directory.substr(found+1, directory.size()-found+1)+"_"+repmethod+"_"+kmerlen+ "kmer_"+transmethod+"_"+ translvlstr + "complvl_" + sns + "_KNN_" + strknn + ".sam";
    
    //Sort reads in a descending order based on their alingment position
    sort(reads->begin(), reads->end(), NGSComp);
    std::ofstream o(dir);
    
    //Generate the sam file
    //Line 1
    o << "@HD" << "\t" << "VN:1.0" << "\t" << "SO:unsorted" << "\n";
    
    //Line 2
    for(std::map<std::string, std::string>::iterator it = reference->begin(); it != reference->end(); it++){
        o << "@SQ" << "\t" << "SN:" << it->first << "\t" << "LN:" << it->second.size() << "\n";
    }
    //Line 3
    o << "@PG" << "\t" << "ID:TimeSeriesAlignmentScript" << "\t" << "PN:TimeSeriesScript" << "\t" << "VN:1"  << "\t" << "CL:Fill Later" << "\n";
    for (int i =0; i < reads->size(); i++)
    {
        
        if( reads->at(i).returnAlPos().size() > 1){
            std::cout << "> " << reads->at(i).returnID() << "\n" ;
        }
        
        for(int j = 0; j < reads->at(i).returnAlPos().size(); j++){
            //tab_1
            o << reads->at(i).returnID() << "\t" ;
            
            //tab_2
            if (!reads->at(i).returnAlDir().at(j))
            {
                o << "0" << "\t";
            }
            else
            {
                o << "16" << "\t";
            }
            
            //tab_3
            //o << refid << "\t";
            o << reads->at(i).returnAlRef().at(j) <<"\t";
            
            //tab_4
            o << reads->at(i).returnAlPos().at(j)+1 << "\t";
            
            //tab_5
            o << "40" << "\t";
            
            //tab_6
            o << reads->at(i).returnCigar().at(j) << "\t";
            
            //tab_7
            o << "*" << "\t";
            
            //tab_8
            o << "0" << "\t";
            
            //tab_9
            o << "0" << "\t";
            
            //tab_10
            // Not save reverse read in sam files
            //o << reads->at(i).getFwRead() << "\t";
            if(!reads->at(i).returnAlDir().at(j))
            {
                o << reads->at(i).returnRead() << "\t";
            }
            else{
                //r.setReverese(reads->at(i).returnRead());
                //o << r.getReverse() << "\t";
                o << reads->at(i).returnRvRead() << "\t";
            }
            
            //tab_11
            o << reads->at(i).returnQual() <<"\n";
        }
    
    }
    
    o.close();
    //end =clock();
    
    //std::cout << "> Creating SAM file elapse time " << (end-start)/CLOCKS_PER_SEC << ".\n";
};
