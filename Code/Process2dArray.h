//
//  Process2dArray.h
//  Al_by_num_v2_array
//
//  Created by Avraam_Tapinos on 19/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#ifndef __Al_by_num_v2_array__Process2dArray__
#define __Al_by_num_v2_array__Process2dArray__

#include <stdio.h>

template <typename T>
class Process2dArray{
private:
    
public:
    //Constractor
    //Create empty ContinuousArray object for class T
    Process2dArray(){};
    
    //Return statement
    //Create and return the 2d array of class T
    T** CreateArray(int N1, int N2){
        T ** array = new T * [N1];
        
        array[0] = new T  [N1*N2];
        
        int i,j;
        
        for( i = 0; i < N1; i++) {
            
            if (i < N1 -1 ) {
                
                array[i+1] = &(array[0][(i+1)*N2]);
                
            }
            
            for( j = 0; j < N2; j++) {
                if (j > 0) array[i][j] = array[i][j-1];
            }
            
        }
        
        return array;
    };
    
    //Create  the 2d array
    void CreateArray(T** &array, int N1, int N2){
        array = new T * [N1];
        
        array[0] = new T  [N1*N2];
        
        int i,j;
        
        for( i = 0; i < N1; i++) {
            
            if (i < N1 -1 ) {
                //
                //    array[0][(i+1)*N2] = &(array[0][0][(i+1)*N2]);
                
                array[i+1] = &(array[0][(i+1)*N2]);
                
            }
            
            for( j = 0; j < N2; j++) {
                if (j > 0) array[i][j] = array[i][j-1];
            }
            
        }
    };


    //Delete the 2d array of class T
    void DeleteArray(T** &array){
        delete[] array[0];
        delete[] array;
    };
};
#endif /* defined(__Al_by_num_v2_array__Process2dArray__) */
