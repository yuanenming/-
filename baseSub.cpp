#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <locale>

#include "baseSub.h"


using namespace std;

//Pull a random base from bases array depenting on their weights
char getBaseFromPool(string * bases, int weights){
    locale loc;    
    srand (time(NULL));
    
    int rnd = rand() % 100;
    int weights_sum = 100;
    
    for( int i = 0; i < (*bases).length(); i++ ) {
        if(rnd <= weights){
            return toupper( (*bases)[i], loc);
        }
        rnd -= weights;
    }        
}

char getTheRespectiveBase( char representation ){
    string as;
    switch(representation){
        //In case of b (not A) the rest have the same occurrence probability (0.33)
        case 'B':
            as = "tcg";
            return getBaseFromPool( &as, 33 );
            
        //In case of d (not C) the rest have the same occurrence probability (0.33)
        case 'D':
            as = "atg";
            return getBaseFromPool( &as, 33);
            
        //In case of h (not G) the rest have the same occurrence probability (0.33)
        case 'H':
            as = "atc";
            return getBaseFromPool( &as, 33);
            
        //In case of k, G and T (Ketones) have the same occurrence probability (0.5)
        case 'K':
            as = "gt";
            return getBaseFromPool( &as, 50);
            
        //In case of m, C and A (with Amino group) have the same occurrence probability (0.5)
        case 'M':
            as = "ca";
            return getBaseFromPool( &as, 50);
            
        //In case of n all bases have the same occurrence probability (0.25)
        case 'N':
            as = "atcg";
            return getBaseFromPool( &as, 25);
            
        //In case of r, A and G (Purines) have the same occurrence probability (0.5)
        case 'R':
            as = "ag";
            return getBaseFromPool( &as, 50);
            
        //In case of s, G and C (Strong interaction) have the same occurrence probability (0.5)
        case 'S':
            as = "cg";
            return getBaseFromPool( &as, 50);
            
        //In case of v (not T) the rest have the same occurrence probability (0.33)
        case 'V':
            as = "agc";
            return getBaseFromPool( &as, 33);
            
        //In case of w, A and T (Weak interaction) have the same occurrence probability (0.5)
        case 'W':
            as = "at";
            return getBaseFromPool( &as, 50);
            
        //In case of y, C and T (Pyrimidines) have the same occurrence probability (0.5)
        case 'Y':
            as = "ct";
            return getBaseFromPool( &as, 50);
            
        
    }
}

std::vector<std::string> cleanData(std::vector<std::string> Set){
    
    for( vector<string>::iterator seq = Set.begin(); seq != Set.end(); ++seq ){
        for ( int i = 0; i < (*seq).length(); i++ ){
            // If the base is not one the simple 4.
            if ( (*seq)[i] != 'A' && (*seq)[i] != 'T' && 
                    (*seq)[i] != 'C' && (*seq)[i] != 'G' ){
                // Find one.
                (*seq)[i] = getTheRespectiveBase((*seq)[i]);
            }
        }
    }
    return Set;
}

