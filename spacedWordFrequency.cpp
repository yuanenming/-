#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <time.h>
#include <ctime>
#include <iomanip>
#include <math.h> 

#include "spacedWordFrequency.h"

using namespace std;

double** spacedWordFrequency(map<int, map<int, map<string,int> > > f) {
    vector<vector<double> > rel_freq_avg;
    string M[256];
    string S[256];
    string E[256];
    ifstream mid ("5mid.txt");
    ifstream start ("5start.txt");
    ifstream end ("5end.txt");
    string line;
    for(int i=0;i<256;i++) {
        getline(mid,line);
        M[i] = line;
        getline(start,line);
        S[i] = line;
        getline(end,line);
        E[i] = line;
    }
    
    int k = 5; //we choose k-mers 
    
    vector<int> tmpV;
    //ofstream off ("SWFDistanceMatrix.txt");    
    map<string, int> tmpM;
    typedef std::map<string, int>::iterator it_type;
    int freq = 0;
    double all5mers;
    
    for(int s=0;s<f.size();s++) { //for every sequence and k=5
        vector<double> tmp;
        double patternsMid[256] = {0}; 
        double patternsStart[256] = {0};  
        double patternsEnd[256] = {0}; 
        all5mers = 0;    
        tmpM = f[s][k];
        for(it_type iterator = tmpM.begin(); iterator != tmpM.end(); iterator++) { // here we input k-tuples for each sequence
                all5mers += iterator->second;
                for(int i=0;i<256;i++) {
                    if(((iterator->first)[0] == M[i][0])&&((iterator->first)[1] == M[i][1])&&((iterator->first)[3] == M[i][3])&&((iterator->first)[4] == M[i][4])) patternsMid[i] += iterator->second;
                    if(((iterator->first)[1] == S[i][1])&&((iterator->first)[2] == S[i][2])&&((iterator->first)[3] == S[i][3])&&((iterator->first)[4] == S[i][4])) patternsStart[i] += iterator->second;
                    if(((iterator->first)[0] == E[i][0])&&((iterator->first)[1] == E[i][1])&&((iterator->first)[2] == E[i][2])&&((iterator->first)[3] == E[i][3])) patternsEnd[i] += iterator->second;
                }
        }
        for(int i=0;i<256;i++) {
            tmp.push_back((patternsMid[i]/all5mers + patternsStart[i]/all5mers + patternsEnd[i]/all5mers)/3);
        }
        rel_freq_avg.push_back(tmp);
    }

    double** distanceMatrix = new double*[f.size()]; 

    for(int i=0;i<f.size();i++) { 
        distanceMatrix[i] = new double[f.size()];
        for(int j=0;j<f.size();j++) { 
            distanceMatrix[i][j] = 0;
        }
    }

    
    for(int i=0;i<f.size();i++) {
        double nom = 0;
        double denoma = 0;
        double denomb = 0;
        for(int j=i+1;j<f.size();j++) {
            for(int x=0;x<256;x++) {
                nom += rel_freq_avg[i][x]*rel_freq_avg[j][x];
                denoma += rel_freq_avg[i][x]*rel_freq_avg[i][x]; 
                denomb += rel_freq_avg[j][x]*rel_freq_avg[j][x]; 
            }
            
            distanceMatrix[i][j] = (1 - nom/(sqrt(denoma*denomb)))/2;
            distanceMatrix[j][i] = (1 - nom/(sqrt(denoma*denomb)))/2;
        }
    }

    mid.close();
    start.close();
    end.close();
    
    return distanceMatrix;
}

