#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "wordContext.h"

using namespace std;


double** wordContext(map<int, map<int, map<string,int> > > f) {
    

    vector<vector<string> > contexts3;
    vector<vector<string> > contexts4;
    vector<vector<string> > contexts5;
    vector<vector<string> > contexts6;
    vector<vector<string> > contexts7;
    vector<vector<string> > contexts8;
    vector<vector<string> > contexts9;
    vector<vector<string> > contexts10;

    map<string, int> tmpM;
    typedef std::map<string, int>::iterator it_type;
    for (int k=3;k<11;k++) { //for every k-mer
        for(int s=0;s<f.size();s++) { //for every sequence (144)
            tmpM = f[s][k];
            vector<string> tmpV;
            for(it_type iterator = tmpM.begin(); iterator != tmpM.end(); iterator++) { // here we input k-tuples for each sequence
                tmpV.reserve(1);
                tmpV.push_back(iterator->first);
            }
            if(k==3) contexts3.push_back(tmpV);
            if(k==4) contexts4.push_back(tmpV);
            if(k==5) contexts5.push_back(tmpV);
            if(k==6) contexts6.push_back(tmpV);
            if(k==7) contexts7.push_back(tmpV);
            if(k==8) contexts8.push_back(tmpV);
            if(k==9) contexts9.push_back(tmpV);
            if(k==10) contexts10.push_back(tmpV);
        }
    }//here ends the proceduce for each k-mer
    double** distanceMatrix = new double*[f.size()]; 

    for(int i=0;i<f.size();i++) { 
        distanceMatrix[i] = new double[f.size()];
        for(int j=0;j<f.size();j++) { 
            distanceMatrix[i][j] = 0;
        }
    }
        
    double distance;
    int sumI = 0;
    int R = 0;
    for (int i = 0;i<f.size();i++) { //example for 5mers
        for(int j=i+1; j<f.size();j++) {
            sumI = 0;
            R = 0;
            for(int y=0;y<contexts5[i].size();y++) {
                for(int z=y;z<contexts5[j].size();z++) {
                    if (contexts5[i][y]!=contexts5[j][z]) sumI++;
                    if ((contexts5[i][y][0] == contexts5[j][z][0])&&(contexts5[i][y][1] == contexts5[j][z][1])&&(contexts5[i][y][3] == contexts5[j][z][3])&&(contexts5[i][y][4] == contexts5[j][z][4])) R++;
                }
            }
            distanceMatrix[i][j] = double(sumI)/double(R);
            distanceMatrix[j][i] = double(sumI)/double(R);
        }
    }

    return distanceMatrix;
}

