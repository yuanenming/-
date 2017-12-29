#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <math.h> 
#include <time.h>
#include <iomanip>

#include "GraphClass.h"

/*EDGE CLASS METHODS*/
Edge::Edge(int org, int dest, double dist)
{
    origin = org;
    destination = dest;
    distance = dist;
}    
int Edge::getOrigin() {return origin;}
int Edge::getDestination() {return destination;}
double Edge::getDistance() {return distance;}
    

/*VERTEX CLASS METHODS*/    
Vertex::Vertex(std::string n, int i )
{
    name = n;
    id = i;
}
std::string Vertex::getName() {return name;} 
int Vertex::getId(){return id;}

/*GRAPH CLASS METHODS*/
void Graph::addVertex(std::string name){
    int id;
    if ( !vertices.empty() ) { id = vertices.size(); }
    else { id = 0;}
    Vertex newVertex(name, id);
    vertices.push_back(newVertex);
}
Vertex Graph::findVertex(std::string name){
    
    for( std::vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); ++i ){
        if ( name.compare( (*i).getName() ) == 0 ){
            return (*i);
        }
    }
}
void Graph::addEdge(std::string name1, std::string name2, double dist){
    
    Vertex node1 = findVertex(name1);
    Vertex node2 = findVertex(name2);
    
    Edge edge(node1.getId(), node2.getId(), dist);
    edges.push_back(edge);
}

bool Graph::destroyEdge(std::string name1, std::string name2){
        
    Vertex node1 = findVertex(name1);
    Vertex node2 = findVertex(name2);
    
    for( std::vector<Edge>::iterator i = edges.begin(); i != edges.end(); ++i ){
        if ( ((*i).getOrigin() == node1.getId() && (*i).getDestination() == node2.getId()) ||
            ((*i).getOrigin() == node2.getId() && (*i).getDestination() == node1.getId()) ){
            //delete the edge
            edges.erase(i);
            return true;
        }
    }
    return false;
}
void Graph::print(){
    std::ofstream of("Graph.txt");
    for( std::vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); ++i ){
        of << (*i).getName() << std::endl;
    }
    for( std::vector<Edge>::iterator i = edges.begin(); i != edges.end(); ++i ){
        of << (*i).getOrigin() << ", "<< (*i).getDestination() <<  " => "<< (*i).getDistance() << std::endl;
    }
    of.close();
}

std::string Graph::getNetwickTreeFormat(){
    return vertices.back().getName();
}

