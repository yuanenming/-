#ifndef GRAPHCLASS_H_
#define GRAPHCLASS_H_

/*
 * Graph Class Header File
 * We use it for constructing the Phylogenetic Tree in Netwick Tree Format.
 * Contains also Edge and Vertex Classes which Graph Class uses.
 * Graph Class:
 *        @attr: 
 *             vertices: The set of vertices for the Graph
 *            edges: The set of edges for the Graph
 *        @methods:
 *            findVertex: Returns Vertex with a specific name
 *            addEdge: Adds an edge to the Graph
 *            addVertex: Adds an vertex to the Graph
 *            destroyEdge: Destroys an edge in the Graph.
 *            getNetwickTreeFormat: Returns Phylogenetic Tree in Netwick Tree Format
 */

class Edge
{
    private:
        int origin;
        int destination;
        double distance;
        
    public:
        Edge(int org, int dest, double dist);
        int getOrigin();
        int getDestination();
        double getDistance();
    
};

class Vertex
{
    private:
        std::string name;
        int id;
    public:
        Vertex(std::string name, int id);
        std::string getName();
        int getId();
};


class Graph
{
    private:
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
    public:
        Graph() {};    
        Vertex findVertex(std::string name);
        void addVertex(std::string name);
        void addEdge(std::string name1, std::string name2, double distance);
        bool destroyEdge(std::string name1, std::string name2);
        
        void print();
        std::string getNetwickTreeFormat();
};

#endif 





