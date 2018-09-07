#ifndef SET_H
#define SET_H

#include <vector>
#include <algorithm>

class Graph
{
public:
    int V;    // No. of vertices
    std::vector<int> *adj; // adj list for each vertice

    void DFSUtil(int v, bool visited[], int j, std::vector<int>* components); // Depth first search

    Graph(int V);   // Constructor
    ~Graph();
    void addEdge(int v, int w);
    std::vector<int>* connectedComponents(int V);// return connected part as vector of vectors
};




#endif
