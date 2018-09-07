#include "set.h"

using namespace std;

vector<int>* Graph::connectedComponents(int V)
{

    vector<int> *components=new vector<int>[V];
    int j=0;

    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for(int v = 0; v < V; v++)
        visited[v] = false;

    for (int v=0; v<V; v++)
    {
        if (visited[v] == false)
        {

            DFSUtil(v, visited, j, components);

            j++;
        }
    }
    delete[] visited;

    for(int cc=0;cc<V;cc++)
    {
        sort(components[cc].begin(), components[cc].end());
    }

    return components;
}


void Graph::DFSUtil(int v, bool visited[], int j, vector<int>* components)
{

    visited[v] = true;
    components[j].push_back(v);

    vector<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            DFSUtil(*i, visited,j,components);
}

Graph::Graph(int V)
{
    this->V = V;
    adj = new vector<int>[V];
}

Graph::~Graph()
{
    delete []adj;
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

