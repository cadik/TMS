#include "Graph.h"

using namespace Cheryl11;

Graph::Graph()
{ }

Graph::~Graph()
{ }

void Graph::addEdge(int c0, int c1)
{
    if ( ! checkEdgeExistance(c0, c1))
    {
        edges.push_back(Edge(c0, c1));
    }
}

bool Graph::checkEdgeExistance(int c0, int c1)
{
    for (const Edge& edge : edges)
    {
        if (edge.c0 == c0) {
            if (edge.c1 == c1) {
                return true;
            }
        }
        else if (edge.c1 == c0) {
            if (edge.c0 == c1) {
                return true;
            }
        }
    }
    return false;
}