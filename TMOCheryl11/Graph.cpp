#include "Graph.h"

using namespace Cheryl11;

Graph::Graph()
{
}

Graph::~Graph()
{
}

void Graph::addEdge(int c0, int c1, float lenght, float colorL2)
{
    if (!checkEdgeExistance(c0, c1))
    {
        edges.push_back(Edge(c0, c1, lenght, colorL2));
    }
}

void Graph::setPsi(int i, double psi)
{
    edges.at(i).psi = psi;
}

void Graph::setScaleFactor(double scale_factor)
{
    scaleFactor = scale_factor;
}

bool Graph::checkEdgeExistance(int c0, int c1)
{
    if (c0 == c1)
    {
        return true;
    }

    for (const Edge &edge : edges)
    {
        if (edge.c0 == c0)
        {
            if (edge.c1 == c1)
            {
                return true;
            }
        }
        else if (edge.c1 == c0)
        {
            if (edge.c0 == c1)
            {
                return true;
            }
        }
    }
    return false;
}