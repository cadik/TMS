#ifndef CHERYL11_GRAPH_H
#define CHERYL11_GRAPH_H

#include <vector>

namespace Cheryl11 {

    class Graph
    {
    public:
        struct Edge {
            Edge(int c_0, int c_1)
            {
                c0 = c_0;
                c1 = c_1;
            }
            
            int c0;
            int c1;
        };
        
        Graph();
        ~Graph();
        
        void addEdge(int c0, int c1);
        
        int getEdgesCount() { return edges.size(); }
        const std::vector<Edge>& getEdges() const { return edges; }
    protected:
        std::vector<Edge> edges;
        
        bool checkEdgeExistance(int c0, int c1);
    };

} // namespace Cheryl11

#endif /* CHERYL11_GRAPH_H */