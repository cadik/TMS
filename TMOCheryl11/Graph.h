#ifndef CHERYL11_GRAPH_H
#define CHERYL11_GRAPH_H

#include <vector>

namespace Cheryl11
{

    class Graph
    {
    public:
        struct Edge
        {
            Edge(int c_0, int c_1, float in_lenght = 0.0, float colorL_2 = 0.0)
            {
                c0 = c_0;
                c1 = c_1;
                lenght = in_lenght;
                colorL2 = colorL_2;
                psi = 0.0;
            }

            int c0;
            int c1;
            float lenght;
            float colorL2;
            double psi;
        };

        Graph();
        ~Graph();

        void addEdge(int c0, int c1, float lenght, float colorL2);
        void setPsi(int i, double psi);
        void setScaleFactor(double scale_factor);

        int getEdgesCount() { return edges.size(); }
        const std::vector<Edge> &getEdges() const { return edges; }
        double getScaleFactor() { return scaleFactor; }

    protected:
        std::vector<Edge> edges;
        double scaleFactor = 0.0;

        bool checkEdgeExistance(int c0, int c1); // also disable edge with themself
    };

} // namespace Cheryl11

#endif /* CHERYL11_GRAPH_H */