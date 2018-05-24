// Copyright (C) 2012-2015 Czech Technical University in Prague - All Rights Reserved

#ifndef GRIDGRAPH_3D_6C_H_
#define GRIDGRAPH_3D_6C_H_

#include <cstdlib>
#include <cstring>
#include <climits>
#include <new>

#ifdef __GNUC__
  #include <stdint.h>
#endif

template
<
  typename type_tcap, // Type used to represent capacities of edges between nodes and terminals.
  typename type_ncap, // Type used to represent capacities of edges between nodes and their neighbors.
  typename type_flow  // Type used to represent the total flow.
>
class GridGraph_3D_6C
{
public:
  GridGraph_3D_6C(int width,int height,int depth);
  ~GridGraph_3D_6C();

  // Returns the index of a node at grid coordinates [x,y,z].
  // The index is used to set the capacities of node's outgoing edges
  // using the set_neighbor_cap function, to set the capacities of edges
  // between node and source/sink terminals using the set_terminal_cap
  // function, and to retrieve the segment to which the node belongs
  // after the maxflow computation using the get_segment function.
  inline int node_id(int x,int y,int z) const;

  // Sets the capacities of the source-->node and node-->sink edges.
  // This function can be called only once per node.
  inline void set_terminal_cap(int node_id,type_tcap cap_source,type_tcap cap_sink);

  // Sets the capacity of the edge between node and its neighbor at [offset_x,offset_y,offset_z].
  // For example, to set the capacity of the edge from node [x,y,z] to node [x+1,y,z], call:
  // graph->set_neighbor_cap(graph->node_id(x,y,z),+1,0,0,edge_capacity);
  inline void set_neighbor_cap(int node_id,int offset_x,int offset_y,int offset_z,type_ncap cap);

  // Alternative way to set the edge capacities is to use the set_caps function which
  // sets capacities of all edges at once based on values from input arrays.
  // Each array has width*height*depth elements, where each element corresponds to one node.
  // For example, cap_lee[x+y*width+z*width*height] is the capacity of the outgoing
  // edge from node [x,y,z] to node [x-1,y,z], and cap_sink[x+y*width+z*width*height]
  // is capacity of edge from node [x,y,z] to sink.
  template<typename type_arg_tcap,typename type_arg_ncap>
  inline void set_caps(const type_arg_tcap* cap_source,
                       const type_arg_tcap* cap_sink,
                       const type_arg_ncap* cap_lee,        // [-1, 0, 0]
                       const type_arg_ncap* cap_gee,        // [+1, 0, 0]
                       const type_arg_ncap* cap_ele,        // [ 0,-1, 0]
                       const type_arg_ncap* cap_ege,        // [ 0,+1, 0]
                       const type_arg_ncap* cap_eel,        // [ 0, 0,-1]
                       const type_arg_ncap* cap_eeg);       // [ 0, 0,+1]

  // Computes the maxflow.
  void compute_maxflow();

  // After the maxflow is computed, this function returns the segment
  // to which the node with index node_id belongs.
  // When node belongs to the source segment, the return value is 0.
  // When node belongs to the sink segment, the return value is 1.
  inline int get_segment(int node_id) const;

  // After the maxflow is computed, this function returns the value of maximum flow.
  inline type_flow get_flow() const;

  // Returns true when memory allocation failed inside the constructor.
  // This can be used to detect allocation failures when exceptions are disabled
  // with the GRIDCUT_NO_EXCEPTIONS macro.
  bool bad_alloc() const;

private:

  template<typename T> struct FQueue
  {
    FQueue()
    {
      buffer = 0;
    }

    void clear()
    {
      ifront = buffer;
      iback = buffer;
    }

    bool empty()
    {
      return ifront==iback;
    }

    void push_back(T item)
    {
      *iback = item;
      iback++;
    }

    T pop_back()
    {
      iback--;
      return *iback;
    }

    T pop_front()
    {
      ifront++;
      return *(ifront-1);
    }

    T front()
    {
      return *ifront;
    }

    T* buffer;
    T* ifront;
    T* iback;
   };

  static const int LABEL_S = 1;
  static const int LABEL_T = 2;
  static const int LABEL_F = 0;

  enum Parent { LEE=0,
                GEE,
                ELE,
                EGE,
                EEL,
                EEG,
                NONE,
                TERMINAL };


  static const unsigned char SISTER_TABLE[6];

  static const int INF_D = INT_MAX;

  union Label_Sat
  {
    struct
    {
    unsigned char label      : 2;
    unsigned char nonsat_lee : 1;
    unsigned char nonsat_gee : 1;
    unsigned char nonsat_ele : 1;
    unsigned char nonsat_ege : 1;
    unsigned char nonsat_eel : 1;
    unsigned char nonsat_eeg : 1;
    };
    unsigned char bits;
  };

  Label_Sat* label_sat;

  unsigned char* parent;

  int* parent_id;

  int* dist;

  type_ncap* rc[6];

  type_tcap* rc_st;

  int* timestamp;

  int* QN;
  int QF;
  int QB;

  inline bool Q_EMPTY()           { return QF==1; }
  inline int  Q_FRONT()           { return QF;    }
  inline void Q_PUSH_BACK1(int v) { if (!QN[v]) { QN[QB]=v; QN[v]=1; QB=v; } }
  inline void Q_PUSH_BACK2(int v) { if (!QN[v]) { QN[QB]=v; QN[v]=1; QB=v; } }
  inline void Q_POP_FRONT()       { const int tmp = QN[QF]; QN[QF]=0; QF=tmp; }

  FQueue<int> orphans;
  FQueue<int> orphans2;

  FQueue<int> free_nodes;

  int TIME;

  type_flow MAXFLOW;

  const int ow;
  const int oh;
  const int od;

  const int W;
  const int H;
  const int D;

  const int WB;
  const int WBHB;

  const int YOFS;
  const int ZOFS;

  int nodeId(unsigned int x,unsigned int y,unsigned int z) const;

  bool grow(int&    vs,
            int&    vt,
            Parent& st,
            const   int YOFS,
            const   int ZOFS);

  type_ncap find_minrf_s(int v,type_ncap minrf) const;

  type_ncap find_minrf_t(int v,type_ncap minrf) const;

  void aug_s(int v,const type_ncap minrf);

  void aug_t(int v,const type_ncap minrf);

  void augment(const int    vs,
               const int    vt,
               const Parent st);

  int find_origin(int v,const int TIME);

  void adopt(const int TIME,const int YOFS,const int ZOFS);

  template<typename T>inline type_ncap mincap(type_ncap a,T b) const
  {
    return a < b ? a : b;
  }

  int next_higher_mul4(int x)
  {
    return ((x-1)/4)*4+4;
  }

  void* align(void* mem,size_t boundary)
  {
    uintptr_t mask = ~(uintptr_t)(boundary - 1);
    void *ptr = (void *)(((uintptr_t)mem+boundary-1) & mask);
    return ptr;
  }

  unsigned char* mem_pool;
  unsigned char* mem_pool_tail;

  void* pool_malloc(size_t size)
  {
    void* ptr = mem_pool_tail;
    mem_pool_tail += size;
    return ptr;
  }

};

#define SISTER(E) (SISTER_TABLE[(E)])

#define LABEL(X)      label_sat[(X)].label
#define RC(E,V)       rc[(E)][(V)]
#define RC_ST(V)      rc_st[(V)]

#define PARENT(X)     parent[(X)]
#define PARENT_ID(X)  parent_id[(X)]

#define TIMESTAMP(V)  timestamp[(V)]
#define DIST(V)       dist[(V)]

#define NONSAT(E,V) ( label_sat[(V)].bits & ((1<<2)<<(E)) )
#define ENSAT(E,V)    label_sat[(V)].bits &= (~((1<<2)<<(E)))
#define DESAT(E,V)    label_sat[(V)].bits |= ((1<<2)<<(E))

#define N_LEE(v) ( ( (  (v)  & 0x00000003) == 0 ) ? (v)-61   : (v)-1  )
#define N_GEE(v) ( ( ((~(v)) & 0x00000003) == 0 ) ? (v)+61   : (v)+1  )
#define N_ELE(v) ( ( (  (v)  & 0x0000000C) == 0 ) ? (v)-YOFS : (v)-4  )
#define N_EGE(v) ( ( ((~(v)) & 0x0000000C) == 0 ) ? (v)+YOFS : (v)+4  )
#define N_EEL(v) ( ( (  (v)  & 0x00000030) == 0 ) ? (v)-ZOFS : (v)-16 )
#define N_EEG(v) ( ( ((~(v)) & 0x00000030) == 0 ) ? (v)+ZOFS : (v)+16 )

template <typename type_tcap,typename type_ncap,typename type_flow>
const unsigned char GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::SISTER_TABLE[6] =
{
  GEE,
  LEE,
  EGE,
  ELE,
  EEG,
  EEL
};

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::nodeId(unsigned int x,unsigned int y,unsigned int z) const
{
  return ( ((x>>2) + ((y>>2)*WB) + ((z>>2)*WBHB)) << 6 ) + ( (x&3) + ((y&3)<<2) + ((z&3)<<4) );
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::grow(int&    vs,
                                                          int&    vt,
                                                          Parent& st,
                                                          const   int YOFS,
                                                          const   int ZOFS)
{

   while(!Q_EMPTY())
  {
    const int P = Q_FRONT();

    const int tree_p = LABEL(P);

    const int N_ID[6] = { N_LEE(P),
                          N_GEE(P),
                          N_ELE(P),
                          N_EGE(P),
                          N_EEL(P),
                          N_EEG(P) };

    if (tree_p==LABEL_S)
    {
      for(int n=0;n<6;n++)
      {
        const int N = N_ID[n];
        if (NONSAT(n,P))
        {
          if (LABEL(N)==LABEL_T) { vs=P; vt=N; st=(Parent)SISTER(n); return true; }
          if (LABEL(N)==LABEL_F) { LABEL(N) = LABEL_S; Q_PUSH_BACK2(N); PARENT(N) = SISTER(n); PARENT_ID(N) = P; }
        }
      }
    }
    else if (tree_p==LABEL_T)
    {
      for(int n=0;n<6;n++)
      {
        const int N = N_ID[n];
        if (NONSAT(SISTER(n),N))
        {
          if (LABEL(N)==LABEL_S ) { vt=P; vs=N; st=(Parent)n; return true; }
          if (LABEL(N)==LABEL_F ) { LABEL(N) = LABEL_T; Q_PUSH_BACK2(N); PARENT(N) = SISTER(n); PARENT_ID(N) = P; }
        }
      }
    }

      Q_POP_FRONT();
  }

  return false;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
type_ncap GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::find_minrf_s(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,RC( SISTER(PARENT(v)) , PARENT_ID(v)));
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,RC_ST(v));

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
type_ncap GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::find_minrf_t(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,RC( PARENT(v) ,v));
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,RC_ST(v));

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::aug_s(int v,const type_ncap minrf)
{
  while(PARENT(v) != TERMINAL)
  {
    RC( SISTER(PARENT(v)) , PARENT_ID(v) ) -= minrf;
    RC( PARENT(v) ,v) += minrf;

    DESAT(PARENT(v),v);

    if(! (RC( SISTER(PARENT(v)) , PARENT_ID(v) )) )
    {
      ENSAT( SISTER(PARENT(v)), PARENT_ID(v) );
      PARENT(v) = NONE;
      orphans.push_back(v);
    }

    v = PARENT_ID(v);
  }

  RC_ST(v) -= minrf;
  if (!RC_ST(v))
  {
    PARENT(v) = NONE;
    orphans.push_back(v);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::aug_t(int v,const type_ncap minrf)
{
  while(PARENT(v) != TERMINAL)
  {
    RC( PARENT(v) ,v) -= minrf;
    RC( SISTER( PARENT(v) ) ,PARENT_ID(v)) += minrf;

    DESAT( SISTER( PARENT(v) ), PARENT_ID(v) );

    if (! (RC(PARENT(v) ,v)))
    {
      ENSAT(PARENT(v),v);

      PARENT(v) = NONE;
      orphans.push_back(v);
    }

    v = PARENT_ID(v);
  }

  RC_ST(v) -= minrf;
  if (!RC_ST(v))
  {
    PARENT(v) = NONE;
    orphans.push_back(v);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::augment(const int    vs,
                                                             const int    vt,
                                                             const Parent st)
{
  type_ncap minrf = RC( SISTER(st) ,vs);

  minrf=find_minrf_s(vs,minrf);
  minrf=find_minrf_t(vt,minrf);

  RC( SISTER(st) ,vs) -= minrf;
  RC( st ,vt) += minrf;

  DESAT(st,vt);

  if (!(RC( SISTER(st) ,vs)))
  {
    ENSAT(SISTER(st),vs);
  }

  aug_s(vs,minrf);
  aug_t(vt,minrf);

  MAXFLOW += minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
int GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::find_origin(int v,const int TIME)
{
  const int start_v = v;
  int d2;
  int d = 1;

  while(1)
  {
    if (TIMESTAMP(v)==TIME)    { goto L1; }
    if (PARENT(v) == NONE)     { return INF_D; }
    if (PARENT(v) == TERMINAL) { goto L2; }

    d++;
    v = PARENT_ID(v);
  }

  L1:
    d = DIST(v)+d;
    v = start_v;
    d2 = d;

    while(TIMESTAMP(v)!=TIME)
    {
      DIST(v) = d;
      d--;
      TIMESTAMP(v)=TIME;
      v = PARENT_ID(v);
    }

    return d2;

  L2:
    v = start_v;
    d2 = d;

    while(PARENT(v)!=TERMINAL)
    {
      DIST(v) = d;
      d--;
      TIMESTAMP(v)=TIME;
      v = PARENT_ID(v);
    }

    TIMESTAMP(v)=TIME;
    DIST(v) = 1;

    return d2;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::adopt(const int TIME,const int YOFS,const int ZOFS)
{
  orphans2.clear();
  free_nodes.clear();

  while(!(orphans.empty()&&orphans2.empty()))
  {
    const int V = (!orphans2.empty()) ? orphans2.pop_front() : orphans.pop_back();

    const int N_ID[6] = { N_LEE(V),
                          N_GEE(V),
                          N_ELE(V),
                          N_EGE(V),
                          N_EEL(V),
                          N_EEG(V) };

    const int tree_p = LABEL(V);

    int min_d = INF_D;
    int best_a = -1;

    if      (tree_p==LABEL_S)
    {
      for(int i=0;i<6;i++)
      {
        const int N = N_ID[i];

        if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S)
        {
          const int d = find_origin(N,TIME);
          if (d<min_d)
          {
            best_a = i;
            min_d=d;
          }
        }
      }
    }
    else if (tree_p==LABEL_T)
    {
      for(int i=0;i<6;i++)
      {
        const int N = N_ID[i];

        if (NONSAT(i,V) && LABEL(N)==LABEL_T)
        {
          const int d = find_origin(N,TIME);
          if (d<min_d)
          {
            best_a = i;
            min_d = d;
          }
        }
      }
    }

    if (best_a!=-1)
    {
      DIST(V) = min_d+1;
      TIMESTAMP(V)=TIME;

      PARENT(V) = best_a;
      PARENT_ID(V) = N_ID[best_a];

      goto next;
    }

    LABEL(V) = LABEL_F;
    free_nodes.push_back(V);

    for(int i=0;i<6;i++)
    {
      const int N = N_ID[i];
      if (LABEL(N)==tree_p && PARENT(N)==SISTER(i))  { PARENT(N) = NONE; orphans2.push_back(N); }
    }

    next:
      ;
  }

  while(!free_nodes.empty())
  {
    const int V = free_nodes.pop_front();

    const int N_ID[6] = { N_LEE(V),
                          N_GEE(V),
                          N_ELE(V),
                          N_EGE(V),
                          N_EEL(V),
                          N_EEG(V) };

    for(int i=0;i<6;i++)
    {
      const int N = N_ID[i];
      if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S) {  Q_PUSH_BACK2(N); }
      if (NONSAT(i,V) && LABEL(N)==LABEL_T) {  Q_PUSH_BACK2(N); }
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::compute_maxflow()
{
  while(!free_nodes.empty())
  {
    const int v = free_nodes.pop_front();

    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[6] = { N_LEE(v),
                          N_GEE(v),
                          N_ELE(v),
                          N_EGE(v),
                          N_EEL(v),
                          N_EEG(v) };

    for(int arc=0;arc<6;arc++)
    {
      const int N = N_ID[arc];

      if (lv==LABEL_S)
      {
        if (NONSAT(arc,v) && LABEL(N)!=lv) {  Q_PUSH_BACK1(v); goto next_node; }
      }
      else
      {
        if (NONSAT(SISTER(arc),N) && LABEL(N)==LABEL_F) { Q_PUSH_BACK1(v); goto next_node; }
      }
    }

    next_node:
      ;
  }

  QF = QN[0];

  while(1)
  {
    int vs;
    int vt;
    Parent st;

    const bool path_found = grow(vs,vt,st,YOFS,ZOFS);

    if (!path_found) break;
    TIME++;

    augment(vs,vt,st);
    adopt(TIME,YOFS,ZOFS);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::node_id(int x,int y,int z) const
{
  return nodeId(x+1,y+1,z+1);
}


template <typename type_tcap,typename type_ncap,typename type_flow>
inline void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::set_terminal_cap(int v,type_tcap cap_s,type_tcap cap_t)
{
  if (cap_s > 0 && cap_t > 0)
  {
    if      (cap_s > cap_t)
    {
      RC_ST(v) = cap_s-cap_t;

      MAXFLOW += cap_t;

      LABEL(v) = LABEL_S;

      PARENT(v) = TERMINAL;

      DIST(v) = 1;

      free_nodes.push_back(v);
    }
    else if (cap_s < cap_t)
    {
      RC_ST(v) = cap_t-cap_s;

      MAXFLOW += cap_s;

      LABEL(v) = LABEL_T;

      PARENT(v) = TERMINAL;

      DIST(v) = 1;

      free_nodes.push_back(v);
    }
    else
    {
      RC_ST(v) = 0;

      MAXFLOW += cap_s;

      PARENT(v) = NONE;
    }
  }
  else
  {
    if (cap_s > 0)
    {
      LABEL(v) = LABEL_S;

      RC_ST(v) = cap_s;

      PARENT(v) = TERMINAL;

      DIST(v) = 1;

      free_nodes.push_back(v);
    }
    else if (cap_t > 0)
    {
      LABEL(v) = LABEL_T;

      RC_ST(v) = cap_t;

      PARENT(v) = TERMINAL;

      DIST(v) = 1;

      free_nodes.push_back(v);
    }
  }
}


template<typename type_tcap,typename type_ncap,typename type_flow>
inline void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::set_neighbor_cap(int node_id,int offset_x,int offset_y,int offset_z,type_ncap cap)
{
  static const unsigned char offsets2arc[3][3][3] = {  { {255,255,255},
                                                         {255,EEL,255},
                                                         {255,255,255} },

                                                       { {255,ELE,255},
                                                         {LEE,255,GEE},
                                                         {255,EGE,255} },

                                                       { {255,255,255},
                                                         {255,EEG,255},
                                                         {255,255,255} }  };

  const int arc = offsets2arc[offset_z+1][offset_y+1][offset_x+1];
  rc[arc][node_id] = cap;

  if (cap!=0) { DESAT(arc,node_id); }
}

template<typename type_tcap,typename type_ncap,typename type_flow> template<typename type_arg_tcap,typename type_arg_ncap>
inline void GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::set_caps(const type_arg_tcap* cap_s,
                                                                     const type_arg_tcap* cap_t,
                                                                     const type_arg_ncap* cap_lee,
                                                                     const type_arg_ncap* cap_gee,
                                                                     const type_arg_ncap* cap_ele,
                                                                     const type_arg_ncap* cap_ege,
                                                                     const type_arg_ncap* cap_eel,
                                                                     const type_arg_ncap* cap_eeg)
{
  for(int xyz=0,z=0;z<od;z++)
  for(int       y=0;y<oh;y++)
  for(int       x=0;x<ow;x++,xyz++)
  {
    const int v = nodeId(x+1,y+1,z+1);

    if (cap_s[xyz] > 0 && cap_t[xyz] > 0)
    {
      if      (cap_s[xyz] > cap_t[xyz])
      {
        RC_ST(v) = cap_s[xyz]-cap_t[xyz];

        MAXFLOW += cap_t[xyz];

        LABEL(v) = LABEL_S;

        PARENT(v) = TERMINAL;

        DIST(v) = 1;
      }
      else if (cap_s[xyz] < cap_t[xyz])
      {
        RC_ST(v) = cap_t[xyz]-cap_s[xyz];

        MAXFLOW += cap_s[xyz];

        LABEL(v) = LABEL_T;

        PARENT(v) = TERMINAL;

        DIST(v) = 1;
      }
      else
      {
        RC_ST(v) = 0;

        MAXFLOW += cap_s[xyz];

        PARENT(v) = NONE;
      }
    }
    else
    {
      if (cap_s[xyz] > 0)
      {
        LABEL(v) = LABEL_S;

        RC_ST(v) = cap_s[xyz];

        PARENT(v) = TERMINAL;

        DIST(v) = 1;
      }
      else if (cap_t[xyz] > 0)
      {
        LABEL(v) = LABEL_T;

        RC_ST(v) = cap_t[xyz];

        PARENT(v) = TERMINAL;

        DIST(v) = 1;
      }
    }

    RC(LEE,v) = cap_lee[xyz]; if (cap_lee[xyz]!=0) { DESAT(LEE,v); }
    RC(GEE,v) = cap_gee[xyz]; if (cap_gee[xyz]!=0) { DESAT(GEE,v); }
    RC(ELE,v) = cap_ele[xyz]; if (cap_ele[xyz]!=0) { DESAT(ELE,v); }
    RC(EGE,v) = cap_ege[xyz]; if (cap_ege[xyz]!=0) { DESAT(EGE,v); }
    RC(EEL,v) = cap_eel[xyz]; if (cap_eel[xyz]!=0) { DESAT(EEL,v); }
    RC(EEG,v) = cap_eeg[xyz]; if (cap_eeg[xyz]!=0) { DESAT(EEG,v); }
  }

  for(int v=0;v<W*H*D;v++)
  {
    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[6] = { N_LEE(v),
                          N_GEE(v),
                          N_ELE(v),
                          N_EGE(v),
                          N_EEL(v),
                          N_EEG(v) };

    for(int arc=0;arc<6;arc++)
    {
      const int N = N_ID[arc];

      if (lv==LABEL_S)
      {
        if (NONSAT(arc,v) && LABEL(N)!=lv) {  Q_PUSH_BACK1(v); goto next_node; }
      }
      else
      {
        if (NONSAT(SISTER(arc),N) && LABEL(N)==LABEL_F) { Q_PUSH_BACK1(v); goto next_node; }
      }
    }

    next_node:
      ;
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::get_segment(int v) const
{
  static const int map_label[3] = { 0, 0, 1 };

  return map_label[LABEL(v)];
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline type_flow GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::get_flow() const
{
  return MAXFLOW;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::bad_alloc() const
{
  return (mem_pool==NULL);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::GridGraph_3D_6C(int w,int h,int d) :
  ow(w),
  oh(h),
  od(d),
  W(next_higher_mul4(w+2)),
  H(next_higher_mul4(h+2)),
  D(next_higher_mul4(d+2)),
  WB(W/4),
  WBHB((W/4)*(H/4)),
  YOFS(WB*64 - 3*4),
  ZOFS(WBHB*64 - 3*(4*4)),
  mem_pool(NULL)
{
  size_t mem_pool_size = (W*H*D*sizeof(Label_Sat)+64+             // label_sat
                          W*H*D*sizeof(unsigned char)+64+         // parent
                          W*H*D*sizeof(int)+64+                   // parent_id
                          6*(W*H*D*sizeof(type_ncap)+64)+         // rc[6]
                          W*H*D*sizeof(type_tcap)+64+             // rc_st
                          W*H*D*sizeof(int)+64+                   // dist
                          W*H*D*sizeof(int)+64+                   // timestamp
                          W*H*D*sizeof(int)+64+                   // orphans
                          W*H*D*sizeof(int)+64+                   // orphans2
                          W*H*D*sizeof(int)+64+                   // free_nodes
                          W*H*D*sizeof(int)+64+                   // QN
                          0);

  mem_pool = (unsigned char*)calloc(mem_pool_size,1);

  if (mem_pool==NULL)
  {
#ifdef GRIDCUT_NO_EXCEPTIONS
    return;
#else    
    throw std::bad_alloc();
#endif    
  }

  mem_pool_tail = mem_pool;

  label_sat = (Label_Sat*)align(pool_malloc(W*H*D*sizeof(Label_Sat)+64),64);

  parent = (unsigned char*)align(pool_malloc(W*H*D*sizeof(unsigned char)+64),64);

  parent_id = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  for(int i=0;i<6;i++)
  {
    rc[i] = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  }

  rc_st = (type_tcap*)align(pool_malloc(W*H*D*sizeof(type_tcap)+64),64);


  dist = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  timestamp = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  memset(parent,NONE,W*H*D);

  QN = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  QF = 0;
  QB = 0;
  QN[0] = 1;

  orphans.buffer = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  orphans.clear();

  orphans2.buffer = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  orphans2.clear();

  free_nodes.buffer = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  free_nodes.clear();

  MAXFLOW = 0;

  TIME = 0;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_3D_6C<type_tcap,type_ncap,type_flow>::~GridGraph_3D_6C()
{
  if (mem_pool!=NULL) { free(mem_pool); }
}

#undef SISTER

#undef LABEL

#undef RC
#undef RC_ST

#undef PARENT
#undef PARENT_ID

#undef TIMESTAMP
#undef DIST

#undef NONSAT
#undef ENSAT
#undef DESAT

#undef N_LEE
#undef N_GEE
#undef N_ELE
#undef N_EGE
#undef N_EEL
#undef N_EEG

#endif
