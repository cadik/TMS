// Copyright (C) 2012-2015 Czech Technical University in Prague - All Rights Reserved

#ifndef GRIDGRAPH_2D_8C_H_
#define GRIDGRAPH_2D_8C_H_

#include <cstdlib>
#include <cstring>
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
class GridGraph_2D_8C
{
public:
  GridGraph_2D_8C(int width,int height);
  ~GridGraph_2D_8C();

  // Returns the index of a node at grid coordinates [x,y].
  // The index is used to set the capacities of node's outgoing edges
  // using the set_neighbor_cap function, to set the capacities of edges
  // between node and source/sink terminals using the set_terminal_cap
  // function, and to retrieve the segment to which the node belongs
  // after the maxflow computation using the get_segment function.
  inline int node_id(int x,int y) const;

  // Sets the capacities of the source-->node and node-->sink edges.
  // This function can be called only once per node.
  inline void set_terminal_cap(int node_id,type_tcap cap_source,type_tcap cap_sink);

  // Sets the capacity of the edge between node and its neighbor at [offset_x,offset_y].
  // For example, to set the capacity of the edge from node [x,y] to node [x+1,y], call:
  // graph->set_neighbor_cap(graph->node_id(x,y),+1,0,edge_capacity);
  inline void set_neighbor_cap(int node_id,int offset_x,int offset_y,type_ncap cap);

  // Alternative way to set the edge capacities is to use the set_caps function which
  // sets capacities of all edges at once based on values from input arrays.
  // Each array has width*height elements, where each element corresponds to one node.
  // For example, cap_le[x+y*width] is the capacity of the outgoing edge from node [x,y]
  // to node [x-1,y], and cap_sink[x+y*width] is capacity of edge from node [x,y] to sink.
  template<typename type_arg_tcap,typename type_arg_ncap>
  inline void set_caps(const type_arg_tcap* cap_source,
                       const type_arg_tcap* cap_sink,
                       const type_arg_ncap* cap_le,         // [-1, 0]
                       const type_arg_ncap* cap_ge,         // [+1, 0]
                       const type_arg_ncap* cap_el,         // [ 0,-1]
                       const type_arg_ncap* cap_eg,         // [ 0,+1]
                       const type_arg_ncap* cap_ll,         // [-1,-1]
                       const type_arg_ncap* cap_gl,         // [+1,-1]
                       const type_arg_ncap* cap_lg,         // [-1,+1]
                       const type_arg_ncap* cap_gg          // [+1,+1]
                       );


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

  enum Parent { LE=0,GE,EL,EG,LL,GL,LG,GG,NONE,TERMINAL };

  static const unsigned char SISTER_TABLE[8];

  unsigned char* label;

  unsigned char* parent;

  int* parent_id;

  type_ncap* rc[8];
  type_ncap* rc_sister[8];

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

  const int W;
  const int H;

  const int WB;

  const int YOFS;

  int nodeId(unsigned int x,unsigned int y) const;

  bool grow(int&    vs,
            int&    vt,
            Parent& st,
            const   int YOFS);

  type_ncap find_minrf_s(int v,type_ncap minrf) const;

  type_ncap find_minrf_t(int v,type_ncap minrf) const;

  void aug_s(int v,const type_ncap minrf);

  void aug_t(int v,const type_ncap minrf);

  void augment(const int    vs,
               const int    vt,
               const Parent st);

  int find_origin(int v,const int TIME);

  void adopt(const int TIME,const int YOFS);

  template<typename T>inline type_ncap mincap(type_ncap a,T b) const
  {
    return a < b ? a : b;
  }

  int next_higher_mul8(int x)
  {
    return ((x-1)/8)*8+8;
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

#define LABEL(X)       label[(X)]
#define RC(E,V)        rc[(E)][(V)]
#define RC_SISTER(E,V) rc_sister[(E)][(V)]
#define RC_ST(V)       rc_st[(V)]

#define PARENT(X)      parent[(X)]
#define PARENT_ID(X)   parent_id[(X)]

#define TIMESTAMP(V)   timestamp[(V)]

#define NONSAT(E,V) ( rc[(E)][(V)] )
#define NONSAT_SISTER(E,V) ( rc_sister[E][V] )
#define ENSAT(E,V)
#define DESAT(E,V)

#define N_LE(v) ( ( (  (v)  & 0x00000007) == 0 ) ? (v)-57   : (v)-1 )
#define N_GE(v) ( ( ((~(v)) & 0x00000007) == 0 ) ? (v)+57   : (v)+1 )
#define N_EL(v) ( ( (  (v)  & 0x00000038) == 0 ) ? (v)-YOFS : (v)-8 )
#define N_EG(v) ( ( ((~(v)) & 0x00000038) == 0 ) ? (v)+YOFS : (v)+8 )
#define N_LL(v) ( (v) + ( ( (  (v)  & 0x00000007) == 0 ) ? -57   : -1  ) + ( ( (  (v)  & 0x00000038) == 0 ) ? -YOFS : -8  ) +  0  )
#define N_GL(v) ( (v) + ( ( ((~(v)) & 0x00000007) == 0 ) ? +57   : +1  ) + ( ( (  (v)  & 0x00000038) == 0 ) ? -YOFS : -8  ) +  0  )
#define N_LG(v) ( (v) + ( ( (  (v)  & 0x00000007) == 0 ) ? -57   : -1  ) + ( ( ((~(v)) & 0x00000038) == 0 ) ? +YOFS : +8  ) +  0  )
#define N_GG(v) ( (v) + ( ( ((~(v)) & 0x00000007) == 0 ) ? +57   : +1  ) + ( ( ((~(v)) & 0x00000038) == 0 ) ? +YOFS : +8  ) +  0  )

template <typename type_tcap,typename type_ncap,typename type_flow>
const unsigned char GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::SISTER_TABLE[8] = {GE,LE,EG,EL,GG,LG,GL,LL};

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::nodeId(unsigned int x,unsigned int y) const
{
  return (((x>>3)+(y>>3)*WB)<<6) + (x&7)+((y&7)<<3);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::grow(int&    vs,
                                                          int&    vt,
                                                          Parent& st,
                                                          const   int YOFS)
{
  while(!Q_EMPTY())
  {
    const int P = Q_FRONT();

    const int tree_p = LABEL(P);

    const int N_ID[8] = { N_LE(P),N_GE(P),N_EL(P),N_EG(P),N_LL(P),N_GL(P),N_LG(P),N_GG(P) };

    if (tree_p==LABEL_S)
    {

      for(int n=0;n<8;n++)
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
      for(int n=0;n<8;n++)
      {
        const int N = N_ID[n];
        if (NONSAT_SISTER(n,N))
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
type_ncap GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::find_minrf_s(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,RC_SISTER( PARENT(v) , PARENT_ID(v) ));
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,RC_ST(v));

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
type_ncap GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::find_minrf_t(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,RC( PARENT(v) , v));
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,RC_ST(v));

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::aug_s(int v,const type_ncap minrf)
{
  while(PARENT(v) != TERMINAL)
  {
    RC_SISTER( PARENT(v) , PARENT_ID(v) ) -= minrf;
    RC( PARENT(v) ,v) += minrf;

    DESAT(PARENT(v),v);

    if(! (RC_SISTER( PARENT(v) , PARENT_ID(v) )) )
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
void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::aug_t(int v,const type_ncap minrf)
{
  while(PARENT(v) != TERMINAL)
  {
    RC( PARENT(v) ,v) -= minrf;
    RC_SISTER(  PARENT(v)  ,PARENT_ID(v)) += minrf;

    DESAT( SISTER( PARENT(v) ), PARENT_ID(v) );

    if (! (RC( PARENT(v) ,v)))
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
void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::augment(const int    vs,
                                                             const int    vt,
                                                             const Parent st)
{
  type_ncap minrf = RC_SISTER( st , vs);

  minrf=find_minrf_s(vs,minrf);
  minrf=find_minrf_t(vt,minrf);

  RC_SISTER( st ,vs) -= minrf;
  RC( st ,vt) += minrf;

  DESAT(st,vt);

  if (!(RC( SISTER(st) , vs)))
  {
    ENSAT(SISTER(st),vs);
  }

  aug_s(vs,minrf);
  aug_t(vt,minrf);

  MAXFLOW += minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
int GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::find_origin(int v,const int TIME)
{
  const int start_v = v;

  while(1)
  {
    if (TIMESTAMP(v)==TIME)    { goto L1; }
    if (PARENT(v) == NONE)     { return NONE; }
    if (PARENT(v) == TERMINAL) { goto L2; }

    v = PARENT_ID(v);
  }

  L1:
    v = start_v;

    while(TIMESTAMP(v)!=TIME)
    {
      TIMESTAMP(v)=TIME;
      v = PARENT_ID(v);
    }

    return TERMINAL;

  L2:
    v = start_v;

    while(PARENT(v)!=TERMINAL)
    {
      TIMESTAMP(v)=TIME;
      v = PARENT_ID(v);
    }

    TIMESTAMP(v)=TIME;

    return TERMINAL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::adopt(const int TIME,const int YOFS)
{
  orphans2.clear();
  free_nodes.clear();

  while(!(orphans.empty()&&orphans2.empty()))
  {
    const int V = (!orphans2.empty()) ? orphans2.pop_front() : orphans.pop_back();

    const int N_ID[8] = { N_LE(V),N_GE(V),N_EL(V),N_EG(V),N_LL(V),N_GL(V),N_LG(V),N_GG(V) };

    const int tree_p = LABEL(V);

    if      (tree_p==LABEL_S)
    {
      for(int i=0;i<8;i++)
      {
        const int N = N_ID[i];

        if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S && find_origin(N,TIME)==TERMINAL)
        {
          TIMESTAMP(V)=TIME; PARENT(V) = i; PARENT_ID(V) = N; goto next;
        }
      }
    }
    else if (tree_p==LABEL_T)
    {
      for(int i=0;i<8;i++)
      {
        const int N = N_ID[i];

        if (NONSAT(i,V) && LABEL(N)==LABEL_T && find_origin(N,TIME)==TERMINAL)
        {
          TIMESTAMP(V)=TIME; PARENT(V) = i; PARENT_ID(V) = N; goto next;
        }
      }
    }

    LABEL(V) = LABEL_F;
    free_nodes.push_back(V);

    for(int i=0;i<8;i++)
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

    const int N_ID[8] = { N_LE(V),N_GE(V),N_EL(V),N_EG(V),N_LL(V),N_GL(V),N_LG(V),N_GG(V) };

    for(int i=0;i<8;i++)
    {
      const int N = N_ID[i];
      if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S) {  Q_PUSH_BACK2(N); }
      if (NONSAT(i,V) && LABEL(N)==LABEL_T) {  Q_PUSH_BACK2(N); }
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::compute_maxflow()
{
  while(!free_nodes.empty())
  {
    const int v = free_nodes.pop_front();

    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[8] = { N_LE(v),N_GE(v),N_EL(v),N_EG(v),N_LL(v),N_GL(v),N_LG(v),N_GG(v) };

    for(int arc=0;arc<8;arc++)
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

    const bool path_found = grow(vs,vt,st,YOFS);

    if (!path_found) break;
    TIME++;

    augment(vs,vt,st);
    adopt(TIME,YOFS);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::node_id(int x,int y) const
{
  return nodeId(x+1,y+1);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::set_terminal_cap(int v,type_tcap cap_s,type_tcap cap_t)
{
  if (cap_s > 0 && cap_t > 0)
  {
    if      (cap_s > cap_t)
    {
      RC_ST(v) = cap_s-cap_t;

      MAXFLOW += cap_t;

      LABEL(v) = LABEL_S;

      PARENT(v) = TERMINAL;

      free_nodes.push_back(v);
    }
    else if (cap_s < cap_t)
    {
      RC_ST(v) = cap_t-cap_s;

      MAXFLOW += cap_s;

      LABEL(v) = LABEL_T;

      PARENT(v) = TERMINAL;

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

      free_nodes.push_back(v);
    }
    else if (cap_t > 0)
    {
      LABEL(v) = LABEL_T;

      RC_ST(v) = cap_t;

      PARENT(v) = TERMINAL;

      free_nodes.push_back(v);
    }
  }
}

template<typename type_tcap,typename type_ncap,typename type_flow>
inline void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::set_neighbor_cap(int node_id,int offset_x,int offset_y,type_ncap cap)
{
  static const unsigned char offsets2arc[3][3] = { { LL,  EL, GL },
                                                   { LE, 255, GE },
                                                   { LG,  EG, GG } };

  const int arc = offsets2arc[offset_y+1][offset_x+1];
  RC(arc,node_id) = cap;
}

template<typename type_tcap,typename type_ncap,typename type_flow> template<typename type_arg_tcap,typename type_arg_ncap>
inline void GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::set_caps(const type_arg_tcap* cap_s,
                                                                     const type_arg_tcap* cap_t,
                                                                     const type_arg_ncap* cap_le,
                                                                     const type_arg_ncap* cap_ge,
                                                                     const type_arg_ncap* cap_el,
                                                                     const type_arg_ncap* cap_eg,
                                                                     const type_arg_ncap* cap_ll,
                                                                     const type_arg_ncap* cap_gl,
                                                                     const type_arg_ncap* cap_lg,
                                                                     const type_arg_ncap* cap_gg)
{
  for(int xy=0,y=0;y<oh;y++)
  for(int      x=0;x<ow;x++,xy++)
  {
    const int v = nodeId(x+1,y+1);

    if (cap_s[xy] > 0 && cap_t[xy] > 0)
    {
      if      (cap_s[xy] > cap_t[xy])
      {
        RC_ST(v) = cap_s[xy]-cap_t[xy];

        MAXFLOW += cap_t[xy];

        LABEL(v) = LABEL_S;

        PARENT(v) = TERMINAL;
      }
      else if (cap_s[xy] < cap_t[xy])
      {
        RC_ST(v) = cap_t[xy]-cap_s[xy];

        MAXFLOW += cap_s[xy];

        LABEL(v) = LABEL_T;

        PARENT(v) = TERMINAL;
      }
      else
      {
        RC_ST(v) = 0;

        MAXFLOW += cap_s[xy];

        PARENT(v) = NONE;
      }
    }
    else
    {
      if (cap_s[xy] > 0)
      {
        LABEL(v) = LABEL_S;

        RC_ST(v) = cap_s[xy];

        PARENT(v) = TERMINAL;
      }
      else if (cap_t[xy] > 0)
      {
        LABEL(v) = LABEL_T;

        RC_ST(v) = cap_t[xy];

        PARENT(v) = TERMINAL;
      }
    }

    RC(LE,v) = cap_le[xy]; if (cap_le[xy]!=0) { DESAT(LE,v); }
    RC(GE,v) = cap_ge[xy]; if (cap_ge[xy]!=0) { DESAT(GE,v); }
    RC(EL,v) = cap_el[xy]; if (cap_el[xy]!=0) { DESAT(EL,v); }
    RC(EG,v) = cap_eg[xy]; if (cap_eg[xy]!=0) { DESAT(EG,v); }
    RC(LL,v) = cap_ll[xy]; if (cap_ll[xy]!=0) { DESAT(LL,v); }
    RC(GL,v) = cap_gl[xy]; if (cap_gl[xy]!=0) { DESAT(GL,v); }
    RC(LG,v) = cap_lg[xy]; if (cap_lg[xy]!=0) { DESAT(LG,v); }
    RC(GG,v) = cap_gg[xy]; if (cap_gg[xy]!=0) { DESAT(GG,v); }
  }

  for(int v=0;v<W*H;v++)
  {
    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[8] = { N_LE(v),N_GE(v),N_EL(v),N_EG(v),N_LL(v),N_GL(v),N_LG(v),N_GG(v) };

    for(int arc=0;arc<8;arc++)
    {
      const int N = N_ID[arc];
      if (lv==LABEL_S)
      {
        if (NONSAT(arc,v) && LABEL(N)!=lv) { Q_PUSH_BACK1(v); goto next_node; }
      }
      else
      {
        if (NONSAT_SISTER(arc,N) && LABEL(N)==LABEL_F) { Q_PUSH_BACK1(v); goto next_node; }
      }
    }

    next_node:
      ;
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::get_segment(int v) const
{
  static const int map_label[3] = { 0, 0, 1 };

  return map_label[LABEL(v)];
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline type_flow GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::get_flow() const
{
  return MAXFLOW;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::bad_alloc() const
{
  return (mem_pool==NULL);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::GridGraph_2D_8C(int w,int h) :
  ow(w),
  oh(h),
  W(next_higher_mul8(w+2)),
  H(next_higher_mul8(h+2)),
  WB(W/8),
  YOFS((WB-1)*64+8),
  mem_pool(NULL)
{
  size_t mem_pool_size = (W*H*sizeof(unsigned char)+64+         // label
                          W*H*sizeof(unsigned char)+64+         // parent
                          W*H*sizeof(int)+64+                   // parent_id
                          8*(W*H*sizeof(type_ncap)+64)+         // rc[8]
                          W*H*sizeof(type_tcap)+64+             // rc_st
                          W*H*sizeof(int)+64+                   // timestamp
                          W*H*sizeof(int)+64+                   // orphans
                          W*H*sizeof(int)+64+                   // orphans2
                          W*H*sizeof(int)+64+                   // free_nodes
                          W*H*sizeof(int)+64+                   // QN
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

  label = (unsigned char*)align(pool_malloc(W*H*sizeof(unsigned char)+64),64);

  parent = (unsigned char*)align(pool_malloc(W*H*sizeof(unsigned char)+64),64);

  parent_id = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  for(int i=0;i<8;i++)
  {
    rc[i] = (type_ncap*)align(pool_malloc(W*H*sizeof(type_ncap)+64),64);
  }

  rc_st = (type_tcap*)align(pool_malloc(W*H*sizeof(type_tcap)+64),64);

  rc_sister[GE] = rc[LE];
  rc_sister[LE] = rc[GE];
  rc_sister[EG] = rc[EL];
  rc_sister[EL] = rc[EG];
  rc_sister[LL] = rc[GG];
  rc_sister[GL] = rc[LG];
  rc_sister[LG] = rc[GL];
  rc_sister[GG] = rc[LL];

  timestamp = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  memset(parent,NONE,W*H);

  QN = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);
  QF = 0;
  QB = 0;
  QN[0] = 1;

  orphans.buffer = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);
  orphans.clear();

  orphans2.buffer = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);
  orphans2.clear();

  free_nodes.buffer = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);
  free_nodes.clear();

  MAXFLOW = 0;

  TIME = 0;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_2D_8C<type_tcap,type_ncap,type_flow>::~GridGraph_2D_8C()
{
  if (mem_pool!=NULL) { free(mem_pool); }
}

#undef SISTER

#undef LABEL

#undef RC
#undef RC_SISTER
#undef RC_ST

#undef PARENT
#undef PARENT_ID

#undef TIMESTAMP

#undef NONSAT
#undef NONSAT_SISTER
#undef ENSAT
#undef DESAT

#undef N_LE
#undef N_GE
#undef N_EL
#undef N_EG
#undef N_LL
#undef N_GL
#undef N_LG
#undef N_GG

#endif
