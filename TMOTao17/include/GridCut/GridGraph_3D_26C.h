// Copyright (C) 2012-2015 Czech Technical University in Prague - All Rights Reserved

#ifndef GRIDGRAPH_3D_26C_H_
#define GRIDGRAPH_3D_26C_H_

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
class GridGraph_3D_26C
{
public:
  GridGraph_3D_26C(int width,int height,int depth);
  ~GridGraph_3D_26C();

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
                       const type_arg_ncap* cap_eeg,        // [ 0, 0,+1]
                       const type_arg_ncap* cap_lle,        // [-1,-1, 0]
                       const type_arg_ncap* cap_gle,        // [+1,-1, 0]
                       const type_arg_ncap* cap_lge,        // [-1,+1, 0]
                       const type_arg_ncap* cap_gge,        // [+1,+1, 0]
                       const type_arg_ncap* cap_ell,        // [ 0,-1,-1]
                       const type_arg_ncap* cap_egl,        // [ 0,+1,-1]
                       const type_arg_ncap* cap_elg,        // [ 0,-1,+1]
                       const type_arg_ncap* cap_egg,        // [ 0,+1,+1]
                       const type_arg_ncap* cap_lel,        // [-1, 0,-1]
                       const type_arg_ncap* cap_leg,        // [-1, 0,+1]
                       const type_arg_ncap* cap_gel,        // [+1, 0,-1]
                       const type_arg_ncap* cap_geg,        // [+1, 0,+1]
                       const type_arg_ncap* cap_lll,        // [-1,-1,-1]
                       const type_arg_ncap* cap_gll,        // [+1,-1,-1]
                       const type_arg_ncap* cap_lgl,        // [-1,+1,-1]
                       const type_arg_ncap* cap_ggl,        // [+1,+1,-1]
                       const type_arg_ncap* cap_llg,        // [-1,-1,+1]
                       const type_arg_ncap* cap_glg,        // [+1,-1,+1]
                       const type_arg_ncap* cap_lgg,        // [-1,+1,+1]
                       const type_arg_ncap* cap_ggg         // [+1,+1,+1]
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

  enum Parent { LEE=0,
                GEE,
                ELE,
                EGE,
                EEL,
                EEG,
                LLE,
                GLE,
                LGE,
                GGE,
                ELL,
                EGL,
                ELG,
                EGG,
                LEL,
                LEG,
                GEL,
                GEG,
                LLL,
                GLL,
                LGL,
                GGL,
                LLG,
                GLG,
                LGG,
                GGG,
                NONE,
                TERMINAL };

  static const unsigned char SISTER_TABLE[26];

  static const int INF_D = INT_MAX;

  union SatFlags
  {
    struct
    {
      unsigned int nonsat_lee : 1;
      unsigned int nonsat_gee : 1;
      unsigned int nonsat_ele : 1;
      unsigned int nonsat_ege : 1;
      unsigned int nonsat_eel : 1;
      unsigned int nonsat_eeg : 1;
      unsigned int nonsat_lle : 1;
      unsigned int nonsat_gle : 1;
      unsigned int nonsat_lge : 1;
      unsigned int nonsat_gge : 1;
      unsigned int nonsat_ell : 1;
      unsigned int nonsat_egl : 1;
      unsigned int nonsat_elg : 1;
      unsigned int nonsat_egg : 1;
      unsigned int nonsat_lel : 1;
      unsigned int nonsat_leg : 1;
      unsigned int nonsat_gel : 1;
      unsigned int nonsat_geg : 1;
      unsigned int nonsat_lll : 1;
      unsigned int nonsat_gll : 1;
      unsigned int nonsat_lgl : 1;
      unsigned int nonsat_ggl : 1;
      unsigned int nonsat_llg : 1;
      unsigned int nonsat_glg : 1;
      unsigned int nonsat_lgg : 1;
      unsigned int nonsat_ggg : 1;
    };
    unsigned int bits;
  };

  unsigned char* label;

  SatFlags* sat_flags;

  unsigned char* parent;

  int* parent_id;

  int* dist;

  type_ncap* rc[26];

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

#define LABEL(X)      label[(X)]
#define RC(E,V)       rc[(E)][(V)]
#define RC_ST(V)      rc_st[(V)]

#define PARENT(X)     parent[(X)]
#define PARENT_ID(X)  parent_id[(X)]

#define TIMESTAMP(V)  timestamp[(V)]
#define DIST(V)       dist[(V)]

#define NONSAT(E,V) ( sat_flags[(V)].bits & (1<<(E)) )
#define ENSAT(E,V)    sat_flags[(V)].bits &= (~(1<<(E)))
#define DESAT(E,V)    sat_flags[(V)].bits |= (1<<(E))

#define N_LEE(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) +  0  +  0  )
#define N_GEE(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) +  0  +  0  )
#define N_ELE(v) ( (v) +  0  + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) +  0  )
#define N_EGE(v) ( (v) +  0  + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) +  0  )
#define N_EEL(v) ( (v) +  0  +  0  + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_EEG(v) ( (v) +  0  +  0  + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_LLE(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) +  0  )
#define N_GLE(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) +  0  )
#define N_LGE(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) +  0  )
#define N_GGE(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) +  0  )
#define N_ELL(v) ( (v) +  0  + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_EGL(v) ( (v) +  0  + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_ELG(v) ( (v) +  0  + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_EGG(v) ( (v) +  0  + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_LEL(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) +  0  + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_LEG(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) +  0  + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_GEL(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) +  0  + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_GEG(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) +  0  + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_LLL(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_GLL(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_LGL(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_GGL(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( (  (v)  & 0x00000030) == 0 ) ? -ZOFS : -16 ) )
#define N_LLG(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_GLG(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( (  (v)  & 0x0000000C) == 0 ) ? -YOFS : -4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_LGG(v) ( (v) + ( ( (  (v)  & 0x00000003) == 0 ) ? -61   : -1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )
#define N_GGG(v) ( (v) + ( ( ((~(v)) & 0x00000003) == 0 ) ? +61   : +1  ) + ( ( ((~(v)) & 0x0000000C) == 0 ) ? +YOFS : +4  ) + ( ( ((~(v)) & 0x00000030) == 0 ) ? +ZOFS : +16 ) )

template <typename type_tcap,typename type_ncap,typename type_flow>
const unsigned char GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::SISTER_TABLE[26] =
{
  GEE,
  LEE,
  EGE,
  ELE,
  EEG,
  EEL,
  GGE,
  LGE,
  GLE,
  LLE,
  EGG,
  ELG,
  EGL,
  ELL,
  GEG,
  GEL,
  LEG,
  LEL,
  GGG,
  LGG,
  GLG,
  LLG,
  GGL,
  LGL,
  GLL,
  LLL
};

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::nodeId(unsigned int x,unsigned int y,unsigned int z) const
{
  return ( ((x>>2) + ((y>>2)*WB) + ((z>>2)*WBHB)) << 6 ) + ( (x&3) + ((y&3)<<2) + ((z&3)<<4) );
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::grow(int&    vs,
                                                          int&    vt,
                                                          Parent& st,
                                                          const   int YOFS,
                                                          const   int ZOFS)
{
 while(!Q_EMPTY())
  {
    const int P = Q_FRONT();

    const int tree_p = LABEL(P);

    const int N_ID[26] = {  N_LEE(P),
                            N_GEE(P),
                            N_ELE(P),
                            N_EGE(P),
                            N_EEL(P),
                            N_EEG(P),
                            N_LLE(P),
                            N_GLE(P),
                            N_LGE(P),
                            N_GGE(P),
                            N_ELL(P),
                            N_EGL(P),
                            N_ELG(P),
                            N_EGG(P),
                            N_LEL(P),
                            N_LEG(P),
                            N_GEL(P),
                            N_GEG(P),
                            N_LLL(P),
                            N_GLL(P),
                            N_LGL(P),
                            N_GGL(P),
                            N_LLG(P),
                            N_GLG(P),
                            N_LGG(P),
                            N_GGG(P)  };

    if (tree_p==LABEL_S)
    {
      for(int n=0;n<26;n++)
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
      for(int n=0;n<26;n++)
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
type_ncap GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::find_minrf_s(int v,type_ncap minrf) const
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
type_ncap GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::find_minrf_t(int v,type_ncap minrf) const
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
void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::aug_s(int v,const type_ncap minrf)
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
void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::aug_t(int v,const type_ncap minrf)
{
  while(PARENT(v) != TERMINAL)
  {
    RC( PARENT(v) , v) -= minrf;
    RC( SISTER( PARENT(v) ) ,PARENT_ID(v)) += minrf;

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
void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::augment(const int    vs,
                                                             const int    vt,
                                                             const Parent st)
{
  type_ncap minrf = RC( SISTER(st) , vs);

  minrf=find_minrf_s(vs,minrf);
  minrf=find_minrf_t(vt,minrf);

  RC( SISTER(st) , vs) -= minrf;
  RC( st , vt) += minrf;

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
int GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::find_origin(int v,const int TIME)
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
void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::adopt(const int TIME,const int YOFS,const int ZOFS)
{
  orphans2.clear();
  free_nodes.clear();

  while(!(orphans.empty()&&orphans2.empty()))
  {
    const int V = (!orphans2.empty()) ? orphans2.pop_front() : orphans.pop_back();

    const int N_ID[26] = {  N_LEE(V),
                            N_GEE(V),
                            N_ELE(V),
                            N_EGE(V),
                            N_EEL(V),
                            N_EEG(V),
                            N_LLE(V),
                            N_GLE(V),
                            N_LGE(V),
                            N_GGE(V),
                            N_ELL(V),
                            N_EGL(V),
                            N_ELG(V),
                            N_EGG(V),
                            N_LEL(V),
                            N_LEG(V),
                            N_GEL(V),
                            N_GEG(V),
                            N_LLL(V),
                            N_GLL(V),
                            N_LGL(V),
                            N_GGL(V),
                            N_LLG(V),
                            N_GLG(V),
                            N_LGG(V),
                            N_GGG(V)  };

    const int tree_p = LABEL(V);

    int min_d = INF_D;
    int best_a = -1;

    if      (tree_p==LABEL_S)
    {
      for(int i=0;i<26;i++)
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
      for(int i=0;i<26;i++)
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

    for(int i=0;i<26;i++)
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

    const int N_ID[26] = {  N_LEE(V),
                            N_GEE(V),
                            N_ELE(V),
                            N_EGE(V),
                            N_EEL(V),
                            N_EEG(V),
                            N_LLE(V),
                            N_GLE(V),
                            N_LGE(V),
                            N_GGE(V),
                            N_ELL(V),
                            N_EGL(V),
                            N_ELG(V),
                            N_EGG(V),
                            N_LEL(V),
                            N_LEG(V),
                            N_GEL(V),
                            N_GEG(V),
                            N_LLL(V),
                            N_GLL(V),
                            N_LGL(V),
                            N_GGL(V),
                            N_LLG(V),
                            N_GLG(V),
                            N_LGG(V),
                            N_GGG(V)  };

    for(int i=0;i<26;i++)
    {
      const int N = N_ID[i];

      if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S) { Q_PUSH_BACK2(N); }
      if (NONSAT(i,V) && LABEL(N)==LABEL_T) { Q_PUSH_BACK2(N); }
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::compute_maxflow()
{
  while(!free_nodes.empty())
  {
    const int v = free_nodes.pop_front();

    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[26] = {  N_LEE(v),
                            N_GEE(v),
                            N_ELE(v),
                            N_EGE(v),
                            N_EEL(v),
                            N_EEG(v),
                            N_LLE(v),
                            N_GLE(v),
                            N_LGE(v),
                            N_GGE(v),
                            N_ELL(v),
                            N_EGL(v),
                            N_ELG(v),
                            N_EGG(v),
                            N_LEL(v),
                            N_LEG(v),
                            N_GEL(v),
                            N_GEG(v),
                            N_LLL(v),
                            N_GLL(v),
                            N_LGL(v),
                            N_GGL(v),
                            N_LLG(v),
                            N_GLG(v),
                            N_LGG(v),
                            N_GGG(v)  };

    for(int arc=0;arc<26;arc++)
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
inline int GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::node_id(int x,int y,int z) const
{
  return nodeId(x+1,y+1,z+1);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::set_terminal_cap(int v,type_tcap cap_s,type_tcap cap_t)
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
inline void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::set_neighbor_cap(int node_id,int offset_x,int offset_y,int offset_z,type_ncap cap)
{
  static const unsigned char offsets2arc[3][3][3] = {  { {LLL,ELL,GLL},
                                                         {LEL,EEL,GEL},
                                                         {LGL,EGL,GGL} },

                                                       { {LLE,ELE,GLE},
                                                         {LEE,255,GEE},
                                                         {LGE,EGE,GGE} },

                                                       { {LLG,ELG,GLG},
                                                         {LEG,EEG,GEG},
                                                         {LGG,EGG,GGG} }  };

  const int arc = offsets2arc[offset_z+1][offset_y+1][offset_x+1];
  rc[arc][node_id] = cap;

  if (cap!=0) { DESAT(arc,node_id); }
}

template<typename type_tcap,typename type_ncap,typename type_flow> template<typename type_arg_tcap,typename type_arg_ncap>
inline void GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::set_caps(const type_arg_tcap* cap_s,
                                                                      const type_arg_tcap* cap_t,
                                                                      const type_arg_ncap* cap_lee,
                                                                      const type_arg_ncap* cap_gee,
                                                                      const type_arg_ncap* cap_ele,
                                                                      const type_arg_ncap* cap_ege,
                                                                      const type_arg_ncap* cap_eel,
                                                                      const type_arg_ncap* cap_eeg,
                                                                      const type_arg_ncap* cap_lle,
                                                                      const type_arg_ncap* cap_gle,
                                                                      const type_arg_ncap* cap_lge,
                                                                      const type_arg_ncap* cap_gge,
                                                                      const type_arg_ncap* cap_ell,
                                                                      const type_arg_ncap* cap_egl,
                                                                      const type_arg_ncap* cap_elg,
                                                                      const type_arg_ncap* cap_egg,
                                                                      const type_arg_ncap* cap_lel,
                                                                      const type_arg_ncap* cap_leg,
                                                                      const type_arg_ncap* cap_gel,
                                                                      const type_arg_ncap* cap_geg,
                                                                      const type_arg_ncap* cap_lll,
                                                                      const type_arg_ncap* cap_gll,
                                                                      const type_arg_ncap* cap_lgl,
                                                                      const type_arg_ncap* cap_ggl,
                                                                      const type_arg_ncap* cap_llg,
                                                                      const type_arg_ncap* cap_glg,
                                                                      const type_arg_ncap* cap_lgg,
                                                                      const type_arg_ncap* cap_ggg
                                                                      )
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
    RC(LLE,v) = cap_lle[xyz]; if (cap_lle[xyz]!=0) { DESAT(LLE,v); }
    RC(GLE,v) = cap_gle[xyz]; if (cap_gle[xyz]!=0) { DESAT(GLE,v); }
    RC(LGE,v) = cap_lge[xyz]; if (cap_lge[xyz]!=0) { DESAT(LGE,v); }
    RC(GGE,v) = cap_gge[xyz]; if (cap_gge[xyz]!=0) { DESAT(GGE,v); }
    RC(ELL,v) = cap_ell[xyz]; if (cap_ell[xyz]!=0) { DESAT(ELL,v); }
    RC(EGL,v) = cap_egl[xyz]; if (cap_egl[xyz]!=0) { DESAT(EGL,v); }
    RC(ELG,v) = cap_elg[xyz]; if (cap_elg[xyz]!=0) { DESAT(ELG,v); }
    RC(EGG,v) = cap_egg[xyz]; if (cap_egg[xyz]!=0) { DESAT(EGG,v); }
    RC(LEL,v) = cap_lel[xyz]; if (cap_lel[xyz]!=0) { DESAT(LEL,v); }
    RC(LEG,v) = cap_leg[xyz]; if (cap_leg[xyz]!=0) { DESAT(LEG,v); }
    RC(GEL,v) = cap_gel[xyz]; if (cap_gel[xyz]!=0) { DESAT(GEL,v); }
    RC(GEG,v) = cap_geg[xyz]; if (cap_geg[xyz]!=0) { DESAT(GEG,v); }
    RC(LLL,v) = cap_lll[xyz]; if (cap_lll[xyz]!=0) { DESAT(LLL,v); }
    RC(GLL,v) = cap_gll[xyz]; if (cap_gll[xyz]!=0) { DESAT(GLL,v); }
    RC(LGL,v) = cap_lgl[xyz]; if (cap_lgl[xyz]!=0) { DESAT(LGL,v); }
    RC(GGL,v) = cap_ggl[xyz]; if (cap_ggl[xyz]!=0) { DESAT(GGL,v); }
    RC(LLG,v) = cap_llg[xyz]; if (cap_llg[xyz]!=0) { DESAT(LLG,v); }
    RC(GLG,v) = cap_glg[xyz]; if (cap_glg[xyz]!=0) { DESAT(GLG,v); }
    RC(LGG,v) = cap_lgg[xyz]; if (cap_lgg[xyz]!=0) { DESAT(LGG,v); }
    RC(GGG,v) = cap_ggg[xyz]; if (cap_ggg[xyz]!=0) { DESAT(GGG,v); }
  }

  for(int v=0;v<W*H*D;v++)
  {
    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int N_ID[26] = {  N_LEE(v),
                            N_GEE(v),
                            N_ELE(v),
                            N_EGE(v),
                            N_EEL(v),
                            N_EEG(v),
                            N_LLE(v),
                            N_GLE(v),
                            N_LGE(v),
                            N_GGE(v),
                            N_ELL(v),
                            N_EGL(v),
                            N_ELG(v),
                            N_EGG(v),
                            N_LEL(v),
                            N_LEG(v),
                            N_GEL(v),
                            N_GEG(v),
                            N_LLL(v),
                            N_GLL(v),
                            N_LGL(v),
                            N_GGL(v),
                            N_LLG(v),
                            N_GLG(v),
                            N_LGG(v),
                            N_GGG(v)  };

    for(int arc=0;arc<26;arc++)
    {
      const int N = N_ID[arc];

      if (lv==LABEL_S)
      {
        if (NONSAT(arc,v) && LABEL(N)!=lv) { Q_PUSH_BACK1(v); goto next_node; }
      }
      else
      {
        if (NONSAT(SISTER(arc),N) && LABEL(N)==LABEL_F) {  Q_PUSH_BACK1(v); goto next_node; }
      }
    }

    next_node:
      ;
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::get_segment(int v) const
{
  static const int map_label[3] = { 0, 0, 1 };

  return map_label[LABEL(v)];
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline type_flow GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::get_flow() const
{
  return MAXFLOW;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::bad_alloc() const
{
  return (mem_pool==NULL);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::GridGraph_3D_26C(int w,int h,int d) :
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
   size_t mem_pool_size = (
                          W*H*D*sizeof(unsigned char)+64+         // label
                          W*H*D*sizeof(SatFlags)+64+              // sat_flags
                          W*H*D*sizeof(unsigned char)+64+         // parent
                          W*H*D*sizeof(int)+64+                   // parent_id
                          26*(W*H*D*sizeof(type_ncap)+64)+        // rc[26]
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

  label = (unsigned char*)align(pool_malloc(W*H*D*sizeof(unsigned char)+64),64);

  sat_flags = (SatFlags*)align(pool_malloc(W*H*D*sizeof(SatFlags)+64),64);

  parent = (unsigned char*)align(pool_malloc(W*H*D*sizeof(unsigned char)+64),64);

  parent_id = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  for(int i=0;i<26;i++)
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
GridGraph_3D_26C<type_tcap,type_ncap,type_flow>::~GridGraph_3D_26C()
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
#undef N_LLE
#undef N_GLE
#undef N_LGE
#undef N_GGE
#undef N_ELL
#undef N_EGL
#undef N_ELG
#undef N_EGG
#undef N_LEL
#undef N_LEG
#undef N_GEL
#undef N_GEG
#undef N_LLL
#undef N_GLL
#undef N_LGL
#undef N_GGL
#undef N_LLG
#undef N_GLG
#undef N_LGG
#undef N_GGG

#endif
