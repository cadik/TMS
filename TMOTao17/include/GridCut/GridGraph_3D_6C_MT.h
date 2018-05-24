// Copyright (C) 2012-2015 Czech Technical University in Prague - All Rights Reserved

#ifndef GRIDGRAPH_3D_6C_MT_H_
#define GRIDGRAPH_3D_6C_MT_H_

#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <new>

#ifdef __GNUC__
  #include <stdint.h>
#endif

#ifdef _WIN32
  extern "C" struct _SECURITY_ATTRIBUTES;
  extern "C" __declspec(dllimport) void* __stdcall
  CreateThread(_SECURITY_ATTRIBUTES*,
               #if defined(_WIN64)
                 unsigned __int64,
               #else
                 unsigned long,
               #endif
               unsigned long (__stdcall *)(void*),void*,unsigned long,unsigned long*);
  extern "C" __declspec(dllimport) unsigned long __stdcall WaitForSingleObject(void*,unsigned long);    
  #ifdef _WIN64
    #include <intrin.h>
    #define InterlockedCompareExchange _InterlockedCompareExchange
    #define InterlockedIncrement _InterlockedIncrement
  #else
    extern "C" __declspec(dllimport) long __stdcall InterlockedCompareExchange(long volatile*,long,long);
    extern "C" __declspec(dllimport) long __stdcall InterlockedIncrement(long volatile*);
  #endif
#else
  #include <pthread.h>
  #ifdef __APPLE__
    #include <libkern/OSAtomic.h>
  #endif
#endif

template
<
  typename type_tcap, // Type used to represent capacities of edges between nodes and terminals.
  typename type_ncap, // Type used to represent capacities of edges between nodes and their neighbors.
  typename type_flow  // Type used to represent the total flow.
>
class GridGraph_3D_6C_MT
{
public:
  GridGraph_3D_6C_MT(int width,int height,int depth,int num_threads,int block_size);
  ~GridGraph_3D_6C_MT();

  // Returns the index of a node at grid coordinates [x,y].
  // The index is used to set the capacities of node's outgoing edges
  // using the set_neighbor_cap function, to set the capacities of edges
  // between node and source/sink terminals using the set_terminal_cap
  // function, and to retrieve the segment to which the node belongs
  // after the maxflow computation using the get_segment function.
  inline int node_id(int x,int y,int z) const;

  // Alternative way to set the edge capacities is to use the set_caps function which
  // sets capacities of all edges at once based on values from input arrays.
  // Each array has width*height elements, where each element corresponds to one node.
  // For example, cap_le[x+y*width] is the capacity of the outgoing edge from node [x,y]
  // to node [x-1,y], and cap_sink[x+y*width] is capacity of edge from node [x,y] to sink.
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

  type_ncap* rc_lee;
  type_ncap* rc_gee;
  type_ncap* rc_ele;
  type_ncap* rc_ege;
  type_ncap* rc_eel;
  type_ncap* rc_eeg;

  type_ncap* rc[6];

  type_tcap* rc_st;

  int* timestamp;

  int* active_next;

  int* orphans1_prev;
  int* orphans2_next;

  int* free_nodes_prev;

  inline void ACTIVE_PUSH_BACK(int V,int& QF,int& QB)
  {
    if (active_next[V]==0) { active_next[QB] = V; QF=QF?QF:V; QB=V; active_next[V]=V; }
  }

  inline void ACTIVE_POP_FRONT(int& QF,int& QB)
  {
    const int T = QF; QF = active_next[QF]; active_next[T]=0; if (T==QF) { QB=0; QF=0; }
  }

  inline void ORPHANS1_PUSH_BACK(int v,int& orphans1_back)
  {
    orphans1_prev[v] = orphans1_back; orphans1_back = v;
  }

  inline int ORPHANS1_POP_BACK(int& orphans1_back)
  {
    const int tmp = orphans1_back; orphans1_back = orphans1_prev[tmp]; orphans1_prev[tmp] = 0; return tmp;
  }

  inline void ORPHANS2_PUSH_BACK(int V,int& o2_front,int& o2_back)
  {
    orphans2_next[o2_back] = V; o2_front=o2_front?o2_front:V; o2_back=V;
  }

  inline int  ORPHANS2_POP_FRONT(int& o2_front,int& o2_back)
  {
    const int T = o2_front; o2_front = orphans2_next[T]; orphans2_next[T]=0; if (!o2_front) { o2_back=0; } return T;
  }

  inline void FREE_NODES_PUSH_BACK(int v,int& free_nodes_back)
  {
    free_nodes_prev[v] = free_nodes_back; free_nodes_back = v;
  }

  inline int FREE_NODES_POP_BACK(int& free_nodes_back)
  {
    const int tmp = free_nodes_back; free_nodes_back = free_nodes_prev[tmp]; free_nodes_prev[tmp] = 0; return tmp;
  }

  const int _w;
  const int _h;
  const int _d;

  const int W;
  const int H;
  const int D;

  const int WB;
  const int WBHB;

  int _YOFS;
  int _ZOFS;

  type_flow MAXFLOW_TOTAL;

  int nodeId(unsigned int x,unsigned int y,unsigned int z) const;

  bool grow(int&         vs,
            int&         vt,
            Parent&      st,
            int& active_front,
            int& active_back,
            const int    YOFS,
            const int    ZOFS);

  type_ncap find_minrf_s(int v,type_ncap minrf) const;

  type_ncap find_minrf_t(int v,type_ncap minrf) const;

  void aug_s(int v,const type_ncap minrf,int& orphans1_back);

  void aug_t(int v,const type_ncap minrf,int& orphans1_back);

  void augment(const int    vs,
               const int    vt,
               const Parent st,
               type_flow&   MAXFLOW,
               int& orphans1_back);

  int find_origin(int v,const int TIME);

  void adopt(int& active_front,
             int& active_back,
             int& orphans1_back,
             const int TIME,
             const int YOFS,
             const int ZOFS);

  void run_BK(int& active_front,
              int& active_back,
              int& orphans1_back,
              int&         TIME,
              type_flow&   MAXFLOW,
              const int    YOFS,
              const int    ZOFS);

  enum SegmentType { SEGMENT_X, SEGMENT_Y, SEGMENT_Z };

  struct Segment
  {
    SegmentType type;
    int length;
    int* node_id;
  };

  struct Boundary
  {
    int block_ids[2];
    Segment segment;

    Boundary() {}
    Boundary(int first_block_id,int second_block_id,SegmentType type,int length,int* node_id)
    {
      block_ids[0] = first_block_id;
      block_ids[1] = second_block_id;

      segment.type = type;
      segment.length = length;
      segment.node_id = node_id;
    }
  };

#ifdef _WIN32
  typedef void* ThreadHandle;

  ThreadHandle threadCreate(unsigned long (__stdcall *start_routine)(void*),void* arg)
  {
    return CreateThread(NULL,0,start_routine,arg,0,NULL);
  }

  void threadJoin(ThreadHandle handle)
  {
    WaitForSingleObject(handle,0xFFFFFFFF);
  }

  #define THREAD_FUNC_WRAPPER_DECL static unsigned long __stdcall thread_func_wrapper(void* arg)

  typedef long Spinlock;

  void spin_init(volatile Spinlock* spinlock)
  {
    *spinlock = 1;
  }

  void spin_lock(volatile Spinlock* spinlock)
  {
    while (InterlockedCompareExchange(spinlock,2,1)==2)
    {
    }
  }

  void spin_unlock(volatile Spinlock* spinlock)
  {
    InterlockedCompareExchange(spinlock,1,2);
  }

  void spin_destroy(volatile Spinlock* spinlock)
  {

  }
#else
  typedef pthread_t ThreadHandle;

  ThreadHandle threadCreate(void*(*start_routine)(void*),void* arg)
  {
    pthread_t thread_handle;
    pthread_attr_t thread_attr;

    pthread_attr_init(&thread_attr);
    pthread_attr_setdetachstate(&thread_attr,PTHREAD_CREATE_JOINABLE);

    pthread_create(&thread_handle,&thread_attr,start_routine,arg);

    return thread_handle;
  }

  void threadJoin(ThreadHandle handle)
  {
    pthread_join(handle,NULL);
  }

  #define THREAD_FUNC_WRAPPER_DECL static void* thread_func_wrapper(void* arg)

  #ifdef __APPLE__
    typedef OSSpinLock Spinlock;

    void spin_init(volatile Spinlock* spinlock)
    {
      *spinlock = OS_SPINLOCK_INIT;
    }

    void spin_lock(volatile Spinlock* spinlock)
    {
      OSSpinLockLock(spinlock);
    }

    void spin_unlock(volatile Spinlock* spinlock)
    {
      OSSpinLockUnlock(spinlock);
    }

    void spin_destroy(volatile Spinlock* spinlock)
    {
    }
  #elif __ANDROID__
    typedef long Spinlock;

    void spin_init(volatile Spinlock* spinlock)
    {
      *spinlock = 1;
    }

    void spin_lock(volatile Spinlock* spinlock)
    {
      while (__sync_val_compare_and_swap(spinlock,1,2)==2)
      {
      }
    }

    void spin_unlock(volatile Spinlock* spinlock)
    {
      __sync_val_compare_and_swap(spinlock,2,1);
    }

    void spin_destroy(volatile Spinlock* spinlock)
    {
    }
  #else
    typedef pthread_spinlock_t Spinlock;

    void spin_init(volatile Spinlock* spinlock)
    {
      pthread_spin_init(spinlock,PTHREAD_PROCESS_PRIVATE);
    }

    void spin_lock(volatile Spinlock* spinlock)
    {
      pthread_spin_lock(spinlock);
    }

    void spin_unlock(volatile Spinlock* spinlock)
    {
      pthread_spin_unlock(spinlock);
    }

    void spin_destroy(volatile Spinlock* spinlock)
    {
      pthread_spin_destroy(spinlock);
    }
  #endif
#endif

  struct Thread
  {
    ThreadHandle handle;

    int active_front;
    int active_back;

    int orphans1_back;

    type_flow MAXFLOW;
    int YOFS;
    int ZOFS;
  };

  std::vector<Boundary> boundary_list;
  std::vector<Boundary> boundary_list_new;

  int BLOCK_WIDTH;
  int BLOCK_HEIGHT;
  int BLOCK_DEPTH;

  int NUM_BLOCKS_X;
  int NUM_BLOCKS_Y;
  int NUM_BLOCKS_Z;

  int NUM_BLOCKS;

  int NUM_X_SEGMENTS;
  int NUM_Y_SEGMENTS;
  int NUM_Z_SEGMENTS;

  bool* block_locked;
  int* block_TIME;
  int next_block_id;

  int NUM_THREADS;

  Thread** threads;

  Spinlock spinlock;

  void init_block(const int    block_id,
                  int& active_front,
                  int& active_back,
                  const int    YOFS,
                  const int    ZOFS);

  void activate_segment(const Segment& segment,
                        int& active_front,
                        int& active_back,
                        const int      YOFS,
                        const int      ZOFS);

  void deactivate_segment(const Segment& segment,
                          const int      YOFS,
                          const int      ZOFS);

  void* thread_func(void* arg);

  int* uf_root;

  void uf_init(int num)
  {
    uf_root = new int[num];
    for(int i=0;i<num;i++) uf_root[i] = i;
  }

  int uf_find(int x)
  {
    int tmp_x = x;

    while(uf_root[x]!=x)
    {
      x = uf_root[x];
    }

    int y = tmp_x;

    while(uf_root[y]!=x)
    {
      const int tmp = uf_root[y];
      uf_root[y] = x;
      y = tmp;
    }

    return x;
  }

  int uf_unify(int x,int y)
  {
    uf_root[x] = uf_find(y);
    uf_root[y] = uf_find(y);

    return uf_root[y];
  }

  void uf_cleanup()
  {
    delete[] uf_root;
  }

  template<typename T>inline type_ncap mincap(type_ncap a,T b) const
  {
    return a < b ? a : b;
  }

  int next_higher_mul4(int x)
  {
    return ((x-1)/4)*4+4;
  }

#ifdef _WIN32
  #define ATOMIC_INCREMENT(X) (InterlockedIncrement((long*)(X)) - 1)
#else
  #define ATOMIC_INCREMENT(X) (__sync_fetch_and_add(X,1))
#endif

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

  template<typename T>inline const T& minimum(const T& a,const T& b) const
  {
    return a < b ? a : b;
  }

  template<typename T>inline const T& maximum(const T& a,const T& b) const
  {
    return a > b ? a : b;
  }

  struct thread_func_wrapper_arg
  {
    void* thisptr;
    void* arg;
  };

  THREAD_FUNC_WRAPPER_DECL
  {
    thread_func_wrapper_arg* wrapper_arg = (thread_func_wrapper_arg*)arg;

    ((GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>*)wrapper_arg->thisptr)->thread_func(wrapper_arg->arg);

    return 0;
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
const unsigned char GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::SISTER_TABLE[6] =
{
  GEE,
  LEE,
  EGE,
  ELE,
  EEG,
  EEL
};

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::nodeId(unsigned int x,unsigned int y,unsigned int z) const
{
  return ( ((x>>2) + ((y>>2)*WB) + ((z>>2)*WBHB)) << 6 ) + ( (x&3) + ((y&3)<<2) + ((z&3)<<4) );
}


template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::grow(int& vs,
                                                             int& vt,
                                                             Parent& st,
                                                             int& active_front,
                                                             int& active_back,
                                                             const int YOFS,
                                                             const int ZOFS)
{
  while(active_front!=0)
  {
    const int P = active_front;

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
          if (LABEL(N)==LABEL_F) { LABEL(N) = LABEL_S; ACTIVE_PUSH_BACK(N,active_front,active_back); PARENT(N) = SISTER(n); PARENT_ID(N) = P; }
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
          if (LABEL(N)==LABEL_F ) { LABEL(N) = LABEL_T; ACTIVE_PUSH_BACK(N,active_front,active_back); PARENT(N) = SISTER(n); PARENT_ID(N) = P; }
        }
      }
    }

    ACTIVE_POP_FRONT(active_front,active_back);
  }

  return false;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
type_ncap GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::find_minrf_s(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,rc[ SISTER(PARENT(v)) ][PARENT_ID(v)]);
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,rc_st[v]);

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
type_ncap GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::find_minrf_t(int v,type_ncap minrf) const
{
  while(PARENT(v) != TERMINAL)
  {
    minrf = mincap(minrf,rc[ PARENT(v) ][v]);
    v = PARENT_ID(v);
  }

  minrf = mincap(minrf,rc_st[v]);

  return minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::aug_s(int v,const type_ncap minrf,int& orphans1_back)
{
  while(PARENT(v) != TERMINAL)
  {
    rc[ SISTER(PARENT(v)) ][ PARENT_ID(v) ] -= minrf;
    rc[ PARENT(v) ][v] += minrf;

    DESAT(PARENT(v),v);

    if(! (rc[ SISTER(PARENT(v)) ][ PARENT_ID(v) ]) )
    {
      ENSAT( SISTER(PARENT(v)), PARENT_ID(v) );
      PARENT(v) = NONE;
      ORPHANS1_PUSH_BACK(v,orphans1_back);
    }

    v = PARENT_ID(v);
  }

  rc_st[v] -= minrf;
  if (!rc_st[v])
  {
    PARENT(v) = NONE;
    ORPHANS1_PUSH_BACK(v,orphans1_back);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::aug_t(int v,const type_ncap minrf,int& orphans1_back)
{
  while(PARENT(v) != TERMINAL)
  {
    rc[ PARENT(v) ][v] -= minrf;
    rc[ SISTER( PARENT(v) ) ][PARENT_ID(v)] += minrf;

    DESAT( SISTER( PARENT(v) ), PARENT_ID(v) );

    if (! (rc[ PARENT(v) ][v]))
    {
      ENSAT(PARENT(v),v);

      PARENT(v) = NONE;
      ORPHANS1_PUSH_BACK(v,orphans1_back);
    }

    v = PARENT_ID(v);
  }

  rc_st[v] -= minrf;
  if (!rc_st[v])
  {
    PARENT(v) = NONE;
    ORPHANS1_PUSH_BACK(v,orphans1_back);
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::augment(const int    vs,
                                                                const int    vt,
                                                                const Parent st,
                                                                type_flow&   MAXFLOW,
                                                                int& orphans1_back)
{
  type_ncap minrf = rc[ SISTER(st) ][vs];

  minrf=find_minrf_s(vs,minrf);
  minrf=find_minrf_t(vt,minrf);

  rc[ SISTER(st) ][vs] -= minrf;
  rc[ st ][vt] += minrf;

  DESAT(st,vt);

  if (!(rc[ SISTER(st) ][vs]))
  {
    ENSAT(SISTER(st),vs);
  }

  aug_s(vs,minrf,orphans1_back);
  aug_t(vt,minrf,orphans1_back);

  MAXFLOW += minrf;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
int GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::find_origin(int v,const int TIME)
{
  const int start_v = v;
  int d2;
  int d = 1;

  while(1)
  {
    if (timestamp[v]==TIME)    { goto L1; }
    if (PARENT(v) == NONE)     { return INF_D; }
    if (PARENT(v) == TERMINAL) { goto L2; }

    d++;
    v = PARENT_ID(v);
  }

  L1:
    d = dist[v]+d;
    v = start_v;
    d2 = d;

    while(timestamp[v]!=TIME)
    {
      dist[v] = d;
      d--;
      timestamp[v]=TIME;
      v = PARENT_ID(v);
    }

    return d2;

  L2:
    v = start_v;
    d2 = d;

    while(PARENT(v)!=TERMINAL)
    {
      dist[v] = d;
      d--;
      timestamp[v]=TIME;
      v = PARENT_ID(v);
    }

    timestamp[v]=TIME;
    dist[v] = 1;

    return d2;
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::adopt(int& active_front,
                                                              int& active_back,
                                                              int& orphans1_back,
                                                              const int TIME,
                                                              const int YOFS,
                                                              const int ZOFS)
{
  int orphans2_front = 0;
  int orphans2_back = 0;

  int free_nodes_back = 0;

  while(!((orphans1_back==0)&&(orphans2_front==0)))
  {
    const int V = (!(orphans2_front==0)) ? ORPHANS2_POP_FRONT(orphans2_front,orphans2_back) : ORPHANS1_POP_BACK(orphans1_back);


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
      dist[V] = min_d+1;
      timestamp[V]=TIME;

      PARENT(V) = best_a;
      PARENT_ID(V) = N_ID[best_a];

      goto next;
    }

    LABEL(V) = LABEL_F;
    FREE_NODES_PUSH_BACK(V,free_nodes_back);

    for(int i=0;i<6;i++)
    {
      const int N = N_ID[i];
      if (LABEL(N)==tree_p && PARENT(N)==SISTER(i)) { PARENT(N) = NONE; ORPHANS2_PUSH_BACK(N,orphans2_front,orphans2_back); }
    }

    next:
      ;
  }

  while(!(free_nodes_back==0))
  {
    const int V = FREE_NODES_POP_BACK(free_nodes_back);

    const int N_ID[6] = { N_LEE(V),
                          N_GEE(V),
                          N_ELE(V),
                          N_EGE(V),
                          N_EEL(V),
                          N_EEG(V) };

    for(int i=0;i<6;i++)
    {
      const int N = N_ID[i];
      if (NONSAT(SISTER(i),N) && LABEL(N)==LABEL_S) { ACTIVE_PUSH_BACK(N,active_front,active_back); }
      if (NONSAT(i,V)         && LABEL(N)==LABEL_T) { ACTIVE_PUSH_BACK(N,active_front,active_back); }
    }
  }
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::
run_BK(int& active_front,
       int& active_back,
       int& orphans1_back,
       int&         TIME,
       type_flow&   MAXFLOW,
       const int    YOFS,
       const int    ZOFS)
{
  while(1)
  {
    int vs;
    int vt;
    Parent st;

    const bool path_found = grow(vs,vt,st,active_front,active_back,YOFS,ZOFS);

    if (!path_found) break;
    TIME++;

    augment(vs,vt,st,MAXFLOW,orphans1_back);

    adopt(active_front,active_back,orphans1_back,TIME,YOFS,ZOFS);
  }
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::
init_block(const int    block_id,
           int& active_front,
           int& active_back,
           const int    YOFS,
           const int    ZOFS)
{
  const int x_min = (block_id % NUM_BLOCKS_X)                                 * BLOCK_WIDTH;
  const int y_min = ((block_id % (NUM_BLOCKS_X*NUM_BLOCKS_Y)) / NUM_BLOCKS_X) * BLOCK_HEIGHT;
  const int z_min = (block_id / (NUM_BLOCKS_X*NUM_BLOCKS_Y))                  * BLOCK_DEPTH;

  const int x_max = minimum(x_min+BLOCK_WIDTH,W);
  const int y_max = minimum(y_min+BLOCK_HEIGHT,H);
  const int z_max = minimum(z_min+BLOCK_DEPTH,D);

  for(int z=z_min;z<z_max;z++)
  for(int y=y_min;y<y_max;y++)
  for(int x=x_min;x<x_max;x++)
  {
    const int v = nodeId(x,y,z);
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
        if (NONSAT(arc,v) && LABEL(N)!=lv) { ACTIVE_PUSH_BACK(v,active_front,active_back);  goto next_node; }
      }
      else
      {
        if (NONSAT(SISTER(arc),N) && LABEL(N)==LABEL_F) { ACTIVE_PUSH_BACK(v,active_front,active_back);  goto next_node; }
      }
    }

    next_node:
      ;
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::
activate_segment(const Segment& segment,
                 int& active_front,
                 int& active_back,
                 const int      YOFS,
                 const int      ZOFS)
{
  if (segment.type==SEGMENT_X)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      if (rc_gee[v]!=0) { DESAT(GEE,v); }
      if (rc_lee[N_GEE(v)]!=0) { DESAT(LEE,N_GEE(v)); }

      if (LABEL(v)==LABEL_S && NONSAT(GEE,v)        && LABEL(N_GEE(v))!=LABEL_S) ACTIVE_PUSH_BACK(v,active_front,active_back);
      if (LABEL(v)==LABEL_T && NONSAT(LEE,N_GEE(v)) && LABEL(N_GEE(v))!=LABEL_T) ACTIVE_PUSH_BACK(v,active_front,active_back);

      if (LABEL(N_GEE(v))==LABEL_S && NONSAT(LEE,N_GEE(v)) && LABEL(v)!=LABEL_S) ACTIVE_PUSH_BACK(N_GEE(v),active_front,active_back);
      if (LABEL(N_GEE(v))==LABEL_T && NONSAT(GEE,v)        && LABEL(v)!=LABEL_T) ACTIVE_PUSH_BACK(N_GEE(v),active_front,active_back);
    }
  }
  else if (segment.type==SEGMENT_Y)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      if (rc_ege[v]!=0) { DESAT(EGE,v); }
      if (rc_ele[N_EGE(v)]!=0) { DESAT(ELE,N_EGE(v)); }

      if (LABEL(v)==LABEL_S && NONSAT(EGE,v)        && LABEL(N_EGE(v))!=LABEL_S) ACTIVE_PUSH_BACK(v,active_front,active_back);
      if (LABEL(v)==LABEL_T && NONSAT(ELE,N_EGE(v)) && LABEL(N_EGE(v))!=LABEL_T) ACTIVE_PUSH_BACK(v,active_front,active_back);

      if (LABEL(N_EGE(v))==LABEL_S && NONSAT(ELE,N_EGE(v)) && LABEL(v)!=LABEL_S) ACTIVE_PUSH_BACK(N_EGE(v),active_front,active_back);
      if (LABEL(N_EGE(v))==LABEL_T && NONSAT(EGE,v)        && LABEL(v)!=LABEL_T) ACTIVE_PUSH_BACK(N_EGE(v),active_front,active_back);
    }
  }
  else if (segment.type==SEGMENT_Z)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      if (rc_eeg[v]!=0) { DESAT(EEG,v); }
      if (rc_eel[N_EEG(v)]!=0) { DESAT(EEL,N_EEG(v)); }

      if (LABEL(v)==LABEL_S && NONSAT(EEG,v)        && LABEL(N_EEG(v))!=LABEL_S) ACTIVE_PUSH_BACK(v,active_front,active_back);
      if (LABEL(v)==LABEL_T && NONSAT(EEL,N_EEG(v)) && LABEL(N_EEG(v))!=LABEL_T) ACTIVE_PUSH_BACK(v,active_front,active_back);

      if (LABEL(N_EEG(v))==LABEL_S && NONSAT(EEL,N_EEG(v)) && LABEL(v)!=LABEL_S) ACTIVE_PUSH_BACK(N_EEG(v),active_front,active_back);
      if (LABEL(N_EEG(v))==LABEL_T && NONSAT(EEG,v)        && LABEL(v)!=LABEL_T) ACTIVE_PUSH_BACK(N_EEG(v),active_front,active_back);
    }
  }
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::
deactivate_segment(const Segment& segment,
                   const int      YOFS,
                   const int      ZOFS)
{
  if (segment.type==SEGMENT_X)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      ENSAT(GEE,v);
      ENSAT(LEE,N_GEE(v));
    }
  }
  else if (segment.type==SEGMENT_Y)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      ENSAT(EGE,v);
      ENSAT(ELE,N_EGE(v));
    }
  }
  else if (segment.type==SEGMENT_Z)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      ENSAT(EEG,v);
      ENSAT(EEL,N_EEG(v));
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void* GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::
thread_func(void* arg)
{
  const int thread_id = (intptr_t)arg;

  Thread& thread = (*threads[thread_id]);

  while(1)
  {
    const int block_id = ATOMIC_INCREMENT(&next_block_id);

    if (block_id>=NUM_BLOCKS) break;

    init_block(block_id,thread.active_front,thread.active_back,thread.YOFS,thread.ZOFS);

    run_BK(thread.active_front,thread.active_back,thread.orphans1_back,block_TIME[block_id],thread.MAXFLOW,thread.YOFS,thread.ZOFS);

    block_locked[block_id] = false;
  }

  std::vector<Segment> segment_list;
  segment_list.reserve(NUM_X_SEGMENTS+NUM_Y_SEGMENTS+NUM_Z_SEGMENTS);

  while(1)
  {
    int block_id = -1;

    spin_lock(&spinlock);

    boundary_list_new.clear();
    segment_list.clear();

    for(unsigned int i=0;i<boundary_list.size();i++)
    {
      const Boundary& boundary = boundary_list[i];

      const int first_block_id = uf_find(boundary.block_ids[0]);
      const int second_block_id = uf_find(boundary.block_ids[1]);

      if (block_locked[first_block_id]==false &&
          block_locked[second_block_id]==false)
      {
        const int TIME = maximum(block_TIME[first_block_id],
                                 block_TIME[second_block_id]);

        block_id = uf_unify(first_block_id,second_block_id);

        block_TIME[block_id] = TIME;
        block_locked[block_id] = true;

        for(unsigned int j=i;j<boundary_list.size();j++)
        {
          const Boundary& boundary = boundary_list[j];

          if (uf_find(boundary.block_ids[0])==block_id &&
              uf_find(boundary.block_ids[1])==block_id)
          {
            segment_list.push_back(boundary.segment);
          }
          else
          {
            boundary_list_new.push_back(boundary);
          }
        }

        swap(boundary_list,boundary_list_new);
        break;
      }
      else
      {
        boundary_list_new.push_back(boundary);
      }
    }

    spin_unlock(&spinlock);

    if (segment_list.empty()) break;

    for(unsigned int i=0;i<segment_list.size();i++)
    {
      activate_segment(segment_list[i],thread.active_front,thread.active_back,thread.YOFS,thread.ZOFS);
    }

    run_BK(thread.active_front,thread.active_back,thread.orphans1_back,block_TIME[block_id],thread.MAXFLOW,thread.YOFS,thread.ZOFS);

    // pthread_mutex_lock(&mutex);
    block_locked[block_id] = false;
    // pthread_mutex_unlock(&mutex);
  }

  return NULL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::compute_maxflow()
{
  for(unsigned int i=0;i<boundary_list.size();i++)
  {
    deactivate_segment(boundary_list[i].segment,_YOFS,_ZOFS);
  }

  spin_init(&spinlock);

  thread_func_wrapper_arg* wrapper_args = new thread_func_wrapper_arg[NUM_THREADS];

  for(int i=0;i<NUM_THREADS;i++)
  {
    wrapper_args[i].thisptr = (void*)this;
    wrapper_args[i].arg     = (void*)((intptr_t)i);

    threads[i]->handle = threadCreate(thread_func_wrapper,(void*)&wrapper_args[i]);
  }

  for(int i=0;i<NUM_THREADS;i++) threadJoin(threads[i]->handle);

  for(int i=0;i<NUM_THREADS;i++) MAXFLOW_TOTAL += threads[i]->MAXFLOW;

  delete[] wrapper_args;

  spin_destroy(&spinlock);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::node_id(int x,int y,int z) const
{
  return nodeId(x+1,y+1,z+1);
}

template<typename type_tcap,typename type_ncap,typename type_flow> template<typename type_arg_tcap,typename type_arg_ncap>
inline void GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::set_caps(const type_arg_tcap* cap_s,
                                                                        const type_arg_tcap* cap_t,
                                                                        const type_arg_ncap* cap_lee,
                                                                        const type_arg_ncap* cap_gee,
                                                                        const type_arg_ncap* cap_ele,
                                                                        const type_arg_ncap* cap_ege,
                                                                        const type_arg_ncap* cap_eel,
                                                                        const type_arg_ncap* cap_eeg)
{
  for(int z=1;z<_d+1;z++)
  for(int y=1;y<_h+1;y++)
  for(int x=1;x<_w+1;x++)
  {
    const int v = nodeId(x,y,z);

    int xyz = (x-1) + (y-1)*_w + (z-1)*_w*_h;

    if (cap_s[xyz] > 0 && cap_t[xyz] > 0)
    {
      if      (cap_s[xyz] > cap_t[xyz])
      {
        rc_st[v] = cap_s[xyz]-cap_t[xyz];

        MAXFLOW_TOTAL += cap_t[xyz];

        LABEL(v) = LABEL_S;

        PARENT(v) = TERMINAL;

        dist[v] = 1;
      }
      else if (cap_s[xyz] < cap_t[xyz])
      {
        rc_st[v] = cap_t[xyz]-cap_s[xyz];

        MAXFLOW_TOTAL += cap_s[xyz];

        LABEL(v) = LABEL_T;

        PARENT(v) = TERMINAL;

        dist[v] = 1;
      }
      else
      {
        rc_st[v] = 0;

        MAXFLOW_TOTAL += cap_s[xyz];

        PARENT(v) = NONE;
      }
    }
    else
    {
      if (cap_s[xyz] > 0)
      {
        LABEL(v) = LABEL_S;

        rc_st[v] = cap_s[xyz];

        PARENT(v) = TERMINAL;

        dist[v] = 1;
      }
      else if (cap_t[xyz] > 0)
      {
        LABEL(v) = LABEL_T;

        rc_st[v] = cap_t[xyz];

        PARENT(v) = TERMINAL;

        dist[v] = 1;
      }
    }

    rc_lee[v] = cap_lee[xyz]; if (cap_lee[xyz]!=0) { DESAT(LEE,v); }
    rc_gee[v] = cap_gee[xyz]; if (cap_gee[xyz]!=0) { DESAT(GEE,v); }
    rc_ele[v] = cap_ele[xyz]; if (cap_ele[xyz]!=0) { DESAT(ELE,v); }
    rc_ege[v] = cap_ege[xyz]; if (cap_ege[xyz]!=0) { DESAT(EGE,v); }
    rc_eel[v] = cap_eel[xyz]; if (cap_eel[xyz]!=0) { DESAT(EEL,v); }
    rc_eeg[v] = cap_eeg[xyz]; if (cap_eeg[xyz]!=0) { DESAT(EEG,v); }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::get_segment(int v) const
{
  static const int map_label[3] = { 0, 0, 1 };

  return map_label[LABEL(v)];
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline type_flow GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::get_flow() const
{
  return MAXFLOW_TOTAL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::bad_alloc() const
{
  return (mem_pool==NULL);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::GridGraph_3D_6C_MT(int w,int h,int d,int num_threads,int block_size) :
  _w(w),
  _h(h),
  _d(d),
  W(next_higher_mul4(w+2)),
  H(next_higher_mul4(h+2)),
  D(next_higher_mul4(d+2)),
  WB(W/4),
  WBHB((W/4)*(H/4)),
  mem_pool(NULL)
{
  const int YOFS = WB*64 - 3*4;
  const int ZOFS = WBHB*64 - 3*(4*4);

  _YOFS = YOFS;
  _ZOFS = ZOFS;

  NUM_THREADS = num_threads;

  BLOCK_WIDTH  = block_size;
  BLOCK_HEIGHT = block_size;
  BLOCK_DEPTH  = block_size;

  NUM_BLOCKS_X = (W+BLOCK_WIDTH-1)/BLOCK_WIDTH;
  NUM_BLOCKS_Y = (H+BLOCK_HEIGHT-1)/BLOCK_HEIGHT;
  NUM_BLOCKS_Z = (D+BLOCK_DEPTH-1)/BLOCK_DEPTH;

  NUM_BLOCKS = NUM_BLOCKS_X * NUM_BLOCKS_Y * NUM_BLOCKS_Z;

  NUM_X_SEGMENTS = (NUM_BLOCKS_X-1) *  NUM_BLOCKS_Y    *  NUM_BLOCKS_Z;
  NUM_Y_SEGMENTS =  NUM_BLOCKS_X    * (NUM_BLOCKS_Y-1) *  NUM_BLOCKS_Z;
  NUM_Z_SEGMENTS =  NUM_BLOCKS_X    *  NUM_BLOCKS_Y    * (NUM_BLOCKS_Z-1);

  mem_pool = (unsigned char*)calloc(W*H*D*sizeof(Label_Sat)+64+
                                    W*H*D*sizeof(unsigned char)+64+
                                    W*H*D*sizeof(int)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_ncap)+64+
                                    W*H*D*sizeof(type_tcap)+64+
                                    W*H*D*sizeof(int)+64+
                                    W*H*D*sizeof(int)+64+

                                    W*H*D*sizeof(int)+64+
                                    W*H*D*sizeof(int)+64+
                                    W*H*D*sizeof(int)+64+
                                    W*H*D*sizeof(int)+64+

                                    NUM_THREADS*(sizeof(Thread)+64)+

                                    NUM_X_SEGMENTS*BLOCK_HEIGHT*BLOCK_DEPTH*sizeof(int)+
                                    NUM_Y_SEGMENTS*BLOCK_WIDTH*BLOCK_DEPTH*sizeof(int)+
                                    NUM_Z_SEGMENTS*BLOCK_WIDTH*BLOCK_HEIGHT*sizeof(int),
                                    1);

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

  rc_lee = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  rc_gee = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  rc_ele = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  rc_ege = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  rc_eel = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);
  rc_eeg = (type_ncap*)align(pool_malloc(W*H*D*sizeof(type_ncap)+64),64);

  rc_st = (type_tcap*)align(pool_malloc(W*H*D*sizeof(type_tcap)+64),64);

  rc[LEE] = rc_lee;
  rc[GEE] = rc_gee;
  rc[ELE] = rc_ele;
  rc[EGE] = rc_ege;
  rc[EEL] = rc_eel;
  rc[EEG] = rc_eeg;

  timestamp = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  dist = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  active_next = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  orphans1_prev = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);
  orphans2_next = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  free_nodes_prev = (int*)align(pool_malloc(W*H*D*sizeof(int)+64),64);

  memset(parent,NONE,W*H*D);

  MAXFLOW_TOTAL = 0;

  boundary_list.clear();

  boundary_list.reserve(NUM_X_SEGMENTS + NUM_Y_SEGMENTS + NUM_Z_SEGMENTS);
  boundary_list_new.reserve(NUM_X_SEGMENTS + NUM_Y_SEGMENTS + NUM_Z_SEGMENTS);

  for(int block_z=0;block_z<NUM_BLOCKS_Z;block_z++)
  for(int block_y=0;block_y<NUM_BLOCKS_Y;block_y++)
  for(int block_x=0;block_x<NUM_BLOCKS_X;block_x++)
  {
    const int x_min = block_x*BLOCK_WIDTH;
    const int y_min = block_y*BLOCK_HEIGHT;
    const int z_min = block_z*BLOCK_DEPTH;

    const int x_max = minimum(x_min+BLOCK_WIDTH, W);
    const int y_max = minimum(y_min+BLOCK_HEIGHT,H);
    const int z_max = minimum(z_min+BLOCK_DEPTH, D);

    if (block_x<NUM_BLOCKS_X-1)
    {
      int* seg_node_id = (int*)pool_malloc(sizeof(int)*BLOCK_HEIGHT*BLOCK_DEPTH);

      for(int i=0, z=z_min; z<z_max; z++)
      for(int      y=y_min; y<y_max; y++, i++)
      {
        const int v = nodeId(x_max-1,y,z);

        seg_node_id[i] = v;
      }

      boundary_list.push_back(Boundary((block_x  ) + block_y*NUM_BLOCKS_X + block_z*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       (block_x+1) + block_y*NUM_BLOCKS_X + block_z*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       SEGMENT_X,(y_max-y_min)*(z_max-z_min),seg_node_id));
    }

    if (block_y<NUM_BLOCKS_Y-1)
    {
      int* seg_node_id = (int*)pool_malloc(sizeof(int)*BLOCK_WIDTH*BLOCK_DEPTH);

      for(int i=0, z=z_min; z<z_max; z++)
      for(int      x=x_min; x<x_max; x++, i++)
      {
        const int v = nodeId(x,y_max-1,z);

        seg_node_id[i] = v;
      }

      boundary_list.push_back(Boundary(block_x + (block_y  )*NUM_BLOCKS_X + block_z*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       block_x + (block_y+1)*NUM_BLOCKS_X + block_z*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       SEGMENT_Y,(x_max-x_min)*(z_max-z_min),seg_node_id));
    }

    if (block_z<NUM_BLOCKS_Z-1)
    {
      int* seg_node_id = (int*)pool_malloc(sizeof(int)*BLOCK_WIDTH*BLOCK_HEIGHT);

      for(int i=0, y=y_min; y<y_max; y++)
      for(int      x=x_min; x<x_max; x++, i++)
      {
        const int v = nodeId(x,y,z_max-1);

        seg_node_id[i] = v;
      }

      boundary_list.push_back(Boundary(block_x + block_y*NUM_BLOCKS_X + (block_z  )*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       block_x + block_y*NUM_BLOCKS_X + (block_z+1)*NUM_BLOCKS_X*NUM_BLOCKS_Y,
                                       SEGMENT_Z,(x_max-x_min)*(y_max-y_min),seg_node_id));
    }
  }

  block_locked = new bool[NUM_BLOCKS];
  block_TIME = new int[NUM_BLOCKS];

  for(int i=0;i<NUM_BLOCKS;i++)
  {
    block_locked[i] = true;
    block_TIME[i] = 0;
  }

  next_block_id = 0;

  uf_init(NUM_BLOCKS);


  threads = new Thread*[NUM_THREADS];

  for(int i=0;i<NUM_THREADS;i++)
  {
    Thread* thread = (Thread*)align(pool_malloc(sizeof(Thread)+64),64);

    thread->active_front = 0;
    thread->active_back = 0;

    thread->orphans1_back = 0;

    thread->MAXFLOW = 0;

    thread->YOFS = YOFS;
    thread->ZOFS = ZOFS;

    threads[i] = thread;
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_3D_6C_MT<type_tcap,type_ncap,type_flow>::~GridGraph_3D_6C_MT()
{
  uf_cleanup();

  delete[] block_locked;
  delete[] block_TIME;
  delete[] threads;

  if (mem_pool!=NULL) { free(mem_pool); }
}

#undef THREAD_FUNC_WRAPPER_DECL
#undef ATOMIC_INCREMENT

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
