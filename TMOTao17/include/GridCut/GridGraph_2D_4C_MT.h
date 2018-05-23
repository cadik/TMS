// Copyright (C) 2012-2015 Czech Technical University in Prague - All Rights Reserved

#ifndef GRIDGRAPH_2D_4C_MT_H_
#define GRIDGRAPH_2D_4C_MT_H_

#include <cstdlib>
#include <cstring>
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
class GridGraph_2D_4C_MT
{
public:
  GridGraph_2D_4C_MT(int width,int height,int num_threads,int block_size);
  ~GridGraph_2D_4C_MT();

  // Returns the index of a node at grid coordinates [x,y].
  // The index is used to set the capacities of node's outgoing edges
  // using the set_neighbor_cap function, to set the capacities of edges
  // between node and source/sink terminals using the set_terminal_cap
  // function, and to retrieve the segment to which the node belongs
  // after the maxflow computation using the get_segment function.
  inline int node_id(int x,int y) const;

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
                       const type_arg_ncap* cap_eg);        // [ 0,+1]


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

  enum Parent { GE=0,EG,EL,LE,NONE, TERMINAL };

  static const unsigned char SISTER_TABLE[4];

  union Label_Sat
  {
    struct
    {
    unsigned char label : 2;
    unsigned char nonsat_le : 1;
    unsigned char nonsat_el : 1;
    unsigned char nonsat_eg : 1;
    unsigned char nonsat_ge : 1;
    unsigned char dummy     : 2;
    };
    unsigned char bits;
  };

  Label_Sat* label_sat;

  type_ncap* rc_ge;
  type_ncap* rc_le;
  type_ncap* rc_eg;
  type_ncap* rc_el;

  type_tcap* rc_st;

  type_ncap* rc[4];

  unsigned char* parent;

  int* parent_id;

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

  const int W;
  const int H;

  const int WB;

  int _YOFS;

  type_flow MAXFLOW_TOTAL;

  int nodeId(unsigned int x,unsigned int y) const;

  bool grow(int&         vs,
            int&         vt,
            Parent&      st,
            int& active_front,
            int& active_back,
            const int    YOFS);

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
             const int YOFS);

  void run_BK(int& active_front,
              int& active_back,
              int& orphans1_back,
              int&         TIME,
              type_flow&   MAXFLOW,
              const int    YOFS);

  enum SegmentType { SEGMENT_HORIZONTAL, SEGMENT_VERTICAL };

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
  };

  std::vector<Boundary> boundary_list;
  std::vector<Boundary> boundary_list_new;

  int BLOCK_WIDTH;
  int BLOCK_HEIGHT;

  int NUM_BLOCKS_X;
  int NUM_BLOCKS_Y;

  int NUM_BLOCKS;

  int NUM_HORIZONTAL_SEGMENTS;
  int NUM_VERTICAL_SEGMENTS;

  bool* block_locked;
  int* block_TIME;
  int next_block_id;

  int NUM_THREADS;

  Thread** threads;

  Spinlock spinlock;

  void init_block(const int    block_id,
                  int& active_front,
                  int& active_back,
                  const int    YOFS);

  void activate_segment(const Segment& segment,
                        int& active_front,
                        int& active_back,
                        const int      YOFS);

  void deactivate_segment(const Segment& segment,
                          const int      YOFS);

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

  int next_higher_mul8(int x)
  {
    return ((x-1)/8)*8+8;
  }

#ifdef _WIN32
  #define ATOMIC_INCREMENT(X) (InterlockedIncrement((long*)(X))-1)
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

    ((GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>*)wrapper_arg->thisptr)->thread_func(wrapper_arg->arg);

    return 0;
  }
};

#define SISTER(E) (SISTER_TABLE[(E)])

#define LABEL(X)       label_sat[(X)].label
#define RC(E,V)        rc[(E)][(V)]
#define RC_SISTER(E,V) rc_sister[(E)][(V)]
#define RC_ST(V)       rc_st[(V)]

#define PARENT(X)      parent[(X)]
#define PARENT_ID(X)   parent_id[(X)]

#define TIMESTAMP(V)   timestamp[(V)]

#define NONSAT(E,V)   ( label_sat[(V)].bits & ((1<<2)<<(E)) )
#define ENSAT(E,V)      label_sat[(V)].bits &= (~((1<<2)<<(E)))
#define DESAT(E,V)      label_sat[(V)].bits |= ((1<<2)<<(E))

#define N_LE(v) ( ( (  (v)  & 0x00000007) == 0 ) ? (v)-57   : (v)-1 )
#define N_GE(v) ( ( ((~(v)) & 0x00000007) == 0 ) ? (v)+57   : (v)+1 )
#define N_EL(v) ( ( (  (v)  & 0x00000038) == 0 ) ? (v)-YOFS : (v)-8 )
#define N_EG(v) ( ( ((~(v)) & 0x00000038) == 0 ) ? (v)+YOFS : (v)+8 )

template <typename type_tcap,typename type_ncap,typename type_flow>
const unsigned char GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::SISTER_TABLE[4] = {LE,EL,EG,GE};

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::nodeId(unsigned int x,unsigned int y) const
{
  return (((x>>3)+(y>>3)*WB)<<6) + (x&7)+((y&7)<<3);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::grow(int&         vs,
                                                             int&         vt,
                                                             Parent&      st,
                                                             int& active_front,
                                                             int& active_back,
                                                             const int    YOFS)
{
  while(active_front!=0)
  {
    const int P = active_front;

    const int tree_p = LABEL(P);

    const int N_ID[4] = { N_GE(P),N_EG(P),N_EL(P),N_LE(P) };

    if (tree_p==LABEL_S)
    {
      for(int n=0;n<4;n++)
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
      for(int n=0;n<4;n++)
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
type_ncap GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::find_minrf_s(int v,type_ncap minrf) const
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
type_ncap GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::find_minrf_t(int v,type_ncap minrf) const
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
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::aug_s(int v,const type_ncap minrf,int& orphans1_back)
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
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::aug_t(int v,const type_ncap minrf,int& orphans1_back)
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
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::augment(const int vs,
                                                                const int vt,
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
int GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::find_origin(int v,const int TIME)
{
  const int start_v = v;

  while(1)
  {
    if (timestamp[v]==TIME)    { goto L1; }
    if (PARENT(v) == NONE)     { return NONE; }
    if (PARENT(v) == TERMINAL) { goto L2; }

    v = PARENT_ID(v);
  }

  L1:
    v = start_v;

    while(timestamp[v]!=TIME)
    {
      timestamp[v]=TIME;
      v = PARENT_ID(v);
    }

    return TERMINAL;

  L2:
    v = start_v;

    while(PARENT(v)!=TERMINAL)
    {
      timestamp[v]=TIME;
      v = PARENT_ID(v);
    }

    timestamp[v]=TIME;

    return TERMINAL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::adopt(int& active_front,
                                                              int& active_back,
                                                              int& orphans1_back,
                                                              const int TIME,
                                                              const int YOFS)
{
  int orphans2_front = 0;
  int orphans2_back = 0;

  int free_nodes_back = 0;

  while(!((orphans1_back==0)&&(orphans2_front==0)))
  {
    const int V = (!(orphans2_front==0)) ? ORPHANS2_POP_FRONT(orphans2_front,orphans2_back) : ORPHANS1_POP_BACK(orphans1_back);

    const int L = N_LE(V);
    const int R = N_GE(V);
    const int T = N_EL(V);
    const int B = N_EG(V);

    const int tree_p = LABEL(V);

    if      (tree_p==LABEL_S)
    {
      if (NONSAT(GE,L) && LABEL(L)==LABEL_S && find_origin(L,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = LE; PARENT_ID(V) = L; goto next; }
      if (NONSAT(LE,R) && LABEL(R)==LABEL_S && find_origin(R,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = GE; PARENT_ID(V) = R; goto next; }
      if (NONSAT(EG,T) && LABEL(T)==LABEL_S && find_origin(T,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = EL; PARENT_ID(V) = T; goto next; }
      if (NONSAT(EL,B) && LABEL(B)==LABEL_S && find_origin(B,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = EG; PARENT_ID(V) = B; goto next; }
    }
    else if (tree_p==LABEL_T)
    {
      if (NONSAT(LE,V) && LABEL(L)==LABEL_T && find_origin(L,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = LE; PARENT_ID(V) = L; goto next; }
      if (NONSAT(GE,V) && LABEL(R)==LABEL_T && find_origin(R,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = GE; PARENT_ID(V) = R; goto next; }
      if (NONSAT(EL,V) && LABEL(T)==LABEL_T && find_origin(T,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = EL; PARENT_ID(V) = T; goto next; }
      if (NONSAT(EG,V) && LABEL(B)==LABEL_T && find_origin(B,TIME)==TERMINAL) { timestamp[V]=TIME; PARENT(V) = EG; PARENT_ID(V) = B; goto next; }
    }

    LABEL(V) = LABEL_F;
    FREE_NODES_PUSH_BACK(V,free_nodes_back);

    if (LABEL(L)==tree_p && PARENT(L)==GE)  { PARENT(L) = NONE; ORPHANS2_PUSH_BACK(L,orphans2_front,orphans2_back); }
    if (LABEL(R)==tree_p && PARENT(R)==LE)  { PARENT(R) = NONE; ORPHANS2_PUSH_BACK(R,orphans2_front,orphans2_back); }
    if (LABEL(T)==tree_p && PARENT(T)==EG)  { PARENT(T) = NONE; ORPHANS2_PUSH_BACK(T,orphans2_front,orphans2_back); }
    if (LABEL(B)==tree_p && PARENT(B)==EL)  { PARENT(B) = NONE; ORPHANS2_PUSH_BACK(B,orphans2_front,orphans2_back); }

    next:
      ;
  }

  while(!(free_nodes_back==0))
  {
    const int V = FREE_NODES_POP_BACK(free_nodes_back);

    const int L = N_LE(V);
    const int R = N_GE(V);
    const int T = N_EL(V);
    const int B = N_EG(V);

    if (NONSAT(GE,L) && LABEL(L)==LABEL_S) { ACTIVE_PUSH_BACK(L,active_front,active_back); }
    if (NONSAT(LE,R) && LABEL(R)==LABEL_S) { ACTIVE_PUSH_BACK(R,active_front,active_back); }
    if (NONSAT(EG,T) && LABEL(T)==LABEL_S) { ACTIVE_PUSH_BACK(T,active_front,active_back); }
    if (NONSAT(EL,B) && LABEL(B)==LABEL_S) { ACTIVE_PUSH_BACK(B,active_front,active_back); }

    if (NONSAT(LE,V) && LABEL(L)==LABEL_T) { ACTIVE_PUSH_BACK(L,active_front,active_back); }
    if (NONSAT(GE,V) && LABEL(R)==LABEL_T) { ACTIVE_PUSH_BACK(R,active_front,active_back); }
    if (NONSAT(EL,V) && LABEL(T)==LABEL_T) { ACTIVE_PUSH_BACK(T,active_front,active_back); }
    if (NONSAT(EG,V) && LABEL(B)==LABEL_T) { ACTIVE_PUSH_BACK(B,active_front,active_back); }
  }
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::
run_BK(int& active_front,
       int& active_back,
       int& orphans1_back,
       int&         TIME,
       type_flow&   MAXFLOW,
       const int    YOFS)
{
  while(1)
  {
    int vs;
    int vt;
    Parent st;

    const bool path_found = grow(vs,vt,st,active_front,active_back,YOFS);

    if (!path_found) break;
    TIME++;

    augment(vs,vt,st,MAXFLOW,orphans1_back);

    adopt(active_front,active_back,orphans1_back,TIME,YOFS);
  }
}


template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::
init_block(const int    block_id,
           int& active_front,
           int& active_back,
           const int    YOFS)
{
  const int x_min = (block_id % NUM_BLOCKS_X)*BLOCK_WIDTH;
  const int y_min = (block_id / NUM_BLOCKS_X)*BLOCK_HEIGHT;

  const int x_max = minimum(x_min+BLOCK_WIDTH,W);
  const int y_max = minimum(y_min+BLOCK_HEIGHT,H);

  for(int y=y_min;y<y_max;y++)
  for(int x=x_min;x<x_max;x++)
  {
    const int v = nodeId(x,y);
    const int lv = LABEL(v);
    if (lv==LABEL_F) continue;

    const int L = N_LE(v);
    const int R = N_GE(v);
    const int T = N_EL(v);
    const int B = N_EG(v);

    if (lv==LABEL_S)
    {
      if (NONSAT(LE,v) && LABEL(L)!=lv) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(GE,v) && LABEL(R)!=lv) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(EL,v) && LABEL(T)!=lv) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(EG,v) && LABEL(B)!=lv) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
    }
    else
    {
      if (NONSAT(GE,L) && LABEL(L)==LABEL_F) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(LE,R) && LABEL(R)==LABEL_F) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(EG,T) && LABEL(T)==LABEL_F) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
      if (NONSAT(EL,B) && LABEL(B)==LABEL_F) { ACTIVE_PUSH_BACK(v,active_front,active_back); continue; }
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::
activate_segment(const Segment& segment,
                 int& active_front,
                 int& active_back,
                 const int      YOFS)
{
  if      (segment.type==SEGMENT_HORIZONTAL)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      if (rc_eg[v]!=0)        { DESAT(EG,v); }
      if (rc_el[N_EG(v)]!= 0) { DESAT(EL,N_EG(v)); }

      if (LABEL(v)==LABEL_S && NONSAT(EG,v)       && LABEL(N_EG(v))!=LABEL_S) ACTIVE_PUSH_BACK(v,active_front,active_back);
      if (LABEL(v)==LABEL_T && NONSAT(EL,N_EG(v)) && LABEL(N_EG(v))!=LABEL_T) ACTIVE_PUSH_BACK(v,active_front,active_back);

      if (LABEL(N_EG(v))==LABEL_S && NONSAT(EL,N_EG(v)) && LABEL(v)!=LABEL_S) ACTIVE_PUSH_BACK(N_EG(v),active_front,active_back);
      if (LABEL(N_EG(v))==LABEL_T && NONSAT(EG,v)       && LABEL(v)!=LABEL_T) ACTIVE_PUSH_BACK(N_EG(v),active_front,active_back);

    }
  }
  else if (segment.type==SEGMENT_VERTICAL)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      if (rc_ge[v]!=0)       { DESAT(GE,v); }
      if (rc_le[N_GE(v)]!=0) { DESAT(LE,N_GE(v)); }

      if (LABEL(v)==LABEL_S && NONSAT(GE,v)       && LABEL(N_GE(v))!=LABEL_S) ACTIVE_PUSH_BACK(v,active_front,active_back);
      if (LABEL(v)==LABEL_T && NONSAT(LE,N_GE(v)) && LABEL(N_GE(v))!=LABEL_T) ACTIVE_PUSH_BACK(v,active_front,active_back);

      if (LABEL(N_GE(v))==LABEL_S && NONSAT(LE,N_GE(v)) && LABEL(v)!=LABEL_S) ACTIVE_PUSH_BACK(N_GE(v),active_front,active_back);
      if (LABEL(N_GE(v))==LABEL_T && NONSAT(GE,v)       && LABEL(v)!=LABEL_T) ACTIVE_PUSH_BACK(N_GE(v),active_front,active_back);
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::
deactivate_segment(const Segment& segment,
                   const int      YOFS)
{
  if      (segment.type==SEGMENT_HORIZONTAL)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      ENSAT(EG,v);
      ENSAT(EL,N_EG(v));
    }
  }
  else if (segment.type==SEGMENT_VERTICAL)
  {
    for(int i=0;i<segment.length;i++)
    {
      const int v = segment.node_id[i];

      ENSAT(GE,v);
      ENSAT(LE,N_GE(v));
    }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void* GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::
thread_func(void* arg)
{
  const int thread_id = (intptr_t)arg;

  Thread& thread = (*threads[thread_id]);

  while(1)
  {
    const int block_id = ATOMIC_INCREMENT(&next_block_id);

    if (block_id>=NUM_BLOCKS) break;

    init_block(block_id,thread.active_front,thread.active_back,thread.YOFS);

    run_BK(thread.active_front,thread.active_back,thread.orphans1_back,block_TIME[block_id],thread.MAXFLOW,thread.YOFS);

    block_locked[block_id] = false;
  }

  std::vector<Segment> segment_list;
  segment_list.reserve(NUM_HORIZONTAL_SEGMENTS+NUM_VERTICAL_SEGMENTS);

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
      activate_segment(segment_list[i],thread.active_front,thread.active_back,thread.YOFS);
    }

    run_BK(thread.active_front,thread.active_back,thread.orphans1_back,block_TIME[block_id],thread.MAXFLOW,thread.YOFS);

    // pthread_mutex_lock(&mutex);
    block_locked[block_id] = false;
    // pthread_mutex_unlock(&mutex);
  }

  return NULL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::compute_maxflow()
{
  for(unsigned int i=0;i<boundary_list.size();i++)
  {
    deactivate_segment(boundary_list[i].segment,_YOFS);
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
inline int GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::node_id(int x,int y) const
{
  return nodeId(x+1,y+1);
}

template<typename type_tcap,typename type_ncap,typename type_flow> template<typename type_arg_tcap,typename type_arg_ncap>
inline void GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::set_caps(const type_arg_tcap* cap_s,
                                                                        const type_arg_tcap* cap_t,
                                                                        const type_arg_ncap* cap_le,
                                                                        const type_arg_ncap* cap_ge,
                                                                        const type_arg_ncap* cap_el,
                                                                        const type_arg_ncap* cap_eg)
{
  for(int y=1;y<_h+1;y++)
  for(int x=1;x<_w+1;x++)
  {
    const int v = nodeId(x,y);

    int xy = (x-1)+(y-1)*_w;

    if (cap_s[xy] > 0 && cap_t[xy] > 0)
    {
      if      (cap_s[xy] > cap_t[xy])
      {
        rc_st[v] = cap_s[xy]-cap_t[xy];

        MAXFLOW_TOTAL += cap_t[xy];

        LABEL(v) = LABEL_S;

        PARENT(v) = TERMINAL;
      }
      else if (cap_s[xy] < cap_t[xy])
      {
        rc_st[v] = cap_t[xy]-cap_s[xy];

        MAXFLOW_TOTAL += cap_s[xy];

        LABEL(v) = LABEL_T;

        PARENT(v) = TERMINAL;
      }
      else
      {
        rc_st[v] = 0;

        MAXFLOW_TOTAL += cap_s[xy];

        PARENT(v) = NONE;
      }
    }
    else
    {
      if (cap_s[xy] > 0)
      {
        LABEL(v) = LABEL_S;

        rc_st[v] = cap_s[xy];

        PARENT(v) = TERMINAL;
      }
      else if (cap_t[xy] > 0)
      {
        LABEL(v) = LABEL_T;

        rc_st[v] = cap_t[xy];

        PARENT(v) = TERMINAL;
      }
    }

    rc_le[v] = cap_le[xy]; if (cap_le[xy]!=0) { DESAT(LE,v); }
    rc_ge[v] = cap_ge[xy]; if (cap_ge[xy]!=0) { DESAT(GE,v); }
    rc_el[v] = cap_el[xy]; if (cap_el[xy]!=0) { DESAT(EL,v); }
    rc_eg[v] = cap_eg[xy]; if (cap_eg[xy]!=0) { DESAT(EG,v); }
  }
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline int GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::get_segment(int v) const
{
  static const int map_label[3] = { 0, 0, 1 };

  return map_label[LABEL(v)];
}

template <typename type_tcap,typename type_ncap,typename type_flow>
inline type_flow GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::get_flow() const
{
  return MAXFLOW_TOTAL;
}

template <typename type_tcap,typename type_ncap,typename type_flow>
bool GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::bad_alloc() const
{
  return (mem_pool==NULL);
}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::GridGraph_2D_4C_MT(int w,int h,int num_threads,int block_size) :
  _w(w),
  _h(h),
  W(next_higher_mul8(w+2)),
  H(next_higher_mul8(h+2)),
  WB(W/8),
  mem_pool(NULL)
{
  const int YOFS = (WB-1)*64+8;

  _YOFS = YOFS;

  NUM_THREADS = num_threads;

  BLOCK_WIDTH = block_size;
  BLOCK_HEIGHT = block_size;

  NUM_BLOCKS_X = (W+BLOCK_WIDTH-1)/BLOCK_WIDTH;
  NUM_BLOCKS_Y = (H+BLOCK_HEIGHT-1)/BLOCK_HEIGHT;

  NUM_BLOCKS = NUM_BLOCKS_X * NUM_BLOCKS_Y;

  NUM_HORIZONTAL_SEGMENTS = (NUM_BLOCKS_X)*(NUM_BLOCKS_Y-1);
  NUM_VERTICAL_SEGMENTS = (NUM_BLOCKS_X-1)*NUM_BLOCKS_Y;

  mem_pool = (unsigned char*)calloc(W*H*sizeof(Label_Sat)+64+
                                    W*H*sizeof(unsigned char)+64+
                                    W*H*sizeof(int)+64+
                                    W*H*sizeof(type_ncap)+64+
                                    W*H*sizeof(type_ncap)+64+
                                    W*H*sizeof(type_ncap)+64+
                                    W*H*sizeof(type_ncap)+64+
                                    W*H*sizeof(type_tcap)+64+
                                    W*H*sizeof(int)+64+
                                    W*H*sizeof(int)+64+
                                    W*H*sizeof(int)+64+
                                    W*H*sizeof(int)+64+
                                    W*H*sizeof(int)+64+
                                    NUM_THREADS*(sizeof(Thread)+64)+
                                    NUM_HORIZONTAL_SEGMENTS*BLOCK_WIDTH*sizeof(int)+
                                    NUM_VERTICAL_SEGMENTS*BLOCK_HEIGHT*sizeof(int),
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

  label_sat = (Label_Sat*)align(pool_malloc(W*H*sizeof(Label_Sat)+64),64);

  parent = (unsigned char*)align(pool_malloc(W*H*sizeof(unsigned char)+64),64);

  parent_id = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  rc_ge = (type_ncap*)align(pool_malloc(W*H*sizeof(type_ncap)+64),64);
  rc_le = (type_ncap*)align(pool_malloc(W*H*sizeof(type_ncap)+64),64);
  rc_eg = (type_ncap*)align(pool_malloc(W*H*sizeof(type_ncap)+64),64);
  rc_el = (type_ncap*)align(pool_malloc(W*H*sizeof(type_ncap)+64),64);

  rc_st = (type_tcap*)align(pool_malloc(W*H*sizeof(type_tcap)+64),64);

  rc[GE] = rc_ge;
  rc[LE] = rc_le;
  rc[EG] = rc_eg;
  rc[EL] = rc_el;

  timestamp = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  active_next = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  orphans1_prev = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);
  orphans2_next = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  free_nodes_prev = (int*)align(pool_malloc(W*H*sizeof(int)+64),64);

  memset(parent,NONE,W*H);

  MAXFLOW_TOTAL = 0;

  boundary_list.clear();

  boundary_list.reserve(NUM_HORIZONTAL_SEGMENTS+NUM_VERTICAL_SEGMENTS);
  boundary_list_new.reserve(NUM_HORIZONTAL_SEGMENTS+NUM_VERTICAL_SEGMENTS);

  for(int block_y=0;block_y<NUM_BLOCKS_Y;block_y++)
  for(int block_x=0;block_x<NUM_BLOCKS_X;block_x++)
  {
    const int x_min = block_x*BLOCK_WIDTH;
    const int y_min = block_y*BLOCK_HEIGHT;

    const int x_max = minimum(x_min+BLOCK_WIDTH,W);
    const int y_max = minimum(y_min+BLOCK_HEIGHT,H);

    if (block_x<NUM_BLOCKS_X-1)
    {
      int* seg_node_id = (int*)pool_malloc(sizeof(int)*BLOCK_HEIGHT);

      for(int y=y_min,i=0;y<y_max;y++,i++)
      {
        const int v = nodeId(x_max-1,y);

        seg_node_id[i] = v;
      }

      boundary_list.push_back(Boundary((block_x  )+block_y*NUM_BLOCKS_X,
                                       (block_x+1)+block_y*NUM_BLOCKS_X,
                                       SEGMENT_VERTICAL,y_max-y_min,seg_node_id));
    }

    if (block_y<NUM_BLOCKS_Y-1)
    {
      int* seg_node_id = (int*)pool_malloc(sizeof(int)*BLOCK_WIDTH);

      for(int x=x_min,i=0;x<x_max;x++,i++)
      {
        const int v = nodeId(x,y_max-1);

        seg_node_id[i] = v;
      }

      boundary_list.push_back(Boundary(block_x+(block_y  )*NUM_BLOCKS_X,
                                       block_x+(block_y+1)*NUM_BLOCKS_X,
                                       SEGMENT_HORIZONTAL,x_max-x_min,seg_node_id));
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

    threads[i] = thread;
  }

}

template <typename type_tcap,typename type_ncap,typename type_flow>
GridGraph_2D_4C_MT<type_tcap,type_ncap,type_flow>::~GridGraph_2D_4C_MT()
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

#undef PARENT
#undef PARENT_ID

#undef NONSAT
#undef ENSAT
#undef DESAT

#undef N_LE
#undef N_GE
#undef N_EL
#undef N_EG

#endif
