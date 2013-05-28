/* Reverse Cuthill-McKee
 *
 * Authors: Raniere Silva <ra092767@ime.unicamp.br>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "memory.h"
#include "solver.h"
#include "Ng-Peyton.h"
#include "string.h"
#include "rcm.h"

/** 
 * \fn create_CC
 * \brief Memory allocation for a CC struct.
 * \param neqns The size of the CC struct members.
 * \param comp The CC struct.
 */
CC* create_CC(int neqns, int ncomp, int *comp){
  int i, j;
  int end;

  CC *comps = (CC *) malloc(sizeof(CC));
  if (comps == NULL) {
    printf("\nCan't allocate memory for CC\n");
    exit(MEMORY_ERROR);
  }
  comps->N = ncomp;
  comps->begin = (int *) calloc((ncomp + 1), sizeof(int));
  if (comps->begin == NULL) {
    printf("\nCan't allocate memory for CC.begin\n");
    exit(MEMORY_ERROR);
  }
  comps->node = (int *) calloc(neqns, sizeof(int));
  if (comps->node == NULL) {
    printf("\nCan't allocate memory for CC.node\n");
    exit(MEMORY_ERROR);
  }

  /* Set up begin and node. */
  end = 0;
  comps->begin[0] = end + 1;
  for (i = 1; i <= ncomp; i++) {
    for (j = 0; j < neqns; j++) {
      if (comp[j] == i){
        comps->node[end] = j + 1;
        end++;
      }
    }
    comps->begin[i] = end + 1;
  }

  comps->comp = (int *) malloc(neqns * sizeof(int));
  memcpy(comps->comp, comp, neqns * sizeof(int));
  return comps;
}

/** 
 * \fn free_CC
 * \brief Free memory of a CC struct.
 * \param comp The CC struct.
 */
void free_CC(CC *comps){
  free(comps->begin);
  comps->begin = NULL;
  free(comps->node);
  comps->node = NULL;
  free(comps->comp);
  comps->comp = NULL;
  free(comps);
  comps = NULL;
}


/** 
 * \fn create_CC_trans
 * \brief Create CC_trans
 * \param neqns The size of the stack of CC_trans.
 */
CC_trans *create_CC_trans(int neqns){
  CC_trans *new = (CC_trans *) malloc(sizeof(CC_trans));
  new->end = -1;
  new->s = (int *) malloc(neqns * sizeof(int));
  if (new->s == NULL) {
    printf("\nCan't allocate memory for CC_trans.s\n");
    exit(MEMORY_ERROR);
  }
  return new;
}


/** 
 * \fn free_CC_trans
 * \brief Free memory of a CC_trans struct.
 * \param trans The CC_trans struct.
 */
void free_CC_trans(CC_trans *trans){
  free(trans->s);
  trans->s = NULL;
  free(trans);
  trans = NULL;
}


/** 
 * \fn is_CC_trans_empty
 * \brief Check if the CC_trans stack is empty.
 * \param trans The CC_trans struct.
 * \return 1 if the stack is empty.
 */
int is_CC_trans_empty(CC_trans *trans){
  if (trans->end == -1)
    return 1;
  else
    return 0;
}


/** 
 * \fn add_CC_trans
 * \brief Add a element in the CC_trans stack.
 * \param trans The CC_trans struct.
 * \param n The element to be add.
 */
void add_CC_trans(CC_trans *trans, int n){
  trans->end++;
  trans->s[trans->end] = n;
}


/** 
 * \fn pop_CC_trans
 * \brief Pop a element of the CC_trans stack.
 * \param trans The CC_trans struct.
 * \return The element poped from the stack.
 */
int pop_CC_trans(CC_trans *trans){
  if (!is_CC_trans_empty(trans)) {
    int p = trans->s[trans->end];
    trans->end--;
    return p;
  }
  else
    abort();
}

/** 
 * \fn Q_enq
 * \brief Check if queue is empty.
 * \param Q The queue.
 */
int Q_is_empty(CMQ *Q){
  return (Q->qh == Q->qt);
}

/** 
 * \fn Q_enq
 * \brief Enqueue operation.
 * \param Q The queue.
 * \param o The node to be add at the tail.
 */
void Q_enq(CMQ *Q, CMN *o){
  Q->Q[Q->qt].r_idx = o->r_idx;
  Q->Q[Q->qt].deg = o->deg;
  Q->Q[Q->qt].dist = o->dist;
  Q->qt = (Q->qt + 1) % (Q->N + 1);
}

/** 
 * \fn Q_deq
 * \brief Dequeue operation.
 * \param Q The queue.
 * \param o The node to store the node from the head.
 */
void Q_deq(CMQ *Q, CMN *o){
  o->r_idx = Q->Q[Q->qh].r_idx;
  o->deg = Q->Q[Q->qh].deg;
  o->dist = Q->Q[Q->qh].dist;
  Q->qh = (Q->qh + 1) % (Q->N + 1);
}

/* Macros for rcm */
#define QUEUE_IS_EMPTY(h, t) (h == t)
#define QUEUE_ENQ(n, q, h, t, v) do{ \
  q[t] = v; \
  t = (t + 1) % n; \
}while(0)
#define QUEUE_DEQ(n, q, h, t, v) do{ \
  v = q[h]; \
  h = (h + 1) % n; \
}while(0)

/** 
 * Cuthill-McKee reordering algorithm.
 * \param neqns Number of rows.
 * \param xadj is a vector of ``neqns + 1`` positions where ``xadj[i]`` inform
 * the ``adjncy`` index where the ``i`` column begin.
 * \param adjncy is a vector of ``xadj[neqns + 1]`` positions  where
 * ``adjncy[i]`` inform a row.
 * \param invp inverse permutation
 * \param perm permutation
 * \param nofsub an upper bound on the number of nonzero subscripts for the
 * compressed storage scheme.
 */
int rcm(int *neqns, int *xadj, int *adjncy, int *invp, int *perm, int *nofsub){
  int c, i, j;  /* Auxiliar variables. */
  int v;  /* vertice initial. */
  int head = 0;  /* head of the queue. */
  int tail = 0;  /* tail of the queue. */
  int neqns_half;
  int *degree;
  int *work = (int *) malloc((xadj[*neqns] - 1) * sizeof(int));
  int *already_visited = (int *) calloc(*neqns, sizeof(int));

  /* Copy adjncy. */
  memcpy(work, adjncy, (xadj[*neqns] - 1) * sizeof(int));

  /* Get the connected components. */
  CC *comp = fcc(*neqns, xadj, work);

  /* The degree of nodes is importante because it is need in RCM and PPNF. */
  degree = calc_degree(*neqns, xadj);
  /* Order the degree will improve the performance. */
  order_by_degree(*neqns, xadj, work, degree);
  
  /* Loop over every connected component. */
  i = 0;
  for (c = 0; c < comp->N; c++) {
    /* Get the first node of the `c` component. */
    v = comp->node[comp->begin[c] - 1];
    /* Use `v` to find a pseudo-peripheral node. */
    v = ppnf(comp->begin[c + 1] - comp->begin[c], xadj, work, degree, v);
    /* Here invp are used as a queue. */
    QUEUE_ENQ(*neqns, invp, head, tail, v);
    already_visited[v - 1] = 1;
    while (!QUEUE_IS_EMPTY(head, tail)) {
      QUEUE_DEQ(*neqns, invp, head, tail, v);
      perm[i] = v;
      i++;
      for (j = xadj[v - 1]; j < xadj[v]; j++) {
        if (!already_visited[work[j - 1] - 1]) {
          QUEUE_ENQ(*neqns, invp, head, tail, work[j - 1]);
          already_visited[work[j - 1] - 1] = 1;
        }
      }
    }
  }
  /* Inverse the permutation order. */
  neqns_half = floor((*neqns) / 2);
  for (i = 0; i < neqns_half; i++) {
    c = perm[i];
    perm[i] = perm[(*neqns) - (i + 1)];
    perm[(*neqns) - (i + 1)] = c;
  }
  /* Set up invp from perm. */
  for (i = 0; i < *neqns; i++) {
    invp[perm[i] - 1] = i + 1;
  }
  *nofsub = *neqns * *neqns;

/*   free(work);
 *   work = NULL;
 *   free(already_visited);
 *   already_visited = NULL;
 */

  return 0;
}

/** 
 * \fn fcc
 * \brief Find the connect component.
 * \param neqns Number of rows.
 * \param xadj list of adj
 * \param adjncy adj
 * \return CC struct;
 */
CC *fcc(int neqns, int *xadj, int *adjncy){
  int i, j, v;
  int comp_number = 0; /* Number of componets. */
  int *comp_info = (int *) calloc(neqns, sizeof(int)); /* Information about the component. */
  CC_trans *trans;

  /* Degenerate case. */
  if (xadj == NULL && adjncy == NULL) {
    comp_number = neqns;
    for (i = 1; i <= neqns; i++)
      comp_info[i - 1] = i;
  }
  else {
    trans = create_CC_trans(neqns);
    while(1) {
      /* Check if travell for all vertices. */
      if (is_CC_trans_empty(trans)) {
        for (i = 0; i < neqns; i++) {
          if (comp_info[i] == 0) {
            add_CC_trans(trans, i);
            comp_number++;
            comp_info[i] = comp_number;
            break;
          }
        }
      }
      if (is_CC_trans_empty(trans))
        break;
      /* End check if travell for all vertices. */
      while (!is_CC_trans_empty(trans)) {
        v = pop_CC_trans(trans);
        /* Check wich node `v` are connect. */
        for (i = xadj[v]; i < xadj[v + 1]; i++) {
          if (comp_info[adjncy[i - 1] - 1] == 0) {
            add_CC_trans(trans, adjncy[i - 1] - 1);
            comp_info[adjncy[i - 1] - 1] = comp_number;
          }
        }
        break;
      }
    }
  }
  return create_CC(neqns, comp_number, comp_info);
}

/** 
 * \fn calc_degree
 * \brief Calculate the degree.
 * \param neqns
 * \param xadj
 * \return
 */
int *calc_degree(int neqns, int *xadj){
  int i;
  int *degree = (int *) malloc(neqns * sizeof(int));

  /* Degenerate case */
  if (xadj == NULL) {
    for (i = 0; i < neqns; i++)
      degree[i] = 0;
  }
  else {
    for (i = 0; i < neqns; i++) {
      degree[i] = xadj[i + 1] - xadj[i];
    }
  }
  return degree;
}

/* Macros for heapsort */
#define PARENT(i) i/2
#define LCHILDREN(i) 2 * i + 1
#define RCHILDREN(i) 2 * i + 2

/** 
 * \fn buid_heap
 * \brief build a max heap (either the keys of parent nodes are always greater
 * than or equal to those of the children and the highest key is in the root
 * node)
 * \param adjncy
 * \param degree
 * \param min
 * \param max
 */
int* build_heap(int *adjncy, int *degree, int min, int max){
  int *heap;
  int heap_last;
  int i, j, aux;

  heap = (int *) calloc(max - min, sizeof(int));
  if (heap  == NULL) {
    printf("\nCan't allocate memory for heap\n");
    exit(MEMORY_ERROR);
  }

  heap_last = -1; // No elements in the heap.
  if (adjncy != NULL)
    // Loop over the elements to be put into the heap.
    for (i = min; i < max; i++) {
      heap_last++;
      j = heap_last;
      heap[heap_last] = adjncy[i];
      /* Correct heap after add new element. */
      while (j != 0) {
        if (degree[heap[j] - 1] > degree[heap[PARENT(j)] - 1]) {
          /* Swap */
          aux = heap[j];
          heap[j] = heap[PARENT(j)];
          heap[PARENT(j)] = aux;
          j = PARENT(j);
        }
        else {
          break;
        }
      }
    }

  return heap;
}

/** 
 * \fn heap_sort
 * \brief Sort part of a vector using heap sort.
 * \param adjncy
 * \param degree
 * \param min
 * \param max
 */
void heap_sort(int *adjncy, int *degree, int min, int max){
  int *heap;
  int i, j, aux;
  int heap_last;

  if (adjncy != NULL) {
    heap = build_heap(adjncy, degree, min, max);

    for (heap_last = max - min - 1; heap_last >= 0; heap_last--) {
      adjncy[min + heap_last] = heap[0];

      /* Correct heap after remove root element. */
      i = 0;
      while (RCHILDREN(i) <= heap_last) {
        if (degree[heap[LCHILDREN(i)] - 1] > degree[heap[RCHILDREN(i)] - 1]) {
          heap[i] = heap[LCHILDREN(i)];
          i = LCHILDREN(i);
        }
        else {
          heap[i] = heap[RCHILDREN(i)];
          i = RCHILDREN(i);
        }
      }
      heap[i] = heap[heap_last];
    }
  }

  free(heap);
}

/** 
 * \fn order_by_degree
 * \brief Order vector base on degree
 * \param neqns
 * \param xadj
 * \param adjncy
 */
void order_by_degree(int neqns, int *xadj, int *adjncy, int *degree){
  int i;
  for (i = 0; i < neqns; i++) {
    /* This is slow. */
    heap_sort(adjncy, degree, xadj[i] - 1, xadj[i + 1] - 1);
  }
}

/** 
 * \fn ppnf
 * \brief Return a pseudo-peripheral node.
 * \param neqns Number of rows.
 * \param xadj list of adj
 * \param adjncy adj
 * \param r The root note.
 * \return A pseudo-peripheral node.
 */
int ppnf(int neqns, int *xadj, int *adjncy, int *degree, int r){
  int it;  /* The iteration counter. */
  int bll;  /* Begin of the higher level. */
  int *L;  /* The level strucuture. */
  int nL, pnL;  /* The number of levels. */
  int x;
  int stop = 0;

  L = (int *) malloc(neqns * sizeof(int));
  if (L == NULL) {
    printf("\nCan't allocate memory for L\n");
    exit(MEMORY_ERROR);
  }

  it = 0;

  /* Build the level strucuture. */
  bll = rls(neqns, xadj, adjncy, r, L, &nL);
  while (!stop && it < MAX_PPNF_LOOP) {
    it++;
    /* Choose a node `x` of minimum degree in the higher level. */
    x = min_degree(neqns, bll, degree, L);
    pnL = nL;
    bll = rls(neqns, xadj, adjncy, x, L, &nL);
    if (nL > pnL) {
      r = x;
    }
    else {
      stop = 1;
    }
  }

//  free(L);
  return x;
}

/** 
 * \fn rls
 * \brief Build the level structure of `r`.
 * \param neqns Number of rows.
 * \param xadj list of adj
 * \param adjncy adj
 * \param r The root of the level strucuture.
 * \param L The level structure.
 * \param nL The number of level.
 */
int rls(int neqns, int *xadj, int *adjncy, int r, int *L, int *nL){
  int i, j, k, itmp; /* Auxiliar variable */
  int p; /* Counter. */
  int blp, bll; /* Position in level structure. */

  blp = 0;
  bll = 1;
  L[0] = r;
  *nL = 0;
  /* Loop to control the elements already in the level structure. */
  for (i = 1; i < neqns; ) {
    (*nL)++;
    /* Look passing through the elements in the previous level. */
    for (j = blp; j < bll; j++) {
      /* Loop passing thorugh the adjacent element of the previous level. */
      for (k = xadj[L[j] - 1]; k < xadj[L[j]]; k++) {
        /* Check if node already in the level structure. */
        for (itmp = 0; itmp < i; itmp++) {
          if (L[itmp] == adjncy[k - 1])
            break;
        }
        /* Add node to level structure if it not present. */
        if (itmp == i) {
          L[i] = adjncy[k - 1];
          i++;
        }
      }
    }
    /* Update position in level. */
    blp = bll;
    bll = i;
  }
  return blp;
}

/** 
 * \fn min_degree
 * \brief Select a node of minimum degree from the higher level of `L`.
 * \param neqns Number of rows.
 * \param bll The begin of the higher level.
 * \param degree The degree of the nodes.
 * \param L The level structure.
 */
int min_degree(int neqns, int bll, int *degree, int *L){
  int i;  /* Auxiliar variable. */
  int x;  /* The node of minimum degree from the higher level of `L`. */

  x = L[bll];  /* The initial node. */
  for (i = bll + 1; i < neqns; i++) {
    if (degree[L[i]] < degree[x])
      x = L[i];
  }
  return x;
}
