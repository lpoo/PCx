/* Reverse Cuthill-McKee
 *
 * Authors: Raniere Silva <ra092767@ime.unicamp.br>.
 */

#include <stdio.h>
#include <stdlib.h>
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
CC* create_CC(int neqns, int ncomp, int *comp) {
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
void free_CC(CC *comps) {
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
CC_trans *create_CC_trans(int neqns) {
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
void free_CC_trans(CC_trans *trans) {
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
int is_CC_trans_empty(CC_trans *trans) {
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
void add_CC_trans(CC_trans *trans, int n) {
  trans->end++;
  trans->s[trans->end] = n;
}


/** 
 * \fn pop_CC_trans
 * \brief Pop a element of the CC_trans stack.
 * \param trans The CC_trans struct.
 * \return The element poped from the stack.
 */
int pop_CC_trans(CC_trans *trans) {
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
int Q_is_empty(CMQ *Q) {
  return (Q->qh == Q->qt);
}

/** 
 * \fn Q_enq
 * \brief Enqueue operation.
 * \param Q The queue.
 * \param o The node to be add at the tail.
 */
void Q_enq(CMQ *Q, CMN *o) {
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
void Q_deq(CMQ *Q, CMN *o) {
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
 */
int rcm(int *neqns, int *xadj, int *adjncy, int *invp, int *perm) {
  int c, i, j;
  int v; /* vertice initial */
  int head = 0; /* head of the queue */
  int tail = 0; /* tail of the queue */
  int *already_visited = (int *) calloc(*neqns, sizeof(int));
  /* Get the connected components. */
  CC *comp = fcc(*neqns, xadj, adjncy);
  order_by_degree(*neqns, xadj, adjncy);
  /* Loop over every connected component. */
  for (c = 0; c < comp->N; c++) {
    /* Here invp are used as a queue. */
    i = 0;
    v = comp->node[comp->begin[c] - 1];
    QUEUE_ENQ(*neqns, invp, head, tail, v);
    already_visited[v - 1] = 1;
    while (!QUEUE_IS_EMPTY(head, tail)) {
      QUEUE_DEQ(*neqns, invp, head, tail, v);
      perm[i] = v;
      i++;
      for (j = xadj[v - 1]; j < xadj[v]; j++) {
        if (!already_visited[adjncy[j - 1] - 1]) {
          QUEUE_ENQ(*neqns, invp, head, tail, adjncy[j - 1]);
          already_visited[adjncy[j - 1] - 1] = 1;
        }
      }
    }
  }
  /* Set up invp from perm. */
  return -1;
}

/** 
 * \fn fcc
 * \brief Find the connect component.
 * \param neqns Number of rows.
 * \param xadj list of adj
 * \param adjncy adj
 * \return CC struct;
 */
CC *fcc(int neqns, int *xadj, int *adjncy) {
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
int *calc_degree(int neqns, int *xadj) {
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


/** 
 * \fn bubble_sort
 * \brief Sort part of a vector using bubble sort.
 * \param adjncy
 * \param degree
 * \param min
 * \param max
 */
void bubble_sort(int *adjncy, int *degree, int min, int max) {
  int i, j, aux;
  if (adjncy != NULL)
    for (i = min; i < max; i++) {
      for (j = i; j < max - 1; j++) {
        if (degree[adjncy[j] - 1] > degree[adjncy[j + 1] - 1]) {
          aux = adjncy[j];
          adjncy[i] = adjncy[j + 1];
          adjncy[j + 1] = aux;
        }
      }
    }
}


/** 
 * \fn order_by_degree
 * \brief Order vector base on degree
 * \param neqns
 * \param xadj
 * \param adjncy
 */
void order_by_degree(int neqns, int *xadj, int *adjncy) {
  int i, j;
  int *degree = calc_degree(neqns, xadj);
  for (i = 0; i < neqns; i++) {
    /* This is slow. */
    bubble_sort(adjncy, degree, xadj[i], xadj[i + 1]);
  }
}

/** 
 * \fn ppnf
 * \brief Return a pseudo-peripheral node.
 * \param neqns Number of rows.
 * \param xadj list of adj
 * \param adjncy adj
 * \return A pseudo-peripheral node.
 */
int ppnf(int neqns, int *xadj, int *adjncy){
}

/*
rls(){
}

adj_list(){
}

adj(){
}

*/


/** 
 * \fn min_degree
 * \brief Select a node of minimum degree from the list `l`.
 * \param A non-empty list of nodes.
 * \param A node of minimum degree from the list `l`.
 */
int min_degree(){
  return 1;
}

