/* Constants */
#define MAX_PPNF_LOOP 100

/* Structs */
typedef struct {
  int end;
  int *s;
} CC_trans;

typedef struct {
  int N;  // Number of connected component in the matrix.
  int *begin;  // Index where the ``i`` component begin.
  int *node;  // Nodes in some component based in ``begin``.
  int *comp;  // Component that a given node belongs.
} CC;

typedef struct {
  int r_idx;  // Matrix row index.
  int deg;  // Node's' degree.
  int dist;  // Minimal distance to root.
} CMN;

typedef struct {
  int N;  // Fixed maximum size of the queue.
  CMN *Q;  // The queue.
  int qh;  // The head of the queue.
  int qt;  // The tail of the queue.
} CMQ;

/* Functions to handle the structs */
CC_trans *create_CC_trans(int neqns);
void free_CC_trans(CC_trans *trans);
int is_CC_trans_empty(CC_trans *trans);
void add_CC_trans(CC_trans *trans, int n);
int pop_CC_trans(CC_trans *trans);
CC* create_CC(int neqns, int ncomp, int *comp);
void free_CC(CC *comps);
int Q_is_empty(CMQ *Q);
void Q_enq(CMQ *Q, CMN *o);
void Q_deq(CMQ *Q, CMN *o);

/* Auxiliar functions */
CC *fcc(int neqns, int *xadj, int *adncy);
int *calc_degree(int neqns, int *xadj);
void bubble_sort(int *adjncy, int *degree, int min, int max);
void order_by_degree(int neqns, int *xadj, int *adjncy, int *degree);
int ppnf(int neqns, int *xadj, int *adjncy, int *degree, int r);
int rls(int neqns, int *xadj, int *adjncy, int r, int *L, int *nL);

/* External functions */
int rcm(int *neqns, int *xadj, int *adjncy, int *invp, int *perm, int *nofsub);
