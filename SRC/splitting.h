#include <math.h>
#define TRUE 1
#define FALSE 0
#define SWAP(X,Y,Z) Z= X; X= Y; Y= Z
#define STILL_ACTIVE 0

#ifdef MAIN /*verifica se MAIN ja foi definida */
int MD, MATCH, MIX, REFAC, MARK, ORDER, LU, INC, TOL, INE, TROCA, LLt;
int iii= 1;
float OPS= 0.0;
#else
extern int MD, MATCH, MIX, REFAC, MARK, ORDER, LU, INC, TOL, INE, TROCA, LLt, iii;
extern float OPS;
#endif
#define one &iii

/* Estruturas utilizadas por Daniele */
/*struct BASE { int m;
              int *col, *row, *diag;
              double *val;
};

typedef struct BASE base; */

/*struct ELEMENTAR { int n, entsE;
                   int *row, *col, *pos;
                   double *val;
};

typedef struct ELEMENTAR elementar; */

struct MATRIX { int m, n;
                unsigned long int nnulos;
                int *row, *col;
                double *val;
              };

typedef struct MATRIX matrix;

struct DECOMP { int scc, cc;
                int *diag, *comp;
                matrix *L;
              };
typedef struct DECOMP decomp;


/*acrescentei*/
typedef struct {
  matrix *At, *ADAt, *L;
  int iter,info;
  }fcc;


typedef struct {
  fcc *ptr;
  int stat, corretor, maxnz, *diag, *pos, *ant;
  double tol, tola, gap, gapa, gap0;
  double *G, *H, *dwork;
  matrix *A, *L;
  decomp *S;
int *compr;
} SplittingType;


void block(int *, int *, matrix *, int *, int *, int *, int);
int blocktri(int, matrix *, matrix *, int *, int *, int *, int *, int *, int);
void degree(matrix *, int *, int *, int *, int);
void gather(int, int *, double *, double *);
void scatter(int, int, int *, int *, double *, double *);

double dnrm2(int *, double *, int *);
double ddot(int *, double *, int *, double *, int *);

void ddott(int *, double *, double *, double *);
void ddivt(int *, double *, double *, double *);
void demy(int *, double *, double *, double *);
void dxmy(int *, double *, double *, double *);
void dxpy(int *, double *, double *, double *);
void dinv(int *, double *, double *);
void daxpyx(int *, double *, double *, int *, double *, int *);
void daxpby(int *, double *, double *, int *, double *, double *, int *);

void iiddownheap(int *, int *, double *, int, int);
void sort(int *, int *, double *, int);
void ddownheap(double *, int, int);
void dsort(double *, int);
void idownheap(int *, int, int);
void isort(int *, int);
void diownheap(double *, int *,int, int);
void sortd(double *, int *, int);
void iidownheap(int *, int *, int, int);
void sorti(int *, int *, int);
void iiidownheap(int *, int *, int *, int, int);
void sortii(int *, int *, int *, int);

void error(int);
void options(void);
void Mz(decomp *, double *, int *, int *, double *, double *, double *);
void ADAtx(matrix *, int, double *, double *, double *, double *);

/* definicao das rotinas usadas para precondicionar com fcc */
unsigned long int calculannulosAAt(matrix *, matrix *);
int calculaestruturaAAt(matrix *, matrix *, matrix *);
int calculaAAt(matrix *A, matrix *At, double *d, matrix *, int *, double *);
int precond(matrix *, matrix *, int, double *);
int solsis(matrix *, matrix *, double *, double *);
int VaRitz(double *,double *, int, int);

/* Rutinas usadas pelo jair*/
void jair1(int);

/* Rotinas usadas por Daniele */
int baseinicial(matrix *, matrix *, int *, int *, int *, int *, int, int *, int *);
//int sistema(int, int, base *, elementar *, double *); 

//int atualizacao(elementar *, matrix *, base *, int, int *, int *, int);

int redundante(int *, int *, matrix *, int *, int *, int, int *, int);


