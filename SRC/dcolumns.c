/*  dense column handling routines for PCx()
 * 
 *  PCx 1.1 11/97
 *
 *  coded by Marc Wenzel, Argonne, Fall 1996.
 *  revised by Steve Wright and Joe Czyzyk, Spring, 1997.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */
         
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "main.h"


/******************************************************************************
contains the functions
******************************************************************************/
/*          
int  LookforDenseColumns(MMTtype *A, int *nonzDense, int *maskDense);

void StripA(MMTtype *A, MMTtype *Afree, int *maskDense, int label);

void StripScale(double *scale, double *scalefree, int *maskDense, 
		int NumCols, int label);

int  PreConjGrad(MMTtype *A, double *scale, FactorType *Factor, 
		 double *residual, double *pcgtol, double rhsn, 
		 int maxit, double *sol);
		 */
/******************************************************************************
  checks for dense columns in the constraint matrix A
	input:  original matrix A in MMTtype
	output: maskDense as 01-vector with the coding 0=sparse, 1=dense;
		nonzDense as the number of nonzeros in the dense part;
		function itself as number of dense columns
*****************************************************************************/ 



int		
LookforDenseColumns(A, nonzDense, maskDense)
     MMTtype	*A;
     int	*nonzDense, *maskDense;
{
  int		i, j, Ndense, Target, *nonzerosCol, LastCol, tnz;
  double	threshold;
  int           m = 0;

  threshold=1.0;
  if (A->NumRows > 500) 
     threshold=0.1;
  if (A->NumRows > 1000) 
     threshold=0.1;
  if (A->NumRows > 2000) 
     threshold=0.05;
  if (A->NumRows > 10000) 
     threshold=0.01;
  if (A->NumRows > 20000) 
     threshold=0.005;
  if (A->NumRows > 50000) /* close your eyes and pray */ 
     threshold=0.002;
  
  tnz = threshold * A->NumRows;

  Target = MIN(A->NumRows / 10, 100);
  
  Ndense = 0; *nonzDense = 0;
  nonzerosCol = NewInt(A->NumCols, "nonzerosCol");
  LastCol = A->NumCols-1;
  for (i = 0; i < A->NumCols; i++) 
    {
      if (A->pEndRow[i] - A->pBeginRow[i] + 1 > tnz) 
	 Ndense++;

      m = MAX(m, A->pEndRow[i] - A->pBeginRow[i] + 1);
    }
     
  if(Ndense == 0) 
    {
      Free((char *) nonzerosCol);
      return 0;
    }

  for (i = 0; i < A->NumCols; i++)
    nonzerosCol[i] = A->pEndRow[i] - A->pBeginRow[i] + 1;
  
  quicksort(nonzerosCol, 0, LastCol);
  
  Target = MIN(Target, Ndense); 
  Target = MIN(Target, A->NumCols-1);
  for(i=0, Ndense=0; i<Target; i++) 
    {
      if(nonzerosCol[i] >= 5*nonzerosCol[i+1])
	{
	  Ndense=i+1; tnz = nonzerosCol[i]; 
	  break;
	}
    }
  Free((char *) nonzerosCol);

  if(Ndense == 0) 
     return 0; 
  
  *nonzDense = 0; Ndense=0;
  for (i = 0; i < A->NumCols; i++)
    {
      if ( (A->pEndRow[i] - A->pBeginRow[i] + 1) >= tnz)
	{
	  *nonzDense += (A->pEndRow[i] - A->pBeginRow[i] + 1);
	  Ndense++; maskDense[i] = 1;
	}
      else  
	 maskDense[i] = 0;
    }
  
  printf("\nNumber of dense columns extracted: %d\n", Ndense);
  
  return Ndense;
}


/***************************************************************************** 
strips off the sparse (for label=1) or dense (for label=0) part of A
	input:  original matrix A in MMTtype;
                maskDense as indicator of the dense columns
                label as integer-indicator described above
	output: Asparse/Adense as the sparse/dense part of A in MMTtype
******************************************************************************/
void		
StripA(A, Afree, maskDense, label)
     MMTtype	*A, *Afree;
     int	*maskDense, label;
{
  int		i, j, k, l;
  /* AT- part in Adense not used */
  /* AT- part in Asparse later fixed */
  l = 0;
  k = 1;
  for (i = 0; i < A->NumCols; i++) 
    {
      if ( maskDense[i] == label )
	{
	  Afree->pBeginRow[l] = k;
	  for (j = A->pBeginRow[i]; j <= A->pEndRow[i]; j++)
	    {
	      Afree->Value[k-1] = A->Value[j-1];
	      Afree->Row[k-1] = A->Row[j-1];
	      k = k + 1;
	    };
	  Afree->pEndRow[l] = k-1;
	  l = l + 1;
		};
    };
  Afree->NumRows = A->NumRows;
  Afree->NumCols = l;
  Afree->Nonzeros = k-1;
}

/***************************************************************************** 
gets the part of scale corresponding to the sparse (label=0) or
dense (label=1) columns
	input:  scale as double-vector;
		maskDense as indicator of the dense columns;
		NumCols as the length of scale
                label as integer-indicator described above
	output: scaleSparse/scaleDense as the part of scale 
		corresponding to sparse/dense columns
******************************************************************************/
void		
StripScale(scale, scalefree, maskDense, NumCols, label)
     double	*scale, *scalefree;
     int	*maskDense, NumCols, label;
{
  int i, j=0;
  
  for (i = 0; i < (NumCols); i++)
    {
      if (maskDense[i] == label)
	{
	  scalefree[j] = scale[i];
	  j = j + 1;
	};
    };
}

/***************************************************************************** 
does an refinement of the solution via preconditioned conjugate gradient
method using the sparse part Asparse*Dsparse*Asparse^t as the preconditioner
uses routine Solve() for the sparse part;
uses routine SparseSaxpyTM(), SparseSaxpyM() for matrix-vector-product
	input:  the constraint matrix A for computing ADAT*p;
		the scaling-vector scale for computing ADAT*p;
		Factor with the cholesky-factor of the sparse part of A, here 
			used as the preconditioner, to be solved in every step;
		the residual residual=ADAT*sol-rhs, here (-1)*residual is used
			as the initializing of pcg;
		pcgtol as the desired accuracy of the solution;
		maxit as the maximum number of iterations allowed;
		sol as the up to here best-known solution of the system
	output: refined Solution of the system as double-vector sol;
		pcgtol as the achieved accuracy of the solution
******************************************************************************/
int
PreConjGrad(A, scale, Factor, residual, pcgtol, rhsn, maxit, sol)
     MMTtype            *A;
     FactorType         *Factor;
     double             *pcgtol, *scale, *residual, *sol, rhsn;
     int                 maxit;
{
  int           i, step;
  double        alpha, betha, rTtimesz, rTtimeszold, normr, 
                *r, *z, *p, *ADATp, *temp;
  int     status, Solve();                                /* from cholNg.c */
  int     SparseSaxpyM();                                /* from wrappers.c */
  int     SparseSaxpyTM();                               /* from wrappers.c */
  double  TwoNorm2();                                    /* from wrappers.c */
  void    PreConjGradCleanup();

  r     = NewDouble(A->NumRows,     "r in PreConjGrad()");
  z     = NewDouble(A->NumRows,     "z in PreConjGrad()");
  p     = NewDouble(A->NumRows,     "p in PreConjGrad()");
  ADATp = NewDouble(A->NumRows, "ADATp in PreConjGrad()");
  temp  = NewDouble(A->NumCols,  "temp in PreConjGrad()");
  step = 0;
  for (i = 0; i < A->NumRows; i++)
    r[i] = -residual[i];
  normr = sqrt( TwoNorm2(r, &(A->NumRows)) );
  /* preconditioned conjugate gradient method; reference:
     Golub, van Loan; Matrix Computations, 2nd edition; page 529. */
  while ( (step < maxit) && (normr > *pcgtol) )
    {
      Solve(Factor, r, z);

      step = step + 1;
      if (step == 1) {
        rTtimesz = 0.0;
        for (i = 0; i < A->NumRows; i++) {
          p[i] = z[i];
          rTtimesz = rTtimesz + r[i]*z[i];
        };
      }
      else {
        rTtimeszold = rTtimesz;
        rTtimesz = 0.0;
        for (i = 0; i < A->NumRows; i++)
          rTtimesz = rTtimesz + r[i]*z[i];
        if(fabs(rTtimeszold) < 1.e-15 * fabs(rTtimesz))
          { 
             PreConjGradCleanup(r, z, p, ADATp, temp); 
             return FACTORIZE_ERROR; 
          }
        betha = rTtimesz / rTtimeszold;
        for (i = 0; i < A->NumRows; i++)
          p[i] = z[i] + betha*p[i];
      };
      for (i = 0; i < A->NumRows; i++)
        ADATp[i] = 0.0;
      for (i = 0; i < A->NumCols; i++)
        temp[i] = 0.0;
      if(SparseSaxpyTM(A, p, temp)) 
          { 
             printf("Error: Returned from SparseSaxpyTM() with error condition\n");
             PreConjGradCleanup(r, z, p, ADATp, temp); 
             return FACTORIZE_ERROR; 
          }
      for (i = 0; i < A->NumCols; i++) temp[i] = scale[i]*temp[i];
      if(SparseSaxpyM(A, temp, ADATp)) {
        printf("Error: Returned from SparseSaxpyM() with error condition\n");
        PreConjGradCleanup(r, z, p, ADATp, temp); 
        return FACTORIZE_ERROR; 
      }
      alpha = 0.0;
      for (i = 0; i < A->NumRows; i++)
        alpha = alpha + p[i]*ADATp[i];
      if(fabs(alpha) < 1.e-15 * fabs(rTtimesz))
        { 
          PreConjGradCleanup(r, z, p, ADATp, temp); 
          return FACTORIZE_ERROR; 
        }
      alpha = rTtimesz / alpha;
      for (i = 0; i < A->NumRows; i++) {
        sol[i] = sol[i] + alpha*p[i];
        r[i] = r[i] - alpha*ADATp[i];
      };
      normr = sqrt( TwoNorm2(r, &(A->NumRows)) );
    };
  printf(" PCG reduced it to %7.1e in %d iterations\n", normr/rhsn,
step);
  if ( step >= maxit )
    printf("     (Max PCG iterations exceeded, use the latest answer)\n");
  *pcgtol = normr;

  return 0;
}

void PreConjGradCleanup(r, z, p, ADATp, temp)
char *r, *z, *p, *ADATp, *temp;
{
  Free((char *)      r);
  Free((char *)      z);
  Free((char *)      p);
  Free((char *)  ADATp);
  Free((char *)   temp);
  return;
}

/*** quicksort ***/

int
quicksort(a, l, r)
     int *a, l, r;
{
  int i, j, v;
  
  if (r > l)
    {
      v = a[r];
      i = l - 1;
      j = r;
      /* for the whole section, put everything less than the pivot (j) 
	 to the left, and greater than the pivot to the right */
      for (;;)
        {
	  while (a[++i] > v) ;
	  while (a[--j] < v) ;
	  if (i >= j) break;
	  swap(a, i, j);
        }
      swap(a, i, r);
      quicksort(a, l, i-1);        /* sort the left side  */
      quicksort(a, i+1, r);        /* sort the right side */
    }
  return 0;
}

int
swap(a, l, r)
     int *a, l, r;
{
  int i;
  
  i = a[l];
  a[l] = a[r];
  a[r] = i;
  return 0;
}
