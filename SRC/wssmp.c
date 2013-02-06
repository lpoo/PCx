/* interface to the sparse Cholesky routines
 *
 * PCx 1.1 11/97
 *
 * Author: Michael Wagner (with hints from A. Gupta)
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */


#include <stdio.h>
#include "main.h"
#include "memory.h"
#include "solver.h"

/*****************************************************************/
/* This is an implementation of solver.h for the IBM wssmp Solver*/
/*****************************************************************/


/*****************************************************************/
/* Allocation and deallocation routines for the FactorType data  */
/* structure                                                     */
/* These are specific to the solver since it might require the   */
/* use of the space in FactorType->ptr (not in this case...      */
/*****************************************************************/


FactorType *
NewFactorType(A, Ndense, NumCols)
     MMTtype *A;
     int Ndense, NumCols;
{
  int             N;
  FactorType     *FactorSpace;
  double          EndUserTime, EndSysTime, StartUserTime, StartSysTime;

  FactorSpace = (FactorType *) Malloc(sizeof(FactorType), "FactorSpace");
  FactorSpace->AAT = (MMTtype *) Malloc(sizeof(MMTtype), "AAT");

  /*   GetTime(&StartUserTime, &StartSysTime); */
  ComputeStructureAAT(A, FactorSpace->AAT);
  /*   GetTime(&EndUserTime, &EndSysTime); */

  N = FactorSpace->AAT->NumCols;

  FactorSpace->N = N;
  FactorSpace->Perm = NewInt(N, "Perm");
  FactorSpace->InvPerm = NewInt(N, "InvPerm");
  FactorSpace->maskDense = NewInt(NumCols, "maskDense");
  FactorSpace->W = NewDouble2(Ndense, N, "W");
  FactorSpace->Ldense = NewDouble2(Ndense, Ndense, "Ldense");
  FactorSpace->SmallDiagonals = NO;
  FactorSpace->FactorizationCode = (char *)Malloc(25*sizeof(char), 
						  "FactorizationCode");
  strcpy(FactorSpace->FactorizationCode,"IBM's WSSMP library");

  return FactorSpace;
}

void
FreeFactorType(Factor)
     FactorType *Factor;
{
   void FreeNgPeytonType();
   
   Free((char *) Factor->Perm);
   Free((char *) Factor->InvPerm);
   Free((char *) Factor->maskDense);
   FreeDouble2(Factor->W);
   FreeDouble2(Factor->Ldense);
   Free((char *) Factor->AAT->pBeginRow);
   Free((char *) Factor->AAT->pEndRow);
   Free((char *) Factor->AAT->Row);
   Free((char *) Factor->AAT->Value);
   Free((char *) Factor->AAT);
   Free((char *) Factor->FactorizationCode);
   
   /* since we don't use Factor->ptr, there's nothing else
      to free here. */
   
   Free((char *) Factor);
}

/****************************************************************/

int             
Order(Factor)
     FactorType       *Factor;
{
   int             dimension, k, i;
   
   /* parameters for WSSMP */ 
   int            *iparm, ldb, nrhs = 1, naux = 0;
   double         *dparm, *avals, *b, *aux, *mrp;
   
   /* local variables */
   int            *pBeginRow, *Row, nonzeros, nnzt;

   /* initialize stuff */
   
   dimension = Factor->AAT->NumCols;
   nnzt = (Factor->AAT->Nonzeros + dimension)/2; /* enough for the triangle */

   pBeginRow = NewInt(dimension+1, "pBeginRow in Order");
   Row = NewInt(nnzt, "Row in Order");
   
   avals = NULL; 
   b = NULL;
   aux = NULL; 
   mrp = NULL;
   ldb = dimension;
   
   iparm = NewInt(64, "iparm in Order()"); 
   dparm = NewDouble(64, "dparm in Order()");
   
   /* extract LOWER TRIANGLE only into local variables */
   
   nonzeros = 1;
   pBeginRow[0] = 1;
   for (i = 0; i < dimension; i++)
      {
	 /* ignore all superdiagnals */
	 for (k = Factor->AAT->pBeginRow[i]; Factor->AAT->Row[k-1] < i+1; k++);
	 
	 for (; k <= Factor->AAT->pEndRow[i]; k++)
	    {
	       Row[nonzeros-1] = Factor->AAT->Row[k-1];
	       nonzeros++;
	    }
	 pBeginRow[i+1] = nonzeros;
      }
   
   /* set up input parameters for WSSMP */
   
   iparm[ 0] = 1; /* choose non-standard options */
   iparm[ 1] = 1; /* number of the starting task */
   iparm[ 2] = 2; /* numer of the ending task (symb. fact) */
   iparm[ 3] = 0; /* use CSR format for matrix storage */
   iparm[ 4] = 1; /* numbering style */
   iparm[ 5] = 0; /* max number of iterative refinement steps */
   iparm[ 6] = 0; /* some more crap about accuracy of it. refinement */
   iparm[ 7] = 0; /* use internal permutation */
   iparm[ 8] = 0; /* rhs is not permuted */
   iparm[ 9] = 0; /* no diagonal scaling */
   iparm[10] = 2; /* handle small pivots by replacing them with DPARM(21) */
   iparm[11] = 0; /* more pivoting stuff */
   iparm[12] = 0; /* output flag, I think */
   iparm[13] = 0; /* permute avals at will and save storage */
   iparm[14] = 0; /* number of columns to be factored before others */
   iparm[15] = 1; /* as recommended in 4.1.4 */
   iparm[16] = 0; /* as recommended in 4.1.4 */
   iparm[17] = 1; /* as recommended in 4.1.4 */
   iparm[18] = 0; /* as recommended in 4.1.4 */
   iparm[19] = 1; /* as recommended in 4.1.4 */
   iparm[20] = 0; /* output: no pivots < tol */
   iparm[21] = 0; /* output referring to iparm(11) */
   iparm[22] = 0; /* output: no double words required for factorization */
   iparm[23] = 0; /* output: nnz in factor */
   iparm[24] = 0; /* condition number estimate is not computed */
   iparm[25] = 0; /* for parallel version only */
   iparm[26] = 0; /* for parallel verions only */
   iparm[27] = 0; /* for parallel verions only */
   iparm[63] = 0; /* output: success of computation */
   
   dparm[ 0] = 1.0; /* unused */
   dparm[ 1] = 0.0; /* output: contains estimate for norm(A) */
   dparm[ 2] = 0.0; /* output: contains estimate for norm(inv(A)) */ 
   dparm[ 3] = 0.0; /* output: contains largest diagonal value in fact. */
   dparm[ 4] = 0.0; /* output: contains smallest diagonal value in fact. */ 
   dparm[ 5] = 0.0; /* stopping criterion (rel. err) for it. refinement */
   dparm[ 6] = 0.0; /* output: contains relative error after it. ref. */ 
   dparm[ 9] = 0.0; /* used by iparm(9) as tolerance for scaling */ 
   dparm[10] = 0.0; /* lower threshold for good diagonal value */ 
   dparm[11] = 0.0; /* another pivoting parameter */ 
   dparm[20] = 1e128; /* BIG number to replace small pivots with */
   dparm[21] = 0.0; /* used by iparm(11), another BIG number */
   dparm[22] = 0.0; /* output: number of flops required in Cholesky */
   dparm[23] = 0.0; /* output: number of flops during solves */
   dparm[63] = 0.0; /* output: value of first pivot < dparm(10) */
   
   
   /* call the solver */
   if (dimension/nonzeros > 100) {
	iparm[16] = dimension + 1;
	iparm[17] = 0;
	iparm[19] = 0;
   }
   
   wssmp(&dimension, pBeginRow, Row, 
	 avals, NULL, Factor->Perm, Factor->InvPerm, b, &ldb, &nrhs,
	 aux, &naux, mrp, iparm, dparm);
   
   /* process output */
   
   if (iparm[63])
      {
	 printf("\nwssmp exited with code %d in routine Order()\n", iparm[63]);
	 return(iparm[63]);
      }

   Factor->NonzerosL = iparm[23];
   
   /* free space */
   
   Free((char *) iparm);
   Free((char *) dparm);
   Free((char *) Row);
   Free((char *) pBeginRow);
   
   return 0;
}

/*****************************************************************/

int
Factorize(Factor, Inputs)
     FactorType *Factor;
     Parameters *Inputs;
{
   int             dimension, k, i, j, unpi;
   
   /* static parameters for WSSMP */
   int             ldb, nrhs = 1, naux = 0;
   int            *iparm, *mrp;
   double         *dparm, *b, *aux, m;
   
   /* local variables */
   int            *pBeginRow, *Row, counter, nnzt;
   double         *Value;

   /* initialize variables */
   
   dimension = Factor->AAT->NumCols;
   nnzt = (Factor->AAT->Nonzeros + dimension)/2; /* enough for the triangle */
   
   b = NULL;
   aux = NULL; 
   /* mrp = NewInt(dimension, "mrp in Factorize()"); */
   ldb = dimension;

   iparm = NewInt(64, "iparm in Factorize()"); 
   dparm = NewDouble(64, "dparm in Factorize()");
   
   pBeginRow = NewInt(dimension+1, "pBeginRow in Factorize");
   Row = NewInt(nnzt, "Row in Factorize");
   Value = NewDouble(nnzt, "Value in Factorize");   
   
   /* extract PERMUTED LOWER TRIANGLE only into local variables */
   
   counter = 1;
   pBeginRow[0] = 1;
   for (i = 0; i < dimension; i++) {	/* assume i is the new index */
	unpi = Factor->Perm[i];		/* unpi is unpermuted index for i */
	for (k = Factor->AAT->pBeginRow[unpi-1]; 
	     k <= Factor->AAT->pEndRow[unpi-1]; k++) {

		j = Factor->InvPerm[Factor->AAT->Row[k-1]-1];
					/* j is an index in column i */
		if (j > i) counter++;	/* only lower tri. j's are useful */
	}				/* can't store j's 'coz unsorted */
	pBeginRow[i+1] = counter;
   }					/* now I know start loc of all cols */
	
   for (i = 0; i < dimension; i++) {	/* assume i is the new index */
	unpi = Factor->Perm[i];		/* unpi is unpermuted index for i */
	for (k = Factor->AAT->pBeginRow[unpi-1]; 
	     k <= Factor->AAT->pEndRow[unpi-1]; k++) {

		j = Factor->InvPerm[Factor->AAT->Row[k-1]-1] - 1;
					/* i is an index in column j */
		if (i >= j) {		/* place element (i,j) in the data- */
		  counter = pBeginRow[j];	/* -structure if i >= j */
		  Row[counter-1] = i + 1;
		  Value[counter-1] = Factor->AAT->Value[k-1];
		  pBeginRow[j] = counter + 1;
		}
	}
   }

   /* restore pBeginRow and calculate maximum diagonal */
   
   m = Value[0];
   for (k = dimension - 1; k > 0; k--) {
      counter = pBeginRow[k-1];
      pBeginRow[k] = counter;
      if (m < Value[counter-1]) m = Value[counter-1];
   }
   pBeginRow[0] = 1;
   
   /* set up input parameters for WSSMP */
   
   iparm[ 0] = 1; /* choose non-standard options */
   iparm[ 1] = 3; /* number of the starting task */
   iparm[ 2] = 3; /* numer of the ending task (num. fact) */
   iparm[ 3] = 0; /* use CSR format for matrix storage */
   iparm[ 4] = 1; /* numbering style */
   iparm[ 5] = 0; /* max number of iterative refinement steps */
   iparm[ 6] = 0; /* some more crap about accuracy of it. refinement */
   iparm[ 7] = 1; /* input permuted matrix */
   iparm[ 8] = 0; /* rhs is not permuted */
   iparm[ 9] = 0; /* no diagonal scaling */
   iparm[10] = 1; /* handle small pivots by replacing them with DPARM[20] */
   iparm[11] = 0; /* more pivoting stuff */
   iparm[12] = 0; /* output flag, I think */
   iparm[13] = 0; /* some crap about internal storage */
   iparm[14] = 0; /* number of columns to be factored before others */
   iparm[15] = 1; /* as recommended in 4.1.4 */
   iparm[16] = 0; /* as recommended in 4.1.4 */
   iparm[17] = 1; /* as recommended in 4.1.4 */
   iparm[18] = 0; /* as recommended in 4.1.4 */
   iparm[19] = 1; /* as recommended in 4.1.4 */
   iparm[20] = 0; /* output: no pivots < tol */
   iparm[21] = 0; /* output referring to iparm(11) */
   iparm[22] = 0; /* output: no double words required for factorization */
   iparm[23] = 0; /* output: nnz in factor */
   iparm[24] = 0; /* condition number estimate is not computed */
   iparm[25] = 0; /* for parallel version only */
   iparm[26] = 0; /* for parallel verions only */
   iparm[27] = 0; /* for parallel verions only */
   iparm[63] = 0; /* output: success of computation */
   
   dparm[ 0] = 1.0; /* unused */
   dparm[ 1] = 0.0; /* output: contains estimate for norm(A) */
   dparm[ 2] = 0.0; /* output: contains estimate for norm(inv(A)) */ 
   dparm[ 3] = 0.0; /* output: contains largest diagonal value in fact. */
   dparm[ 4] = 0.0; /* output: contains smallest diagonal value in fact. */ 
   dparm[ 5] = 0.0; /* stopping criterion (rel. err) for it. refinement */
   dparm[ 6] = 0.0; /* output: contains relative error after it. ref. */ 
   dparm[ 9] = 1e-30; /* used by iparm(9) as tolerance for scaling */ 
   dparm[10] = m*(1e-30); /* lower threshold for good diagonal value */ 
   dparm[11] = 0.0; /* another pivoting parameter */ 
   dparm[20] = 1.0e128; /* BIG number to replace small pivots with */
   dparm[21] = 0.0; /* used by iparm(11), another BIG number */
   dparm[22] = 0.0; /* output: number of flops required in Cholesky */
   dparm[23] = 0.0; /* output: number of flops during solves */
   dparm[63] = 0.0; /* output: value of first pivot < dparm(10) */
   
   
   /* call factorization */

   /* wsmemrep(); */
   wssmp(&dimension, pBeginRow, Row, Value, NULL,
	 Factor->Perm, Factor->InvPerm, b, &ldb, &nrhs,
	 aux, &naux, NULL, iparm, dparm);

   /* process output */
  
   /* wsmemrep(); */
   if (iparm[63])
      {
	 printf("\nWSSMP returned with error %d in Factorize()\n", iparm[63]);
	 return(iparm[63]);
      }
   
   
   if (iparm[20])
      {
	 printf("found %d tiny diagonals; replaced with inf\n",  iparm[20]);
	 Factor->SmallDiagonals = YES;
      }
   else
      Factor->SmallDiagonals = NO;
   
   /* free up space */
   
   Free((char *) iparm);
   Free((char *) dparm);
   Free((char *) Row);
   Free((char *) pBeginRow);
   Free((char *) Value);
/* Free((char *) mrp); */
}

/*****************************************************************/

int
Solve(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   
   /* stuff for WSSMP */
   int            *iparm, ldb = 1, nrhs = 1, naux = 0;
   double         *dparm, *aux, *mrp, m;
   
   /* local variables */
   int            counter, dimension, k, i;
      
   /* set up and initialize */
   
   dimension = Factor->AAT->NumCols;
   aux = NULL; 
   mrp = NULL;
   ldb = dimension;

   
   iparm = NewInt(64, "iparm in Solve()"); 
   dparm = NewDouble(64, "dparm in Solve()");
      
   /* copy rhs -> Solution, since WSSMP modifies rhs-input */
      for (k=0; k < dimension; k++)
      Solution[k] = rhs[k];
   
   /* set up input parameters for WSSMP */
   
   iparm[ 0] = 1; /* choose non-standard options */
   iparm[ 1] = 4; /* number of the starting task */
   iparm[ 2] = 4; /* numer of the ending task (symb. fact) */
   iparm[ 3] = 0; /* use CSR format for matrix storage */
   iparm[ 4] = 1; /* numbering style */
   iparm[ 5] = 0; /* max number of iterative refinement steps */
   iparm[ 6] = 0; /* some more crap about accuracy of it. refinement */
   iparm[ 7] = 0; /* use given permutation */
   iparm[ 8] = 0; /* rhs is not permuted */
   iparm[ 9] = 0; /* no diagonal scaling */
   iparm[10] = 2; /* handle small pivots by replacing them with DPARM(21) */
   iparm[11] = 0; /* more pivoting stuff */
   iparm[12] = 0; /* output flag, I think */
   iparm[13] = 0; /* some crap about internal storage */
   iparm[14] = 0; /* number of columns to be factored before others */
   iparm[15] = 1; /* as recommended in 4.1.4 */
   iparm[16] = 0; /* as recommended in 4.1.4 */
   iparm[17] = 1; /* as recommended in 4.1.4 */
   iparm[18] = 0; /* as recommended in 4.1.4 */
   iparm[19] = 1; /* as recommended in 4.1.4 */
   iparm[20] = 0; /* output: no pivots < tol */
   iparm[21] = 0; /* output referring to iparm(11) */
   iparm[22] = 0; /* output: no double words required for factorization */
   iparm[23] = 0; /* output: nnz in factor */
   iparm[24] = 0; /* condition number estimate is not computed */
   iparm[25] = 0; /* for parallel version only */
   iparm[26] = 0; /* for parallel verions only */
   iparm[27] = 0; /* for parallel verions only */
   iparm[28] = 0; /* for parallel version only */
   iparm[29] = 0; /* perform both forward and backward solve */
   iparm[63] = 0; /* output: success of computation */
   
   dparm[ 0] = 1.0; /* unused */
   dparm[ 1] = 0.0; /* output: contains estimate for norm(A) */
   dparm[ 2] = 0.0; /* output: contains estimate for norm(inv(A)) */ 
   dparm[ 3] = 0.0; /* output: contains largest diagonal value in fact. */
   dparm[ 4] = 0.0; /* output: contains smallest diagonal value in fact. */ 
   dparm[ 5] = 0.0; /* stopping criterion (rel. err) for it. refinement */
   dparm[ 6] = 0.0; /* output: contains relative error after it. ref. */ 
   dparm[ 9] = 0.0; /* used by iparm(9) as tolerance for scaling */ 
   dparm[10] = 1e-30*m; /* lower threshold for good diagonal value */ 
   dparm[11] = 0.0; /* another pivoting parameter */ 
   dparm[20] = 1e128; /* BIG number to replace small pivots with */
   dparm[21] = 0.0; /* used by iparm(11), another BIG number */
   dparm[22] = 0.0; /* output: number of flops required in Cholesky */
   dparm[23] = 0.0; /* output: number of flops during solves */
   dparm[63] = 0.0; /* output: value of first pivot < dparm(10) */
   
   
   /* call solver */
   
   wssmp(&dimension, NULL, NULL, NULL, NULL,
	 Factor->Perm, Factor->InvPerm, Solution, &ldb, &nrhs,
	 aux, &naux, mrp, iparm, dparm);
   
   if (iparm[63])
      {
	 printf("\nWSSMP returned with error %d in Solve()\n", iparm[63]);
	 return(iparm[63]);
      }
   
   /* free up space */
   
   Free((char *) iparm);
   Free((char *) dparm);
   
   return 0;
}

/**********************************************************************/

int
SolveForward(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   
   /* stuff for WSSMP */
   int            *iparm, ldb = 1, nrhs = 1, naux = 0;
   double         *dparm, *aux, *mrp, m;
   
   /* local variables */
   int            counter, dimension, k, i;
      
   /* set up and initialize */
   
   dimension = Factor->AAT->NumCols;
   aux = NULL; 
   mrp = NULL;
   ldb = dimension;

   
   iparm = NewInt(64, "iparm in Solve()"); 
   dparm = NewDouble(64, "dparm in Solve()");
      
   /* copy rhs -> Solution, since WSSMP modifies rhs-input */
      for (k=0; k < dimension; k++)
      Solution[k] = rhs[k];
   
   /* set up input parameters for WSSMP */
   
   iparm[ 0] = 1; /* choose non-standard options */
   iparm[ 1] = 4; /* number of the starting task */
   iparm[ 2] = 4; /* numer of the ending task (symb. fact) */
   iparm[ 3] = 0; /* use CSR format for matrix storage */
   iparm[ 4] = 1; /* numbering style */
   iparm[ 5] = 0; /* max number of iterative refinement steps */
   iparm[ 6] = 0; /* some more crap about accuracy of it. refinement */
   iparm[ 7] = 0; /* use given permutation */
   iparm[ 8] = 0; /* rhs is not permuted */
   iparm[ 9] = 0; /* no diagonal scaling */
   iparm[10] = 2; /* handle small pivots by replacing them with DPARM(21) */
   iparm[11] = 0; /* more pivoting stuff */
   iparm[12] = 0; /* output flag, I think */
   iparm[13] = 0; /* some crap about internal storage */
   iparm[14] = 0; /* number of columns to be factored before others */
   iparm[15] = 1; /* as recommended in 4.1.4 */
   iparm[16] = 0; /* as recommended in 4.1.4 */
   iparm[17] = 1; /* as recommended in 4.1.4 */
   iparm[18] = 0; /* as recommended in 4.1.4 */
   iparm[19] = 1; /* as recommended in 4.1.4 */
   iparm[20] = 0; /* output: no pivots < tol */
   iparm[21] = 0; /* output referring to iparm(11) */
   iparm[22] = 0; /* output: no double words required for factorization */
   iparm[23] = 0; /* output: nnz in factor */
   iparm[24] = 0; /* condition number estimate is not computed */
   iparm[25] = 0; /* for parallel version only */
   iparm[26] = 0; /* for parallel verions only */
   iparm[27] = 0; /* for parallel verions only */
   iparm[28] = 0; /* for parallel version only */
   iparm[29] = 1; /* perform forward solve only! */
   iparm[63] = 0; /* output: success of computation */
   
   dparm[ 0] = 1.0; /* unused */
   dparm[ 1] = 0.0; /* output: contains estimate for norm(A) */
   dparm[ 2] = 0.0; /* output: contains estimate for norm(inv(A)) */ 
   dparm[ 3] = 0.0; /* output: contains largest diagonal value in fact. */
   dparm[ 4] = 0.0; /* output: contains smallest diagonal value in fact. */ 
   dparm[ 5] = 0.0; /* stopping criterion (rel. err) for it. refinement */
   dparm[ 6] = 0.0; /* output: contains relative error after it. ref. */ 
   dparm[ 9] = 0.0; /* used by iparm(9) as tolerance for scaling */ 
   dparm[10] = 1e-30*m; /* lower threshold for good diagonal value */ 
   dparm[11] = 0.0; /* another pivoting parameter */ 
   dparm[20] = 1e128; /* BIG number to replace small pivots with */
   dparm[21] = 0.0; /* used by iparm(11), another BIG number */
   dparm[22] = 0.0; /* output: number of flops required in Cholesky */
   dparm[23] = 0.0; /* output: number of flops during solves */
   dparm[63] = 0.0; /* output: value of first pivot < dparm(10) */
   
   
   /* call solver */
   
   wssmp(&dimension, NULL, NULL, NULL, NULL,
	 Factor->Perm, Factor->InvPerm, Solution, &ldb, &nrhs,
	 aux, &naux, mrp, iparm, dparm);
   
   if (iparm[63])
      {
	 printf("\nWSSMP returned with error %d in solveForward()\n", iparm[63]);
	 return(iparm[63]);
      }
   
   /* free up space */
   
   Free((char *) iparm);
   Free((char *) dparm);
   
   return 0;
}
/**********************************************************************/


int
SolveBackward(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   
   /* stuff for WSSMP */
   int            *iparm, ldb = 1, nrhs = 1, naux = 0;
   double         *dparm, *aux, *mrp, m;
   
   /* local variables */
   int            counter, dimension, k, i;
      
   /* set up and initialize */
   
   dimension = Factor->AAT->NumCols;
   aux = NULL; 
   mrp = NULL;
   ldb = dimension;

   
   iparm = NewInt(64, "iparm in Solve()"); 
   dparm = NewDouble(64, "dparm in Solve()");
      
   /* copy rhs -> Solution, since WSSMP modifies rhs-input */
      for (k=0; k < dimension; k++)
      Solution[k] = rhs[k];
   
   /* set up input parameters for WSSMP */
   
   iparm[ 0] = 1; /* choose non-standard options */
   iparm[ 1] = 4; /* number of the starting task */
   iparm[ 2] = 4; /* numer of the ending task (symb. fact) */
   iparm[ 3] = 0; /* use CSR format for matrix storage */
   iparm[ 4] = 1; /* numbering style */
   iparm[ 5] = 0; /* max number of iterative refinement steps */
   iparm[ 6] = 0; /* some more crap about accuracy of it. refinement */
   iparm[ 7] = 0; /* use given permutation */
   iparm[ 8] = 0; /* rhs is not permuted */
   iparm[ 9] = 0; /* no diagonal scaling */
   iparm[10] = 2; /* handle small pivots by replacing them with DPARM(21) */
   iparm[11] = 0; /* more pivoting stuff */
   iparm[12] = 0; /* output flag, I think */
   iparm[13] = 0; /* some crap about internal storage */
   iparm[14] = 0; /* number of columns to be factored before others */
   iparm[15] = 1; /* as recommended in 4.1.4 */
   iparm[16] = 0; /* as recommended in 4.1.4 */
   iparm[17] = 1; /* as recommended in 4.1.4 */
   iparm[18] = 0; /* as recommended in 4.1.4 */
   iparm[19] = 1; /* as recommended in 4.1.4 */
   iparm[20] = 0; /* output: no pivots < tol */
   iparm[21] = 0; /* output referring to iparm(11) */
   iparm[22] = 0; /* output: no double words required for factorization */
   iparm[23] = 0; /* output: nnz in factor */
   iparm[24] = 0; /* condition number estimate is not computed */
   iparm[25] = 0; /* for parallel version only */
   iparm[26] = 0; /* for parallel verions only */
   iparm[27] = 0; /* for parallel verions only */
   iparm[28] = 0; /* for parallel version only */
   iparm[29] = 2; /* perform backward solve only! */
   iparm[63] = 0; /* output: success of computation */
   
   dparm[ 0] = 1.0; /* unused */
   dparm[ 1] = 0.0; /* output: contains estimate for norm(A) */
   dparm[ 2] = 0.0; /* output: contains estimate for norm(inv(A)) */ 
   dparm[ 3] = 0.0; /* output: contains largest diagonal value in fact. */
   dparm[ 4] = 0.0; /* output: contains smallest diagonal value in fact. */ 
   dparm[ 5] = 0.0; /* stopping criterion (rel. err) for it. refinement */
   dparm[ 6] = 0.0; /* output: contains relative error after it. ref. */ 
   dparm[ 9] = 0.0; /* used by iparm(9) as tolerance for scaling */ 
   dparm[10] = 1e-30*m; /* lower threshold for good diagonal value */ 
   dparm[11] = 0.0; /* another pivoting parameter */ 
   dparm[20] = 1e128; /* BIG number to replace small pivots with */
   dparm[21] = 0.0; /* used by iparm(11), another BIG number */
   dparm[22] = 0.0; /* output: number of flops required in Cholesky */
   dparm[23] = 0.0; /* output: number of flops during solves */
   dparm[63] = 0.0; /* output: value of first pivot < dparm(10) */
   
   
   /* call solver */
   
   wssmp(&dimension, NULL, NULL, NULL, NULL,
	 Factor->Perm, Factor->InvPerm, Solution, &ldb, &nrhs,
	 aux, &naux, mrp, iparm, dparm);
   
   if (iparm[63])
      {
	 printf("\nWSSMP returned with error %d in SolveBackward()\n", iparm[63]);
	 return(iparm[63]);
      }
   
   /* free up space */
   
   Free((char *) iparm);
   Free((char *) dparm);
   
   return 0;
}

/****************************************************************************** 
computes W = (L^-1 * P * Adense) and 
	 the (dense) cholesky-factor of scaleDense^{-1} + W^T W
	 in the lower part of Ldense
uses routine SolveForward() for the sparse part
	input:  Adense as the dense part of A;
		Factor as the cholesky-factor of Asparse;
		scale as the diagonal scaling-matrix (double-vector);
		NumCols as the number of columns of the matrix A
	output: members W and Ldense of Factor
******************************************************************************/
int		
ComputeWandLdense(Adense, Factor, scale, NumCols)
     MMTtype	 *Adense;
     FactorType *Factor;
     double	 *scale;
     int	 NumCols;
{
   int		i, j, k;
   double      *temprhs, *scaleDense, temp, sqrt();
   int          status;
   void         StripScale();
   
   temprhs = NewDouble(Adense->NumRows, "temprhs in ComputeWandLdense");
   scaleDense = NewDouble(Factor->Ndense, 
			  "scaleDense in ComputeWandLdense");
   StripScale(scale, scaleDense, Factor->maskDense, NumCols, 1);
   
   for (i = 0; i < Factor->Ndense; i++) 
      {
	 for (k=0; k < (Adense->NumRows); k++)  
	    temprhs[k] = 0.0;
	 
	 for (j = ((Adense->pBeginRow)[i])-1; 
	      j <= ((Adense->pEndRow)[i])-1; j++)
	    temprhs[((Adense->Row)[j])-1] = (Adense->Value)[j];

	 status = SolveForward(Factor, temprhs, (Factor->W)[i]);
      };
   
   for (i=0; i < (Factor->Ndense); i++)
      for(j=0; j <= i; j++) 
	 {
	    temp=0.0;
	    for(k=0; k < (Adense->NumRows); k++) 
	       temp += (Factor->W)[i][k] * (Factor->W)[j][k];
	    (Factor->Ldense)[i][j] = temp;
	    (Factor->Ldense)[j][i] = temp;
	 } 
   
   for (i=0; i < (Factor->Ndense); i++) 
      (Factor->Ldense)[i][i] += 1.0 / scaleDense[i];
   
   for (i = 0; i < (Factor->Ndense); i++)
      for (k = i; k < (Factor->Ndense); k++) 
	 {
	    temp = (Factor->Ldense)[i][k];
	    for (j = (i-1); j >= 0; j--)  
	       temp -= ((Factor->Ldense)[k][j])*((Factor->Ldense)[i][j]);
	    if ( i==k ) 
	       {
		  if ( temp <= 0.0 )	
		     return FACTORIZE_ERROR;
		  else  
		     scaleDense[i] = 1.0/sqrt(temp);
	       } 
	    else
	       (Factor->Ldense)[k][i] = temp*scaleDense[i];
	 };

   for (i = 0; i < (Factor->Ndense); i++) 
      (Factor->Ldense)[i][i] = 1.0/scaleDense[i];

   Free((char *) temprhs);
   Free((char *) scaleDense);
   return 0;
}

/***************************************************************************** 
solves the psd-equation (A*scale*A^t)*Solution=rhs via Sherman-Morrison
uses routine Solve() for the sparse part
	input:  Factor with ready-to-use members L, Ldense, W;
		rhs as double-vector
	output: Solution of the system as double-vector
******************************************************************************/
int		
EnhancedSolve(Factor, rhs, Solution)
     FactorType *Factor;
     double	 *rhs, *Solution;
{
   int		i, j, k;
   double      *temprhs, *tempd, temp;
   int          status, SolveForward(), SolveBackward();
   
   /* start with the forward substitution */
   
   temprhs = NewDouble(Factor->AAT->NumCols, "temprhs in EnhancedSolve");
   SolveForward(Factor, rhs, temprhs);
   
   if ( (Factor->Ndense) > 0 ) 
      {	
	 
	 tempd = NewDouble(Factor->Ndense, "tempdense in EnhancedSolve");
	 /* compute W^T temprhs */
	 for(i=0; i < (Factor->Ndense); i++) 
	    {
	       temp = 0.0; 
	       for(k=0; k< (Factor->AAT->NumCols); k++)
		  temp += (Factor->W)[i][k] * temprhs[k];
	       tempd[i] = temp;
	    }

	 /* back- and forward-substitution with Ldense */
	 
	 for (i = 0; i < Factor->Ndense; i++) 
	    {
	       temp = 0.0;
	       for (j = 0; j < i; j++)
		  temp = temp + (Factor->Ldense)[i][j] * tempd[j];
	       tempd[i] = (tempd[i] - temp) / (Factor->Ldense)[i][i];
	    };
	 
	 for (i = (Factor->Ndense)-1; i>=0; i--) 
	    {
	       temp = 0.0;
	       for (j = i+1; j < (Factor->Ndense); j++)
		  temp = temp + (Factor->Ldense)[j][i] * tempd[j];
	       tempd[i] = (tempd[i] - temp) / (Factor->Ldense)[i][i];
	    };

	 /* modify temprhs by adding W * tempd */
	 
	 for(i=0; i<(Factor->Ndense); i++)
	    for(k=0; k < (Factor->AAT->NumCols); k++)
	       temprhs[k] -= (Factor->W)[i][k] * tempd[i];
	       
      } /* end if ( (Factor->Ndense) > 0 ) */
   
   /* finally, do the back substitution */
   
   SolveBackward(Factor, temprhs, Solution);
   
   Free((char *) temprhs);
   return 0;	
}



