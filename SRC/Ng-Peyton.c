/* interface to the Ng/Peyton sparse Cholesky routines
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include "main.h"
#include "memory.h"
#include "solver.h"
#include "Ng-Peyton.h"
#include "rcm.h"
#include "string.h"

/*****************************************************************/
/* This is an implementation of solver.h for the Ng-Peyton Solver*/
/*****************************************************************/

#ifndef rs6000
#ifndef hpux
#define ordmmd    ordmmd_
#define sfinit    sfinit_
#define bfinit    bfinit_
#define blklvl    blklvl_
#define blkslv    blkslv_
#define blkslf    blkslf_
#define blkslb    blkslb_
#define symfct    symfct_
#define inpnv     inpnv_
#define ordnat    ordnat_
#endif
#endif

/*****************************************************************/
/* Allocation and deallocation routines for the FactorType data  */
/* structure                                                     */
/* These are specific to the solver since it might require the   */
/* use of the space in FactorType->ptr                           */
/*****************************************************************/

FactorType *
NewFactorType(A, Ndense, NumCols)
     MMTtype *A;
     int      Ndense, NumCols;
{
  int             ComputeStructureAAT(), N;
  FactorType     *FactorSpace;
  double          EndUserTime, EndSysTime, StartUserTime, StartSysTime;

  NgPeytonType   *NgPeyton, *NewNgPeytonType();

  FactorSpace = (FactorType *) Malloc(sizeof(FactorType), "FactorSpace");
  FactorSpace->AAT = (MMTtype *) Malloc(sizeof(MMTtype), "AAT");

  ComputeStructureAAT(A, FactorSpace->AAT);

  N = FactorSpace->AAT->NumCols;

  FactorSpace->N = N;
  FactorSpace->Perm = NewInt(N, "Perm");
  FactorSpace->InvPerm = NewInt(N, "InvPerm");
  FactorSpace->maskDense = NewInt(NumCols, "maskDense");
  FactorSpace->W = NewDouble2(Ndense, N, "W");
  FactorSpace->Ldense = NewDouble2(Ndense, Ndense, "Ldense");
  FactorSpace->SmallDiagonals = NO;

  NgPeyton = NewNgPeytonType(N);
  FactorSpace->ptr = NgPeyton;
  
  FactorSpace->FactorizationCode = (char *) Malloc(50*sizeof(char), 
						  "FactorizationCode");
  strcpy(FactorSpace->FactorizationCode,"Ng-Peyton sparse Cholesky library");

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
   
   /* THIS IS THE ONLY SOLVER SPECIFIC LINE */
   FreeNgPeytonType(Factor->ptr);
   
   Free((char *) Factor);
}

/*****************************************************************/
/* Allocation and deallocation routines for NgPeyton data        */
/*****************************************************************/


NgPeytonType *
NewNgPeytonType(N)
     int N;
{
   NgPeytonType *NgPeytonSpace;

   NgPeytonSpace = (NgPeytonType *) Malloc(sizeof(NgPeytonType), "NgPeyton");

   NgPeytonSpace->SuperPartitioning = NewInt(N + 1, "SuperPartitioning");
   NgPeytonSpace->mapColumnToSupernode = NewInt(N, "mapColumnToSupernode");
   NgPeytonSpace->pSuperNodeCols = NewInt(N + 1, "pSuperNodeCols");
   NgPeytonSpace->SuperNodeRows = NULL;
   NgPeytonSpace->pBeginRowL = NewInt(N + 1, "pBeginRowL");



   /* space for L will be allocated after the ordering phase */

   return NgPeytonSpace;
}


void 
FreeNgPeytonType(NgPeyton)
     NgPeytonType *NgPeyton;
{
   Free((char *) NgPeyton->SuperPartitioning);
   Free((char *) NgPeyton->mapColumnToSupernode);
   Free((char *) NgPeyton->pSuperNodeCols);
   Free((char *) NgPeyton->SuperNodeRows);
   Free((char *) NgPeyton->pBeginRowL);
   Free((char *) NgPeyton->L);  
   Free((char *) NgPeyton);
}

/*****************************************************************/
/* Ordering auxiliar functions                                   */
/*****************************************************************/

int ordinc(int dimension, int *InvPerm, int *Perm){
    int i;
    for (i = 1; i <= dimension; i++) {
        InvPerm[i - 1] = i;
        Perm[i - 1] = i;
    }
    return 0;
}

/*****************************************************************/
/* Ordering and symbolic factorization routine                   */
/*****************************************************************/

int             
Order(Factor, OrderAlg)
     FactorType     *Factor;
     int             OrderAlg;
{
   int       WorkSize, *Work, dimension, col, entry, flag;
   int      *ColumnCount, offset, row, NumCols, nonzeros;
   int      *TempBeginRow, *TempRow;
   long int  temp;

   NgPeytonType   *NgPeyton;

  /* Given A A^T (formed elsewhere), fill in pSuperNodeCols, SuperNodeRows
   * and compute multiple minimum degree ordering */

   dimension = Factor->AAT->NumCols;
   nonzeros = Factor->AAT->Nonzeros - dimension;
   
   NgPeyton = (NgPeytonType *) Factor->ptr;

   WorkSize = 4 * dimension;
   Work = NewInt(WorkSize, "Work in Order()");
   
   /* Compute sparse structure without diagonal elements */
   
   TempBeginRow = NewInt(dimension + 1, "TempBeginRow in Order()");
   TempRow = NewInt(nonzeros, "TempRow in Order()");
   
   NumCols = Factor->AAT->NumCols;
   
   for (col = 0; col < NumCols; col++)
      TempBeginRow[col] = Factor->AAT->pBeginRow[col] - col;
   
   TempBeginRow[NumCols] = Factor->AAT->pEndRow[dimension - 1] + 1 - NumCols;
   
   offset = 0;
   for (col = 0; col < NumCols; col++) 
      {
        if (Factor->AAT->pEndRow[col] -
            Factor->AAT->pBeginRow[col] + 1 == 0)
           {
              /* empty column in AAT */
              printf("There is an empty column in the matrix A A^T.\n");
              printf("This will cause an error in the Ng-Peyton matrix");
              printf(" ordering\n");
              printf("routine.  You can avoid this problem by selecting");
              printf(" the\n");
              printf("presolve option in the specifications file.\n");
              return FACTORIZE_ERROR;
           }
        for (entry = Factor->AAT->pBeginRow[col] - 1;
             entry <= Factor->AAT->pEndRow[col] - 1; entry++)
           {
              row = Factor->AAT->Row[entry];
              if (row != col + 1)
                 TempRow[entry - offset] = Factor->AAT->Row[entry];
              else
                 offset++;
           }
      }

  /* Fill in pSuperNodeCols and SuperNodeRows (as temp space) since this data
   * structure will be destroyed by ordmmd.
   *
   * (pSuperNodeCols, SuperNodeRows) - The Adjacency Structure.
   * */

   for (col = 0; col < NumCols + 1; col++)
      NgPeyton->pSuperNodeCols[col] = TempBeginRow[col];
   
   NgPeyton->SuperNodeRows = NewInt(nonzeros, "NgPeyton->SuperNodeRows");
   
   for (entry = 0; entry < nonzeros; entry++)
      NgPeyton->SuperNodeRows[entry] = TempRow[entry];
   
   /* Select and apply order algorithm. */

   switch (OrderAlg) {
       case 0:
           /* call natural ordering */

           flag = ordinc(dimension, Factor->InvPerm, Factor->Perm);

           break;
       case 1:
           /* call multiple minimum degree routine */

           ordmmd(&dimension, NgPeyton->pSuperNodeCols,
                  NgPeyton->SuperNodeRows,
                  Factor->InvPerm, Factor->Perm,
                  &WorkSize, Work, &(NgPeyton->NumCompressedCols), &flag);

           break;
       case 2:
           /* call the Reverse Cuthill-McKee routine */

           flag = rcm(&dimension, NgPeyton->pSuperNodeCols,
                  NgPeyton->SuperNodeRows,
                  Factor->InvPerm, Factor->Perm, &(NgPeyton->NumCompressedCols));

           break;
   }
   if (flag)
      printf("ordmmd error flag = %d\n", flag);

   if (flag == -1)
      {
         printf("Size of work array, %d, is larger than WorkSize\n", WorkSize);
         printf("in Order().\n");
         return FACTORIZE_ERROR;
      }

   /* Symbolic Factorization Initialization: Compute supernode partition and
    * storage requirements */
   
   WorkSize = 7 * dimension + 3;
   Work = (int *) Realloc(Work, WorkSize * sizeof(int), "Work in Order()");
   ColumnCount = NewInt(dimension, "ColumnCount");
   
   sfinit(&dimension, &nonzeros,
         TempBeginRow, TempRow,
         Factor->Perm, Factor->InvPerm, ColumnCount,
         &(Factor->NonzerosL), &(NgPeyton->NumCompressedCols),
         &(NgPeyton->NumSuperNodes), NgPeyton->mapColumnToSupernode,
         NgPeyton->SuperPartitioning,
         &WorkSize, Work, &flag);

   if (flag)
      printf("sfinit error flag = %d\n", flag);
   
   if (flag == -1) 
      {
        printf("Size of Work array, %d, is larger than WorkSize\n", WorkSize);
        printf("in Order().\n");
        return FACTORIZE_ERROR;
      }
   
   printf("Cholesky factor will have density %8.5f\n", 
                  (( 2.0 * Factor->NonzerosL - dimension ) / 
                  dimension) / dimension);

   /* allocate memory for Cholesky factor here */
   NgPeyton->L = NewDouble(Factor->NonzerosL, "L");

 
   NgPeyton->SuperNodeRows =
      (int *) Realloc(NgPeyton->SuperNodeRows,
                     Factor->NonzerosL * sizeof(int), "SuperNodeRows");
   
   /* Perform supernodal symbolic factorization */
   
   WorkSize = NgPeyton->NumSuperNodes + 2 * dimension + 1;
   Work = (int *) Realloc(Work, WorkSize * sizeof(int), "Work");
   
   symfct(&dimension, &nonzeros,
         TempBeginRow, TempRow,
         Factor->Perm, Factor->InvPerm,
         ColumnCount,
         &(NgPeyton->NumSuperNodes), NgPeyton->SuperPartitioning,
         NgPeyton->mapColumnToSupernode, &(NgPeyton->NumCompressedCols),
         NgPeyton->pSuperNodeCols, NgPeyton->SuperNodeRows,
         NgPeyton->pBeginRowL, &WorkSize, Work, &flag);
   
   if (flag)
      printf("symfct error flag = %d\n", flag);

   if (flag == -1) 
      {
        printf("Size of Work array, %d, is larger than WorkSize\n", WorkSize);
        printf("in Order().\n");
        return FACTORIZE_ERROR;
      }
   if (flag == -2) 
      {
        printf("Inconsistency in the input to symfct in Order().\n");
        return FACTORIZE_ERROR;
      }
   Free((char *) ColumnCount);
   Free((char *) Work);
   Free((char *) TempBeginRow);
   Free((char *) TempRow);

   return 0;
}

/*****************************************************************/
/*****************************************************************/

int             
EnterNumbers(Factor)
     FactorType     *Factor;
{
   int            *Work, WorkSize;
   NgPeytonType   *NgPeyton;

   NgPeyton = (NgPeytonType *) Factor->ptr;
 
   WorkSize = 2 * (Factor->AAT->NumCols) + (NgPeyton->NumSuperNodes) + 1;
   Work = NewInt(WorkSize, "Work in EnterNumbers()");
  
   inpnv(&(Factor->AAT->NumCols),
	 Factor->AAT->pBeginRow,
	 Factor->AAT->Row,
	 Factor->AAT->Value,
	 Factor->Perm,
	 Factor->InvPerm,
	 &(NgPeyton->NumSuperNodes),
	 NgPeyton->SuperPartitioning,
	 NgPeyton->pSuperNodeCols,
	 NgPeyton->SuperNodeRows,
	 NgPeyton->pBeginRowL,
	 NgPeyton->L, Work);
   
   Free((char *) Work);
   return 0;
}

/*****************************************************************/
/* Compute Cholesky factor and do the singularity handling       */
/*****************************************************************/

int             
Factorize(Factor, Inputs)
     FactorType     *Factor;
     Parameters     *Inputs;
{
   int             WorkSize, *Work, TmpSize, *Split, UnrollingLevel, CacheSize;
   int             ErrorFlag;
   double         *Tmp;
   NgPeytonType   *NgPeyton;

   EnterNumbers(Factor);
   
   NgPeyton = (NgPeytonType *) Factor->ptr;

   CacheSize = Inputs->CacheSize;
   
   Split = NewInt(Factor->AAT->NumCols, "Split");

   Factor->SmallDiagonals = NO;
   bfinit(&(Factor->AAT->NumCols),
	  &(NgPeyton->NumSuperNodes),
	  NgPeyton->SuperPartitioning,
	  NgPeyton->mapColumnToSupernode,
	  NgPeyton->pSuperNodeCols,
	  NgPeyton->SuperNodeRows, &CacheSize, &TmpSize, Split);
   
   /* allocate the amount of memory returned as TmpSize */
   
   Tmp = NewDouble(TmpSize, "Tmp");
   
   UnrollingLevel = Inputs->UnrollingLevel;
   WorkSize = 2 * (Factor->AAT->NumCols) + 2 * (NgPeyton->NumSuperNodes);
   Work = NewInt(WorkSize, "Work");
   
   blklvl(&(Factor->AAT->NumCols),
	  &(NgPeyton->NumSuperNodes),
	  NgPeyton->SuperPartitioning,
	  NgPeyton->mapColumnToSupernode,
	  Split,
	  NgPeyton->pSuperNodeCols,
	  NgPeyton->SuperNodeRows,
	  NgPeyton->pBeginRowL,
	  NgPeyton->L, &WorkSize, Work,
	  &TmpSize, Tmp, &ErrorFlag, &UnrollingLevel);
   
   switch (ErrorFlag) 
      {
      case -1:
	 Factor->SmallDiagonals = YES;
	 break;
      case -2:
	 printf("Insufficient work storage (temp) in Factorize().\n");
	 break;
      case -3:
	 printf("Insufficient work storage (iwork) in Factorize().\n");
	 break;
      }
   
   /* lstats_(&(Factor->NumSuperNodes), Factor->SuperPartitioning,
    * Factor->pSuperNodeCols, Factor->SuperNodeRows, Factor->pBeginRowL,
    * &TmpSize, &six); */
   
   Free((char *) Tmp);
   Free((char *) Work);
   Free((char *) Split);
   
   if (ErrorFlag == -2 || ErrorFlag == -3) 
      {
	 printf("blkLVL error flag = %d\n", ErrorFlag);
	 return FACTORIZE_ERROR;
      } 
   else 
      return 0;
}

/*****************************************************************/
/*****************************************************************/

int             
Solve(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   int             i;
   double         *TempRHS;
   NgPeytonType   *NgPeyton;
   
   TempRHS = NewDouble(Factor->AAT->NumCols, "TempRHS in Solve()");
      
   NgPeyton = (NgPeytonType *) Factor->ptr;

   for (i = 0; i < Factor->AAT->NumCols; i++)
      TempRHS[i] = rhs[Factor->Perm[i] - 1];
   
   blkslv(&(NgPeyton->NumSuperNodes),
	  NgPeyton->SuperPartitioning,
	  NgPeyton->pSuperNodeCols,
	  NgPeyton->SuperNodeRows,
	  NgPeyton->pBeginRowL,
	  NgPeyton->L,
	  TempRHS);
   
   for (i = 0; i < Factor->AAT->NumCols; i++)
      Solution[i] = TempRHS[Factor->InvPerm[i] - 1];
   
   Free((char *) TempRHS);
   return 0;
}


/*****************************************************************/
/* Does a forward triangular substitution only;                  */
/* Skips the final inverse permutation step.                     */
/*****************************************************************/

int             
SolveForward(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   int             i;
   double         *TempRHS;
   NgPeytonType   *NgPeyton;

   TempRHS = NewDouble(Factor->AAT->NumCols, "TempRHS in SolveForward()");
      
   NgPeyton = (NgPeytonType *) Factor->ptr;

   for (i = 0; i < Factor->AAT->NumCols; i++)
      TempRHS[i] = rhs[Factor->Perm[i] - 1];
   
   blkslf(&(NgPeyton->NumSuperNodes),
	  NgPeyton->SuperPartitioning,
	  NgPeyton->pSuperNodeCols,
	  NgPeyton->SuperNodeRows,
	  NgPeyton->pBeginRowL,
	  NgPeyton->L,
	  TempRHS);
   
   for (i = 0; i < Factor->AAT->NumCols; i++)
      Solution[i] = TempRHS[i];
   
   Free((char *) TempRHS);
   return 0;
}


/*****************************************************************/
/* Does a back triangular substitution only;                     */
/* Skips the initial permutation step, but applies the final     */
/* inverse permutation step.                                     */
/*****************************************************************/

int             
SolveBackward(Factor, rhs, Solution)
     FactorType     *Factor;
     double         *rhs, *Solution;
{
   int             i;
   double         *TempRHS;
   NgPeytonType   *NgPeyton;

   TempRHS = NewDouble(Factor->AAT->NumCols, "TempRHS in SolveBackward()");
      
   NgPeyton = (NgPeytonType *) Factor->ptr;

   for (i = 0; i < Factor->AAT->NumCols; i++)
      TempRHS[i] = rhs[i];
   
   blkslb(&(NgPeyton->NumSuperNodes),
	  NgPeyton->SuperPartitioning,
	  NgPeyton->pSuperNodeCols,
	  NgPeyton->SuperNodeRows,
	  NgPeyton->pBeginRowL,
	  NgPeyton->L,
	  TempRHS);
   
   for (i = 0; i < Factor->AAT->NumCols; i++)
      Solution[i] = TempRHS[Factor->InvPerm[i] - 1];
   
   Free((char *) TempRHS);
   return 0;
}

/****************************************************************************** 
computes W = (L^-1 * P * Adense) and 
	 the (dense) cholesky-factor of scaleDense^{-1} + W^T W
	 in the lower part of Ldense
uses routine SolveForward() for the sparse part
uses routine SparseSaxpyTM() for matrix-vector-product
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
   int          status, SolveForward();	   /* from solverfiles */
   int     SparseSaxpyTM();        	   /* from wrappers.c */
   void    StripScale();
   
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
         Free((char*) tempd);
      } /* end if ( (Factor->Ndense) > 0 ) */
   
   /* finally, do the back substitution */
   
   SolveBackward(Factor, temprhs, Solution);
   
   Free((char *) temprhs);
   return 0;	
}
