/* auxiliary routines for PCx() 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 *
 * Modified 2/27/01 to fix max steplength adjustment in ComputeStepFactor()
 *
 */ 

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "memory.h"


/*****************************************************************/
/* this file contains the following functions:                   */
/*****************************************************************/

/*
int ComputeObjInf(MMTtype *A, double *b, double *c, double *upbound, 
	      Iterate *Current, double *xibound, double *xib, 
	      double *xis, double *PriInf, double *DualInf, 
	      double *mu, double *primal_objective, 
	      double *dual_objective, double *phi, double *rmu);

int RecomputeDualVariables(LPtype *LP, solution *Solution);

int Trivial_No_Rows(LPtype *LP, solution *Solution);


int StoreHistory(solution *Solution, int Iteration, double primal_objective,
	     double dual_objective, double PriInf, double DualInf,
	     double mu,  int NumCorrections, double phi,
	     double rhs2norm, double cost2norm, int NumCols, 
	     int NumBounds);

int ComputeStepFactor(Iterate *Current, Iterate *Predictor, 
		  double *alpha_P, double *alpha_D, 
		  Parameters *Inputs);

int FreePCxMemory(double *xib, double *xibound, double *xis, 
	      double *scale, double *scaleSparse, 
	      Iterate *Current, Iterate *Predictor, 
	      Iterate *Corrector, Iterate *Copy, 
	      int *iupbound, double *upbound, 
	      double *b, MMTtype *A, MMTtype *Adense, 
	      MMTtype *Asparse, FactorType *Factor, int *maskDense);

int MaxGondzioCorrections(FactorType *Factor, double FactorTime, 
			  double SolveTime);
 
			  */
/*****************************************************************/
/* compute objectives and infeasibilities at the current iterate */
/*****************************************************************/

int
ComputeObjInf(A, b, c, upbound, Current, xibound, xib, xis,
	      PriInf, DualInf, mu, primal_objective, dual_objective, phi, rmu)
  MMTtype        *A;
  double         *b, *c, *upbound, *xibound, *xib, *xis, 
                 *PriInf, *DualInf, *mu, *primal_objective, *dual_objective;
  Iterate        *Current;
  double         *phi, *rmu;
{
  int             NumRows, NumCols, NumBounds, i, j, 
                  SparseSaxpy(), SparseSaxpyT(), 
                  SparseSaxpyM(), SparseSaxpyTM();
  double          TwoNorm2();
  double          b_norm2, c_norm2, up_norm2, max1bu, max1c, max1buc;
  double          xib_norm2, xibound_norm2, xis_norm2;

  NumRows = A->NumRows;
  NumCols = A->NumCols;
  NumBounds = Current->NumBounds;

  /* compute the residuals of the constraints              */

  /* primal constraints first        Ax-b                  */

  for (i = 0; i < NumRows; i++)
    xib[i] = -b[i];
  if (SparseSaxpyM(A, Current->x, xib)) 
     {
	printf("Error: Returned from SparseSaxpyM() with error condition\n");
	return FACTORIZE_ERROR;
     }
  /* primal bound x + w - u                                */

  for (i = 0; i < NumBounds; i++) 
     {
	j = Current->BoundIndex[i] - 1;
	xibound[i] = Current->x[j] + Current->w[i] - upbound[i];
     }

  /* Dual constraints Next        A^T pi + s - r -c        */

  for (i = 0; i < NumCols; i++)
    xis[i] = -c[i];

  if (SparseSaxpyTM(A, Current->pi, xis)) 
     {
	printf("Error: Returned from SparseSaxpyTM() with error condition\n");
	return FACTORIZE_ERROR;
     }
  for (i = 0; i < NumCols; i++)
    xis[i] += Current->s[i];

  for (i = 0; i < NumBounds; i++) 
     {
	j = Current->BoundIndex[i] - 1;
	xis[j] -= Current->r[i];
     }

  /* Compute the norms of residual vectors                            */

  xib_norm2     = TwoNorm2(xib, &NumRows);
  xibound_norm2 = TwoNorm2(xibound, &NumBounds);
  xis_norm2     = TwoNorm2(xis, &NumCols);

  b_norm2       = TwoNorm2(b, &NumRows);
  c_norm2       = TwoNorm2(c, &NumCols);
  up_norm2      = TwoNorm2(upbound, &NumBounds);

  max1bu  = MAX(1.0, sqrt(b_norm2 + up_norm2));
  max1c   = MAX(1.0, sqrt(c_norm2));
  max1buc = MAX(max1bu, max1c);

  *PriInf = sqrt(xib_norm2 + xibound_norm2) / max1bu;
  *DualInf = sqrt(xis_norm2) /  max1c;

  /* Compute gap, and primal and dual objectives                      */

  *mu = 0.0;
  for (i = 0; i < NumCols; i++)
    *mu += Current->x[i] * Current->s[i];
  for (i = 0; i < NumBounds; i++)
    *mu += Current->w[i] * Current->r[i];
  *mu = *mu / (NumCols + NumBounds);

  *primal_objective = 0.0;
  for (i = 0; i < NumCols; i++)
    *primal_objective += c[i] * Current->x[i];

  *dual_objective = 0.0;
  for (i = 0; i < NumRows; i++)
    *dual_objective += b[i] * Current->pi[i];
  for (i = 0; i < NumBounds; i++)
    *dual_objective -= Current->r[i] * upbound[i];

  /* compute ratio of infeasibility to mu */
  *rmu = MAX(*PriInf, *DualInf) / *mu;

  *phi = *PriInf + *DualInf + 
          fabs(*primal_objective - *dual_objective) / max1buc;
  return 0;
}

  /*******************************************************************/
  /* Do one more matrix-vector multiplication to recover the dual    */
  /* variables.                                                      */
  /*******************************************************************/

int
RecomputeDualVariables(LP, Solution)
     LPtype         *LP;
     solution       *Solution;
{
   
   
   int             col, NumCols, *VarType, SparseSaxpyT();
   double         *DualLower, *DualUpper;
   
   /* Transfer to local pointers                                      */
   
   NumCols = LP->Cols;
   VarType = LP->VarType;
   DualLower = Solution->DualLower;
   DualUpper = Solution->DualUpper;
   
   /* Find c - A^T pi, store in DualLower                             */
   
   for (col = 0; col < NumCols; col++)
      DualLower[col] = 0.0;
   SparseSaxpyT(LP->A, Solution->pi, Solution->DualLower);
   
   /* By checking VarType and the sign of (c - A^T pi), figure out *
    * whether it's the upper-bound or lower-bound dual variable that *
    * needs to be set to a positive value */
   
   for (col = 0; col < NumCols; col++)
      DualLower[col] = LP->c[col] - DualLower[col];
   for (col = 0; col < NumCols; col++) {
      if (VarType[col] == NORMAL) 
	 {
	    if (DualLower[col] < 0.0) 
	       DualLower[col] = 0.0;
	 } 
      else if (VarType[col] == UPPER) 
	 {
	    if (DualLower[col] >= 0.0)
	       DualUpper[col] = 0.0;
	    else if (DualLower[col] < 0.0) 
	       {
		  DualUpper[col] = -DualLower[col];
		  DualLower[col] = 0.0;
	       }
	 }
      else if (VarType[col] == FREE) 
	 DualLower[col] = 0.0;
      else 
	 printf(" What are we doing here in PCx? \n");

  }
  return 0;
}

int 
Trivial_No_Rows(LP, Solution)
     LPtype         *LP;
     solution       *Solution;
{
   int             i;
   
   for (i = 0; i < LP->Cols; i++) 
      {
	 if (LP->c[i] >= 0.0 && LP->VarType[i] != FREE) 
	    {
	       Solution->x[i] = 0.0;
	       Solution->DualLower[i] = LP->c[i];
	    } 
	 else if (LP->c[i] < 0.0 && LP->VarType[i] == UPPER) 
	    {
	       Solution->x[i] = LP->UpBound[i];
	       Solution->DualUpper[i] = -LP->c[i];
	    } 
	 else 
	    {
	       printf(" Unbounded objective detected in case \
                        of null A matrix\n");
	       Solution->Status = INFEASIBLE_SOL;
	       return 0;
	    }
      }
   Solution->Status = OPTIMAL_SOL;
   return 0;
}

int
StoreHistory(Solution, Iteration,
	     primal_objective, dual_objective,
	     PriInf, DualInf, mu, NumCorrections, phi)
     solution       *Solution;
     int             Iteration, NumCorrections;
     double          primal_objective, dual_objective, 
	             PriInf, DualInf, mu, phi;
{
   
   /* Store data for each iteration in the "history" data structure */
   
   Solution->IterationHistory[Iteration].PrimalObjective =
      primal_objective;
   Solution->IterationHistory[Iteration].DualObjective =
      dual_objective;
   Solution->IterationHistory[Iteration].PriInf =  PriInf;
   Solution->IterationHistory[Iteration].DualInf = DualInf;
   if(mu>0) 
     Solution->IterationHistory[Iteration].logmu = log10(mu);
   else
     Solution->IterationHistory[Iteration].logmu = -100.0;
   Solution->IterationHistory[Iteration].NumCorrections = NumCorrections;
   Solution->IterationHistory[Iteration].phi = phi;
   return 0;
}

int
ComputeStepFactor(Current, Predictor, alpha_P, alpha_D, Inputs)
     Iterate  *Current, *Predictor;
     double   *alpha_P, *alpha_D;
     Parameters     *Inputs;
{
   
   double   gamma_f, gamma_a;
   double   mufull, PrimalFactor, DualFactor;
   double   MaxStepPrimal, MaxStepDual;
   int      BlockPrimal,   BlockPrimalBnd, BlockDual, BlockDualBnd;
   int      NumCols, NumBounds, i;
   
   NumCols = Current->NumCols;  NumBounds = Current->NumBounds;
   
   MaxStepPrimal = 1.0;
   BlockPrimal    = -1;
   BlockPrimalBnd =  0;
   
   for (i = 0; i < Current->NumCols; i++)
      if (Current->x[i] - Predictor->x[i] < 0.0)
	 if (Current->x[i] - MaxStepPrimal * Predictor->x[i] < 0.0) 
	    {
	       MaxStepPrimal = Current->x[i] / Predictor->x[i];
	       BlockPrimal = i;
	    }
   
   for (i = 0; i < Current->NumBounds; i++)
      if (Current->w[i] - Predictor->w[i] < 0.0)
	 if (Current->w[i] - MaxStepPrimal * Predictor->w[i] < 0.0) 
	    {
	       MaxStepPrimal = Current->w[i] / Predictor->w[i];
	       BlockPrimal    = i;
	       BlockPrimalBnd = 1;  /* blocking variable is bound type */
	    }
   
   MaxStepDual = 1.0;
   BlockDual    = -1;
   BlockDualBnd =  0;
   
   for (i = 0; i < Current->NumCols; i++)
      if (Current->s[i] - Predictor->s[i] < 0.0)
	 if (Current->s[i] - MaxStepDual * Predictor->s[i] < 0.0) 
	    {
	       MaxStepDual = Current->s[i] / Predictor->s[i];
	       BlockDual = i;
	    }
   
   for (i = 0; i < Current->NumBounds; i++)
      if (Current->r[i] - Predictor->r[i] < 0.0)
	 if (Current->r[i] - MaxStepDual * Predictor->r[i] < 0.0) 
	    {
	       MaxStepDual = Current->r[i] / Predictor->r[i];
	       BlockDual    = i;
	       BlockDualBnd = 1;  /* blocking variable is bound type */
	    }
   
   /* compute scaled complementarity at full step */
   
   gamma_f = Inputs->AlphaScale;
   gamma_a = 1.0 / (1.0 - gamma_f);
   
   mufull = 0.0;
   for (i = 0; i < NumCols; i++)
      mufull += (Current->x[i] - MaxStepPrimal * Predictor->x[i]) *
	 (Current->s[i] - MaxStepDual * Predictor->s[i]);
   
   for (i = 0; i < NumBounds; i++)
      mufull += (Current->w[i] - MaxStepPrimal * Predictor->w[i]) *
	 (Current->r[i] - MaxStepDual * Predictor->r[i]);
   
   mufull /= (NumCols + NumBounds);
   mufull /= gamma_a;
   
   /* perform Mehrotra's adaptive step size procedure */
   
   /* primal step length */
   
   if (BlockPrimal == -1)
      *alpha_P = 1.0;
   else {
      if (BlockPrimalBnd) 
	 {
	    PrimalFactor = (Current->w[BlockPrimal] -
			    mufull / (Current->r[BlockPrimal] - 
			    MaxStepDual * Predictor->r[BlockPrimal])) /
	       Predictor->w[BlockPrimal];
	 }
      else 
	 {
	    PrimalFactor = (Current->x[BlockPrimal] -
			    mufull / (Current->s[BlockPrimal] - 
			      MaxStepDual * Predictor->s[BlockPrimal])) /
	       Predictor->x[BlockPrimal];
	 }
      
      PrimalFactor = MAX(PrimalFactor, gamma_f * MaxStepPrimal);
      PrimalFactor = MIN(PrimalFactor, 1.0);
      
      *alpha_P = PrimalFactor;
   }
   
   /*   dual step length */
   
   if (BlockDual == -1)
      *alpha_D = 1.0;
   else 
      {
	 if (BlockDualBnd) 
	    {
	       DualFactor = (Current->r[BlockDual] -
			     mufull / (Current->w[BlockDual] - 
			       MaxStepPrimal * Predictor->w[BlockDual])) /
		  Predictor->r[BlockDual];
	    } 
	 else 
	    {
	       DualFactor = (Current->s[BlockDual] -
			     mufull / (Current->x[BlockDual] - 
			       MaxStepPrimal * Predictor->x[BlockDual])) /
		  Predictor->s[BlockDual];
	    }
	 
	 DualFactor = MAX(DualFactor, gamma_f * MaxStepDual);
	 DualFactor = MIN(DualFactor, 1.0);
	 
	 *alpha_D = DualFactor;
      }
   return 0;
}

int
FreePCxMemory(xib, xibound, xis, scale, scaleSparse, Current, Predictor, 
		Corrector, Copy, iupbound, upbound, b, 
		A, Adense, Asparse, Factor, maskDense)
     double         *xib, *xibound, *xis, *scale, *scaleSparse, *upbound, *b;
     MMTtype        *A, *Adense, *Asparse;
     int            *iupbound, *maskDense;
     Iterate        *Current, *Predictor, *Corrector, *Copy;
     FactorType     *Factor;
     
{
   Free((char *) xib);
   Free((char *) xibound);
   Free((char *) xis);
   Free((char *) scale);
   
   Free((char *) Predictor->s);
   Free((char *) Predictor->x);
   Free((char *) Predictor->w);
   Free((char *) Predictor->r);
   Free((char *) Predictor->pi);
   /* the next is the same as iupbound!! */
   /* Free((char *) Predictor->BoundIndex); */
   Free((char *) Predictor);
   
   Free((char *) Corrector->s);
   Free((char *) Corrector->x);
   Free((char *) Corrector->w);
   Free((char *) Corrector->r);
   Free((char *) Corrector->pi);
   /* the next is the same as iupbound!! */
   /* Free((char *) Corrector->BoundIndex); */ 
   Free((char *) Corrector);
   
   Free((char *) iupbound);
   /* We need not free maskDense as it is the same as Factor->maskDense */
   Free((char *) upbound);
   
   Free((char *) Current->s);
   Free((char *) Current->x);
   Free((char *) Current->w);
   Free((char *) Current->r);
   Free((char *) Current->pi);
   /* the next is the same as iupbound!! */
   /* Free((char *) Current->BoundIndex); */
   Free((char *) Current);
   Free((char *) b);
   
   if (Copy != NULL) 
      {
	 Free((char *) Copy->s);
	 Free((char *) Copy->x);
	 Free((char *) Copy->w);
	 Free((char *) Copy->r);
	 Free((char *) Copy->pi);
   /* the next is the same as iupbound!! */
   /* Free((char *) Copy->BoundIndex); */
	 Free((char *) Copy);
      }
   Free((char *) A->pBeginRow);
   Free((char *) A->pBeginRowT);
   Free((char *) A->pEndRow);
   Free((char *) A->pEndRowT);
   Free((char *) A->Row);
   Free((char *) A->RowT);
   Free((char *) A->Value);
   Free((char *) A->ValueT);
   Free((char *) A);
   
   if ( Factor->Ndense > 0 )
      {
	 Free((char *) scaleSparse);
	 
	 Free((char *) Adense->pBeginRow);                           
	 Free((char *) Adense->pBeginRowT);                          
	 Free((char *) Adense->pEndRow);                             
	 Free((char *) Adense->pEndRowT);                            
	 Free((char *) Adense->Row);                                 
	 Free((char *) Adense->RowT);                                
	 Free((char *) Adense->Value);                               
	 Free((char *) Adense->ValueT);                               
	 Free((char *) Adense);                              
	 
	 Free((char *) Asparse->pBeginRow);
	 Free((char *) Asparse->pBeginRowT);
	 Free((char *) Asparse->pEndRow);
	 Free((char *) Asparse->pEndRowT);
	 Free((char *) Asparse->Row);
	 Free((char *) Asparse->RowT);
	 Free((char *) Asparse->Value);
	 Free((char *) Asparse->ValueT);
	 Free((char *) Asparse);
      };
   
   FreeFactorType(Factor);

   return 0;
}

int 
MaxGondzioCorrections(Factor, FactorTime, SolveTime)
     FactorType *Factor;
     double      FactorTime, SolveTime;
{
   double Ratio;
   int i, MaxCorrections;
   
   if (SolveTime == 0.0)
      Ratio = 0.0;
   else
      Ratio = FactorTime / SolveTime;
   
   if (Ratio > 50.0)
      MaxCorrections = (int) floor(Ratio / 50.0) + 2;
   else if (Ratio > 30.0)
      MaxCorrections = 2;
   else if (Ratio > 10.0)
      MaxCorrections = 1;
   else
      MaxCorrections = 0;
   
   if (MaxCorrections > 10)
      MaxCorrections = 10;
   
   return MaxCorrections;
}




