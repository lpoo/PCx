/* PCx() 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "memory.h"

/********************************************************************
 *                                                                  *
 * PCx() solves linear programming problems by using Mehrotra's     *
 * predictor-corrector algorithm, a primal-dual interior-point      *
 * method. It assumes that the linear program is stored in a data   *
 * structure of type "LPtype". The formulation used in the "LPtype" *
 * structure is as follows:                                         *
 *                                                                  *
 *    Primal:    min  c^T x                                         *
 * (pi)          s.t. A x   = b          x - primal variable        *
 * (s, r)               0 <= x <= u      "upper"  variables         *
 * (s)                  0 <= x           "normal" variables         *
 *                           x free      "free"   variables         *
 *                                                                  *
 *    Dual:      max  b^T pi - r^T u        pi - dual variable      *
 * ("upper"  x)  s.t. A^T pi + s - r = c     r - dual bound slack   *
 * ("normal" x)  s.t. A^T pi + s     = c     s - dual slack         *
 * ("free"   x)  s.t. A^T pi         = c                            *
 *                    pi free                                       *
 *                    (r,s) >= 0                                    *
 *                                                                  *
 *                                                                  *
 * The input structure "Inputs" contains various algorithmic        *
 * parameters. If PCx() is invoked through the main program,        *
 * default values are assigned to these parameters automatically.   *
 * The user can change them by establishing a "specifications"      *
 * file. See the User Manual for details.                           *
 *                                                                  *
 *    This version of PCx assumes that there are NO FREE VARIABLES, *
 *    since the Cholesky-based algorithm cannot handle              *
 *    them. If PCx() is called from main(), this issue does not     *
 *    arise, since main() calls SplitFreeVars(), since the free     *
 *    variables and split and recombined transparently to the user. *
 *                                                                  *
 ********************************************************************/

int             
PCx(LP, Solution, Inputs)
     LPtype         *LP;
     solution       *Solution;
     Parameters     *Inputs;
{

   FactorType     *Factor, *NewFactorType();
   
   double         *b, *xib, *xis, *xibound, *scale, *scaleSparse,
                   ComputeCentering();
   
   Iterate        *NewIterate(), *Current = NULL,	/* Current Solution */
                  *Predictor = NULL,	/* Predictor step search direction */
                  *Corrector = NULL,	/* Corrector step search direction */
                  *Copy = NULL;	        /* Copy of current solution */

   int             algorithm_loop_flag,	/* While ON, continue iterating */
                   Iteration,	/* current iteration number     */
                   NumRows, NumCols, NumBounds, status, digits, 
                   i, irow, im, ip, nonzDense, Ndense, *maskDense;

   double          rhs2norm, cost2norm, PriInf, DualInf, TwoNorm2(), 
                   mu, sigma,
                   primal_objective, dual_objective, 
                   alpha_P, alpha_D, scaletemp,
                   PriFeasTol, DualFeasTol, OptTol;

   MMTtype        *A, *Asparse, *Adense, *NewMMTtype(), *copyMMT();
  
   int             Factorize(), Permute(), InversePermute(), 
                   InitialPoint(), Trivial_No_Rows(), ComputeWandLdense(),
                   LookforDenseColumns();
   void            StoreHistory(), StripA(), StripScale();

   int            *iupbound, k;	/* this storage is needed because of
				 * incompatibilities between the structs
				 * defined in main.h and the conventions used
				 * in this code */
   double         *upbound;
   double          phi, *min_phi, rmu_ratio, rmu_ratio0;
   int             BestIteration=0;  /* best solution, in Copy */

   double          SolveADATTime=0.0, FactorizationTime=0.0;
   double          InitTime=0.0;
   double          FormADATtime=0.0, PredictorTime=0.0;
   double          CorrectorTime=0.0, LoopTime=0.0;
   double          time1 = 0.0, time2 = 0.0;
 
   /* needed for Gondzio corrections */
   double          FactorTime, SolveTime;
   double          UserTime, SysTime, OldUserTime, OldSysTime;
   double          UserTime2, SysTime2, OldUserTime2, OldSysTime2;
   int             MaxCorrections, NumCorrections=0, MaxGondzioCorrections();

   
   GetTime(&OldUserTime, &OldSysTime);
   
   /* Move the dimensions to local storage */
   NumRows = LP->Rows;
   NumCols = LP->Cols;
   NumBounds = LP->NumberBounds;

   /* Trap a trivial case up front */
   if (NumRows == 0) 
      {
	 Trivial_No_Rows(LP, Solution); 
	 return 0;
      }

  /*******************************************************************/
  /* INITIALIZE - assign space, put the parameters in their place    */
  /*******************************************************************/

  /* Copy the matrix and rhs into more congenial data structures */

   A = copyMMT(LP->A, NumRows, NumCols);
   b = NewDouble(NumRows, "b");

   for (i = 0; i < NumRows; i++)
      b[i] = LP->b[i];

   iupbound = NewInt(NumBounds, "iupbound");
   upbound  = NewDouble(NumBounds, "upbound");
   for (i = 0; i < NumBounds; i++) 
      {
	 k = LP->BoundIndex[i];
	 iupbound[i] = k + 1;
	 upbound[i] = LP->UpBound[k];
      }

   PriFeasTol = Inputs->PriFeasTol;
   DualFeasTol = Inputs->DualFeasTol;
   OptTol = Inputs->OptTol;
   
   Solution->Status = UNKNOWN_SOL;    /* default */
   Solution->Iterations = 0;
   
   xib = NewDouble(NumRows, "xib");
   xibound = NewDouble(NumBounds, "xibound");
   xis = NewDouble(NumCols, "xis");
   scale = NewDouble(NumCols, "scale");
   
   Current = NewIterate(NumRows, NumCols, NumBounds);
   Current->BoundIndex = iupbound;
   Copy    = NewIterate(NumRows, NumCols, NumBounds);
   Predictor = NewIterate(NumRows, NumCols, NumBounds);
   Predictor->BoundIndex = iupbound;
   Corrector = NewIterate(NumRows, NumCols, NumBounds);
   Corrector->BoundIndex = iupbound;
   
   min_phi = NewDouble(Inputs->IterationLimit+1, "min_phi");
   
   /* Compute norms for convergence tests */
   
   rhs2norm  = sqrt(TwoNorm2(b, &NumRows) + TwoNorm2(upbound, &NumBounds));
   cost2norm = sqrt(TwoNorm2(LP->c, &NumCols));
   
   /* flags that control the main iterative loop */
   algorithm_loop_flag = ON;
   Iteration = 0;
   
  /*******************************************************************/
  /* DETERMINE STRUCTURE OF THE L FACTOR of the ADA^T system that    */
  /* is formulated and solved at every iteration, and find the       */
  /* optimal pivot ordering                                          */
  /*******************************************************************/

   maskDense = NewInt(A->NumCols, "maskDense");
   Ndense = LookforDenseColumns(A, &nonzDense, maskDense);
   if ( Ndense > 0 )
      {
	 Asparse = NewMMTtype(A->NumRows, A->NumCols - Ndense, 
			      A->Nonzeros - nonzDense);
	 StripA(A, Asparse, maskDense, 0);
	 scaleSparse = NewDouble(NumCols-Ndense, "scaleSparse");
	 Adense = NewMMTtype(A->NumRows, Ndense, nonzDense);
	 StripA(A, Adense, maskDense, 1);
	 Factor = NewFactorType(Asparse, Ndense, A->NumCols);
	 /* Also compute A^T, as A is not an argument of NewFactorType */
	 A->ValueT = NewDouble(A->Nonzeros, "A->ValueT");
	 status = TransposeSparseRealMatrix
	    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
	     A->ValueT, A->pBeginRowT, A->pEndRowT, A->RowT,
	     &(A->NumRows), &(A->NumCols));
	 if (status) 
	    printf("Error: TransposeSparseRealMatrix = %d\n", status);
      }
   else
      {
	 Asparse = Adense = 0;
	 scaleSparse = 0;
	 Factor = NewFactorType(A, Ndense, A->NumCols);
      }
   Factor->Ndense = Ndense;
   Free((char*) Factor->maskDense);
   Factor->maskDense = maskDense;
      
   status = Order(Factor, Inputs->OrderAlg);
   if (status != 0) 
      {
	 printf("Error return from Order() routine.\n");
	 printf("Aborting.\n");
	 return status;
      }
   
   strcpy(Solution->FactorizationCode, Factor->FactorizationCode);
   
   Solution->FactorizationHistory->Nonzeros = Factor->NonzerosL;
   Solution->FactorizationHistory->Density =
      (2.0 * Factor->NonzerosL - NumRows) / (NumRows * NumRows);
   Solution->FactorizationHistory->NumDenseCols = Ndense;
   
   status=InitialPoint(A, Asparse, Adense, Factor, b, LP->c, upbound, scale,
		       Current, Inputs, &FactorTime, &SolveTime);
   if (status) 
      {
	 printf(" Error: Returned from InitialPoint() with status %d\n", 
		status);
	 return status;
      }
   
   if (Inputs->HOCorrections) 
      {
	 if(Inputs->MaxCorrections == 0) 
	    {
	       MaxCorrections = MaxGondzioCorrections(Factor, 
						      FactorTime, SolveTime);
	       Inputs->MaxCorrections = MaxCorrections;
	    }
	 else 
	    MaxCorrections = Inputs->MaxCorrections;
		if (Inputs->Diagnostics > 0)
	 printf("\nMaximum Gondzio corrections = %d\n", MaxCorrections);  
      }
   

   if (Inputs->ReportingLevel > 1)
     if(Inputs->HOCorrections && MaxCorrections > 0) 
       {
	 printf("\nIter    Primal       Dual      (PriInf  DualInf)  ");
	 printf("log(mu) dgts corr  Merit\n");
       } 
     else 
       {
	 printf("\nIter    Primal       Dual      (PriInf  DualInf)  ");
	 printf("log(mu) dgts   Merit\n");
       }
   
   GetTime(&UserTime, &SysTime);

   InitTime = UserTime - OldUserTime + SysTime - OldSysTime;

   OldUserTime = UserTime;
   OldSysTime = SysTime;

   /*******************************************************************/
   /* START OF MAIN LOOP                                              */
   /*******************************************************************/
   
   do 
      {

	 /* Shift the positive and negative parts of the split free 
	  * variables, if necessary, to keep the two of them from blowing up */
	 
	 ShiftSplitVariables(LP, Current);
	 ComputeObjInf(A, b, LP->c, upbound, Current,
		       xibound, xib, xis,
		       &PriInf, &DualInf,
		       &mu, &primal_objective, &dual_objective, 
		       &phi, &rmu_ratio);
	 
    /*******************************************************************/
    /* Print a progress report - one line of output per iteration      */
    /*******************************************************************/
	 
	 if(Iteration==0) 
	    digits=0;
	 else 
	    {
	       digits = floor(-log10(fabs(primal_objective - dual_objective) /
				     (1.0 + fabs(primal_objective))));
	       if (digits < 0)
		  digits = 0;
	    }
	 
	 /* Keep track of smallest phi encountered here or earlier */
	 
	 if (Iteration == 0) 
	    {
	       min_phi[Iteration] = phi;
	       rmu_ratio0 = rmu_ratio;
	    } 
	 else 
	    {
	       if (phi < min_phi[Iteration-1]) 
		  {
		     min_phi[Iteration] = phi;
		     /* This is the best point so far.  Keep a copy in Copy */
		     BestIteration = Iteration;
		     for (i = 0; i < NumCols; i++) 
			{
			   Copy->x[i] = Current->x[i]; 
			   Copy->s[i] = Current->s[i];
			}
		     for (i = 0; i < NumRows; i++)
			Copy->pi[i] = Current->pi[i];
		     for (i = 0; i < NumBounds; i++) 
			{
			   Copy->w[i] = Current->w[i]; 
			   Copy->r[i] = Current->r[i];
			}
		  } 
	       else 
		  min_phi[Iteration] = min_phi[Iteration-1];
	    }
	 
	 if (Inputs->ReportingLevel > 1)
	   if(Inputs->HOCorrections && MaxCorrections>0)
	     printf("%3d  %11.4e  %11.4e  (%7.1e %7.1e)  %6.2f   %2d   %2d  %7.1e\n",
		    Iteration, primal_objective, dual_objective,
		    PriInf, DualInf, log10(mu), digits, NumCorrections, phi);
	   else
	     printf("%3d  %11.4e  %11.4e  (%7.1e %7.1e)  %6.2f   %2d   %7.1e\n",
		    Iteration, primal_objective, dual_objective,
		    PriInf, DualInf, log10(mu), digits, phi);
	 
	 fflush(stdout);
	 
      /* Save info about the current iterate, for possible later review */
	 
	 StoreHistory(Solution, Iteration, primal_objective, dual_objective,
		      PriInf, DualInf, mu, NumCorrections, phi);
	 
    /*******************************************************************/
    /* TERMINATION TESTS                                               */
    /*******************************************************************/

    /* first, check for OPTIMAL termination */

	 if (PriInf < PriFeasTol && DualInf < DualFeasTol &&
	     /*        (fabs(primal_objective - dual_objective) / 
		      (1.0 + (fabs(primal_objective))) < OptTol)) { */
	     mu / (1.0 + (fabs(primal_objective))) < OptTol) 
	    {
	       Solution->Status = OPTIMAL_SOL;
	       if (Inputs->ReportingLevel > 0)
		 printf("\n--termination with OPTIMAL status\n");
	       
	       algorithm_loop_flag = OFF;
	       break;
	    }
	 
    /* If we are much bigger than best phi so far, declare INFEASIBLE */

	 if (phi > MAX(1.e5 * min_phi[Iteration], 1.0e-8)) 
	    {
	
	      if (Inputs->ReportingLevel > 0)
		Solution->Status = INFEASIBLE_SOL;
	       printf("\n--termination with INFEASIBLE status\n");
	       algorithm_loop_flag = OFF;
	       break;
	    }    
	 
    /* If the ratio of infeasibility to mu is getting out of hand,
       declare UNKNOWN */

	 if((PriInf >= PriFeasTol || DualInf >= DualFeasTol) &&
	    (rmu_ratio / rmu_ratio0) >= 1.0e6) 
	    {
	       Solution->Status = UNKNOWN_SOL_INF_MU;
	       
	       if (Inputs->ReportingLevel > 0)
		 {
		   printf("\n--termination with UNKNOWN status");
		   printf(" (due to large residual/mu ratio)\n");
		   /* printf("(current rmu) / (initial rmu) = %e\n", 
		      rmu_ratio / rmu_ratio0); */
		 }

	       algorithm_loop_flag = OFF;
	       break;
	    }
	 
    /*  If progress in reducing phi is slow, declare UNKNOWN */

	 if((Iteration >= 30) && 
	    (min_phi[Iteration] >= 0.5*min_phi[Iteration-30])) 
	    {
	       Solution->Status = UNKNOWN_SOL_PHI_SLOW;
	       if (Inputs->ReportingLevel > 0)
		 {
		   printf("\n--termination with UNKNOWN status ");
		   printf(" (due to slow merit improvement)\n");
		 }
	       algorithm_loop_flag = OFF;
	       break;
	    }

    /* If we have maxed out on iterations, declare SUBOPTIMAL */
	
	 if (Iteration == Inputs->IterationLimit) 
	    {
	       Solution->Status = SUBOPTIMAL_SOL;
	       
	       if (Inputs->ReportingLevel > 0)
		 printf("Performed Maximum number of iterations.\n");
	       
	       algorithm_loop_flag = OFF;
	    }

    /*******************************************************************/
    /* FORM ADA^T AND COMPUTE ITS CHOLESKY FACTORIZATION               */
    /*******************************************************************/

	 GetTime(&OldUserTime2, &OldSysTime2);

	 /* compute scaling matrix Theta^-1 = (X^{-1}S + W^{-1}R)^{-1} */

	 for (i = 0; i < NumCols; i++)
	    scale[i] = Current->s[i] / Current->x[i];
	 
	 for (i = 0; i < NumBounds; i++) 
	    {
	       k = iupbound[i] - 1;
	       scale[k] += Current->r[i] / Current->w[i];
	    }
	 
	 for (i = 0; i < NumCols; i++)
	    scale[i] = 1.0 / scale[i];

	 /* compute coefficient matrix */
	 
	 if ( Ndense > 0 )
	    {
	       StripScale(scale, scaleSparse, Factor->maskDense, 
			  NumCols, 0);
	       ComputeADAT(Asparse, scaleSparse, Factor->AAT);
	    }
	 else   
	    ComputeADAT(A, scale, Factor->AAT);
	 
	 GetTime(&UserTime2, &SysTime2);

	 FormADATtime += UserTime2 + SysTime2 - OldUserTime2 - OldSysTime2;


	 /* perform the numerical factorization */

	 GetTime(&OldUserTime2, &OldSysTime2);

	 status = Factorize(Factor, Inputs);

	 GetTime(&UserTime2, &SysTime2);

	 FactorizationTime += UserTime2 + SysTime2 - OldUserTime2 - 
	                                                    OldSysTime2;
	 
	 if ( status ) 
	    {
	       printf("Error: returned from Factorize() with status %d\n", 
		      status);
	       return FACTORIZE_ERROR;
	    }
	 if ( Factor->Ndense > 0 ) 
	    {
	       status = ComputeWandLdense(Adense, Factor, scale, NumCols);
	       if (status) 
		  {
		     printf("Error: returned from ComputeWandLdense() with");
		     printf(" status %d\n", status);
		     return FACTORIZE_ERROR;
		  }
	    }

    /*******************************************************************/
    /* COMPUTE PREDICTOR STEP                                          */
    /*******************************************************************/

	 GetTime(&OldUserTime2, &OldSysTime2);

	 ComputePredictor(A, Factor, scale, Current,xis, xibound, xib,
			  upbound, PriFeasTol, Predictor, Inputs, &time2);
	 	
	 GetTime(&UserTime2, &SysTime2);

	 PredictorTime += UserTime2 + SysTime2 - OldUserTime2 - OldSysTime2;

	 SolveADATTime += time2;
	
    /*******************************************************************/
    /* COMPUTE CENTERING PARAMETER                                     */
    /*******************************************************************/

	 sigma = ComputeCentering(Current, Predictor, mu, &alpha_P, &alpha_D, 
				  Inputs);

    /*******************************************************************/
    /* COMPUTE CORRECTOR STEP                                          */
    /*******************************************************************/

	 GetTime(&OldUserTime2, &OldSysTime2);

	 ComputeCorrector(A, Factor, scale, Current, xis, upbound, 
			  sigma, mu, PriFeasTol, Predictor, Corrector, 
			  Inputs, &time2);
	 

	 SolveADATTime += time2;

    /*******************************************************************/
    /* ADD PREDICTOR AND CORRECTOR (result in Predictor)               */
    /*******************************************************************/
	 
	 for (i = 0; i < NumCols; i++) 
	    {
	       Predictor->x[i] += Corrector->x[i];
	       Predictor->s[i] += Corrector->s[i];
	    }
	 for (i = 0; i < NumBounds; i++) 
	    {
	       Predictor->w[i] += Corrector->w[i];
	       Predictor->r[i] += Corrector->r[i];
	    }
	 for (i = 0; i < NumRows; i++)
	    Predictor->pi[i] += Corrector->pi[i];
	 
    /*******************************************************************/
    /* COMPUTE GONDZIO CORRECTIONS                                     */
    /*******************************************************************/
	 
	 if (Inputs->HOCorrections && MaxCorrections>0)
	    NumCorrections = ComputeGondzioCorrections(A, Factor, 
						       scale, Current, 
						       sigma, mu, 
						       Predictor, Corrector, 
						       Inputs, MaxCorrections);
	 
	 GetTime(&UserTime2, &SysTime2);

	 CorrectorTime += UserTime2 + SysTime2 - OldUserTime2 - OldSysTime2;

    /*******************************************************************/
    /* COMPUTE STEP LENGTH  - Use Mehrotra's adaptive heuristic        */
    /*******************************************************************/

	 ComputeStepFactor(Current, Predictor, &alpha_P, &alpha_D, Inputs);
	 
    /*******************************************************************/
    /* TAKE THE STEP -- UPDATE THE CURRENT SOLUTION VECTOR             */
    /*******************************************************************/

	 for (i = 0; i < NumCols; i++) 
	    {
	       Current->x[i] -= alpha_P * Predictor->x[i];
	       Current->s[i] -= alpha_D * Predictor->s[i];
	    }
	 
	 for (i = 0; i < NumRows; i++)
	    Current->pi[i] -= alpha_D * Predictor->pi[i];
	 
	 for (i = 0; i < NumBounds; i++) 
	    {
	       Current->w[i] -= alpha_P * Predictor->w[i];
	       Current->r[i] -= alpha_D * Predictor->r[i];
	    }
	 
	 Iteration++;
	 Solution->Iterations = Iteration;
	 
    /*******************************************************************/
    /* END OF MAIN LOOP                                                */
    /*******************************************************************/

      } while (algorithm_loop_flag);

   GetTime(&UserTime, &SysTime);

   LoopTime = UserTime - OldUserTime + SysTime - OldSysTime;

  /*******************************************************************/
  /* Copy best iteration to Current                                  */
  /*******************************************************************/
  
   if (Solution->Status != OPTIMAL_SOL) 
      {
	/* move solution in Copy to Current */

	if (Inputs->ReportingLevel > 0)
	  printf("Restoring best stored solution (iteration %d).\n", 
		 BestIteration);
	 Solution->RestoredIteration = BestIteration;
	 for (i = 0; i < NumCols; i++) 
	    {
	       Current->x[i] = Copy->x[i]; 
	       Current->s[i] = Copy->s[i];
	    }
	 for (i = 0; i < NumRows; i++)
	    Current->pi[i] = Copy->pi[i];
	 for (i = 0; i < NumBounds; i++) 
	    {
	       Current->w[i] = Copy->w[i]; 
	       Current->r[i] = Copy->r[i];
	    }
      }
   Free((char*) min_phi);

  /*******************************************************************/
  /* TRANSFER RESULTS TO "Solution"                                  */
  /*******************************************************************/
   Solution->Iterations = Iteration;
   
  /* Compute objective and infeasibilities at this final point */

   ComputeObjInf(A, b, LP->c, upbound, Current, xibound, xib, xis,
		 &PriInf, &DualInf, &mu, &primal_objective, &dual_objective, 
		 &phi, &rmu_ratio);
   
   Solution->PrimalObjective = primal_objective;
   Solution->DualObjective = dual_objective;

   
   for (i = 0; i < NumBounds; i++) 
      {
	 k = LP->BoundIndex[i];
	 Solution->DualUpper[k] = Current->r[i];
      }
   for (i = 0; i < NumCols; i++) 
      {
	 Solution->x[i] = Current->x[i];
	 Solution->DualLower[i] = Current->s[i];
      }
   for (i = 0; i < NumRows; i++)
      Solution->pi[i] = Current->pi[i];
   
   /* Recompute DualUpper and DualLower, based on Lagrange multiplier pi */


   RecomputeDualVariables(LP, Solution);

  /*******************************************************************/
  /* recycle!                                                        */
  /*******************************************************************/

   FreePCxMemory(xib, xibound, xis, scale, scaleSparse, Current, 
		 Predictor, Corrector, Copy, iupbound, upbound, 
		 b, A, Adense, Asparse, Factor, maskDense);

   /* Add factorization and iteration times */
   Solution->FactorizationTime = FactorizationTime;
   Solution->SolveADATTime     = SolveADATTime;
   Solution->LoopTime          = LoopTime;
   Solution->PredictorTime     = PredictorTime;
   Solution->CorrectorTime     = CorrectorTime;
   Solution->FormADATtime      = FormADATtime;
   Solution->InitTime          = InitTime;

   return 0;			/* normal return */
}

/*****************************************************************/
/* Mehrotra's heuristic to compute the starting point            */
/*****************************************************************/

int             
InitialPoint(A, Asparse, Adense, Factor, b, c, 
	     upbound, scale, Current, Inputs,
	     FactorTime, SolveTime)
     MMTtype        *A, *Asparse, *Adense;
     FactorType     *Factor;
     double         *b, *c, *upbound, *scale;
     Iterate        *Current;
     Parameters     *Inputs;
     double         *FactorTime, *SolveTime;
{
   int             i, irow, status, Factorize(), EnhancedSolve(), 
                   SparseSaxpyM(), SparseSaxpyTM();
   void            StripScale();
   int             NumRows, NumCols, NumBounds, *iupbound;
   double          one = 1.0, *tempR1, *tempR2, *tempC, *scaleSparse, 
                   delta_primal, delta_dual,
                   xTs, wTr, sr_sum, xw_sum, *pi, *s, *r, *x, *w;

   double          BeginUserTime, EndUserTime, BeginSysTime, EndSysTime;
   double          time1;


   NumRows = Current->NumRows;
   NumCols = Current->NumCols;
   NumBounds = Current->NumBounds;
   
   pi = Current->pi;
   s = Current->s;
   r = Current->r;
   x = Current->x;
   w = Current->w;
   iupbound = Current->BoundIndex;
   
   tempR1 = NewDouble(NumRows, "tempR1");
   tempR2 = NewDouble(NumRows, "tempR2");
   tempC = NewDouble(NumCols, "tempC");
   
/******************************************************************************
 * pi = (AA^T)^-1 Ac;   s = 0.5 (c - A^T pi);  r = -s;                        *
 * x = u/2 - A^T (AA^T)^-1 (A d/2 - b); w = u - x                             *
 *                                                                            *
 * delta_primal = max(-1.5 min{x_i}, -1.5 min{w_i}, 0.0);                     *
 * delta_dual = max(-1.5 min{s_i}, -1.5 min{r_i}, 0.0);                       *
 *                                                                            *
 * xTs = (x + delta_primal e)^T (s + delta_primal e)                          *
 * wTr = (w + delta_primal e)^T (r + delta_dual e)                            *
 *                                                                            *
 *                                           xTs + wTr                        *
 * Delta_Primal = delta_primal + .5 ------------------------------------      *
 *                                  sum (s_i + delta_dual + r_i + delta_dual) *
 *                                                                            *
 *                                        xTs + wTr                           *
 * Delta_Dual = delta_dual + .5 -----------------------------------------     *
 *                              sum (x_i + delta_primal + w_i + delta_primal) *
 *                                                                            *
 * pi = pi; x = x + Delta_Primal e; w = w + Delta_Primal e;                   *
 *          s = s + Delta_Dual e; r = r + Delta_Dual e;                       *
******************************************************************************/

   for (i = 0; i < NumCols; i++)
      scale[i] = 1.0;
   
   GetTime(&BeginUserTime, &BeginSysTime);
   
   if ( Factor->Ndense > 0 )
      {
	 scaleSparse = NewDouble(NumCols-Factor->Ndense, 
  				 "scaleSparse in InitialPoint");
	 StripScale(scale, scaleSparse, Factor->maskDense, NumCols, 0);
	 ComputeADAT(Asparse, scaleSparse, Factor->AAT);
	 Free((char *) scaleSparse);
      }
   else  
      ComputeADAT(A, scale, Factor->AAT);
   
   status = Factorize(Factor, Inputs);
   if ( status ) 
      {
	 printf("Error: Returned from Factorize() with status %d\n",status);
	 return FACTORIZE_ERROR;
      }
   if (Factor->Ndense > 0) 
      {
	 status = ComputeWandLdense(Adense, Factor, scale, NumCols);
	 if (status) 
	    {
	       printf("Error: returned from ComputeWandLdense() ");
	       printf("with status %d\n", status);
	       return FACTORIZE_ERROR;
	    };
      };
   
   GetTime(&EndUserTime, &EndSysTime);
   *FactorTime = EndUserTime - BeginUserTime;
   
   /* compute pi */
   
   for (i = 0; i < NumRows; i++)
      tempR1[i] = 0.0;
   if (SparseSaxpyM(A, c, tempR1)) 
      {
	 printf("Error: Returned from SparseSaxpyM() with error condition\n");
	 return FACTORIZE_ERROR;
      }
   
   GetTime(&BeginUserTime, &BeginSysTime);
   
   EnhancedSolve(Factor, tempR1, pi);
   
   /* compute s and r */
   
   for (i = 0; i < NumCols; i++)
      tempC[i] = 0.0;
   if (SparseSaxpyTM(A, pi, tempC)) 
      {
	 printf("Error: Returned from SparseSaxpyTM() with error condition\n");
	 return FACTORIZE_ERROR;
      }
   for (i = 0; i < NumCols; i++)
      s[i] = 0.5 * (c[i] - tempC[i]);
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = iupbound[i] - 1;
	 r[i] = -s[irow];
      }
   
   /* compute x */
   
   for (i = 0; i < NumRows; i++)
      tempR2[i] = -b[i];
   for (i = 0; i < NumCols; i++)
      tempC[i] = 0.0;
   
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = iupbound[i] - 1;
	 tempC[irow] = upbound[i] / 2.0;
      }
   
   GetTime(&EndUserTime, &EndSysTime);
   *SolveTime = EndUserTime - BeginUserTime;
   
   if (SparseSaxpyM(A, tempC, tempR2)) 
      {
	 printf("Error: Returned from SparseSaxpyM() with error condition\n");
	 return FACTORIZE_ERROR;
      }
   EnhancedSolve(Factor, tempR2, tempR1);
   
   for (i = 0; i < NumCols; i++)
      x[i] = 0.0;
   
   if (SparseSaxpyTM(A, tempR1, x)) 
      {
	 printf("Error: Returned from SparseSaxpyTM() with error condition\n");
	 return FACTORIZE_ERROR;
      }
   
   for (i = 0; i < NumCols; i++)
      x[i] = -x[i];
   
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = iupbound[i] - 1;
	 x[irow] += upbound[i] / 2.0;
     }
   
   /* compute w */
   
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = iupbound[i] - 1;
	 w[i] = upbound[i] - x[irow];
      }

   
   /* compute delta_primal and delta_dual */
   
   delta_primal = delta_dual = 1.0e30;	/* initialize with large number */
   for (i = 0; i < NumCols; i++) 
      {
	 delta_primal = MIN(delta_primal, x[i]);
	 delta_dual = MIN(delta_dual, s[i]);
      }
   
   for (i = 0; i < NumBounds; i++) 
      {
	 delta_primal = MIN(delta_primal, w[i]);
	 delta_dual   = MIN(delta_dual, r[i]);
      }
   
   delta_primal = MAX(-1.5 * delta_primal, 0.01);
   delta_dual   = MAX(-1.5 * delta_dual, 0.01);
   
   /* compute Delta_primal and Delta_dual */
   
   xTs = wTr = sr_sum = xw_sum = 0.0;
   for (i = 0; i < NumCols; i++) 
      {
	 xTs += (x[i] + delta_primal) * (s[i] + delta_dual);
	 sr_sum += s[i];
	 xw_sum += x[i];
      }
   
   for (i = 0; i < NumBounds; i++)
      {
	 wTr += (w[i] + delta_primal) * (r[i] + delta_dual);
	 sr_sum += r[i];
	 xw_sum += w[i];
      }
   
   sr_sum += (NumCols + NumBounds) * delta_dual;
   xw_sum += (NumCols + NumBounds) * delta_primal;
   
   if (sr_sum == 0.0)
      sr_sum = 1.0;
   
   if (xw_sum == 0.0)
      xw_sum = 1.0;
   
   delta_dual += 0.5 * (xTs + wTr) / xw_sum;
   delta_primal += 0.5 * (xTs + wTr) / sr_sum;
   
   /* update x and s */
   
   for (i = 0; i < NumCols; i++) 
      {
	 x[i] += delta_primal;
	 s[i] += delta_dual;
      }
   
   for (i = 0; i < NumBounds; i++) 
      {
	 r[i] += delta_dual;
	 w[i] += delta_primal;
      }
   

   Free((char *) tempR1);
   Free((char *) tempR2);
   Free((char *) tempC);
   return 0;
}
