/* solve for the predictor and corrector steps 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <math.h>
#include <stdio.h>
#include "main.h"
#include "memory.h"

ComputePredictor(A, Factor, scale, Current,
		 xis, xibound, xib,
		 upbound, PriFeasTol, Predictor, Inputs, time)
     MMTtype        *A;
     
     FactorType     *Factor;
     
     double         *scale, *upbound;
     double         *xis, *xibound, *xib, PriFeasTol;
     Iterate        *Current, *Predictor;
     Parameters     *Inputs;
     double         *time;
{
   int             SolveAugmented();
   int             NumRows, NumCols, NumBounds, irow, i, row, col, entry;
   double         *rhs_col, *rhs_row, *x, *s, *pi, *r, *w, 
                  *dx, *ds, *dpi, *dr, *dw;
   
  /*******************************************************************/
  /* Initialize and transfer to local pointers                       */
  /*******************************************************************/

   x = Current->x;     dx = Predictor->x;
   s = Current->s;     ds = Predictor->s;
   w = Current->w;     dw = Predictor->w;
   r = Current->r;     dr = Predictor->r;
   pi = Current->pi;  dpi = Predictor->pi;
   
   NumRows = A->NumRows;
   NumCols = A->NumCols;
   NumBounds = Current->NumBounds;
   
   rhs_col = NewDouble(NumCols, "rhs_column");
   rhs_row = NewDouble(NumRows, "rhs_row");
   
  /*******************************************************************/
  /* rhs = [xis - s + r - W^-1R xiw    xib]                          */
  /*******************************************************************/
   
   for (i = 0; i < NumCols; i++)
      rhs_col[i] = (xis[i] - s[i]);
   
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = Current->BoundIndex[i] - 1;
	 rhs_col[irow] += (r[i] - xibound[i] * r[i] / w[i]);
      }
   
   for (i = 0; i < NumRows; i++)
      rhs_row[i] = xib[i];
   
  /*******************************************************************/
  /* Solve for predictor (affine-scaling) step                       */
  /*******************************************************************/

   /* find dpi and dx by solving the augmented system */
   
   SolveAugmented(A, Factor, rhs_col, rhs_row, dx, dpi, scale,
		  PriFeasTol, NumRows, NumCols, Inputs, time);
   
   /* recover dw */
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = Current->BoundIndex[i] - 1;
	 dw[i] = xibound[i] - dx[irow];
      }
   
   /* recover ds */
   for (i = 0; i < NumCols; i++)
      ds[i] = s[i] * (1.0 - (dx[i] / x[i]));
   
   /* recover dr */
   for (i = 0; i < NumBounds; i++)
      dr[i] = r[i] * (1.0 - (dw[i] / w[i]));
   
   Free((char *) rhs_col);
   Free((char *) rhs_row);
   return 0;
}

/******************************************************************/

ComputeCorrector(A, Factor, scale, Current, xis, upbound, sigma, mu, 
		 PriFeasTol, Predictor, Corrector, Inputs, time)
     MMTtype        *A;
     
     FactorType     *Factor;
     
     double         *scale;
     double         *xis, *upbound;
     double          sigma, mu, PriFeasTol;
     Iterate        *Current, *Predictor, *Corrector;
     Parameters     *Inputs;
     double         *time;
{
   int             SolveAugmented();
   int             NumRows, NumCols, NumBounds, i, irow;
   double         *rhs_col, *rhs_row;
   
   double         *x, *s, *pi, *r, *w;
   double         *dx, *ds, *dpi, *dr, *dw;
   double         *dx2, *ds2, *dpi2, *dr2, *dw2;
   
  /*******************************************************************/
  /* Initialize and transfer to local pointers                       */
  /*******************************************************************/

   x = Current->x;      dx = Predictor->x;        dx2 = Corrector->x;
   s = Current->s;      ds = Predictor->s;        ds2 = Corrector->s;
   w = Current->w;      dw = Predictor->w;        dw2 = Corrector->w;
   r = Current->r;      dr = Predictor->r;        dr2 = Corrector->r;
   pi = Current->pi;   dpi = Predictor->pi;      dpi2 = Corrector->pi;
   
   NumRows = A->NumRows;
   NumCols = A->NumCols;
   NumBounds = Current->NumBounds;
   
   rhs_col = NewDouble(NumCols, "rhs_col");
   rhs_row = NewDouble(NumRows, "rhs_row");

  /*******************************************************************/
  /* compute the right-hand side for the centering direction, which  */
  /* combines a centering component with a second-order correction   */
  /* component.                                                      */
  /*******************************************************************/

   for (i = 0; i < NumCols; i++)
      rhs_col[i] = (sigma * mu - dx[i] * ds[i]) / x[i];
   
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = Current->BoundIndex[i] - 1;
	 rhs_col[irow] += -(sigma * mu - dw[i] * dr[i]) / w[i];
      }
   
   for (i = 0; i < NumRows; i++)
      rhs_row[i] = 0.0;
   
  /*******************************************************************/
  /* Solve for corrector step                                        */
  /*******************************************************************/

  /* find dpi2 and dx2 by solving the augmented system */

   SolveAugmented(A, Factor, rhs_col, rhs_row, dx2, dpi2, scale,
		  PriFeasTol, NumRows, NumCols, Inputs, time);

  /* recover dw2 */
   for (i = 0; i < NumBounds; i++) 
      {
	 irow = Current->BoundIndex[i] - 1;
	 dw2[i] = -dx2[irow];
      }
   
  /* recover ds2 */
   for (i = 0; i < NumCols; i++)
      ds2[i] = s[i] * (-(dx2[i] / x[i])) - (sigma * mu - dx[i] * ds[i]) / x[i];

  /* recover dr2 */
   for (i = 0; i < NumBounds; i++)
      dr2[i] = r[i] * (-(dw2[i] / w[i])) - (sigma * mu - dw[i] * dr[i]) / w[i];

   Free((char *) rhs_col);
   Free((char *) rhs_row);
   return 0;
}

/******************************************************************/

/* The Gondzio correction code is based on a paper by
 * Jacek Gondzio: "Multiple Centrality Corrections in a Primal-Dual
 * Method for Linear Programming," Computational Optimization and 
 * Applications 6 (1996) 137-156, and also by the heuristics in 
 * Gondzio's HOPDM code, version 2.13.
 */

ComputeGondzioCorrections (A, Factor, scale, Current, sigma, mu, 
			   Predictor, Corrector, Inputs, MaxCorrections)

     MMTtype        *A;
     FactorType     *Factor;
     double         *scale;
     double          sigma, mu;
     Iterate        *Current, *Predictor, *Corrector;
     Parameters     *Inputs;
     int             MaxCorrections;
{
   int             SolveAugmented();
   int             NumRows, NumCols, NumBounds, i, irow;
   int             CorrectionNumber;
   int             ContinueCorrecting;
   double         *rhs_col, *rhs_row;
   
   double         *x, *s, *pi, *r, *w;
   double         *dx, *ds, *dpi, *dr, *dw;
   double         *dx2, *ds2, *dpi2, *dr2, *dw2;
   double         *comp_xs, *comp_rw, sigmamu;
   double          munew,signew;
   
   double          Step0_p, Step0_d, Step1_p, Step1_d;
   double          mu_min, mu_max;   /* define the Gondzio centering box */
   
   double          beta_min = 0.1, beta_max = 10.0;  /* Gondzio algorithmic 
							parameters */
   double          StepFactor0 = 0.08, StepFactor1 = 1.08;

   double         *time; /* this is ignored here */

  /*******************************************************************/
  /* Initialize and transfer to local pointers                       */
  /*******************************************************************/
  
   x = Current->x;    dx = Predictor->x;    dx2 = Corrector->x;
   s = Current->s;    ds = Predictor->s;    ds2 = Corrector->s;
   w = Current->w;    dw = Predictor->w;    dw2 = Corrector->w;
   r = Current->r;    dr = Predictor->r;    dr2 = Corrector->r;
   pi = Current->pi; dpi = Predictor->pi;  dpi2 = Corrector->pi;
   
   NumRows = A->NumRows;
   NumCols = A->NumCols;
   NumBounds = Current->NumBounds;
   
   rhs_col = NewDouble(NumCols, "rhs_col");
   rhs_row = NewDouble(NumRows, "rhs_row");
   
   comp_xs = NewDouble(NumCols, "comp_xs");
   comp_rw = NewDouble(NumBounds, "comp_rw");
   
   CorrectionNumber   = 0;
   ContinueCorrecting = 1;
   
  /*******************************************************************/
  /* Compute step length for new point                               */
  /*  Given Mehrotra point, compute step length                      */
  /*******************************************************************/

   StepToBoundary(Current, Predictor, NULL, &Step0_p, &Step0_d);

   while (ContinueCorrecting && (CorrectionNumber < MaxCorrections)) 
      {
	 
	/* use these max steplengths to compute the target munew */
	 munew=0.0;
	 for (i = 0; i < NumCols; i++) 
	    munew += (x[i] - Step0_p * dx[i]) * (s[i] - Step0_d * ds[i]);
	 
	 for (i = 0; i < NumBounds; i++) 
	    munew += (w[i] - Step0_p * dw[i]) * (r[i] - Step0_d * dr[i]);
	 
	 munew /= (NumCols + NumBounds);
	 signew = pow(munew / mu, Inputs->CenterExponent); 
	 
	 /* recompute sigmamu, mu_min, mu_max */
	 sigmamu = signew * munew;
	 mu_min  = sigmamu * beta_min;
	 mu_max  = sigmamu * beta_max;
	 
	 /* set target step lengths */
	 Step1_p = MIN(StepFactor1 * Step0_p + StepFactor0, 1.0);
	 Step1_d = MIN(StepFactor1 * Step0_d + StepFactor0, 1.0);
	
    /*******************************************************************/
    /* Compute Xs and Wr                                               */
    /* Limit elements of Xs and Wr to bounds                           */
    /*******************************************************************/

#if 0
	 /* the following implements the heuristics that appear in
            Gondzio's code, which differ significantly from those in
            his paper. */

	 for (i = 0; i < NumCols; i++) 
	    {
	       comp_xs[i] = (x[i] - Step1_p * dx[i]) * 
		  (s[i] - Step1_d * ds[i]);
	       
	       if (comp_xs[i] < mu_min)
		  comp_xs[i] = sigmamu - comp_xs[i]; 
	       else if (comp_xs[i] > mu_max)
		  comp_xs[i] = -5*sigmamu;
	       else
		  comp_xs[i] = 0.0;
	    }
	 
	 for (i = 0; i < NumBounds; i++) 
	    {
	       comp_rw[i] = (w[i] - Step1_p * dw[i]) * 
		  (r[i] - Step1_d * dr[i]);
	       
	       if (comp_rw[i] < mu_min)
		  comp_rw[i] = sigmamu - comp_rw[i]; 
	       else if (comp_rw[i] > mu_max)
		  comp_rw[i] = -5*sigmamu;
	       else
		  comp_rw[i] = 0.0;
	    }
	
#endif

	 /* This code implements Gondzio's heuristics as described in
            his paper. They are more satisfying than the techniques
            above, and/but the practical performance is not very
            different. */

 	 for (i = 0; i < NumCols; i++) 
	    {
	       comp_xs[i] = (x[i] - Step1_p * dx[i]) * 
		  (s[i] - Step1_d * ds[i]);
	       
	       if (comp_xs[i] < mu_min)
		  comp_xs[i] = mu_min - comp_xs[i]; 
	       else if (comp_xs[i] > mu_max)
		  comp_xs[i] = mu_max - comp_xs[i];
	       else
		  comp_xs[i] = 0.0;

	       if(comp_xs[i] < -mu_max) 
		 comp_xs[i] = -mu_max;
	    }
	 
	 for (i = 0; i < NumBounds; i++) 
	    {
	       comp_rw[i] = (w[i] - Step1_p * dw[i]) * 
		  (r[i] - Step1_d * dr[i]);
	       
	       if (comp_rw[i] < mu_min)
		  comp_rw[i] = mu_min - comp_rw[i]; 
	       else if (comp_rw[i] > mu_max)
		  comp_rw[i] = mu_max - comp_rw[i];
	       else
		  comp_rw[i] = 0.0;

	       if(comp_rw[i] < -mu_max) 
		 comp_rw[i] = -mu_max;
	    }

    /*******************************************************************/
    /* Compute the right-hand side for the Gondzio corrections         */
    /*    (rhs is negative of what is in paper)                        */
    /*******************************************************************/
    
	 for (i = 0; i < NumCols; i++)
	    rhs_col[i] = (comp_xs[i] / x[i]);
	 
	 for (i = 0; i < NumBounds; i++) 
	    {
	       irow = Current->BoundIndex[i] - 1;
	       rhs_col[irow] -= (comp_rw[i] / w[i]);
	    }
	 
	 for (i = 0; i < NumRows; i++)
	    rhs_row[i] = 0.0;
    
    /*******************************************************************/
    /* Solve for correction step                                       */
    /*******************************************************************/
    
	 /* find dpi2 and dx2 by solving the augmented system */

	 /* no timing here so far */
	 time = (double *) Malloc(sizeof(double), "time in Gondzio_corr");
  
	 SolveAugmented(A, Factor, rhs_col, rhs_row, dx2, dpi2, scale,
			1.0e-8, NumRows, NumCols, Inputs, time);

	 Free((char *) time);
	 
	 /* recover dw2 */
	 for (i = 0; i < NumBounds; i++) 
	    {
	       irow = Current->BoundIndex[i] - 1;
	       dw2[i] = -dx2[irow];
	    }
	 
	 /* recover ds2 */
	 for (i = 0; i < NumCols; i++)
	    ds2[i] = (-comp_xs[i] - s[i] * dx2[i]) / x[i];
	 
	 /* recover dr2 */
	 for (i = 0; i < NumBounds; i++)
	    dr2[i] = (-comp_rw[i] - r[i] * dw2[i]) / w[i];

    /*******************************************************************/
    /* Compute step length for predictor + new corrector step          */
    /*******************************************************************/

	 StepToBoundary(Current, Predictor, Corrector, &Step1_p, &Step1_d);

    /* accept the corrected step even if the improvement is minimal */
	 if((Step1_p >= MIN(1.01 * Step0_p, 1.0)) &&
	    (Step1_d >= MIN(1.01 * Step0_d, 1.0))) 
	    {
	       
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
	       
	       CorrectionNumber++;
	       ContinueCorrecting = 1;
	       
	       Step0_p = Step1_p; 
	       Step0_d = Step1_d; 
	    } 
	 else 
	    {
	       /* stop iterating */
	       ContinueCorrecting = 0;
	    }
	 
      } /* end while */
   
   Free((char *) rhs_col);
   Free((char *) rhs_row);
   
   Free((char *) comp_xs);
   Free((char *) comp_rw);
   
   return CorrectionNumber;
}

/******************************************************************/

int 
StepToBoundary(Current, Predictor, Corrector, alpha_P, alpha_D)
     Iterate        *Current, *Predictor, *Corrector;
     double         *alpha_P, *alpha_D;
{
   int    i;
   double step;
   
  /*******************************************************************/
  /* Find the distance to boundary, i.e. the largest alpha that does */
  /* not violate nonnegativity of the x, s, w, and r components      */
  /*******************************************************************/

  if (Corrector == NULL) 
     {

	/* find step size for predictor correction */
	
	*alpha_P = 1.0;
	for (i = 0; i < Current->NumCols; i++)
	   if (Current->x[i] - Predictor->x[i] < 0.0)
	      if (Current->x[i] - *alpha_P * Predictor->x[i] < 0.0) 
		 *alpha_P = Current->x[i] / Predictor->x[i];
	
	
	for (i = 0; i < Current->NumBounds; i++)
	   if (Current->w[i] - Predictor->w[i] < 0.0)
	      if (Current->w[i] - *alpha_P * Predictor->w[i] < 0.0) 
		 *alpha_P = Current->w[i] / Predictor->w[i];
	
	
	*alpha_D = 1.0;
	for (i = 0; i < Current->NumCols; i++)
	   if (Current->s[i] - Predictor->s[i] < 0.0)
	      if (Current->s[i] - *alpha_D * Predictor->s[i] < 0.0)
		 *alpha_D = Current->s[i] / Predictor->s[i];
	
	
	for (i = 0; i < Current->NumBounds; i++)
	   if (Current->r[i] - Predictor->r[i] < 0.0)
	      if (Current->r[i] - *alpha_D * Predictor->r[i] < 0.0) 
		 *alpha_D = Current->r[i] / Predictor->r[i];
	
     } 
  else 
     { 

    /* Find step size for predictor + corrector direction */
    /* This is for Gondzio corrections                    */

    *alpha_P = 1.0;
    for (i = 0; i < Current->NumCols; i++) 
       {
	  step = Predictor->x[i] + Corrector->x[i];
	  if (Current->x[i] - step < 0.0)
	     if (Current->x[i] - *alpha_P * step < 0.0) 
		*alpha_P = Current->x[i] / step;
	  
       }    
    
    for (i = 0; i < Current->NumBounds; i++) 
       {
	  step = Predictor->w[i] + Corrector->w[i];
	  if (Current->w[i] - step < 0.0)
	     if (Current->w[i] - *alpha_P * step < 0.0) 
		*alpha_P = Current->w[i] / step;
       }

    *alpha_D = 1.0;
    for (i = 0; i < Current->NumCols; i++) 
       {
	  step = Predictor->s[i] + Corrector->s[i];
	  if (Current->s[i] - step < 0.0)
	     if (Current->s[i] - *alpha_D * step < 0.0) 
		*alpha_D = Current->s[i] / step;
	  
       }
    
    for (i = 0; i < Current->NumBounds; i++) 
       {
	  step = Predictor->r[i] + Corrector->r[i];
	  if (Current->r[i] - step < 0.0)
	     if (Current->r[i] - *alpha_D * step < 0.0) 
		*alpha_D = Current->r[i] / step;
       }
    
     } /* end if Corrector == NULL */

  return(0);
}

/******************************************************************/

double 
ComputeCentering(Current, Predictor, mu, alpha_P, alpha_D, Inputs)
     Iterate        *Current, *Predictor;
     double         mu, *alpha_P, *alpha_D;
     Parameters     *Inputs;
{
   int             i, NumBounds, NumCols;
   double          mdg, sigma;
   double         *x, *s, *w, *r, *dx, *ds, *dw, *dr;
   
  /*******************************************************************/
  /* Transfer to local pointers                                      */
  /*******************************************************************/

   NumCols = Current->NumCols;
   NumBounds = Current->NumBounds;
   x = Current->x;     dx = Predictor->x;
   s = Current->s;     ds = Predictor->s;
   w = Current->w;     dw = Predictor->w;
   r = Current->r;     dr = Predictor->r;

  /*******************************************************************/
  /* Find the distance to boundary, i.e. the largest alpha that does */
  /* not violate nonnegativity of the x, s, w, and r components      */
  /*******************************************************************/

   StepToBoundary(Current, Predictor, NULL, alpha_P, alpha_D);

  /*******************************************************************/
  /* Use the largest step in a heuristic to compute the              */
  /* centering parameter mu. (We use Mehrotra's original heuristic)  */
  /*******************************************************************/

   mdg = 0.0;
   for (i = 0; i < NumCols; i++)
      mdg += (x[i] - *alpha_P * dx[i]) * (s[i] - *alpha_D * ds[i]);
   for (i = 0; i < NumBounds; i++)
      mdg += (w[i] - *alpha_P * dw[i]) * (r[i] - *alpha_D * dr[i]);
   mdg /= (NumCols + NumBounds);
   
   sigma = pow(mdg / mu, Inputs->CenterExponent); 
   
   return sigma;
}

  int
SolveADAT(A, Factor, rhs, sol, scale, PriFeasTol, NumRows,
          NumCols, Inputs)
     MMTtype        *A;
     FactorType     *Factor;
     double         *rhs, *sol, *scale, PriFeasTol;
     int            *NumRows, *NumCols;
     Parameters     *Inputs;
{
   int             status, i, EnhancedSolve(), PreConjGrad(), density;
   double         *residual, *correction = NULL, residual2norm,
                   rhs2norm, firstnorm, pcgtol, TwoNorm2();

   residual = NewDouble(*NumRows, "residual");
   correction = NewDouble(*NumRows, "correction");
   EnhancedSolve(Factor, rhs, sol);

  /* check the accuracy of the solutions */
   if (FindResidualSolutionofNormalEquations
       (A->Value, A->Row, A->pBeginRow, A->pEndRow, scale,
        sol, rhs, residual, NumRows, NumCols))
      {
         printf("\nError: PCx: Error detected during \
                  FindResidualSolutionofNormalEquations.\n");
         Free((char *) correction); Free((char *) residual);
         return FACTORIZE_ERROR;
      }
   residual2norm = sqrt(TwoNorm2(residual, NumRows));
   rhs2norm      = sqrt(TwoNorm2(rhs, NumRows));

   if ( residual2norm > rhs2norm*PriFeasTol ) {
      printf("     relative resid = %7.1e;", residual2norm / rhs2norm);

      if ( (Factor->SmallDiagonals == NO) &&
           (Factor->Ndense > 0 ||
            (Factor->Ndense == 0 && Inputs->Refinement)))
        {
            density = (Factor->Ndense > 0)? 10*Factor->Ndense: 10;
            firstnorm = residual2norm;
            for (i = 0; i < *NumRows; i++)
               correction[i] = sol[i];

            pcgtol = rhs2norm*PriFeasTol;
            status = PreConjGrad(A, scale, Factor, residual,
                                 &pcgtol, rhs2norm, density, sol);
            if (status)
               {
                  printf("\nError detected in PreConjGrad routine.\n");
                  Free((char *) correction); Free((char *) residual);
                  return FACTORIZE_ERROR;
               };

            if ( pcgtol > firstnorm )
               {
                  printf("     Restoring original solution since \
                           error increased.\n");
                  for (i = 0; i < *NumRows; i++)
                     sol[i] = correction[i];
               }
         }
      else
         {
            printf("\n");
         }
   }
   Free((char *) correction); Free((char *) residual);
   return 0;
}


/******************************************************************/

int  
SolveAugmented(A, Factor, rhs_col, rhs_row, sol_col, sol_row, scale,
	                       PriFeasTol, NumRows, NumCols, Inputs, time)
     MMTtype        *A;
     FactorType     *Factor;
     double         *rhs_col, *rhs_row, *sol_col, *sol_row;
     double         *scale, PriFeasTol;
     int             NumRows, NumCols;
     Parameters     *Inputs;
     double         *time;
{
   int             i, SparseSaxpyM(), SparseSaxpyTM(), SolveADAT();
   double         *temp_col, *temp_row;
   double          UserTime, SysTime, OldUserTime, OldSysTime;

   temp_col = NewDouble(NumCols, "temp_col");
   temp_row = NewDouble(NumRows, "temp_row");
   
   for (i = 0; i < NumCols; i++)
      temp_col[i] = rhs_col[i] * scale[i];
   
   for (i = 0; i < NumRows; i++)
      temp_row[i] = rhs_row[i];
   
   if (SparseSaxpyM(A, temp_col, temp_row)) 
      {
	 printf("Error: PCx: Error detected during SparseSaxpyM().\n");
	 return FACTORIZE_ERROR;
      }
   /* Perform solve for search direction and check accuracy. */
   
   
   GetTime(&OldUserTime, &OldSysTime);
   SolveADAT(A, Factor, temp_row, sol_row, scale,
	     PriFeasTol, &NumRows, &NumCols, Inputs);
   GetTime(&UserTime, &SysTime);
   
   *time = UserTime - OldUserTime + SysTime - OldSysTime;


   for(i=0; i<NumCols; i++)
      sol_col[i] = -rhs_col[i];
   
   if (SparseSaxpyTM(A, sol_row, sol_col)) 
      {
	 printf("Error: PCx: Error detected during SparseSaxpyT().\n");
	 return FACTORIZE_ERROR;
      }
   for(i=0; i<NumCols; i++)
      sol_col[i] *= scale[i];
   Free((char *) temp_col);
   Free((char *) temp_row);
   return 0;
}

/*****************************************************************/

int             
ComputeStructureAAT(A, AAT)
  MMTtype        *A, *AAT;
{
  int             NumRowsA, NumColsA, NumEntsA;
  int             flag_maximum_nonzeros, max_nonzeros, status;

  NumRowsA = A->NumRows;
  NumColsA = A->NumCols;
  NumEntsA = A->Nonzeros;

  /* Also compute A^T */

  A->ValueT = NewDouble(NumEntsA, "A->ValueT");

  status = TransposeSparseRealMatrix
    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
     A->ValueT, A->pBeginRowT, A->pEndRowT, A->RowT,
     &NumRowsA, &NumColsA);

  if (status)
    printf("Error: TransposeSparseRealMatrix = %d\n", status);

  status = FindNonzeroNormalEquations
    (A->pBeginRow, A->pEndRow, A->Row, &NumRowsA, &NumColsA, &NumEntsA,
     &(AAT->Nonzeros));

  if (status)
    printf("Error: FindNonzeroNormalEquations = %d\n", status);

  /* allocate memory for AAT */

  AAT->NumRows = NumRowsA;
  AAT->NumCols = NumRowsA;

  AAT->pBeginRow = NewInt(NumRowsA + 1, "AAT->pBeginRow");
  AAT->pEndRow   = NewInt(NumRowsA + 1, "AAT->pEndRow");
  AAT->Row       = NewInt(AAT->Nonzeros, "AAT->Row");
  AAT->Value     = NewDouble(AAT->Nonzeros, "AAT->Value");

  AAT->pBeginRowT = NULL;
  AAT->pEndRowT   = NULL;
  AAT->RowT       = NULL;
  AAT->ValueT     = NULL;

  status = FindStructureNormalEquations	/* returns structure of whole AAT */
    (A->pBeginRow, A->pEndRow, A->Row,
     AAT->pBeginRow, AAT->pEndRow, AAT->Row,
     &NumRowsA, &NumColsA, &NumEntsA,
     &(AAT->Nonzeros), &flag_maximum_nonzeros);

  if (flag_maximum_nonzeros != 0) {
    printf("Error in ComputeStructureAAT().\
            Too many nonzeros in Normal Equations.\n");
    return FACTORIZE_ERROR;
  }
  if (status)
    printf("Error: FindStructureNormalEquations = %d\n", status);

  status = SortColumnRealSparseMatrix
    (AAT->Value, AAT->Row, AAT->pBeginRow, AAT->pEndRow,
     &NumRowsA, &NumRowsA);

  if (status)
    printf("Error: SortColumnRealSparseMatrix = %d\n", status);

  /* extend AAT->pBeginRow by 1 */

  AAT->pBeginRow[NumRowsA] = AAT->pEndRow[NumRowsA - 1] + 1;

  return 0;
}

/*****************************************************************/

int             
ComputeADAT(A, scale, AAT)
  MMTtype        *A, *AAT;
  double         *scale;
{
  int             NumRowsA, NumColsA, entryAAT;
  int             ent1, ent2, row1, row2, col1, col2, base, top1, 
                  top2, found_symmetric, temp2, i, j, k;
  double          sum, *Temp, *Temp2;

  NumRowsA = A->NumRows;
  NumColsA = A->NumCols;

  Temp  = NewDouble(NumColsA, "Temp in ComputeADAT()");
  Temp2 = NewDouble(NumRowsA, "Temp2 in ComputeADAT()");

  for(i=0; i<NumRowsA; i++) Temp2[i] = 0.0;

  for (col1 = 0; col1 < NumRowsA; col1++) {	/* for each column of A^T */

    /* copy (scaled) elements of A^T column to Temp */

    base = A->pBeginRowT[col1] - 1;
    for (ent1 = A->pBeginRowT[col1] - 1; ent1 <= A->pEndRowT[col1] - 1; ent1++) {
      row1 = A->RowT[ent1] - 1;
      Temp[row1] = A->ValueT[ent1] * scale[row1];
    }

    /* work down column col1 of A^T, looking at one nonzero entry at a time */
    for(i=A->pBeginRowT[col1]-1; i <= A->pEndRowT[col1] - 1; i++) {

      /* what row of A^T is this element in? */
      row1 = A->RowT[i]-1; 

      /* search across row1 of A^T and identify other columns with
         nonzeros in this row */
      for(j=A->pBeginRow[row1]-1; j <= A->pEndRow[row1]-1; j++) {

	/* what column of A^T is this element is? */
        k = A->Row[j]-1;

	/* OK, update element k of Temp2 */
        Temp2[k] += Temp[row1] * A->Value[j];
      }
    }

    /* finish by storing the elements in AAT and zeroing out the Temp2
       elements as appropriate */

    for (entryAAT = AAT->pBeginRow[col1] - 1;
	 entryAAT <= AAT->pEndRow[col1] - 1; entryAAT++) {
      col2 = AAT->Row[entryAAT] - 1;
      AAT->Value[entryAAT] = Temp2[col2];
      Temp2[col2] = 0.0;
    }
  }

  Free((char *) Temp); 
  Free((char *) Temp2);
  return 0;
}

/*****************************************************************/


