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

double CurtisReidMeasure();

/*  These routines are for scaling the linear program to reduce the
 *  ratio of the largest and smallest elements in each column and each
 *  row.
 */

int ScaleLP(LP, Inputs)
  LPtype     *LP;
  Parameters *Inputs;
{
  int           Rows, Cols;
  sparseMatrix *A, *AT;
  double       *b, *c;
  double       *UpBound;
  int           NumberBounds;
  int          *BoundIndex;
  double       *RowScale, *ColScale;  /* local pointers */
  double        ScaleFactor;
  int           row, col, ent, i, NonzerosA;
  int           passes, CurtisReidScaling();

  /* Make local copies */

  Rows = LP->Rows;
  Cols = LP->Cols;

  A  = &(LP->A);             b = LP->b;
  AT = &(LP->Atranspose);    c = LP->c;

  UpBound      = LP->UpBound;     
  BoundIndex   = LP->BoundIndex;
  NumberBounds = LP->NumberBounds;

  LP->NumScale = 1;

  LP->ColScale = NewDouble(Cols, "LP->ColScale");
  LP->RowScale = NewDouble(Rows, "LP->RowScale");

  for (col = 0; col < Cols; col++)
    LP->ColScale[col] = 1.0;

  for (row = 0; row < Rows; row++)
    LP->RowScale[row] = 1.0;

  ColScale = NewDouble(Cols, "ColScale");
  RowScale = NewDouble(Rows, "RowScale");

  /* allocate temporary memory */

  ScaleFactor = CurtisReidMeasure(A, Rows, Cols);

  if (Inputs->Diagnostics)
    printf("\nBefore Scaling: ScaleFactor = %.1f\n", ScaleFactor);
  /* Call Curtis-Reid routine */

  NonzerosA= A->pEndRow[Cols-1] - A->pBeginRow[0];
  if (ScaleFactor < 0.1 * NonzerosA) 
    goto cleanup; /* and return */

  passes = CurtisReidScaling(A, AT, Rows, Cols, RowScale, ColScale);

  /* scale columns and rows of A, entries of c */
    
  for (col = 0; col < Cols; col++) {

    c[col] *= ColScale[col];

    for (ent = A->pBeginRow[col]-1; ent <= A->pEndRow[col]-1; ent++) {
      row = A->Row[ent]-1;
      A->Value[ent] *= (ColScale[col] * RowScale[row]);
    }
  }    

  ScaleFactor = CurtisReidMeasure(A, Rows, Cols);
  
  if (Inputs->Diagnostics)
    printf("After  Scaling: ScaleFactor = %.1f   (%d %s)\n\n", 
	   ScaleFactor, passes, (passes == 1)? "pass" : "passes");

  /* For each entry, scale b by rowscale */
    
  for (row = 0; row < Rows; row++)
    b[row] *= RowScale[row];

    /* scale bounds */

  for (ent = 0; ent < NumberBounds; ent++) {
    col = BoundIndex[ent];
    UpBound[col] /= ColScale[col];
  }

  /* Do the same for Atranspose */
    
  for (col = 0; col < Rows; col++) {
    for (ent = AT->pBeginRow[col]-1; ent <= AT->pEndRow[col]-1; ent++) {
      row = AT->Row[ent] - 1;
      AT->Value[ent] *= (ColScale[row] * RowScale[col]);
    }
  }    
    
  for (col = 0; col < Cols; col++)
    LP->ColScale[col] *= ColScale[col];

  for (row = 0; row < Rows; row++)
    LP->RowScale[row] *= RowScale[row];
cleanup:
  Free((char *) RowScale);
  Free((char *) ColScale);
  return 0;
}
/****************************************************************/

double CRM2(A, Rows, Cols, RowScale, ColScale)
sparseMatrix  *A;
int            Rows, Cols;
double        *RowScale, *ColScale;
{
  /*  Compute the scale factor by the formula
   *
   *  ScaleFactor = SUM (log |Aij| - RowScale[i] - ColScale[j]) ^ 2
   */

  double ScaleFactor;
  int    row, col, ent;
  double value, logvalue;

  ScaleFactor = 0.0;

  for (col = 0; col < Cols; col++)
    for (ent = A->pBeginRow[col] - 1; ent <= A->pEndRow[col] - 1; ent++) {
      row = A->Row[ent] - 1;
      value = fabs(A->Value[ent]);
      if (value > 0.0) {
	logvalue = log(value) - RowScale[row] - ColScale[col];
	ScaleFactor += (logvalue * logvalue);
      }
    }
      
  return(ScaleFactor);
}

/****************************************************************/

double CurtisReidMeasure(A, Rows, Cols)
sparseMatrix  *A;
int            Rows, Cols;
{
  /*  Compute the scale factor by the formula
   *
   *  ScaleFactor = SUM (log |Aij|) ^ 2
   */

  double ScaleFactor;
  int    row, col, ent;
  double value, logvalue;

  ScaleFactor = 0.0;

  for (col = 0; col < Cols; col++)
    for (ent = A->pBeginRow[col] - 1; ent <= A->pEndRow[col] - 1; ent++) {
      row = A->Row[ent] - 1;
      value = fabs(A->Value[ent]);
      if (value > 0.0) {
	logvalue = log(value);
	ScaleFactor += (logvalue * logvalue);
      }
    }
      
  return(ScaleFactor);
}

/****************************************************************/

int CurtisReidScaling(A, AT, Rows, Cols, RowScale, ColScale)
sparseMatrix  *A, *AT;
int            Rows, Cols;
double        *RowScale, *ColScale;

/************************************************************
 * Implement Curtis-Reid scaling based on the paper
 * "On the Automatic Scaling of Matrices for Gaussian
 * Elimination," Journal of the Institute of Mathematics and
 * Its Applications (1972) 10, 118-124.
 *
 * Solve the system | M   E | (r)   (sigma)
 *                  |       | ( ) = (     )
 *                  | E^T N | (c)   ( tau )
 *
 * by the conjugate gradient method (clever recurrences).
 *
 * E is the matrix A with all elements = 1
 *
 * M is diagonal matrix of row    counts (RowCount)
 * N is diagonal matrix of column counts (ColCount)
 *
 * sigma is the vector of row    logarithm sums (RowSum)
 * tau   is the vector of column logarithm sums (ColSum)
 *
 * r, c are returned as row and column scalings (RowScale, ColScale)
 ************************************************************/

{
  int     pass, MatrixCount;
  int     row, col, ent;
  double  StopTolerance;
  double *residual_even, *residual_odd;
  double  sk,   qk,   ek;
  double  skm1, qkm1, ekm1;    /* previous iteration */
  double        qkm2, ekm2;    /* 2 iterations ago   */
  double *RowScalem2, *ColScalem2;

  double *RowSum, *ColSum, value, logvalue, check, error;
  int    *RowCount, *ColCount;
  double  CRM2();

  /*  Allocate temporary memory  and 
   *  find RowSum and ColSum measures 
   */

  RowSum   = NewDouble(Rows, "RowSum");
  ColSum   = NewDouble(Cols, "ColSum");

  RowCount = NewInt(Rows, "RowCount");
  ColCount = NewInt(Cols, "ColCount");

  for (row = 0; row < Rows; row++) {
    RowSum[row]   = 0.0;
    RowCount[row] = 0;
  }

  for (col = 0; col < Cols; col++) {
    ColSum[col]   = 0.0;
    ColCount[col] = 0;
  }

  for (col = 0; col < Cols; col++) {

    for (ent = A->pBeginRow[col]-1; ent <= A->pEndRow[col]-1; ent++) {
      row   = A->Row[ent]-1;
      value = fabs(A->Value[ent]);
      
      if (value > 0.0) {
	logvalue = log(value);
	ColSum[col] += logvalue;
	RowSum[row] += logvalue;
	ColCount[col]++;
	RowCount[row]++;
      }
    } /* end for ent */
  }   /* end for col */

  /* Allocate memory for the scaling method */

  residual_even = NewDouble(Cols, "residual_even");
  residual_odd  = NewDouble(Rows, "residual_odd");

  for (col = 0; col < Cols; col++) residual_even[col] = 0.0;
  for (row = 0; row < Rows; row++) residual_odd[row]  = 0.0;

  RowScalem2 = NewDouble(Rows, "RowScalem2");
  ColScalem2 = NewDouble(Cols, "ColScalem2");
  
  /* Initialize 
   * RowScale = RowCount^-1 RowSum 
   * ColScale = 0.0
   * residual = ColSum - E^T RowCount^-1 RowSum
   */

  MatrixCount = A->pEndRow[Cols-1] - A->pBeginRow[0];

  StopTolerance = 1.0e-2 * (double) MatrixCount;

  for (row = 0; row < Rows; row++) {
    RowScale[row]   = RowSum[row] / (double) RowCount[row];
    RowScalem2[row] = RowScale[row];
  }

  for (col = 0; col < Cols; col++) {
    ColScale[col]   = 0.0;
    ColScalem2[col] = 0.0;
  }

  /* compute initial residual */

  for (col = 0; col < Cols; col++)
    residual_even[col] = ColSum[col];

  for (col = 0; col < Rows; col++)
    for (ent = AT->pBeginRow[col] - 1; ent <= AT->pEndRow[col] - 1; ent++) {
      row = AT->Row[ent] - 1;
      residual_even[row] -= RowSum[col] / (double) RowCount[col];
    }

  /* compute sk */
  skm1 = 0.0;
  sk   = 0.0;
  for (col = 0; col < Cols; col++)
    sk += (residual_even[col] * residual_even[col] / (double) ColCount[col]);

  pass = 0;
  qk   = 1.0;  qkm1 = 0.0;  qkm2 = 0.0;
  ek   = 0.0;  ekm1 = 0.0;  ekm2 = 0.0;

  while (sk > StopTolerance) {

    /* Given the values of residual and sk, construct
     *  ColScale (when pass is even)
     *  RowScale (when pass is odd)
     */

    switch (pass % 2) {
    case 0: /* pass is even */

      /* construct RowScale[pass+1] */

      if (pass != 0) {  /* RowScale[1] = RowScale[0] */

	for (row = 0; row < Rows; row++)
	  RowScalem2[row] = RowScale[row];

	for (row = 0; row < Rows; row++)
	  RowScale[row] *= (1.0 + ek * ekm1 / (qk * qkm1));

	for (row = 0; row < Rows; row++)
	  RowScale[row] += 
	    (residual_odd[row] / (qk * qkm1 * (double) RowCount[row]) -
	    RowScalem2[row] * ek * ekm1 / (qk * qkm1));

      }
      break;

    case 1: /* pass is odd */

      /* construct ColScale[pass+1] */

      for (col = 0; col < Cols; col++)
	ColScalem2[col] = ColScale[col];

      for (col = 0; col < Cols; col++)
	ColScale[col] *= (1.0 + ek * ekm1 / (qk * qkm1));

      for (col = 0; col < Cols; col++)
	ColScale[col] += 
	  (residual_even[col] / ((double) ColCount[col] * qk * qkm1) -
	  ColScalem2[col] * ek * ekm1 / (qk * qkm1));

      break;
    } /* end switch */
  
    /* update residual and sk (pass + 1) */

    switch (pass % 2) {
    case 0:  /* even */
      /* residual */
      for (row = 0; row < Rows; row++)
	residual_odd[row] *= ek;

      for (col = 0; col < Cols; col++) 
	for (ent = A->pBeginRow[col] - 1; ent <= A->pEndRow[col] - 1; ent++) {
	  row = A->Row[ent] - 1;
	  residual_odd[row] += (residual_even[col] / (double) ColCount[col]);
	}

      for (row = 0; row < Rows; row++)
	residual_odd[row] *= (-1.0 / qk);

      /* sk */

      skm1 = sk;

      sk = 0.0;
      for (row = 0; row < Rows; row++)
	sk += (residual_odd[row] * residual_odd[row] / (double) RowCount[row]);

      break;

    case 1: /* odd */
      /* residual */
      for (col = 0; col < Cols; col++)
	residual_even[col] *= ek;

      for (col = 0; col < Rows; col++) 
	for (ent = AT->pBeginRow[col] - 1; ent <= AT->pEndRow[col] - 1; ent++) {
	  row = AT->Row[ent] - 1;
	  residual_even[row] += (residual_odd[col] / (double) RowCount[col]);
	}

      for (col = 0; col < Cols; col++)
	residual_even[col] *= (-1.0 / qk);

      /* sk */

      skm1 = sk;

      sk = 0.0;
      for (col = 0; col < Cols; col++)
	sk += (residual_even[col] * residual_even[col] / (double) ColCount[col]);

      break;

    } /* end switch */

    /* compute ek and qk */

    ekm2 = ekm1;
    ekm1 = ek;
    ek   = qk * sk / skm1;

    qkm2 = qkm1;
    qkm1 = qk;
    qk   = 1.0 - ek;

    /*
    printf("sk = %g  skm1 = %g   Objective = %g\n", 
       sk, skm1, CRM2(A, Rows, Cols, RowScale, ColScale));
    printf("ek = %g  ekm1 = %g  ekm2 = %g\n", ek, ekm1, ekm2);
    printf("qk = %g  qkm1 = %g  qkm2 = %g\n", qk, qkm1, qkm2);
    */

    pass++;
  } /* end while */

  /* Synchronize the RowScale and ColScale vectors */

  if (pass > 0)
    {
      switch (pass % 2) {
      case 0: /* pass is even, compute RowScale */
	
	for (row = 0; row < Rows; row++)
	  RowScale[row] *= (1.0 + ek * ekm1 / qkm1);
	
	for (row = 0; row < Rows; row++)
	  RowScale[row] += (residual_odd[row] / (qkm1 * (double) RowCount[row]) -
			    RowScalem2[row] * ek * ekm1 / qkm1);
        
	break;
	
      case 1: /* pass is odd, compute ColScale */
	
	for (col = 0; col < Cols; col++)
	  ColScale[col] *= (1.0 + ek * ekm1 / qkm1);
	
	for (col = 0; col < Cols; col++)
	  ColScale[col] += (residual_even[col] / ((double) ColCount[col] * qkm1) -
			    ColScalem2[col] * ek * ekm1 / qkm1);
	
	break;
     
 }
    }

  /* CHECK */

  /* M RowScale + E ColScale = RowSum */

  error = 0.0;

  for (row = 0; row < Rows; row++) {

    check = (double) RowCount[row] * RowScale[row];

    for (ent = AT->pBeginRow[row]-1; ent <= AT->pEndRow[row]-1; ent++) {
      col = AT->Row[ent] - 1;
      check += ColScale[col];
    }
    check -= RowSum[row];
    error += (check * check);
  }

  /* E^T RowScale + N ColScale = ColSum */

  for (col = 0; col < Cols; col++) {

    check = (double) ColCount[col] * ColScale[col];

    for (ent = A->pBeginRow[col]-1; ent <= A->pEndRow[col]-1; ent++) {
      row = A->Row[ent]-1;
      check += RowScale[row];
    }
    check -= ColSum[col];
    error += (check * check);
  }

  /* Convert to scaling factors */

  for (col = 0; col < Cols; col++) {
    ColScale[col] = exp(-ColScale[col]);
    if (ColScale[col] < 1.0e-8) ColScale[col] = 1.e-8;
  }

  for (row = 0; row < Rows; row++) {
    RowScale[row] = exp(-RowScale[row]);
    if (RowScale[row] < 1.0e-8) RowScale[row] = 1.0e-8;
  }

  /* free temporary memory */

  Free((char *) residual_even);  Free((char *) residual_odd);
  Free((char *) RowScalem2);     Free((char *) ColScalem2);
  Free((char *) RowSum);         Free((char *) ColSum);
  Free((char *) RowCount);       Free((char *) ColCount);

  return pass;
}

/****************************************************************/

int UnscaleLP(LP, Solution)
  LPtype    *LP;
  solution  *Solution;
{
  int           Rows, Cols;
  sparseMatrix *A, *AT;
  double       *b, *c;
  double       *UpBound;
  int           NumberBounds;
  int          *BoundIndex;
  
  int           row, col, ent, i;

  /* Make local copies */

  Rows = LP->Rows;
  Cols = LP->Cols;

  A  = &(LP->A);             b = LP->b;
  AT = &(LP->Atranspose);    c = LP->c;

  UpBound      = LP->UpBound;     
  BoundIndex   = LP->BoundIndex;
  NumberBounds = LP->NumberBounds;

  /* undo the scalings */

  for (col = 0; col < Cols; col++) {
    c[col] /= LP->ColScale[col];
    Solution->x[col]         *= LP->ColScale[col];
    Solution->DualLower[col] /= LP->ColScale[col];
    
    for (ent = A->pBeginRow[col]-1; ent <= A->pEndRow[col]-1; ent++) {
      row = A->Row[ent]-1;
      A->Value[ent] /= (LP->ColScale[col] * LP->RowScale[row]);
    }
  }    
  
  for (row = 0; row < Rows; row++) {
    b[row] /= LP->RowScale[row];
    Solution->pi[row] *= LP->RowScale[row];
  }
  
  /* scale bounds */
  
  for (ent = 0; ent < NumberBounds; ent++) {
    col = BoundIndex[ent];
    UpBound[col] *= LP->ColScale[col];
    Solution->DualUpper[col] /= LP->ColScale[col];
  }
  
  /* Do the same for Atranspose */
  
  for (col = 0; col < Rows; col++) {
    for (ent = AT->pBeginRow[col]-1;
       ent <= AT->pEndRow[col]-1 && (row = AT->Row[ent] - 1) < Cols;
       ent++) {
      AT->Value[ent] /= (LP->ColScale[row] * LP->RowScale[col]);
    }
  }    
  return 0;
}  /* end for LP->NumScale loop */
