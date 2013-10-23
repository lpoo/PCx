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

TestGondzio
  (A, Factor, scale, Current, Predictor, Corrector, comp_xs, comp_rw)
  MMTtype        *A;

  NgPeyton       *Factor;

  double         *scale;
  Iterate        *Current, *Predictor, *Corrector;
  double         *comp_xs, *comp_rw;
{
  int             SolveAugmented();
  int             NumRows, NumCols, NumBounds, i, irow;
  int             CorrectionNumber;
  int             ContinueCorrecting;
  double         *tmp_col, *tmp_row;

  double         *x, *s, *pi, *r, *w;
  double         *dx, *ds, *dpi, *dr, *dw;
  double         *dx2, *ds2, *dpi2, *dr2, *dw2;

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

  tmp_row = NewDouble(NumRows, "tmp_row");
  tmp_col = NewDouble(NumCols, "tmp_row");

  /* First equation */

  RealSparseMatrixVectorProduct
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, dx2, tmp_row,
     &(A->NumRows), &(A->NumCols));

  printf("First Equation:\n");
  for (i = 0; i < NumRows; i++)
    if (fabs(tmp_row[i]) > 1.0e-8)
      printf("tmp_row[%d] = %f\n", i, tmp_row[i]);

  /* Second equation */

  RealSparseMatrixTransposeVectorProduct
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, dpi2, tmp_col,
     &(A->NumRows), &(A->NumCols));

  for (i = 0; i < NumCols; i++)
    tmp_col[i] += ds2[i];

  for (i = 0; i < NumBounds; i++) {
    irow = Current->BoundIndex[i] - 1;
    tmp_col[irow] -= dr2[i];
  }

  printf("Second Equation:\n");
  for (i = 0; i < NumCols; i++)
    if (fabs(tmp_col[i]) > 1.0e-8)
      printf("tmp_col[%d] = %f\n", i, tmp_col[i]);

  /* third equation */

  printf("Third Equation:\n");
  for (i = 0; i < NumBounds; i++) {
    irow = Current->BoundIndex[i] - 1;
    if (fabs(dx2[irow] + dw2[i]) > 1.0e-8)
      printf("dx2[%d] = %f    dw2 = %f\n", i, dx2[irow], dw2[i]);
  }

  /* fourth equation */

  printf("Fourth Equation:\n");

  for (i = 0; i < NumCols; i++)
    tmp_col[i] = s[i] * dx2[i] + x[i] * ds2[i];

  for (i = 0; i < NumCols; i++)
    if (tmp_col[i] + comp_xs[i] > 1.0e-8)
      printf("%d: tmp_col = %f  comp_xs = %f\n", i, tmp_col[i], comp_xs[i]);

  /* fifth equation */

  printf("Fifth equation:\n");

  for (i = 0; i < NumBounds; i++)
    tmp_col[i] = r[i] * dw2[i] + w[i] * dr2[i];

  for (i = 0; i < NumBounds; i++)
    if (tmp_col[i] - comp_rw[i] > 1.0e-8)
      printf("%d: tmp_col = %f  comp_rw = %f\n", tmp_col[i], comp_rw[i]);

  Free((char *) tmp_col);
  Free((char *) tmp_row);

  return 0;
}
