/* wrapper routines for basic sparse linear algebra 
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


/* copies a sparseMatrix type into an MMTtype, but doesn't fill in the three
 * index arrays for the transpose structure (this is done in a later routine
 * ) */

MMTtype        *copyMMT(A, NumRows, NumCols)
  sparseMatrix    A;
  int             NumRows, NumCols;
{
  int             i, Nonzeros;
  MMTtype        *Anew, *NewMMTtype();

  Nonzeros = A.pEndRow[NumCols - 1];
  Anew = NewMMTtype(NumRows, NumCols, Nonzeros);
  for (i = 0; i < NumCols; i++) {
    Anew->pBeginRow[i] = A.pBeginRow[i];
    Anew->pEndRow[i] = A.pEndRow[i];
  }
  for (i = 0; i < Nonzeros; i++) {
    Anew->Row[i] = A.Row[i];
    Anew->Value[i] = A.Value[i];
  }
  return Anew;
}

/* wrapper for RealSparseMatrixVectorProductPlusx() */
int             Axplusb(A, x, b)
  MMTtype        *A;
  double         *x, *b;

{
  return RealSparseMatrixVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, x, b,
     &(A->NumRows), &(A->NumCols));
}


/* wrapper for RealSparseMatrixVectorProduct */
int             Ax(A, x, b)
  MMTtype        *A;
  double         *x, *b;

{
  return RealSparseMatrixVectorProduct
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, x, b,
     &(A->NumRows), &(A->NumCols));
}


/* The following two routines are for A of type sparseMatrix.  These are
 * usually called with LP->A as the argument, where LP is an LPtype data
 * structure.  */

int             SparseSaxpy(A, x, y)
  sparseMatrix    A;
  double         *x, *y;

{
  int             m, n;

  m = A.NumRows;
  n = A.NumCols;

  return RealSparseMatrixVectorProductPlusx
    (A.Value, A.pBeginRow, A.pEndRow, A.Row,
     x, y, &m, &n);
}


int             SparseSaxpyT(A, x, y)
  sparseMatrix    A;
  double         *x, *y;

{
  int             m, n;

  m = A.NumRows;
  n = A.NumCols;

  return RealSparseMatrixTransposeVectorProductPlusx
    (A.Value, A.pBeginRow, A.pEndRow, A.Row,
     x, y, &m, &n);
}


/* The following two routines are for A of type *MMTtype.  */

int             SparseSaxpyM(A, x, y)
  MMTtype        *A;
  double         *x, *y;

{
  int             m, n;

  m = A->NumRows;
  n = A->NumCols;

  return RealSparseMatrixVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
     x, y, &m, &n);
}


int             SparseSaxpyTM(A, x, y)
  MMTtype        *A;
  double         *x, *y;

{
  int             m, n;

  m = A->NumRows;
  n = A->NumCols;

  return RealSparseMatrixTransposeVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
     x, y, &m, &n);
}


double          TwoNorm2(x, n)
  double         *x;
  int            *n;

{
  double          temp;

  if (*n <= 0)
    return 0.0;
  NormTwoSquareRealDenseVector(x, n, &temp);
  return temp;
}

