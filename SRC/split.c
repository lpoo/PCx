/* split free variables 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "main.h"

int SplitFreeVars(LP)
  LPtype         *LP;
{
  int             col, NextCol, FreeCounter, NextEnt, entry, NumberFree,
                  NewCols, NewEnts, MoreEnts, NextBound, NewBounds;
  LPtype         *MakeTranspose();

  /* count number of free variables */

  for (col = 0, NumberFree = 0, MoreEnts = 0; col < LP->Cols; col++)
    if (LP->VarType[col] == FREE) {
      NumberFree++;
      MoreEnts += (LP->A.pEndRow[col] - LP->A.pBeginRow[col] + 1);
    }
  /* Reallocate memory to allow for larger vectors, etc. */

  NewEnts = LP->Ents + MoreEnts;
  NewCols = LP->Cols + NumberFree;

  if (LP->FreePlus) Free((char *) LP->FreePlus);
  LP->FreePlus = NewInt(NumberFree, "LP->FreePlus");
  if (LP->FreeMinus) Free((char *) LP->FreeMinus);
  LP->FreeMinus = NewInt(NumberFree, "LP->FreeMinus");

  LP->c = (double *)
    Realloc(LP->c, NewCols * sizeof(double), "LP->c");
  LP->VarType = (int *)
    Realloc(LP->VarType, NewCols * sizeof(int), "LP->VarType");
  LP->UpBound = (double *)
    Realloc(LP->UpBound, NewCols * sizeof(double), "LP->UpBound");

  LP->A.pBeginRow = (int *)
    Realloc(LP->A.pBeginRow, NewCols * sizeof(int), "LP->A.pBeginRow");
  LP->A.pEndRow = (int *)
    Realloc(LP->A.pEndRow, NewCols * sizeof(int), "LP->A.pEndRow");
  LP->A.Row = (int *)
    Realloc(LP->A.Row, NewEnts * sizeof(int), "LP->A.Row");
  LP->A.Value = (double *)
    Realloc(LP->A.Value, NewEnts * sizeof(double), "LP->A.Value");

  LP->A.NumCols = NewCols; LP->A.Nonzeros = NewEnts;

  /* LP->A.transpose already has the right number of columns, so just
     allocate extra space for its additional nonzeros */

  LP->Atranspose.Row = (int *)
    Realloc(LP->Atranspose.Row, NewEnts * sizeof(int), 
                                     "LP->Atranspose.Row");
  LP->Atranspose.Value = (double *)
    Realloc(LP->Atranspose.Value, NewEnts * sizeof(double), 
                                     "LP->Atranspose.Value");
 
  LP->Atranspose.NumRows = NewCols; LP->Atranspose.Nonzeros = NewEnts;

  NextCol = LP->Cols; NextEnt = LP->Ents;
  for (col = 0, FreeCounter = 0; col < LP->Cols; col++)
    if (LP->VarType[col] == FREE) {

      /* reset this column to NORMAL */
      LP->VarType[col] = NORMAL;
      LP->FreePlus[FreeCounter] = col;

      /* make additional variable such that x = x^+ - x^- */
      LP->VarType[NextCol] = NORMAL;	/* add new split variable at end */
      LP->FreeMinus[FreeCounter] = NextCol;

      /* copy constraint matrix elements from original variable - flip signs */
      LP->A.pBeginRow[NextCol] = NextEnt + 1;
      for (entry = LP->A.pBeginRow[col] - 1; entry <= LP->A.pEndRow[col] - 1; entry++) {
	LP->A.Row[NextEnt] = LP->A.Row[entry];
	LP->A.Value[NextEnt] = -LP->A.Value[entry];
	NextEnt++;
      }
      LP->A.pEndRow[NextCol] = NextEnt;
      LP->c[NextCol] = -LP->c[col];
      NextCol++;
      FreeCounter++;
    }
  LP->Cols = NewCols;
  LP->A.NumCols = NewCols;
  LP->Ents = NewEnts;
  LP->A.Nonzeros = NewEnts;
  LP->NumberFree = 0;
  LP->NumberSplit = NumberFree;

  if (NextEnt != NewEnts) {
    printf("Error splitting free variables: ");
    printf("NextEnt (%d) != NewEnts (%d)\n", NextEnt, NewEnts);
  }
  if (NextCol != NewCols) {
    printf("Error splitting free variables: ");
    printf("NextCol (%d) != NewCols (%d)\n", NextCol, NewCols);
  }
  if (FreeCounter != NumberFree) {
    printf("Error splitting free variables: ");
    printf("FreeCounter (%d) != NumberFree (%d)\n", FreeCounter, NumberFree);
  }
  LP = MakeTranspose(LP);
  return 0;
}

int UnSplitFreeVars(LP, Solution)
  LPtype         *LP;
  solution       *Solution;
{
  int             i, plus, minus;

  for (i = 0; i < LP->NumberSplit; i++) {
    plus = LP->FreePlus[i];
    minus = LP->FreeMinus[i];
    Solution->x[plus] = Solution->x[plus] - Solution->x[minus];
    Solution->DualLower[plus] = 0.0;
    LP->VarType[plus] = FREE;
  }
  LP->Cols -= LP->NumberSplit;
  LP->A.NumCols -= LP->NumberSplit;
  Solution->Columns -= LP->NumberSplit;
  return 0;
}

/* Implement the split-shift heuristic for dealing with free vars */

int ShiftSplitVariables(LP, Current)
  LPtype         *LP;
  Iterate        *Current;
{
  int             i, im, ip;
  int             NumCols, NumBounds;
  double          mutemp=0.0, splitshift, xmax, xmin, 
                  goal, avx, sxmax=0.0, smin=10.e0, smax=100.e0;

  NumCols = LP->Cols;
  NumBounds = LP->NumberBounds;

  mutemp = 0.0; avx = 0.0;
  for (i = 0; i < NumCols; i++) {
    mutemp += Current->x[i] * Current->s[i];
    avx += Current->x[i];
  }
  for (i = 0; i < NumBounds; i++) {
    mutemp += Current->w[i] * Current->r[i];
    avx += Current->w[i];
  }

  /* Shifts each split-variable pair so that the smaller of the two takes on
   * the value MID(smin, larger-smaller, smax).  Adjusts the s components
   * accordingly.  */

  for (i = 0; i < LP->NumberSplit; i++) {
    im = LP->FreeMinus[i];
    ip = LP->FreePlus[i];    
    mutemp -= ( Current->x[ip] * Current->s[ip] + 
                Current->x[im] * Current->s[im]);
    avx -= (Current->x[ip] + Current->x[im]);
    sxmax = MAX(Current->x[ip], Current->x[im]);
  }

  mutemp /= (NumCols + NumBounds - 2*LP->NumberSplit);
  avx    /= (NumCols + NumBounds - 2*LP->NumberSplit);

  for (i = 0; i < LP->NumberSplit; i++) {
    im = LP->FreeMinus[i];
    ip = LP->FreePlus[i];    

    xmax = MAX(Current->x[ip], Current->x[im]);
    xmin = MIN(Current->x[ip], Current->x[im]);

    /* shift so that the geometric mean of the shifted var equals avx */
    splitshift = (xmax-xmin)*(xmax-xmin) + 4.0*avx*avx;
    splitshift = sqrt(splitshift) - (xmax-xmin);
    splitshift/= 2.0;

    if (splitshift > smax) splitshift = smax;
    if (splitshift < smin) splitshift = smin;

    splitshift  = xmin - splitshift;
    Current->x[ip] -= splitshift;
    Current->x[im] -= splitshift;
    }

  /* finally, adjust duals so that the average pairwise product is
     achieved by the split vars */

  goal = 1.e-2;
  for (i = 0; i < LP->NumberSplit; i++) {
    im = LP->FreeMinus[i];
    ip = LP->FreePlus[i];    

    if (Current->s[im] * Current->x[im] < goal*mutemp)
       Current->s[im] =  goal * mutemp / Current->x[im]; 
    if (Current->s[ip] * Current->x[ip] < goal*mutemp)
       Current->s[ip] =  goal * mutemp / Current->x[ip];
  }

  return 0;
}



