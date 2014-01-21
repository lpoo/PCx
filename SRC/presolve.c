/* Preprocessor
 *
 * PCx 1.1 11/97
 * 
 * Author: Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "main.h"
#include "pre.h"

/* Preprocess an LP supplied in the form
 * 
 * min c.x ,   Ax=b,
 * 
 * where the components of x fall into one of three classes:
 * 
 * FREE NORMAL  0 <= x_i UPPER   0 <= x_i <= u_i
 * 
 * Stores and return a list of changes, so that various components of the 
 * primal and dual solution can be recovered later.
 * 
 * On entry: LP points to a linear program, in the PCx data structure;
 * pReducedLP, pRecord are ignored;
 * 
 * On return: LP is unchanged; pReducedLP points to the reduced LP; pRecord
 * points to a record of all the changes (needed for postprocessing);
 * Preprocess() returns the number of Passes, or one of the three flags
 * INFEASIBLE, UNBOUNDED, or PREPROCESSING_ERROR.
 * 
 * SJW December 94. */

SparseVector   *NewSparseVector(nonzeros)
  int           nonzeros;
{
  SparseVector *sparseVector = (SparseVector *)
	  Malloc(sizeof(SparseVector), "sparseVector");
  sparseVector->size = nonzeros;
  sparseVector->Index = NewInt(nonzeros, "sparseVector->Index");
  sparseVector->Value = NewDouble(nonzeros, "sparseVector->Value");
  return sparseVector;
}

FreeSparseVector(sparseVector)
  SparseVector		*sparseVector;
{
  Free((char*) sparseVector->Value);
  Free((char*) sparseVector->Index);
  Free((char*) sparseVector);
  return 0;
}

int             Preprocess(LP, pReducedLP, pRecord, Inputs)
  LPtype         *LP;
  LPtype        **pReducedLP;
  ChangeStack   **pRecord; 
  Parameters     *Inputs;

{
  int             i, j, k, inx, inxrc, inxt, inxk, size, sizeUnit, rows,
                  cols, top, nonzeros, nonzerosc, singleton, infeasible = 0,
                  unbounded = 0, preprocessing_error = 0, Pass = 0, ReducedOnThisPass=0,
                  implied_bound_inf, still_alive, *nonzerosRow, *nonzerosCol;
  double          xfix, *cNew, *bNew, implied_bound, bLower, bUpper, bLower_inf,
                  bUpper_inf;
  LPtype         *repackLP(), *ReducedLP;
  ChangeStack    *NewRecord(), *Record;
  int             ResizeRecord();


  /* declare the preprocessing component routines */
  int             ObviouslyInfeasible(), DuplicateCols(), nDuplicateCols, 
                  DuplicateRows(), nDuplicateRows;
  void            FindEmptyRows(), FindEmptyCols(), 
                  FindFixedVariables(), FindSingletonRows(),
                  FindColumnSingletons(), FindForcedRows(),
                  DeleteChangeStack();

  /* Allocate a new Record */
  *pRecord = NewRecord(LP);
  Record = *pRecord;

  rows = LP->Rows;
  cols = LP->Cols;
  for (i = 0; i < cols; i++)
    Record->VarType[i] = LP->VarType[i];

  if (Inputs->Diagnostics > 0)
    printf("\nBefore Presolving:  %d rows, %d columns\n",
	   rows, cols);

  /* Initialize the modified cost and rhs vectors */
  cNew = NewDouble(cols, "cNew");
  bNew = NewDouble(rows, "bNew");

  for (i = 0; i < rows; i++)
    bNew[i] = LP->b[i];
  for (i = 0; i < cols; i++)
    cNew[i] = LP->c[i];

  /* Initialize  the nonzero counters */
  nonzerosCol = NewInt(cols, "nonzerosCol");
  nonzerosRow = NewInt(rows, "nonzerosRow");

  for (i = 0; i < rows; i++) {
    nonzerosRow[i] =
      LP->Atranspose.pEndRow[i] - LP->Atranspose.pBeginRow[i] + 1;
  }

  for (i = 0; i < cols; i++) {
    nonzerosCol[i] = LP->A.pEndRow[i] - LP->A.pBeginRow[i] + 1;
  }

  /* EACH PASS STARTS HERE. Continue to loop until no more reductions are
   * made on this pass. */

  for (Pass = 0; Pass == 0 || ReducedOnThisPass >= 1;) {

    Pass++;
    ReducedOnThisPass = 0;
    top = Record->Top;

    /* Look for DUPLICATE COLUMNS */

    nDuplicateCols = DuplicateCols(LP, Record, cNew, nonzerosRow, nonzerosCol,
				   &top, Pass, &ReducedOnThisPass);

    /* it's possible that unboundedness was detected here - exit if so */
    if (nDuplicateCols < 0) {
      printf(" ERROR in Preprocess(): Uncoundedness detected in DuplicateCols()\n");
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return UNBOUNDED;
    }

    /* Look for DUPLICATE ROWS */

    nDuplicateRows = DuplicateRows(LP, Record, bNew, nonzerosRow, nonzerosCol,
				   &top, Pass, &ReducedOnThisPass);

    /* Look for INFEASIBILITY: lower bounds > upper bounds */

    if (infeasible = ObviouslyInfeasible(LP, bNew, nonzerosRow, Record)) {
      printf(" ERROR in Preprocess(): detected %d infeasible bounds\n", infeasible);
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return INFEASIBLE;
    }
    /* look for EMPTY OR ZERO rows */
    FindEmptyRows(LP, Record, nonzerosRow, &top, Pass, &ReducedOnThisPass);

    /* look for EMPTY OR ZERO columns. if we find one, can immediately fix
     * the corresponding x component at its upper bound if c<0, at its lower
     * bound if c>0, at an arbitrary value if c=0. declare the problem
     * unbounded if these strategies give an unbounded objective.  */

    FindEmptyCols(LP, Record, nonzerosCol, cNew, &unbounded,
		  &top, Pass, &ReducedOnThisPass);

    if (unbounded) {
      printf("ERROR in Preprocess(): Unbounded primal objective\n");
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return UNBOUNDED;
    }

    /* Look for FIXED VARIABLES. Don't do any more checking for consistency
     * of upper and lower bounds; this was done in the loop immediately
     * above. */

    FindFixedVariables(LP, Record, nonzerosRow, nonzerosCol,
		       cNew, bNew, &top, Pass, &ReducedOnThisPass);

    /* Look for SINGLETON ROWS. These are really the same as fixed variables,
     * except that we get to delete a row too. */

    FindSingletonRows(LP, Record, nonzerosRow, nonzerosCol,
		      cNew, bNew, &preprocessing_error, &infeasible,
		      &top, Pass, &ReducedOnThisPass);

    if (preprocessing_error) {
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return PREPROCESSING_ERROR;
    }
    if (infeasible) {
      printf(" ERROR in Preprocess(): detected infeasibility\n", infeasible);
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return INFEASIBLE;
    }

    /* look for FORCED ROWS. These are constraints which imply that all
     * variables involved must be at one of their bounds. e.g. x1 + x2 - x3
     * =2,   x1<=1, x2<=1, x3>=0, forces each variable in the row to be at
     * its bound */

    FindForcedRows(LP, Record, nonzerosRow, nonzerosCol,
		   cNew, bNew, &top, Pass, &ReducedOnThisPass);

    /* look for FREE COLUMN SINGLETONS and IMPLIED FREE COLUMN SINGLETONS. We
     * get to express a variable in terms of other variables, so delete a row
     * and column and transform the objective funtion */

    FindColumnSingletons(LP, Record, nonzerosRow, nonzerosCol,
			 cNew, bNew, &preprocessing_error, &top,
			 Pass, &ReducedOnThisPass);
    if (preprocessing_error) {
      Free((char *) cNew);        Free((char *) bNew);
      Free((char *) nonzerosCol); Free((char *) nonzerosRow);
      DeleteChangeStack(Record);
      return PREPROCESSING_ERROR;
    }

    /* End of this pass */
#ifdef PREPROCESS_VERBOSE
    printf("**** %d reductions on pass %d\n", ReducedOnThisPass, Pass);
#endif
  }
  Pass--;

  if (Pass != 0) {
    /* repack the reduced problem */
    *pReducedLP = repackLP(LP, Record->RowMask, Record->ColumnMask,
			   Record->VarType, cNew, bNew, Record->cshift);
    ReducedLP = *pReducedLP;
  } else
    *pReducedLP = LP;


  if (Inputs->Diagnostics > 0)
    printf("After  Presolving:  %d rows, %d columns (%d %s)\n\n",
	   (*pReducedLP)->Rows, (*pReducedLP)->Cols, Pass,
	   (Pass == 1)? "pass" : "passes");

  Record->Passes = Pass;
  Record->ReducedRows = (*pReducedLP)->Rows;
  Record->ReducedColumns = (*pReducedLP)->Cols;

  Free((char *) cNew);        Free((char *) bNew);
  Free((char *) nonzerosCol); Free((char *) nonzerosRow);
  return Pass;
}

/* repacks a processed LP, removing flagged rows and columns and replacing
 * the c, b, and cshift terms as appropriate.  SJW 12/94.  */

LPtype         *repackLP(LP, RowMask, ColumnMask, VarType, cNew, bNew, cshift)
  LPtype         *LP;
  int            *RowMask, *ColumnMask, *VarType;
  double         *cNew, *bNew, cshift;
{
  LPtype         *LPnew, *MakeTranspose();
  int             i, j, k, count, colcount, inx;
  int             rows, cols, ents;
  int            *RowMap, *ColumnMap;
  int             numberBounds, numberFree;

  LPnew = (LPtype *) Malloc(sizeof(LPtype), "LPnew");

  /* count up the new rows, cols, ents, and find the mappings between old and
   * new index sets */

  RowMap = NewInt(LP->Rows + 1, "RowMap");
  ColumnMap = NewInt(LP->Cols + 1, "ColumnMap");

  for (i = 1, rows = 0; i <= LP->Rows; i++)
    if (RowMask[i] == STILL_ACTIVE)
      RowMap[i] = ++rows;
    else
      RowMap[i] = 0;
  for (i = 1, cols = 0; i <= LP->Cols; i++)
    if (ColumnMask[i] == STILL_ACTIVE)
      ColumnMap[i] = ++cols;
    else
      ColumnMap[i] = 0;

  ents = 0;
  for (i = 0; i < LP->Cols; i++)
    if (ColumnMask[i + 1] == STILL_ACTIVE)
      for (j = LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++) {
	inx = LP->A.Row[j - 1];
	if (RowMask[inx] == STILL_ACTIVE)
	  ents++;
      }

  /* now can allocate space for the new LP */

  LPnew->Rows = rows;
  LPnew->Cols = cols;
  LPnew->Ents = ents;
  LPnew->A.NumRows = rows;
  LPnew->Atranspose.NumRows = cols;
  LPnew->A.NumCols = cols;
  LPnew->Atranspose.NumCols = rows;
  LPnew->A.Nonzeros = ents;
  LPnew->Atranspose.Nonzeros = ents;
  LPnew->cshift = cshift;

  LPnew->FreePlus = NULL;
  LPnew->FreeMinus = NULL;
  
  LPnew->RowScale = NULL;
  LPnew->ColScale = NULL;

  LPnew->b = NewDouble(rows, "LPnew->b");
  LPnew->c = NewDouble(cols, "LPnew->c");

  LPnew->VarType = NewInt(cols, "LPnew->VarType");
  LPnew->UpBound = NewDouble(cols, "LPnew->UpBound");

  LPnew->A.pBeginRow = NewInt(cols, "LPnew->A.pBeginRow");
  LPnew->A.pEndRow = NewInt(cols, "LPnew->A.pEndRow");
  LPnew->A.Row = NewInt(ents, "LPnew->A.Row");
  LPnew->A.Value = NewDouble(ents, "LPnew->A.Value");

  LPnew->Atranspose.pBeginRow =
    NewInt(rows, "LPnew->Atranspose.pBeginRow");
  LPnew->Atranspose.pEndRow =
    NewInt(rows, "LPnew->Atranspose.pEndRow");
  LPnew->Atranspose.Row =
    NewInt(ents, "LPnew->Atranspose.Row");
  LPnew->Atranspose.Value =
    NewDouble(ents, "LPnew->Atranspose.Value");

  /* and shift the data */
  numberBounds = 0;
  numberFree = 0;
  for (i = 0, count = 0; i < LP->Cols; i++)
    if (ColumnMask[i + 1] == STILL_ACTIVE) {
      LPnew->c[count] = cNew[i];
      LPnew->UpBound[count] = LP->UpBound[i];
      LPnew->VarType[count] = VarType[i];
      if (VarType[i] == FREE)
	numberFree++;
      if (VarType[i] == UPPER)
	numberBounds++;
      count++;
    }
  /* fill out array of free and bounded indices */
  LPnew->NumberBounds = numberBounds;
  LPnew->NumberFree = numberFree;
  LPnew->BoundIndex = NewInt(numberBounds, "LPnew->BoundIndex");
  LPnew->FreeIndex = NewInt(numberFree, "LPnew->FreeIndex");
  for (i = 0, numberBounds = 0, numberFree = 0; i < LPnew->Cols; i++) {
    if (LPnew->VarType[i] == FREE) {
      LPnew->FreeIndex[numberFree] = i;
      numberFree++;
    } else if (LPnew->VarType[i] == UPPER) {
      LPnew->BoundIndex[numberBounds] = i;
      numberBounds++;
    }
  }

  for (i = 0, count = 0; i < LP->Rows; i++)
    if (RowMask[i + 1] == STILL_ACTIVE) {
      LPnew->b[count] = bNew[i];
      count++;
    }
  for (i = 0, colcount = 0, count = 0; i < LP->Cols; i++)
    if (ColumnMask[i + 1] == STILL_ACTIVE) {
      /* add this column to the reduced LP */
      LPnew->A.pBeginRow[colcount] = count + 1;
      for (j = LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++) {
	inx = LP->A.Row[j - 1];
	if (RowMask[inx] == STILL_ACTIVE) {
	  /* add this element to the reduced matrix */
	  LPnew->A.Row[count] = RowMap[inx];
	  LPnew->A.Value[count] = LP->A.Value[j - 1];
	  count++;
	}
      }
      LPnew->A.pEndRow[colcount] = count;
      colcount++;
    }
  /* Now call the utility routine to fill in Atranspose */
  LPnew = MakeTranspose(LPnew);

  /* all done! */

  Free((char *) RowMap);
  Free((char *) ColumnMap);
  return (LPnew);
}

/* Check for INFEASIBILITY */

int             ObviouslyInfeasible(LP, bNew, nonzerosRow, Record)
  LPtype         *LP;
  double         *bNew;
  int            *nonzerosRow;  
  ChangeStack    *Record;
{
  int             infeasible, i, cols, rows;

  infeasible = 0;
  cols = LP->Cols;
  rows = LP->Rows;

  /* check for an infeasible bound range */
  for (i = 0; i < cols; i++)
    if (Record->ColumnMask[i+1] == STILL_ACTIVE) {
      if (Record->VarType[i] == UPPER) {
	if (LP->UpBound[i] < -TINY) {
	  printf("UpBound[%d] = %e < 0.0\n", i, LP->UpBound[i]);
	  infeasible++;
	}
      }
    }

  /* check for a zero row with a nonzero right-hand side */
  for (i = 0; i < rows; i++)
    if (Record->RowMask[i + 1] == STILL_ACTIVE) {
      if (nonzerosRow[i] == 0 && bNew[i] != 0) {
         printf("Zero row %d has rhs %e != 0.0\n", i, bNew[i]);
         infeasible++;
      }
    }

  return infeasible;
}


/* Look for EMPTY OR ZERO rows */

void            FindEmptyRows(LP, Record, nonzerosRow, top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  int            *nonzerosRow, *top, Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, inx, nonzeros;
  ZeroRow        *pZeroRow;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;

  for (i = 0; i < rows; i++)
    if (Record->RowMask[i + 1] == STILL_ACTIVE) {
      if (nonzerosRow[i] == 0)
	Record->RowMask[i + 1] = Pass;

      /* make a Row deletion record, push it on the stack */

      if (Record->RowMask[i + 1] == Pass) {
#ifdef PREPROCESS_VERBOSE
	printf("Empty Row %d\n", i + 1);
#endif
	(*ReducedOnThisPass)++;
	Record->StackOfChanges[*top] = (SingleChange *)
	  Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	Record->StackOfChanges[*top]->ChangeType = ZERO_ROW;
	pZeroRow = (ZeroRow *) Malloc(sizeof(ZeroRow), "pZeroRow");
	pZeroRow->RowIndex = i + 1;
	Record->StackOfChanges[*top]->pZeroRow = pZeroRow;

	/* increment top-of-stack pointer, resizing the stack if necessary */
	(*top)++;
	Record->Top = *top;
	if (Record->Top >= Record->Size)
	  ResizeRecord(Record);
      }
    }
  return;
}

/* Look for EMPTY OR ZERO columns. if we find one, can immediately fix the
 * corresponding x component at its upper bound if c<0, at its lower bound if
 * c>0, at an arbitrary value if c=0. declare the problem unbounded if these
 * strategies give an unbounded objective.  */

void            FindEmptyCols(LP, Record, nonzerosCol, cNew, unbounded,
			                      top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  double         *cNew;
  int            *nonzerosCol, *unbounded, *top, Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, inx, nonzeros;
  double          xfix;
  ZeroColumn     *pZeroColumn;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;

  for (i = 0; i < cols; i++)
    if (Record->ColumnMask[i + 1] == STILL_ACTIVE) {
      if (nonzerosCol[i] == 0)
	Record->ColumnMask[i + 1] = Pass;

      /* make a Column deletion record, push it on the stack */

      if (Record->ColumnMask[i + 1] == Pass) {
#ifdef PREPROCESS_VERBOSE
	printf("Empty Column %d\n", i + 1);
#endif
	(*ReducedOnThisPass)++;
	Record->StackOfChanges[*top] = (SingleChange *)
	  Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	Record->StackOfChanges[*top]->ChangeType = ZERO_COLUMN;
	pZeroColumn = (ZeroColumn *)
	  Malloc(sizeof(ZeroColumn), "pZeroColumn");
	pZeroColumn->ColumnIndex = i + 1;

	/* fix the value of x at its appropriate bound, or detect
	 * unboundedness */
	if (cNew[i] < 0) {
	  if (Record->VarType[i] == NORMAL ||
	      Record->VarType[i] == FREE) {
	    printf("Problem is unbounded; column %d cost is %f\n",
		   i + 1, cNew[i]);
	    (*unbounded)++;
	    return;
	  } else if (Record->VarType[i] == UPPER) {
	    xfix = LP->UpBound[i];
	  } else
	    printf(" Unknown VarType %d %d\n", i + 1, Record->VarType[i]);

	} else if (cNew[i] > 0) {
	  if (Record->VarType[i] == FREE) {
	    printf("Problem is unbounded; column %d cost is %f\n",
		   i + 1, cNew[i]);
	    (*unbounded)++;
	    return;
	  } else if (Record->VarType[i] == NORMAL ||
		     Record->VarType[i] == UPPER)
	    xfix = 0.0;
	  else
	    printf(" Unknown VarType %d %d\n", i + 1, Record->VarType[i]);
	} else if (cNew[i] == 0) {
	  /* set it to zero; anything feasible would do */
	  xfix = 0.0;
	}
	pZeroColumn->Value = xfix;
	pZeroColumn->cshiftChange = xfix * cNew[i];
	pZeroColumn->cElement = cNew[i];
	Record->StackOfChanges[*top]->pZeroColumn = pZeroColumn;

	/* keep track of the change to the cshift constant term */
	Record->cshift += pZeroColumn->cshiftChange;

	/* increment top-of-stack pointer, resizing the stack if necessary */
	(*top)++;
	Record->Top = *top;
	if (Record->Top >= Record->Size)
	  ResizeRecord(Record);
      }
    }
  return;
}


/* Look for FIXED VARIABLES */
void            FindFixedVariables(LP, Record, nonzerosRow, nonzerosCol,
		                   cNew, bNew, top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  double         *cNew, *bNew;
  int            *nonzerosRow, *nonzerosCol, *top, Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, inx, nonzeros;
  double          xfix;
  FixedVariable  *pFixedVariable;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;

  for (i = 0; i < cols; i++)
    if (Record->ColumnMask[i + 1] == STILL_ACTIVE) {
      if (Record->VarType[i] == UPPER) {
	if (LP->UpBound[i] == 0.0) {
	  xfix = 0.0;
	  Record->ColumnMask[i + 1] = Pass;
	}
      }
      /* if this column was flagged for deletion, create a "FixedVariable"
       * record of the process and push it on the stack */

      if (Record->ColumnMask[i + 1] == Pass) {
	(*ReducedOnThisPass)++;
	Record->StackOfChanges[*top] = (SingleChange *)
	  Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	Record->StackOfChanges[*top]->ChangeType = FIXED_VARIABLE;

	pFixedVariable = (FixedVariable *)
	  Malloc(sizeof(FixedVariable), "pFixedVariable");
	pFixedVariable->ColumnIndex = i + 1;
	pFixedVariable->Value = xfix;
	pFixedVariable->cshiftChange = xfix * cNew[i];
	pFixedVariable->cElement = cNew[i];
#ifdef PREPROCESS_VERBOSE
	printf("Fixed Variable %d cshift %10.2e\n",
	       i + 1, pFixedVariable->cshiftChange);
#endif

	/* store the deleted column. Actually, need only store the elements
	 * that are STILL_ACTIVE at this point, but for simplicity store them
	 * all. (The ones that we don't want can be masked out during the
	 * postsolve.) */

	nonzeros = LP->A.pEndRow[i] - LP->A.pBeginRow[i] + 1;
	pFixedVariable->DeletedColumn = NewSparseVector(nonzeros);
	for (j = 0; j < nonzeros; j++) {
	  pFixedVariable->DeletedColumn->Index[j] =
	    LP->A.Row[LP->A.pBeginRow[i] + j - 1];
	  pFixedVariable->DeletedColumn->Value[j] =
	    LP->A.Value[LP->A.pBeginRow[i] + j - 1];
	}
	Record->StackOfChanges[*top]->pFixedVariable = pFixedVariable;

	/* Also, modify the rhs vector */
	for (j = 0; j < nonzeros; j++) {
	  inx = pFixedVariable->DeletedColumn->Index[j];
	  if (Record->RowMask[inx] == STILL_ACTIVE) {
	    bNew[inx - 1] -= xfix * pFixedVariable->DeletedColumn->Value[j];
	    nonzerosRow[inx - 1] -= 1;
	  }
	}
	/* and keep track of the change to the cshift constant term */
	Record->cshift += pFixedVariable->cshiftChange;

	/* increment top-of-stack pointer, resizing the stack if necessary */
	(*top)++;
	Record->Top = *top;
	if (Record->Top >= Record->Size)
	  ResizeRecord(Record);
      }
    }
  return;
}



/* Look for SINGLETON ROWS. */

void            FindSingletonRows(LP, Record, nonzerosRow, nonzerosCol,
	                   cNew, bNew, preprocessing_error, infeasible, top,
			   Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  double         *cNew, *bNew;
  int            *nonzerosRow, *nonzerosCol, *preprocessing_error, *infeasible, *top,
                  Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, inx, inxt, inxrc, nonzeros, singleton;
  double          xfix;
  SingletonRow   *pSingletonRow;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;
  (*preprocessing_error) = 0;
  (*infeasible) = 0;

  for (i = 0; i < rows; i++)
    if (Record->RowMask[i + 1] == STILL_ACTIVE) {

      /* is there just one (active) element in this row? */
      singleton = 0;
      for (j = LP->Atranspose.pBeginRow[i]; j <= LP->Atranspose.pEndRow[i]; j++) {
	inxt = LP->Atranspose.Row[j - 1];
	if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	  singleton++;
	  inxrc = inxt;
	  inx = j;
	}
      }

      if (singleton == 1) {

	/* Is the single element zero? (this case should have been trapped
	 * above by the zero-row loop) */
	if (LP->Atranspose.Value[inx - 1] == 0.0) {
	  printf("(Found a zero row in the singleton rows section %d)\n", i + 1);
	  nonzerosRow[i] = 0;
	} else {

	  /* flag row and column for deletion */
	  (*ReducedOnThisPass)++;
	  Record->ColumnMask[inxrc] = Pass;
	  Record->RowMask[i + 1] = Pass;

	  /* calculate the fixed value and check that it's within bounds */
	  xfix = bNew[i] / LP->Atranspose.Value[inx - 1];

/* SJW FIXED A HUGE BUG HERE 11/3/03 - uncovered by Alper Yildirim's student */
/* the row index "i" was being used instead of the column index inxrc
   in the lines below  - disgraceful!! */

	  (*infeasible) = 0;
	  if (Record->VarType[inxrc-1] == NORMAL) {
	    if (xfix < -TINY) {
	      (*infeasible) += 1;
	      printf(" ERROR: detected inconsistency in RowSingleton\n row %d column %d xfix %e (component should be nonnegative)\n",
		   i + 1, inxrc, xfix);
	    }
	  } else if (Record->VarType[inxrc-1] == UPPER) {
	    if (xfix < -TINY || xfix > LP->UpBound[inxrc-1]) {
	      (*infeasible)++;
	      printf(" ERROR: detected inconsistency in RowSingleton\n row %d column %d xfix %e (component should be in [0, %e]\n",
		     i + 1, inxrc, xfix, LP->UpBound[inxrc-1]);
	    }
	      
	  } 

	  /* if infeasible, return a null pointer */
	  if (*infeasible) 
	    return;

	  Record->StackOfChanges[*top] = (SingleChange *)
	    Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	  Record->StackOfChanges[*top]->ChangeType = SINGLETON_ROW;

	  pSingletonRow = (SingletonRow *)
	    Malloc(sizeof(SingletonRow), "pSingletonRow");
	  pSingletonRow->ColumnIndex = inxrc;
	  pSingletonRow->RowIndex = i + 1;
	  pSingletonRow->Value = xfix;
	  pSingletonRow->cshiftChange = xfix * cNew[inxrc - 1];
	  pSingletonRow->cElement = cNew[inxrc - 1];
	  pSingletonRow->bElement = bNew[i];
	  pSingletonRow->aElement = LP->Atranspose.Value[inx - 1];
#ifdef PREPROCESS_VERBOSE
	  printf("Row Singleton %d Col %d x %e cshift %e\n",
		 i + 1, inxrc, xfix, pSingletonRow->cshiftChange);
#endif

	  /* store the deleted column. Actually, need only store the elements
	   * that are STILL_ACTIVE at this point, but for simplicity store
	   * them all. (The ones that we don't want can be masked out during
	   * the postsolve.) */

	  nonzeros = LP->A.pEndRow[inxrc - 1] - LP->A.pBeginRow[inxrc - 1] + 1;
	  pSingletonRow->DeletedColumn = NewSparseVector(nonzeros);
	  for (j = 0; j < nonzeros; j++) {
	    inxt = LP->A.pBeginRow[inxrc - 1] + j - 1;
	    pSingletonRow->DeletedColumn->Index[j] =
	      LP->A.Row[inxt];
	    pSingletonRow->DeletedColumn->Value[j] =
	      LP->A.Value[inxt];
	  }
	  Record->StackOfChanges[*top]->pSingletonRow = pSingletonRow;

	  /* Also, modify the rhs vector */
	  for (j = 0; j < nonzeros; j++) {
	    inxt = pSingletonRow->DeletedColumn->Index[j];
	    if (Record->RowMask[inxt] == STILL_ACTIVE) {
	      bNew[inxt - 1] -= xfix * pSingletonRow->DeletedColumn->Value[j];
	      nonzerosRow[inxt - 1] -= 1;
	    }
	  }

	  /* and keep track of the change to the cshift constant term */
	  Record->cshift += pSingletonRow->cshiftChange;

	  /* increment top-of-stack pointer, resizing the stack if necessary */
	  (*top)++;
	  Record->Top = *top;
	  if (Record->Top >= Record->Size)
	    ResizeRecord(Record);
	}
      }
    }
  return;
}

/* look for FREE COLUMN SINGLETONS and IMPLIED FREE COLUMN SINGLETONS. We get
 * to express a variable in terms of other variables, so delete a row and
 * column and transform the objective funtion */

void            FindColumnSingletons(LP, Record, nonzerosRow, nonzerosCol,
				     cNew, bNew, preprocessing_error, top, 
				     Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  double         *cNew, *bNew;
  int            *nonzerosRow, *nonzerosCol, *preprocessing_error, *top,
                  Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, inx, inxt, inxrc, nonzeros, singleton;
  int             implied_bound_inf, still_alive;
  double          implied_bound;
  double          xfix;
  FreeColumnSingleton *pFreeColumnSingleton;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;
  (*preprocessing_error) = FALSE;

  for (i = 0; i < cols; i++)
    if (Record->ColumnMask[i + 1] == STILL_ACTIVE) {

      /* is there just one (active) element in this column? */
      singleton = 0;
      for (j = LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++) {
	inxt = LP->A.Row[j - 1];
	if (Record->RowMask[inxt] == STILL_ACTIVE) {
	  singleton++;
	  inxrc = inxt;
	  inx = j;
	}
      }

      /* Is this column a singleton, and is it free? */
      if (singleton == 1) {

	/* Is the single element zero? if so, set nonzerosCol to zero to trap
	 * it on next pass */
	if (LP->A.Value[inx - 1] == 0.0) {
	  printf("(Found a zero column in the column singleton section %d)\n", i + 1);
	  nonzerosCol[i] = 0;
	} else {

	  /* Is the singleton column free? If so flag it for deletion */

	  if (Record->VarType[i] == FREE) {
	    /* flag row and column for deletion */
	    (*ReducedOnThisPass)++;
	    Record->ColumnMask[i + 1] = Pass;
	    Record->RowMask[inxrc] = Pass;
#ifdef PREPROCESS_VERBOSE
	    printf("Free Column Singleton %d Row %d ", i + 1, inxrc);
#endif
	  }
	  /* OK, it's not free, but now check if it's *implied* free */
	  /* If it has a lower bound, test to see if this bound is redundant. */
	  else {
	    still_alive = TRUE;
	    if (LP->A.Value[inx - 1] > 0.0) {
	      implied_bound_inf = FALSE;
	      implied_bound = bNew[inxrc - 1];
	      for (j = LP->Atranspose.pBeginRow[inxrc - 1] - 1;
		   j <= LP->Atranspose.pEndRow[inxrc - 1] - 1; j++) {
		inxt = LP->Atranspose.Row[j];
		if (inxt != i + 1 && Record->ColumnMask[inxt] == STILL_ACTIVE) {
		  /* printf("subtract effect of col %d A elt is %f\n",
		   * inxt,LP->Atranspose.Value[j]); */
		  if (LP->Atranspose.Value[j] > 0.0) {
		    if (Record->VarType[inxt - 1] != UPPER)
		      implied_bound_inf = TRUE;
		    else
		      implied_bound -= LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
		  } else if (LP->Atranspose.Value[j] < 0.0) {
		    if (Record->VarType[inxt - 1] == FREE)
		      implied_bound_inf = TRUE;
		    /* (for non-free variables, the lower bound is zero, so
		     * we subtract nothing from the running tally) */
		  }
		}
	      }
	      if (implied_bound_inf)
		still_alive = FALSE;
	      else {
		implied_bound /= LP->A.Value[inx - 1];
		if (implied_bound < 0.0)
		  still_alive = FALSE;
	      }

	      /* same case, for negative coefficient term */
	    } else if (LP->A.Value[inx - 1] < 0.0) {
	      implied_bound_inf = FALSE;
	      implied_bound = bNew[inxrc - 1];
	      for (j = LP->Atranspose.pBeginRow[inxrc - 1] - 1;
		   j <= LP->Atranspose.pEndRow[inxrc - 1] - 1; j++) {
		inxt = LP->Atranspose.Row[j];
		if (inxt != i + 1 && Record->ColumnMask[inxt] == STILL_ACTIVE) {
		  if (LP->Atranspose.Value[j] < 0.0) {
		    if (Record->VarType[inxt - 1] != UPPER)
		      implied_bound_inf = TRUE;
		    else
		      implied_bound -= LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
		  } else if (LP->Atranspose.Value[j] > 0.0) {
		    if (Record->VarType[inxt - 1] == FREE)
		      implied_bound_inf = TRUE;
		    /* (for non-free variables, the lower bound is zero, so
		     * we subtract nothing from the running tally) */
		  }
		}
	      }
	      if (implied_bound_inf)
		still_alive = FALSE;
	      else {
		implied_bound /= LP->A.Value[inx - 1];
		if (implied_bound < 0.0)
		  still_alive = FALSE;
	      }
	    }
	    /* if we're still alive at this point, and if the variable has an
	     * upper bound, check to see if it's redundant */
	    if (Record->VarType[i] == UPPER && still_alive
		&& LP->A.Value[inx - 1] > 0.0) {
	      implied_bound_inf = FALSE;
	      implied_bound = bNew[inxrc - 1];
	      for (j = LP->Atranspose.pBeginRow[inxrc - 1] - 1;
		   j <= LP->Atranspose.pEndRow[inxrc - 1] - 1; j++) {
		inxt = LP->Atranspose.Row[j];
		if (inxt != i + 1 && Record->ColumnMask[inxt] == STILL_ACTIVE) {
		  if (LP->Atranspose.Value[j] < 0.0) {
		    if (Record->VarType[inxt - 1] != UPPER)
		      implied_bound_inf = TRUE;
		    else
		      implied_bound -= LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
		  } else if (LP->Atranspose.Value[j] > 0.0) {
		    if (Record->VarType[inxt - 1] == FREE)
		      implied_bound_inf = TRUE;
		    /* (for non-free variables, the lower bound is zero, so
		     * we subtract nothing from the running tally) */
		  }
		}
	      }
	      if (implied_bound_inf)
		still_alive = FALSE;
	      else {
		implied_bound /= LP->A.Value[inx - 1];
		if (implied_bound > LP->UpBound[i])
		  still_alive = FALSE;
	      }

	      /* same case, for negative coefficient term */
	    } else if (Record->VarType[i] == UPPER && still_alive
		       && LP->A.Value[inx - 1] < 0.0) {
	      implied_bound_inf = FALSE;
	      implied_bound = bNew[inxrc - 1];
	      for (j = LP->Atranspose.pBeginRow[inxrc - 1] - 1;
		   j <= LP->Atranspose.pEndRow[inxrc - 1] - 1; j++) {
		inxt = LP->Atranspose.Row[j];
		if (inxt != i + 1 && Record->ColumnMask[inxt] == STILL_ACTIVE) {
		  if (LP->Atranspose.Value[j] > 0.0) {
		    if (Record->VarType[inxt - 1] != UPPER)
		      implied_bound_inf = TRUE;
		    else
		      implied_bound -= LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
		  } else if (LP->Atranspose.Value[j] < 0.0) {
		    if (Record->VarType[inxt - 1] == FREE)
		      implied_bound_inf = TRUE;
		    /* (for non-free variables, the lower bound is zero, so
		     * we subtract nothing from the running tally) */
		  }
		}
	      }
	      if (implied_bound_inf)
		still_alive = FALSE;
	      else {
		implied_bound /= LP->A.Value[inx - 1];
		if (implied_bound > LP->UpBound[i])
		  still_alive = FALSE;
	      }
	    }
	    /* If we're still alive at this point, the column is implied
	     * free, so flag it for deletion. */

	    if (still_alive) {
	      (*ReducedOnThisPass)++;
	      Record->ColumnMask[i + 1] = Pass;
	      nonzerosCol[i] = 0;
	      Record->RowMask[inxrc] = Pass;
	      nonzerosRow[inxrc - 1] = 0;
#ifdef PREPROCESS_VERBOSE
	      printf("Implied Free Column Singleton %d Row %d ", i + 1, inxrc);
#endif
	    }
	  }

	  /* If this column is flagged for deletion, create a
	   * "FreeColumnSingleton" record of the deletion and push it on the
	   * stack */

	  if (Record->ColumnMask[i + 1] == Pass) {

	    Record->StackOfChanges[*top] = (SingleChange *)
	      Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	    Record->StackOfChanges[*top]->ChangeType = FREE_COLUMN_SINGLETON;

	    pFreeColumnSingleton = (FreeColumnSingleton *)
	      Malloc(sizeof(FreeColumnSingleton), "pFreeColumnSingleton");

	    pFreeColumnSingleton->ColumnIndex = i + 1;
	    pFreeColumnSingleton->RowIndex = inxrc;
	    pFreeColumnSingleton->cshiftChange =
	      bNew[inxrc - 1] * cNew[i] / LP->A.Value[inx - 1];
	    pFreeColumnSingleton->cElement = cNew[i];
	    pFreeColumnSingleton->bElement = bNew[inxrc - 1];
	    pFreeColumnSingleton->aElement = LP->A.Value[inx - 1];
#ifdef PREPROCESS_VERBOSE
	    printf("cshift %e\n", pFreeColumnSingleton->cshiftChange);
#endif

	    /* store the deleted row. Actually, need only store the elements
	     * that are STILL_ACTIVE at this point, but for simplicity store
	     * them all. (The ones that we don't want can be masked out
	     * during the postsolve.)  */

	    nonzeros = LP->Atranspose.pEndRow[inxrc - 1] -
	      LP->Atranspose.pBeginRow[inxrc - 1] + 1;
	    pFreeColumnSingleton->DeletedRow = NewSparseVector(nonzeros);
	    for (j = 0; j < nonzeros; j++) {
	      inxt = LP->Atranspose.pBeginRow[inxrc - 1] + j - 1;
	      pFreeColumnSingleton->DeletedRow->Index[j] =
		LP->Atranspose.Row[inxt];
	      pFreeColumnSingleton->DeletedRow->Value[j] =
		LP->Atranspose.Value[inxt];
	    }
	    Record->StackOfChanges[*top]->pFreeColumnSingleton =
	      pFreeColumnSingleton;

	    /* Modify the cost vector */
	    for (j = 0; j < nonzeros; j++) {
	      inxt = pFreeColumnSingleton->DeletedRow->Index[j];
	      if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
		cNew[inxt - 1] -= cNew[i] *
		  pFreeColumnSingleton->DeletedRow->Value[j] /
		  pFreeColumnSingleton->aElement;
		nonzerosCol[inxt - 1] -= 1;
	      }
	    }
	    Record->cshift += pFreeColumnSingleton->cshiftChange;

	    /* increment top-of-stack pointer, resizing the stack if
	     * necessary */
	    (*top)++;
	    Record->Top = (*top);
	    if (Record->Top >= Record->Size)
	      ResizeRecord(Record);
	  }
	}
      }
    }
  return;
}


/* look for FORCED ROWS. These are constraints which imply that all variables
 * involved must be at one of their bounds. e.g. x1 + x2 - x3 =2,   x1<=1,
 * x2<=1, x3>=0, forces each variable in the row to be at its bound */

void            FindForcedRows(LP, Record, nonzerosRow, nonzerosCol,
		                   cNew, bNew, top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  double         *cNew, *bNew;
  int            *nonzerosRow, *nonzerosCol, *top, Pass, *ReducedOnThisPass;
{
  int             rows, cols, i, j, k, inx, inxt, inxk, inxrc, nonzeros,
                  nonzerosc;
  int             bLower_inf, bUpper_inf;
  double          bLower, bUpper;
  ForcedRow      *pForcedRow;
  int             ResizeRecord();

  cols = LP->Cols;
  rows = LP->Rows;

  for (i = 0; i < rows; i++)
    if (Record->RowMask[i + 1] == STILL_ACTIVE) {

      /* find upper and lower bounds for the row sum, or flag if infinite. */

      bLower = 0.0;
      bUpper = 0.0;
      bLower_inf = FALSE;
      bUpper_inf = FALSE;
      for (j = LP->Atranspose.pBeginRow[i] - 1, nonzeros = 0;
	   j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	inxt = LP->Atranspose.Row[j];
	if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	  /* calculate the effect of this entry on the upper and lower bounds */
	  nonzeros++;
	  if (LP->Atranspose.Value[j] > 0.0) {
	    if (Record->VarType[inxt - 1] != UPPER)
	      bUpper_inf = TRUE;
	    else
	      bUpper += LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
	    if (Record->VarType[inxt - 1] == FREE)
	      bLower_inf = TRUE;
	  } else if (LP->Atranspose.Value[j] < 0.0) {
	    if (Record->VarType[inxt - 1] != UPPER)
	      bLower_inf = TRUE;
	    else
	      bLower += LP->Atranspose.Value[j] * LP->UpBound[inxt - 1];
	    if (Record->VarType[inxt - 1] == FREE)
	      bUpper_inf = TRUE;
	  }
	}
      }

      /* if there are zero or one active nonzeros in this row, ignore it - it
       * will be dealt with as an empty row or a row column singleton on the
       * next pass */

      /* if the upper or lower bound matches the RHS, we have a forced row */
      if (nonzeros > 1 && !(bUpper_inf) && bUpper == bNew[i]) {
#ifdef PREPROCESS_VERBOSE
	printf("Forced Row %d upper, %d variables ", i + 1, nonzeros);
#endif

	/* row forced at the UPPER limit. Make a change record and push it on
	 * the stack. */

	(*ReducedOnThisPass)++;
	Record->RowMask[i + 1] = Pass;
	nonzerosRow[i] = 0;

	Record->StackOfChanges[*top] = (SingleChange *)
	  Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	Record->StackOfChanges[*top]->ChangeType = FORCED_ROW;

	/* basic assignment of the ForcedRow data structure */
	pForcedRow = (ForcedRow *) Malloc(sizeof(ForcedRow), "pForcedRow");
	pForcedRow->RowIndex = i + 1;
	pForcedRow->cshiftChange = 0.0;
	pForcedRow->bElement = bNew[i];
	pForcedRow->LowerUpper = RANGEUPPER;

	/* copy the deleted row and corresponding cost elements into the
	 * ForcedRow data structure */

	pForcedRow->DeletedRow = NewSparseVector(nonzeros);
	pForcedRow->cElements = NewSparseVector(nonzeros);
	for (j = LP->Atranspose.pBeginRow[i] - 1, k = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    pForcedRow->DeletedRow->Index[k] = inxt;
	    pForcedRow->DeletedRow->Value[k] = LP->Atranspose.Value[j];
	    pForcedRow->cElements->Index[k] = inxt;
	    pForcedRow->cElements->Value[k] = cNew[inxt - 1];
	    k++;
	  }
	}

	/* assign values to all the eliminated variables, keeping track of
	 * how these assignments affect the cost (cshift) */

	pForcedRow->Values = NewSparseVector(nonzeros);
	for (j = LP->Atranspose.pBeginRow[i] - 1, k = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    pForcedRow->Values->Index[k] = inxt;
	    if (LP->Atranspose.Value[j] > 0.0)
	      pForcedRow->Values->Value[k] = LP->UpBound[inxt - 1];
	    else
	      pForcedRow->Values->Value[k] = 0.0;
	    pForcedRow->cshiftChange +=
	      cNew[inxt - 1] * pForcedRow->Values->Value[k];
	    k++;
	  }
	}
	Record->cshift += pForcedRow->cshiftChange;
#ifdef PREPROCESS_VERBOSE
	printf("cshift %e\n", pForcedRow->cshiftChange);
#endif

	for (k = 0; k < nonzeros; k++)
	  inxt = pForcedRow->DeletedRow->Index[k];

	/* add all the deleted columns to the record (yuk!). Modify the
	 * right-hand side while we're at it */

	pForcedRow->DeletedColumns = (SparseVector **)
	  Malloc(nonzeros * sizeof(SparseVector *), "DeletedRow");
	for (j = LP->Atranspose.pBeginRow[i] - 1, inx = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    Record->ColumnMask[inxt] = Pass;
	    nonzerosc = LP->A.pEndRow[inxt - 1] - LP->A.pBeginRow[inxt - 1] + 1;
	    pForcedRow->DeletedColumns[inx] = NewSparseVector(nonzerosc);
	    for (k = LP->A.pBeginRow[inxt - 1] - 1, inxk = 0;
		 k <= LP->A.pEndRow[inxt - 1] - 1; k++) {
	      pForcedRow->DeletedColumns[inx]->Index[inxk] = LP->A.Row[k];
	      pForcedRow->DeletedColumns[inx]->Value[inxk] = LP->A.Value[k];
	      if (Record->RowMask[LP->A.Row[k]] == STILL_ACTIVE) {
		bNew[LP->A.Row[k] - 1] -=
		  LP->A.Value[k] * pForcedRow->Values->Value[inx];
		nonzerosRow[LP->A.Row[k] - 1] -= 1;
	      }
	      inxk++;
	    }
	    inx++;
	  }
	}

	/* finally, eliminate the forced variables */
	/* for(j=LP->Atranspose.pBeginRow[i]-1, inx=0;
	 * j<=LP->Atranspose.pEndRow[i]-1; j++) { inxt=LP->Atranspose.Row[j];
	 * if(Record->ColumnMask[inxt] == STILL_ACTIVE) {
	 * Record->ColumnMask[inxt]=Pass; } } */

	/* put the pointer to the ForcedRow data structure into the record */

	Record->StackOfChanges[*top]->pForcedRow = pForcedRow;

	/* increment top-of-stack pointer, resizing the stack if necessary */
	(*top)++;
	Record->Top = *top;
	if (Record->Top >= Record->Size)
	  ResizeRecord(Record);

      } else if (nonzeros > 1 && !(bLower_inf) && bLower == bNew[i]) {
#ifdef PREPROCESS_VERBOSE
	printf("Forced Row %d lower, %d variables ", i + 1, nonzeros);
#endif

	/* row forced at the LOWER limit. Make a change record and push it on
	 * the stack. */

	(*ReducedOnThisPass)++;
	Record->RowMask[i + 1] = Pass;

	Record->StackOfChanges[*top] = (SingleChange *)
	  Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	Record->StackOfChanges[*top]->ChangeType = FORCED_ROW;

	/* basic assignment of the ForcedRow data structure */
	pForcedRow = (ForcedRow *) Malloc(sizeof(ForcedRow), "pForcedRow");
	pForcedRow->RowIndex = i + 1;
	pForcedRow->cshiftChange = 0.0;
	pForcedRow->bElement = bNew[i];
	pForcedRow->LowerUpper = RANGELOWER;

	/* copy the deleted row and corresponding cost elements into the
	 * ForcedRow data structure */

	pForcedRow->DeletedRow = NewSparseVector(nonzeros);
	pForcedRow->cElements = NewSparseVector(nonzeros);
	for (j = LP->Atranspose.pBeginRow[i] - 1, k = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    pForcedRow->DeletedRow->Index[k] = inxt;
	    pForcedRow->DeletedRow->Value[k] = LP->Atranspose.Value[j];
	    pForcedRow->cElements->Index[k] = inxt;
	    pForcedRow->cElements->Value[k] = cNew[inxt - 1];
	    k++;
	  }
	}

	/* assign values to all the eliminated variables, keeping track of
	 * how these assignments affect the cost (cshift) */

	pForcedRow->Values = NewSparseVector(nonzeros);
	for (j = LP->Atranspose.pBeginRow[i] - 1, k = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    pForcedRow->Values->Index[k] = LP->Atranspose.Row[j];
	    if (LP->Atranspose.Value[j] > 0.0)
	      pForcedRow->Values->Value[k] = 0.0;
	    else
	      pForcedRow->Values->Value[k] = LP->UpBound[inxt - 1];
	    pForcedRow->cshiftChange +=
	      cNew[inxt - 1] * pForcedRow->Values->Value[k];
	    k++;
	  }
	}
	Record->cshift += pForcedRow->cshiftChange;
#ifdef PREPROCESS_VERBOSE
	printf("cshift %e\n", pForcedRow->cshiftChange);
#endif

	/* add all the deleted columns to the record (yuk!). Modify the
	 * right-hand side while we're at it */

	pForcedRow->DeletedColumns = (SparseVector **)
	  Malloc(nonzeros * sizeof(SparseVector *), "DeletedRow");
	for (j = LP->Atranspose.pBeginRow[i] - 1, inx = 0;
	     j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	  inxt = LP->Atranspose.Row[j];
	  if (Record->ColumnMask[inxt] == STILL_ACTIVE) {
	    Record->ColumnMask[inxt] = Pass;
	    nonzerosc = LP->A.pEndRow[inxt - 1] - LP->A.pBeginRow[inxt - 1] + 1;
	    pForcedRow->DeletedColumns[inx] = NewSparseVector(nonzerosc);
	    for (k = LP->A.pBeginRow[inxt - 1] - 1, inxk = 0;
		 k <= LP->A.pEndRow[inxt - 1] - 1; k++) {
	      pForcedRow->DeletedColumns[inx]->Index[inxk] = LP->A.Row[k];
	      pForcedRow->DeletedColumns[inx]->Value[inxk] = LP->A.Value[k];
	      if (Record->RowMask[LP->A.Row[k]] == STILL_ACTIVE) {
		bNew[LP->A.Row[k] - 1] -=
		  LP->A.Value[k] * pForcedRow->Values->Value[inx];
		nonzerosRow[LP->A.Row[k] - 1] -= 1;
	      }
	      inxk++;
	    }
	    inx++;
	  }
	}

	/* finally, eliminate the forced variables */
	/* for(j=LP->Atranspose.pBeginRow[i]-1, inx=0;
	 * j<=LP->Atranspose.pEndRow[i]-1; j++) { inxt=LP->Atranspose.Row[j];
	 * if(Record->ColumnMask[inxt] == STILL_ACTIVE)
	 * Record->ColumnMask[inxt] = Pass; } */

	/* put the pointer to the ForcedRow data structure into the record */
	Record->StackOfChanges[*top]->pForcedRow = pForcedRow;

	/* increment top-of-stack pointer, resizing the stack if necessary */
	(*top)++;
	Record->Top = *top;
	if (Record->Top >= Record->Size)
	  ResizeRecord(Record);
      }
    }
  return;
}


/* Combine Duplicate Columns */

int             DuplicateCols(LP, Record, c, nonzerosRow, nonzerosCol,
			                      top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  int            *nonzerosRow, *nonzerosCol, *top, Pass, *ReducedOnThisPass;
  double         *c;
{
  int             i, j, jj, k, entry1, entry2, NumInCol, inx1, inx2, ptr1,
                  ptr2, Iisbigger, counter = 0, same = 0;
  DuplicateCol   *pDuplicateCol;
  int             ResizeRecord();
  double          ratioIJ, costdiff;
  int             deleteI = FALSE, deleteJ = FALSE;
  int             between(), SkipToNextActive(), ShortestRow, nShortestRow, RowIndex,
                  NumInRow;

  for (i = 0; i < LP->Cols; i++)
    /* Compare NORMAL column i with other columns j to see if these two
     * columns are multiples of each other. */

    if (Record->ColumnMask[i + 1] == STILL_ACTIVE &&
	nonzerosCol[i] > 1 &&
	Record->VarType[i] == NORMAL) {

      deleteI = FALSE;

      /* look at the nonzeros in this column and find the row index of the
       * shortest row */
      nShortestRow = LP->Cols;
      ShortestRow = 0;
      for (j = LP->A.pBeginRow[i] - 1; j <= LP->A.pEndRow[i] - 1; j++) {
	RowIndex = LP->A.Row[j];
	if (Record->RowMask[RowIndex] == STILL_ACTIVE) {
	  NumInRow = nonzerosRow[RowIndex - 1];
	  if (NumInRow < nShortestRow) {
	    nShortestRow = NumInRow;
	    ShortestRow = RowIndex - 1;
	  }
	}
      }

      /* if there is only one element in the shortest row, then this column
       * is structurally unique */
      if (nShortestRow <= 1)
	same = FALSE;
      else {

	/* loop over the columns represented in the shortest row */

	for (jj = LP->Atranspose.pBeginRow[ShortestRow] - 1;
	  jj <= LP->Atranspose.pEndRow[ShortestRow] - 1 && !deleteI; jj++) {
	  j = LP->Atranspose.Row[jj] - 1;
	  if (j >= LP->Cols)
	    printf("Whats this horseshit j=%d cols=%d\n", j, LP->Cols);
	  if (j > i &&
	      Record->ColumnMask[j + 1] == STILL_ACTIVE &&
	      nonzerosCol[j] == nonzerosCol[i] &&
	      Record->VarType[j] == NORMAL) {

	    deleteJ = FALSE;
	    same = TRUE;
	    NumInCol = nonzerosCol[i];

	    ptr1 = LP->A.pBeginRow[i] - 1;
	    ptr2 = LP->A.pBeginRow[j] - 1;
	    for (k = 0; k < NumInCol && same; k++) {
	      if (SkipToNextActive(LP->A, i, j, &ptr1, &ptr2, Record->RowMask) == 1)
		same = FALSE;
	      else if (LP->A.Row[ptr1] != LP->A.Row[ptr2])
		same = FALSE;
	      else {
		ptr1++;
		ptr2++;
	      }
	    }

	    /* now test to see if one column is a multiple of the other */
	    if (same) {
	      ptr1 = LP->A.pBeginRow[i] - 1;
	      ptr2 = LP->A.pBeginRow[j] - 1;
	      if (SkipToNextActive(LP->A, i, j, &ptr1, &ptr2, Record->RowMask) == 1)
		same = FALSE;
	    }
	    if (same) {

	      if (fabs(LP->A.Value[ptr1]) > fabs(LP->A.Value[ptr2])) {

		Iisbigger = TRUE;
		ratioIJ = LP->A.Value[ptr2] / LP->A.Value[ptr1];

		ptr1++;
		ptr2++;
		for (k = 1; k < NumInCol && same; k++) {
		  if (SkipToNextActive(LP->A, i, j, &ptr1, &ptr2,
				       Record->RowMask) == 1)
		    same = FALSE;
		  else if (LP->A.Value[ptr2] != ratioIJ * LP->A.Value[ptr1])
		    same = FALSE;
		  else {
		    ptr1++;
		    ptr2++;
		  }
		}

		/* calculate the costdiff */
		if (same) {
		  costdiff = c[j] - ratioIJ * c[i];
		  /* trap the infeasible case */
		  if (ratioIJ < 0.0 && costdiff < 0.0) {
		    printf(" ERROR detected in DuplicateCols\n");
		    printf(" Primal problem is unbounded below i=%d j=%d\n", i, j);
		    return UNBOUNDED;
		  } else if (ratioIJ < 0.0 && costdiff >= 0.0) {
		    deleteI = FALSE;
		    deleteJ = FALSE;
		    same = FALSE;
		  } else if (ratioIJ >= 0.0 && costdiff < 0.0) {
		    deleteI = TRUE;
		    deleteJ = FALSE;
		  } else if (ratioIJ >= 0.0 && costdiff >= 0.0) {
		    deleteI = FALSE;
		    deleteJ = TRUE;
		  }
		}
	      } else {

		Iisbigger = FALSE;
		ratioIJ = LP->A.Value[ptr1] / LP->A.Value[ptr2];

		ptr1++;
		ptr2++;
		for (k = 1; k < NumInCol && same; k++) {
		  if (SkipToNextActive(LP->A, i, j, &ptr1, &ptr2,
				       Record->RowMask) == 1)
		    same = FALSE;
		  else if (LP->A.Value[ptr1] != ratioIJ * LP->A.Value[ptr2])
		    same = FALSE;
		  else {
		    ptr1++;
		    ptr2++;
		  }
		}

		/* calculate the costdiff */
		if (same) {
		  costdiff = c[j] - c[i] / ratioIJ;
		  /* trap the infeasible case */
		  if (ratioIJ < 0.0 && costdiff < 0.0) {
		    printf(" ERROR detected in DuplicateCols\n");
		    printf(" Primal problem is unbounded below i=%d j=%d\n", i, j);
		    return UNBOUNDED;
		  } else if (ratioIJ < 0.0 && costdiff >= 0.0) {
		    deleteI = FALSE;
		    deleteJ = FALSE;
		    same = FALSE;
		  } else if (ratioIJ >= 0.0 && costdiff < 0.0) {
		    deleteI = TRUE;
		    deleteJ = FALSE;
		  } else if (ratioIJ >= 0.0 && costdiff >= 0.0) {
		    deleteI = FALSE;
		    deleteJ = TRUE;
		  }
		}
	      }
	    }
	    if (same) {
	      counter++;
#ifdef PREPROCESS_VERBOSE
	      if (costdiff != 0)
		printf(" Dup Cols (different costs) %d %d (%d elements) %f\n",
		       i, j, nonzerosCol[i], ratioIJ);
	      else
		printf(" Dup Cols %d %d (%d elements) %f\n",
		       i, j, nonzerosCol[i], ratioIJ);
#endif
	      (*ReducedOnThisPass)++;
	      Record->StackOfChanges[*top] = (SingleChange *)
		Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	      Record->StackOfChanges[*top]->ChangeType = DUPLICATE_COLUMN;
	      pDuplicateCol = (DuplicateCol *)
		Malloc(sizeof(DuplicateCol), "pDuplicateCol");

	      if (deleteJ) {
		pDuplicateCol->Split1 = i + 1;
		pDuplicateCol->Split2 = j + 1;
	      } else if (deleteI) {
		pDuplicateCol->Split1 = j + 1;
		pDuplicateCol->Split2 = i + 1;
	      } else {
		printf(" what are we doing here???\n");
	      }

	      /* VarType become FREE if ratioIJ is negative; stays NORMAL
	       * otherwise. */
	      if (ratioIJ < 0.0) {
		pDuplicateCol->Type = FREE;
		if (deleteJ)
		  Record->VarType[i] = FREE;
		else
		  Record->VarType[j] = FREE;
	      } else {
		pDuplicateCol->Type = NORMAL;
	      }

	      /* compute the multiplier for the deleted column */
	      if ((deleteI && Iisbigger) || (deleteJ && !Iisbigger)) {
		pDuplicateCol->multiplier = 1.0 / ratioIJ;
		pDuplicateCol->costdiff =
		  c[pDuplicateCol->Split2 - 1] -
		  pDuplicateCol->multiplier * c[pDuplicateCol->Split1 - 1];
	      } else {
		pDuplicateCol->multiplier = ratioIJ;
		pDuplicateCol->costdiff =
		  c[pDuplicateCol->Split2 - 1] -
		  pDuplicateCol->multiplier * c[pDuplicateCol->Split1 - 1];
	      }

	      Record->StackOfChanges[*top]->pDuplicateCol = pDuplicateCol;

	      /* increment top-of-stack pointer, resizing the stack if
	       * necessary */
	      (*top)++;
	      Record->Top = *top;
	      if (Record->Top >= Record->Size)
		ResizeRecord(Record);

	      for (k = LP->A.pBeginRow[j] - 1; k <= LP->A.pEndRow[j] - 1; k++) {
		entry1 = LP->A.Row[k];
		if (Record->RowMask[entry1] == STILL_ACTIVE)
		  nonzerosRow[entry1 - 1] -= 1;
	      }
	      if (deleteJ)
		Record->ColumnMask[j + 1] = Pass;
	      else if (deleteI)
		Record->ColumnMask[i + 1] = Pass;
	    }
	  }			/* end if j == NORMAL */
	}			/* end for jj=... */
      }				/* end if(nShortestRow <= 1) */
    }				/* end if i == NORMAL */
  return (counter);		/* return number of split variables */
}

/* Combine Duplicate Rows */

int             DuplicateRows(LP, Record, b, nonzerosRow, nonzerosCol,
			                      top, Pass, ReducedOnThisPass)
  LPtype         *LP;
  ChangeStack    *Record;
  int            *nonzerosRow, *nonzerosCol, *top, Pass, *ReducedOnThisPass;
  double         *b;
{
  int             i, j, jj, k, entry1, entry2, NumInRow, inx1, inx2, ptr1,
                  ptr2, inx0, counter = 0, same = 0;
  DuplicateRow   *pDuplicateRow;
  int             ResizeRecord();
  double          ratioIJ;
  int             deleteI = 0, deleteJ = 0;
  int             SkipToNextActive(), ColIndex, NumInCol, ShortestCol,
                  nShortestCol;

  for (i = 0; i < LP->Rows; i++)
    /* Compare NORMAL column i with other columns j to see if these two
     * columns are multiples of each other. */

    if (Record->RowMask[i + 1] == STILL_ACTIVE &&
	nonzerosRow[i] > 1) {

      deleteI = FALSE;

      /* look at the nonzeros in this row and find the col index of the
       * shortest col */
      nShortestCol = LP->Rows;
      ShortestCol = 0;
      for (j = LP->Atranspose.pBeginRow[i] - 1; j <= LP->Atranspose.pEndRow[i] - 1; j++) {
	ColIndex = LP->Atranspose.Row[j];
	if (Record->ColumnMask[ColIndex] == STILL_ACTIVE) {
	  NumInCol = nonzerosCol[ColIndex - 1];
	  if (NumInCol < nShortestCol) {
	    nShortestCol = NumInCol;
	    ShortestCol = ColIndex - 1;
	  }
	}
      }

      /* if there is only one element in the shortest col, then this row is
       * structurally unique */
      if (nShortestCol <= 1)
	same = FALSE;
      else {

	/* loop over the columns represented in the shortest row that come
	 * after this column */

	for (jj = LP->A.pBeginRow[ShortestCol] - 1;
	     jj <= LP->A.pEndRow[ShortestCol] - 1 && !deleteI; jj++) {
	  j = LP->A.Row[jj] - 1;
	  if (j > i && Record->RowMask[j + 1] == STILL_ACTIVE &&
	      nonzerosRow[j] == nonzerosRow[i]) {

	    deleteJ = FALSE;
	    same = TRUE;
	    NumInRow = nonzerosRow[i];

	    ptr1 = LP->Atranspose.pBeginRow[i] - 1;
	    ptr2 = LP->Atranspose.pBeginRow[j] - 1;

	    for (k = 0; k < NumInRow && same; k++) {
	      if (SkipToNextActive(LP->Atranspose, i, j, &ptr1, &ptr2,
				   Record->ColumnMask) == 1)
		same = FALSE;
	      else if (LP->Atranspose.Row[ptr1] != LP->Atranspose.Row[ptr2])
		same = FALSE;
	      else {
		ptr1++;
		ptr2++;
	      }
	    }

	    /* now test to see if one row is the multiple of the other */
	    if (same) {
	      ptr1 = LP->Atranspose.pBeginRow[i] - 1;
	      ptr2 = LP->Atranspose.pBeginRow[j] - 1;
	      if (SkipToNextActive(LP->Atranspose, i, j, &ptr1, &ptr2,
				   Record->ColumnMask) == 1)
		same = FALSE;
	    }
	    if (same) {

	      if (fabs(LP->Atranspose.Value[ptr1]) >
		  fabs(LP->Atranspose.Value[ptr2])) {

		ratioIJ = LP->Atranspose.Value[ptr2] /
		  LP->Atranspose.Value[ptr1];
		ptr1++;
		ptr2++;
		for (k = 1; k < NumInRow && same; k++) {
		  if (SkipToNextActive(LP->Atranspose, i, j, &ptr1, &ptr2,
				       Record->ColumnMask) == 1)
		    same = FALSE;
		  else if (LP->Atranspose.Value[ptr2] !=
			   ratioIJ * LP->Atranspose.Value[ptr1])
		    same = FALSE;
		  else {
		    ptr1++;
		    ptr2++;
		  }
		}

		if (same && b[j] != ratioIJ * b[i]) {
		  same = FALSE;
		  printf(" Infeasible duplicate rows %d %d. Ents: %d Ratio: %e\n",
			 i, j, NumInRow, ratioIJ);
		}
		if (same) {
		  deleteI = FALSE;
		  deleteJ = TRUE;
		}
	      } else {

		ratioIJ = LP->Atranspose.Value[ptr1] /
		  LP->Atranspose.Value[ptr2];
		ptr1++;
		ptr2++;
		for (k = 1; k < NumInRow && same; k++) {
		  if (SkipToNextActive(LP->Atranspose, i, j, &ptr1, &ptr2,
				       Record->ColumnMask) == 1)
		    same = FALSE;
		  else if (LP->Atranspose.Value[ptr1] !=
			   ratioIJ * LP->Atranspose.Value[ptr2])
		    same = FALSE;
		  else {
		    ptr1++;
		    ptr2++;
		  }
		}

		if (same && b[i] != ratioIJ * b[j]) {
		  same = FALSE;
		  printf(" Infeasible duplicate rows %d %d. Ents: %d Ratio: %e\n",
			 i, j, NumInRow, ratioIJ);
		}
		if (same) {
		  deleteI = FALSE;
		  deleteJ = TRUE;
		}
	      }
	    }
	    if (same) {
	      counter++;
#ifdef PREPROCESS_VERBOSE
	      printf(" Duplicate Rows %d %d %f\n", i, j, ratioIJ);
#endif
	      (*ReducedOnThisPass)++;
	      Record->StackOfChanges[*top] = (SingleChange *)
		Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
	      Record->StackOfChanges[*top]->ChangeType = DUPLICATE_ROW;
	      pDuplicateRow = (DuplicateRow *)
		Malloc(sizeof(DuplicateRow), "pDuplicateRow");

	      if (deleteJ) {
		pDuplicateRow->Split1 = i + 1;
		pDuplicateRow->Split2 = j + 1;
	      } else if (deleteI) {
		pDuplicateRow->Split1 = j + 1;
		pDuplicateRow->Split2 = i + 1;
	      } else {
		printf(" what are we doing here???\n");
	      }

	      Record->StackOfChanges[*top]->pDuplicateRow = pDuplicateRow;

	      /* increment top-of-stack pointer, resizing the stack if
	       * necessary */
	      (*top)++;
	      Record->Top = *top;
	      if (Record->Top >= Record->Size)
		ResizeRecord(Record);

	      for (k = LP->Atranspose.pBeginRow[j] - 1;
		   k <= LP->Atranspose.pEndRow[j] - 1; k++) {
		entry1 = LP->Atranspose.Row[k];
		if (Record->ColumnMask[entry1] == STILL_ACTIVE)
		  nonzerosCol[entry1 - 1] -= 1;
	      }
	      if (deleteJ)
		Record->RowMask[j + 1] = Pass;
	      else
		Record->RowMask[i + 1] = Pass;
	    }
	  }			/* end (if j == STILL_ACTIVE) */
	}			/* end(for jj=... */
      }				/* end if(nShortestCol <= 1) */
    }				/* end (if i == STILL_ACTIVE) */
  return (counter);		/* return number of split variables */
}

/* returns TRUE if x lies (nonstrictly) between a and b, FALSE otherwise */

int             between(x, a, b)
  double          x, a, b;
{
  double          tmp;

  tmp = (x - a) * (x - b);
  if (tmp <= 0.0)
    return (TRUE);
  else
    return (FALSE);
}

/* Skips to next active indices in a pair of rows of A. returns "1" if
 * there's a problem, "0" otherwise */

int             SkipToNextActive(A, i, j, ptr1, ptr2, Mask)
  sparseMatrix    A;
  int             i, j, *ptr1, *ptr2, *Mask;
{
  int             inx1, inx2, errflag = 1, ok = 0;

  if ((*ptr1) < A.pBeginRow[i] - 1) {
    printf("Error in SkipToNextActive i=%d j=%d\n", i, j);
    return (errflag);
  }
  if ((*ptr2) < A.pBeginRow[j] - 1) {
    printf("Error in SkipToNextActive i=%d j=%d\n", i, j);
    return (errflag);
  }
  inx1 = A.Row[*ptr1];
  for (; (*ptr1) < A.pEndRow[i] && Mask[inx1] != STILL_ACTIVE;) {
    (*ptr1)++;
    inx1 = A.Row[*ptr1];
  }
  if ((*ptr1) == A.pEndRow[i]) {
    printf("SkipToNextActive run over end of row %d\n", i);
    return (errflag);
  }
  inx2 = A.Row[*ptr2];
  for (; (*ptr2) < A.pEndRow[j] && Mask[inx2] != STILL_ACTIVE;) {
    (*ptr2)++;
    inx2 = A.Row[*ptr2];
  }
  if ((*ptr2) == A.pEndRow[j]) {
    printf("SkipToNextActive run over end of row %d\n", j);
    return (errflag);
  }
  return (ok);
}


/* Routines to manipulate the data structure "Record" that keeps
 * track of preprocessing changes.  SJW 2/13/95 */

/* Allocate a new record for a given LP */

ChangeStack    *NewRecord(LP)
  LPtype         *LP;
{
  int             rows, cols, sizeUnit;
  ChangeStack    *Record;

  /* How big do we make the Record of changes? */
  sizeUnit = STACK_SIZE_FACTOR * MAX(LP->Rows, LP->Cols);
  sizeUnit = MAX(sizeUnit, 100);

  rows = LP->Rows;
  cols = LP->Cols;

  /* Allocate space for the record */
  Record = (ChangeStack *) Malloc(sizeof(ChangeStack), "Record");

  /* initialize the Record with the stuff we know about */
  Record->Rows = rows;
  Record->Columns = cols;
  Record->Top = 0;
  Record->Passes = 0;
  Record->ReducedRows = rows;
  Record->ReducedColumns = cols;
  Record->cshift = 0.0;
  Record->Size = sizeUnit;
  Record->SizeUnit = sizeUnit;
  Record->Passes = 0;
  Record->ReducedRows = 0;
  Record->ReducedColumns = 0;
  
  /* NewInt initializes the mask vectors to zeros, meaning that the default
   * status for each row and column is "STILL_ACTIVE" */

  Record->RowMask = NewInt(rows + 1, "Record->RowMask");
  Record->ColumnMask = NewInt(cols + 1, "Record->ColumnMask");
  Record->VarType = NewInt(cols, "Record->VarType");

  /* Allocate space for the stack of changes */
  Record->StackOfChanges = (SingleChange **)
    Malloc(sizeUnit * sizeof(SingleChange *), "Record->StackOfChanges");
  return (Record);
}

void            DeleteChangeStack(Record)
  ChangeStack    *Record;
{
  int             i;

  for (i = 0; i < Record->Top; i++) {
    SingleChange *top = Record->StackOfChanges[i];
    switch(top->ChangeType) {
    case ZERO_ROW:
      Free((char *) top->pZeroRow);
      break;
    case ZERO_COLUMN:
      Free((char *) top->pZeroColumn);
      break;
    case FIXED_VARIABLE:
      FreeSparseVector(top->pFixedVariable->DeletedColumn);
      Free((char *) top->pFixedVariable);
      break;
    case SINGLETON_ROW:
      FreeSparseVector(top->pSingletonRow->DeletedColumn);
      Free((char *) top->pSingletonRow);
      break;
    case FREE_COLUMN_SINGLETON:
      FreeSparseVector(top->pFreeColumnSingleton->DeletedRow);
      Free((char *) top->pFreeColumnSingleton);
      break;
    case FORCED_ROW:
      {
        int i, nonzeros = top->pForcedRow->Values->size;
        SparseVector **DeletedColumns = top->pForcedRow->DeletedColumns;
        for (i = 0; i < nonzeros; ++i) {
          if (DeletedColumns[i])
            FreeSparseVector(DeletedColumns[i]);
        }
        Free((char *) DeletedColumns);
        FreeSparseVector(top->pForcedRow->DeletedRow);
        FreeSparseVector(top->pForcedRow->cElements);
        FreeSparseVector(top->pForcedRow->Values);
        Free((char *) top->pForcedRow);
        break;
      }
    case DUPLICATE_ROW:
      Free((char *) top->pDuplicateRow);
      break;
    case DUPLICATE_COLUMN:
      Free((char *) top->pDuplicateCol);
      break;
    }
    Free((char *) top);
  }
  Free((char *) Record->StackOfChanges);
  Free((char *) Record->RowMask);
  Free((char *) Record->ColumnMask);
  Free((char *) Record->VarType);
  Free((char *) Record);
  return;
}

/* Add space to the Record if we run out */
int             ResizeRecord(Record)
  ChangeStack    *Record;
{
  int             size;

  size = Record->Size + Record->SizeUnit;
  Record->StackOfChanges =
    (SingleChange **) Realloc(Record->StackOfChanges,
		   size * sizeof(SingleChange *), "Record->StackOfChanges");
  Record->Size = size;
  return 0;
}


/* Postprocessing routines SJW  12/95 */


/* Postprocess reads through the record of changes made by Preprocess, in
 * reverse order, and "undoes" them.
 * 
 * On entry: LP points to a linear program, in the PCx data structure; pRecord
 * points to a record of all the changes; pSolution points to a solution of
 * the reduced problem.
 * 
 * On return: pSolution points to a solution of the original problem. The 
 * solution of the reduced problem is deleted. The record  pointed to by
 * pRecord is also deleted and the space is freed.  SJW 6/24/97.
 * 
 * SJW 12/95 */

int             Postprocess(LP, pRecord, pSolution)
  LPtype         *LP;
  ChangeStack   **pRecord;
  solution      **pSolution;
{
  ChangeStack    *Record;
  solution       *Solution;
  solution       *FullSolution;
  int             i, j, k, l, inx, inxr, inxc, inxt, rows, columns, stackTop,
                  col, piUpper_inf, piLower_inf, err;
  double          temp, piUpper, piLower, TOLERANCE = 1.e-8;

  /* some pointers that make the code a little easier to read */
  DuplicateCol   *pDuplicateCol;
  DuplicateRow   *pDuplicateRow;
  ZeroRow        *pZeroRow;
  ZeroColumn     *pZeroColumn;
  FixedVariable  *pFixedVariable;
  SingletonRow   *pSingletonRow;
  FreeColumnSingleton *pFreeColumnSingleton;
  ForcedRow      *pForcedRow;

  /* Set up "Solution" as a pointer to the solution of the reduced problem.
   * FullSolution will gradually be massaged into the solution of the
   * original problem */

  Record = *pRecord;
  Solution = *pSolution;

  /* if no preprocessing was done, return immediately without changing a thing */
  if(Record->Top <= 0) {
    Solution->PriorRows = Solution->Rows;
    Solution->PriorColumns = Solution->Columns;
    Solution->ReducedRows = Solution->Rows;
    Solution->ReducedColumns = Solution->Columns;
    Solution->Passes = 0;
    DeleteChangeStack(Record);
    return 0;
  }

  rows = Record->Rows;
  columns = Record->Columns;
  FullSolution = (solution *) Malloc(sizeof(solution), "FullSolution");
  FullSolution->Rows = rows;
  FullSolution->Columns = columns;
  FullSolution->PriorRows = Record->Rows;
  FullSolution->PriorColumns = Record->Columns;
  FullSolution->ReducedRows = Record->ReducedRows;
  FullSolution->ReducedColumns = Record->ReducedColumns;
  FullSolution->Passes =  Record->Passes;
  FullSolution->x = NewDouble(columns, "FullSolution->x");
  FullSolution->pi = NewDouble(rows, "FullSolution->pi");
  FullSolution->DualUpper = NewDouble(columns, "FullSolution->DualUpper");
  FullSolution->DualLower = NewDouble(columns, "FullSolution->DualLower");

  FullSolution->Activity  = Solution->Activity;

  FullSolution->ReadTime = Solution->ReadTime;
  FullSolution->PreprocessTime = Solution->PreprocessTime;
  FullSolution->SolutionTime = Solution->SolutionTime;
  FullSolution->FactorizationTime = Solution->FactorizationTime;
  FullSolution->SolveADATTime = Solution->SolveADATTime;
  FullSolution->InitTime = Solution->InitTime;
  FullSolution->PredictorTime = Solution->PredictorTime;
  FullSolution->LoopTime = Solution->LoopTime;
  FullSolution->CorrectorTime = Solution->CorrectorTime;
  FullSolution->FormADATtime = Solution->FormADATtime;

  FullSolution->FactorizationCode = (char *) Malloc((strlen(Solution->FactorizationCode)+1) * sizeof(char), "FullSolution->FactorizationCode");
  strcpy(FullSolution->FactorizationCode, Solution->FactorizationCode);

  FullSolution->Iterations = Solution->Iterations;
  FullSolution->RestoredIteration = Solution->RestoredIteration;

  FullSolution->Status = Solution->Status;
  FullSolution->IterationHistory = (IterationRecord *)
    Malloc((Solution->Iterations + 1) * sizeof(IterationRecord),
	   "Solution->IterationHistory");

  for (i = 0; i <= FullSolution->Iterations; i++) {
    FullSolution->IterationHistory[i].PrimalObjective =
      Solution->IterationHistory[i].PrimalObjective;
    FullSolution->IterationHistory[i].DualObjective =
      Solution->IterationHistory[i].DualObjective;
    FullSolution->IterationHistory[i].PriInf =
      Solution->IterationHistory[i].PriInf;
    FullSolution->IterationHistory[i].DualInf =
      Solution->IterationHistory[i].DualInf;
    FullSolution->IterationHistory[i].logmu =
      Solution->IterationHistory[i].logmu;
    FullSolution->IterationHistory[i].NumCorrections =
      Solution->IterationHistory[i].NumCorrections;
    FullSolution->IterationHistory[i].phi =
      Solution->IterationHistory[i].phi;
  }

  FullSolution->Factorizations = Solution->Factorizations;
  FullSolution->FactorizationHistory = (FactorizationRecord *)
    Malloc(sizeof(FactorizationRecord), "Solution->FactorizationHistory");

  FullSolution->PrimalObjective = Solution->PrimalObjective + Record->cshift;
  FullSolution->DualObjective = Solution->DualObjective + Record->cshift;

  FullSolution->PrimalInfeasibility = Solution->PrimalInfeasibility;
  FullSolution->DualInfeasibility = Solution->DualInfeasibility;

  FullSolution->FactorizationHistory->Nonzeros =
    Solution->FactorizationHistory->Nonzeros;
  FullSolution->FactorizationHistory->Density =
    Solution->FactorizationHistory->Density;
  FullSolution->FactorizationHistory->Operations =
    Solution->FactorizationHistory->Operations;
  FullSolution->FactorizationHistory->NumDenseCols =
    Solution->FactorizationHistory->NumDenseCols;

  /* transfer reduced solution components into their rightful places in the
   * full solution vector. (Ignore DualUpper, DualLower; which are recovered
   * at the end) */

  for (i = 0, inx = 0; i < columns; i++)
    if (Record->ColumnMask[i + 1] == STILL_ACTIVE) {
      FullSolution->x[i] = Solution->x[inx];
      inx++;
    }
  if (inx != Solution->Columns) {
    printf(" reduced primal dimension %d contradicts ColumnMask %d\n",
	   inx, Solution->Columns);
    FreeSolution(FullSolution); 
    DeleteChangeStack(Record);
    return PREPROCESSING_ERROR;
  }
  for (i = 0, inx = 0; i < rows; i++)
    if (Record->RowMask[i + 1] == STILL_ACTIVE)
      FullSolution->pi[i] = Solution->pi[inx++];
  if (inx != Solution->Rows) {
    printf(" reduced dual dimension %d contradicts RowMask %d\n",
	   inx, Solution->Rows);
    FreeSolution(FullSolution); 
    DeleteChangeStack(Record);
    return PREPROCESSING_ERROR;
  }
  /* pop the stack an element at a time, adding components to the primal and
   * dual each time. */

  for (stackTop = Record->Top - 1; stackTop >= 0; stackTop--)
    switch (Record->StackOfChanges[stackTop]->ChangeType) {

    case DUPLICATE_COLUMN:
      pDuplicateCol =
	Record->StackOfChanges[stackTop]->pDuplicateCol;

      if (pDuplicateCol->Type == FREE) {

	if (FullSolution->x[pDuplicateCol->Split1 - 1] < 0.0) {
	  FullSolution->x[pDuplicateCol->Split2 - 1] =
	    FullSolution->x[pDuplicateCol->Split1 - 1] /
	    pDuplicateCol->multiplier;
	  FullSolution->x[pDuplicateCol->Split1 - 1] = 0.0;
	} else {
	  FullSolution->x[pDuplicateCol->Split2 - 1] = 0.0;
	}

	Record->VarType[pDuplicateCol->Split1 - 1] = NORMAL;
	Record->VarType[pDuplicateCol->Split2 - 1] = NORMAL;

      } else if (pDuplicateCol->Type == NORMAL) {

	if (pDuplicateCol->costdiff >= 0.0) {
	  FullSolution->x[pDuplicateCol->Split2 - 1] = 0.0;
	} else {
	  FullSolution->x[pDuplicateCol->Split2 - 1] =
	    FullSolution->x[pDuplicateCol->Split1 - 1] /
	    pDuplicateCol->multiplier;
	  FullSolution->x[pDuplicateCol->Split1 - 1] = 0.0;
	}

	Record->VarType[pDuplicateCol->Split1 - 1] = NORMAL;
	Record->VarType[pDuplicateCol->Split2 - 1] = NORMAL;
      }
      Record->ColumnMask[pDuplicateCol->Split2] = STILL_ACTIVE;
      break;

    case DUPLICATE_ROW:
      pDuplicateRow =
	Record->StackOfChanges[stackTop]->pDuplicateRow;
      FullSolution->pi[pDuplicateRow->Split2 - 1] = 0.0;
      Record->RowMask[pDuplicateRow->Split2] = STILL_ACTIVE;
      break;

    case ZERO_ROW:
      pZeroRow =
	Record->StackOfChanges[stackTop]->pZeroRow;
      FullSolution->pi[pZeroRow->RowIndex - 1] = 0.0;

      Record->RowMask[pZeroRow->RowIndex] = STILL_ACTIVE;
      break;

    case ZERO_COLUMN:
      pZeroColumn =
	Record->StackOfChanges[stackTop]->pZeroColumn;
      inx = pZeroColumn->ColumnIndex;
      FullSolution->x[inx - 1] = pZeroColumn->Value;

      Record->ColumnMask[inx] = STILL_ACTIVE;
      break;

    case FIXED_VARIABLE:
      pFixedVariable =
	Record->StackOfChanges[stackTop]->pFixedVariable;
      inx = pFixedVariable->ColumnIndex;
      FullSolution->x[inx - 1] = pFixedVariable->Value;

      temp = pFixedVariable->cElement;
      for (j = 0; j < pFixedVariable->DeletedColumn->size; j++) {
	inxt = pFixedVariable->DeletedColumn->Index[j];
	if (Record->RowMask[inxt] == STILL_ACTIVE)
	  temp -= pFixedVariable->DeletedColumn->Value[j] * FullSolution->pi[inxt - 1];
      }

      Record->ColumnMask[inx] = STILL_ACTIVE;
      break;

    case SINGLETON_ROW:
      pSingletonRow =
	Record->StackOfChanges[stackTop]->pSingletonRow;
      inx = pSingletonRow->ColumnIndex;
      inxr = pSingletonRow->RowIndex;
      FullSolution->x[inx - 1] = pSingletonRow->Value;

      temp = pSingletonRow->cElement;
      for (j = 0; j < pSingletonRow->DeletedColumn->size; j++) {
	inxt = pSingletonRow->DeletedColumn->Index[j];
	if (Record->RowMask[inxt] == STILL_ACTIVE && inxt != inxr)
	  temp -= pSingletonRow->DeletedColumn->Value[j] * FullSolution->pi[inxt - 1];
      }
      temp /= pSingletonRow->aElement;
      FullSolution->pi[inxr - 1] = temp;
      Record->RowMask[inxr] = STILL_ACTIVE;
      Record->ColumnMask[inx] = STILL_ACTIVE;
      break;

    case FREE_COLUMN_SINGLETON:
      pFreeColumnSingleton =
	Record->StackOfChanges[stackTop]->pFreeColumnSingleton;
      inx = pFreeColumnSingleton->ColumnIndex;
      inxr = pFreeColumnSingleton->RowIndex;

      temp = pFreeColumnSingleton->bElement;
      for (j = 0; j < pFreeColumnSingleton->DeletedRow->size; j++) {
	inxt = pFreeColumnSingleton->DeletedRow->Index[j];
	if (Record->ColumnMask[inxt] == STILL_ACTIVE && inxt != inx)
	  temp -= pFreeColumnSingleton->DeletedRow->Value[j] * FullSolution->x[inxt - 1];
      }
      temp /= pFreeColumnSingleton->aElement;
      FullSolution->x[inx - 1] = temp;
      FullSolution->pi[inxr - 1] =
	pFreeColumnSingleton->cElement / pFreeColumnSingleton->aElement;

      Record->RowMask[inxr] = STILL_ACTIVE;
      Record->ColumnMask[inx] = STILL_ACTIVE;
      break;

    case FORCED_ROW:
      pForcedRow = Record->StackOfChanges[stackTop]->pForcedRow;

      /* work through the deleted columns, one by one, assigning the primal
       * values and searching for the correct value of pi */
      inxr = pForcedRow->RowIndex;
      piUpper_inf = TRUE;
      piLower_inf = TRUE;
      for (k = 0; k < pForcedRow->Values->size; k++) {
	j = pForcedRow->Values->Index[k];
	if (Record->ColumnMask[j] == STILL_ACTIVE) {
	  printf("Error in post-process of forced row %d at column %d\n",
		 inxr, j);
	}
	FullSolution->x[j - 1] = pForcedRow->Values->Value[k];

	/* find the bound on pi that's implied by this column */
	temp = pForcedRow->cElements->Value[k];
	for (l = 0; l < pForcedRow->DeletedColumns[k]->size; l++) {
	  inx = pForcedRow->DeletedColumns[k]->Index[l];
	  if (Record->RowMask[inx] == STILL_ACTIVE)
	    temp -=
	      pForcedRow->DeletedColumns[k]->Value[l] * FullSolution->pi[inx - 1];
	}
	temp /= pForcedRow->DeletedRow->Value[k];

	/* does this column give an upper or lower bound on pi? */
	if (pForcedRow->LowerUpper == RANGELOWER) {
	  /* upper bound */
	  if (piUpper_inf) {
	    piUpper_inf = FALSE;
	    piUpper = temp;
	  } else if (temp < piUpper)
	    piUpper = temp;
	} else if (pForcedRow->LowerUpper == RANGEUPPER) {
	  /* lower bound */
	  if (piLower_inf) {
	    piLower_inf = FALSE;
	    piLower = temp;
	  } else if (temp > piLower)
	    piLower = temp;
	}
      }

      /* assign the value of pi */
      if (piUpper_inf && piLower_inf) {
	printf("Error in post-process of forced row %d. Inf range for pi\n",
	       inxr);
      } else if (piUpper_inf && !piLower_inf) {
	FullSolution->pi[inxr - 1] = piLower;
      } else if (!piUpper_inf && piLower_inf) {
	FullSolution->pi[inxr - 1] = piUpper;
      } else if (!piUpper_inf && !piLower_inf) {
	if (piUpper < piLower - TOLERANCE) {
	  printf(" Error in Postprocess(): Forced Row %d.\n", inxr);
	  printf(" pi bounds overlap. Lower %f Upper %f\n", piLower, piUpper);
          FreeSolution(FullSolution); 
          DeleteChangeStack(Record);
          return PREPROCESSING_ERROR;
	}
	/* arbitrarily, take Mr. In-Between. */
	FullSolution->pi[inxr - 1] = (piUpper + piLower) / 2.0;
      }
      /* Now assign duals for the bounds */
      for (k = 0; k < pForcedRow->Values->size; k++) {
	j = pForcedRow->Values->Index[k];
	temp = pForcedRow->cElements->Value[k];
	for (l = 0; l < pForcedRow->DeletedColumns[k]->size; l++) {
	  inx = pForcedRow->DeletedColumns[k]->Index[l];
	  if (Record->RowMask[inx] == STILL_ACTIVE)
	    temp -=
	      pForcedRow->DeletedColumns[k]->Value[l] * FullSolution->pi[inx - 1];
	}
	temp -= pForcedRow->DeletedRow->Value[k] * FullSolution->pi[inxr - 1];

	if (pForcedRow->LowerUpper == RANGELOWER &&
	    pForcedRow->DeletedRow->Value[k] > 0.0) {
	  /* variable j forced to lower bound */
	  if (temp < (-TOLERANCE)) {
	    printf(" Error in Postprocess(): Forced Row %d:\n", inxr);
	    printf(" Negative value %g for DualLower component %d.\n", temp, j);
            FreeSolution(FullSolution); 
            DeleteChangeStack(Record);
            return PREPROCESSING_ERROR;
	  }
	} else if (pForcedRow->LowerUpper == RANGEUPPER &&
		   pForcedRow->DeletedRow->Value[k] < 0.0) {
	  /* variable j forced to lower bound */
	  if (temp < (-TOLERANCE)) {
	    printf(" Error in Postprocess(): Forced row %d:\n", inxr);
	    printf(" Negative value %g for DualLower component %d.\n", temp, j);
            FreeSolution(FullSolution); 
            DeleteChangeStack(Record);
            return PREPROCESSING_ERROR;
	  }
	} else if (pForcedRow->LowerUpper == RANGELOWER &&
		   pForcedRow->DeletedRow->Value[k] < 0.0) {
	  /* variable j forced to upper bound */
	  if (temp > (TOLERANCE)) {
	    printf(" Error in Postprocess(): Forced row %d:\n", inxr);
	    printf(" Negative value %g for DualUpper component %d.\n", -temp, j);
            FreeSolution(FullSolution); 
            DeleteChangeStack(Record);
            return PREPROCESSING_ERROR;
	  }
	} else if (pForcedRow->LowerUpper == RANGEUPPER &&
		   pForcedRow->DeletedRow->Value[k] > 0.0) {
	  /* variable j forced to upper bound */
	  if (temp > (TOLERANCE)) {
	    printf(" Error in Postprocess: Forced row %d:\n", inxr);
	    printf(" Negative value %g for DualUpper component %d.\n", -temp, j);
            FreeSolution(FullSolution); 
            DeleteChangeStack(Record);
            return PREPROCESSING_ERROR;
	  }
	}
      }

      /* finally, activate the row and all the columns */
      Record->RowMask[inxr] = STILL_ACTIVE;
      for (k = 0; k < pForcedRow->Values->size; k++) {
	j = pForcedRow->Values->Index[k];
	Record->ColumnMask[j] = STILL_ACTIVE;
      }
      break;

    default:
      printf(" Change type not recognized %d\n",
	     Record->StackOfChanges[stackTop]->ChangeType);
      printf(" Stacktop = %d\n", stackTop);
    }

  /* Recompute DualUpper and DualLower, based on pi */
  RecomputeDualVariables(LP, FullSolution);

  /* check that all rows and columns are active again */
  err = 0;
  for (i = 1; i <= rows; i++)
    if (Record->RowMask[i] != STILL_ACTIVE) {
      printf("Postprocess: Row %d has not been freed\n", i);
      err--;
    }
  for (i = 1; i <= columns; i++)
    if (Record->ColumnMask[i] != STILL_ACTIVE) {
      printf("Postprocess: Column %d has not been freed\n", i);
      err--;
    }
  if(err<0) { 
    FreeSolution(FullSolution); 
    DeleteChangeStack(Record);
    return PREPROCESSING_ERROR;
  }
  
  /* free memory for old Solution (*pSolution) */
  FreeSolution(Solution); 
  
  /* free memory for Record (no longer needed) */
  DeleteChangeStack(Record);

  *pSolution = FullSolution;
  return 0;
}


