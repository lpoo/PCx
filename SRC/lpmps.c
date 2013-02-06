/* LPtype <--> MPStype
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

/* Routine to expand the LP from the MPS file into the format handled by the
 * LP code (and the preprocessor), namely,
 * 
 * min c.x ,   Ax=b,   0 <= x <= u.
 * 
 * The components of x are allowed to be one of three types: NORMAL, for which
 * x >= 0; UPPER,    for which u >= x >= 0; FREE,     where x is
 * unconstrained. We usually have to do some shifting and sign flipping to
 * squeeze everything into these three categories. We save a record of these
 * changes, which will be inverted as the very last step of the process.
 * 
 * It's assumed that we keep the original MPS structure lying around somewhere,
 * so the row and column names, original ranges, etc, don't have to be
 * duplicated.  */


/*****************************************************************************/
/*  contents:                                                                */
/*****************************************************************************/
/*
LPtype  *Convert_MPS_LP(MPStype *MPS, MPSchanges **pChanges);

solution *MPSsolution(LPtype *LP, MPStype *MPS, solution *LPsolution, 
		      MPSchanges *Changes, Parameters *Inputs);

LPtype   *MakeTranspose(LPtype *LP);
    
int DeleteLP(LPtype *LP);
*/
/*****************************************************************************/

LPtype  *Convert_MPS_LP(MPS, pChanges)
     MPStype        *MPS;
     MPSchanges    **pChanges;
{
   int             i, j, k, rows, cols, sigcols, ents, sigents, extras,
                   inxrow, entry;
   int            *RowMap, count;
   LPtype         *LP, *MakeTranspose(), *NewLP();
   MPSchanges     *NewChanges(), *Changes;
   
   /* allocate space for the record of Changes and for the LP data structure */
   
   *pChanges = NewChanges(MPS->NumRows, MPS->NumCols);
   Changes = *pChanges;
   
   sigcols = cols = MPS->NumCols;
   
   /* count up the number of slacks to add */
   
   for (i = 0, extras = 0; i < MPS->NumRows; i++)
      if (toupper(MPS->RowType[i]) == 'L' ||
	  toupper(MPS->RowType[i]) == 'G' ||
	  (toupper(MPS->RowType[i]) == 'E' && MPS->Ranges[i] != 0.0))
	 extras++;
   
   /* count number of free rows and construct RowMap */
   
   RowMap = NewInt(MPS->NumRows, "RowMap");
   
   for (i = 0, rows = 0; i < MPS->NumRows; i++)
      if (toupper(MPS->RowType[i]) != 'N')
	 RowMap[i] = rows++;
   
   /* count elements not in free rows and not zero.  "sigents" is the number
    * of significant entries, which doesn't include the nonzeros due to
    * slacks, etc.  */
   
   for (i = 0, ents = 0; i < MPS->NumCols; i++)
      for (k = MPS->A.pBeginRow[i] - 1;
	   k <= MPS->A.pEndRow[i] - 1; k++) 
	 {
	    j = MPS->A.Row[k] - 1;
	    if (toupper(MPS->RowType[j]) != 'N' && MPS->A.Value[k] != 0.0)
	       ents++;
	 }
  sigents = ents;
  
  /* now we can assign the remaining space in LP */
  
  cols += extras;
  ents += extras;
  LP = NewLP(rows, cols, ents);
  
  /* copy across the data from MPS, but ignore the zeros in the coefficient
   * matrix.  */
  
  count = 0;
  for (i = 0; i < sigcols; i++) 
     {
	LP->c[i] = MPS->c[i];
	LP->A.pBeginRow[i] = count + 1;
	
	for (entry = MPS->A.pBeginRow[i] - 1;
	     entry <= MPS->A.pEndRow[i] - 1; entry++) 
	   {
	      j = MPS->A.Row[entry] - 1;
	      if (toupper(MPS->RowType[j]) != 'N' && 
		  MPS->A.Value[entry] != 0.0) 
		 {
		    LP->A.Row[count] = RowMap[j] + 1;
		    LP->A.Value[count] = MPS->A.Value[entry];
		    count++;		/* added new element to array */
		 }
	   }
	LP->A.pEndRow[i] = count;
     }
  
  if (count != sigents)
     printf("Error counting elements of matrix in Convert_MPS_LP.\n");
  
  LP->cshift = MPS->cshift;
  
  for (i = 0; i < MPS->NumRows; i++)
     if (toupper(MPS->RowType[i]) != 'N')
	LP->b[RowMap[i]] = MPS->b[i];
  
  /* Process the existing variables, shifting and flipping where necessary */
  
  for (i = 0; i < sigcols; i++) 
     {
	if (MPS->BoundType[i] == PINFTY) 
	   {	/* easy - don't change a thing */
	      LP->VarType[i] = NORMAL;
	      
	   } 
	else if (MPS->BoundType[i] == FREE) 
	   {	/* easy too! */
	      LP->VarType[i] = FREE;
	      LP->FreeIndex[LP->NumberFree] = i;
	      LP->NumberFree++;
	   } 
	else if (MPS->BoundType[i] == UPPER) 
	   {	/* also pretty easy */
	      LP->UpBound[i] = MPS->UpBound[i];
	      LP->VarType[i] = UPPER;
	      LP->BoundIndex[LP->NumberBounds] = i;
	      LP->NumberBounds++;
	   } 
	else if (MPS->BoundType[i] == LOWER) 
	   {	/* shift the sucker */
	      Changes->VarShifts[i] = MPS->LowBound[i];
	      LP->cshift += LP->c[i] * MPS->LowBound[i];
	      for (j = LP->A.pBeginRow[i] - 1; 
		   j <= LP->A.pEndRow[i] - 1; j++) 
		 {
		    inxrow = LP->A.Row[j] - 1;
		    LP->b[inxrow] -= LP->A.Value[j] * MPS->LowBound[i];
		 }
	      LP->VarType[i] = NORMAL;
	   } 
	else if (MPS->BoundType[i] == UPPERLOWER) 
	   {	/* shift the sucker */
	      Changes->VarShifts[i] = MPS->LowBound[i];
	      LP->cshift += LP->c[i] * MPS->LowBound[i];
	      for (j = LP->A.pBeginRow[i] - 1; 
		   j <= LP->A.pEndRow[i] - 1; j++) 
		 {
		    inxrow = LP->A.Row[j] - 1;
		    LP->b[inxrow] -= LP->A.Value[j] * MPS->LowBound[i];
		 }
	      LP->UpBound[i] = MPS->UpBound[i] - MPS->LowBound[i];
	      LP->VarType[i] = UPPER;
	      LP->BoundIndex[LP->NumberBounds] = i;
	      LP->NumberBounds++;
	      
	   } 
	else if (MPS->BoundType[i] == FIX) 
	   {	/* shift the sucker */
	      Changes->VarShifts[i] = MPS->LowBound[i];
	      LP->cshift += LP->c[i] * MPS->LowBound[i];
	      for (j = LP->A.pBeginRow[i] - 1; 
		   j <= LP->A.pEndRow[i] - 1; j++) 
		 {
		    inxrow = LP->A.Row[j] - 1;
		    LP->b[inxrow] -= LP->A.Value[j] * MPS->LowBound[i];
		 }
	      LP->UpBound[i] = 0.0;
	      LP->VarType[i] = UPPER;
	      LP->BoundIndex[LP->NumberBounds] = i;
	      LP->NumberBounds++;
	   } 

/* changed 2/8/01 to handle the case in which the LP->UpBound[i] contains 
   a nonzero, which is meant to be an upper bound on x[i], i.e.
   -inf <= x[i] <= MPS->UpBound[i].
   Need to convert this to a NORMAL variable by flipping and shifting */

	else if (MPS->BoundType[i] == MINFTY) 
	   {	/* flip the sign */
	      Changes->SignChanges[i] = 1;
	      for (j = LP->A.pBeginRow[i] - 1; 
		   j <= LP->A.pEndRow[i] - 1; j++) 
		 {
		    LP->A.Value[j] = -LP->A.Value[j];
		 }
	      LP->c[i] = -LP->c[i];
                /* shift the sucker too, 
                   its lower bound is -MPS->UpBound[i] */ 
              if (MPS->UpBound[i] != 0.0) {
                Changes->VarShifts[i] = -MPS->UpBound[i]; 
                LP->cshift -= LP->c[i] * MPS->UpBound[i]; 
                for (j = LP->A.pBeginRow[i] - 1;
                     j <= LP->A.pEndRow[i] - 1; j++)
                   {
                    inxrow = LP->A.Row[j] - 1;
                    LP->b[inxrow] += LP->A.Value[j] * MPS->UpBound[i];
                   }   
              }
	      LP->VarType[i] = NORMAL;
	   }
     }

  /* now add the slacks */
  
  for (i = 0; i < MPS->NumRows; i++) 
     {
	if (toupper(MPS->RowType[i]) == 'L') 
	   {
	      if (MPS->Ranges[i] == 0.0) 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = 1.0;
		    LP->VarType[sigcols] = NORMAL;
		 } 
	      else 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = 1.0;
		    LP->VarType[sigcols] = UPPER;
		    LP->UpBound[sigcols] = fabs(MPS->Ranges[i]);
		    LP->BoundIndex[LP->NumberBounds] = sigcols;
		    LP->NumberBounds++;
		 }
	      sigcols++;
	      sigents++;
	   } 
	else if (toupper(MPS->RowType[i]) == 'G') 
	   {
	      if (MPS->Ranges[i] == 0.0) 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = -1.0;
		    LP->VarType[sigcols] = NORMAL;
		 } 
	      else 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = -1.0;
		    LP->VarType[sigcols] = UPPER;
		    LP->UpBound[sigcols] = fabs(MPS->Ranges[i]);
		    LP->BoundIndex[LP->NumberBounds] = sigcols;
		    LP->NumberBounds++;
		 }
	      sigcols++;
	      sigents++;
	   } 
	else if ((toupper(MPS->RowType[i]) == 'E') 
		 && (MPS->Ranges[i] != 0.0)) 
	   {
	      if (MPS->Ranges[i] > 0.0) 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = -1.0;
		    LP->VarType[sigcols] = UPPER;
		    LP->UpBound[sigcols] = MPS->Ranges[i];
		 } 
	      else 
		 {
		    LP->c[sigcols] = 0.0;
		    LP->A.pBeginRow[sigcols] = sigents + 1;
		    LP->A.pEndRow[sigcols] = sigents + 1;
		    LP->A.Row[sigents] = RowMap[i] + 1;
		    LP->A.Value[sigents] = 1.0;
		    LP->VarType[sigcols] = UPPER;
		    LP->UpBound[sigcols] = -MPS->Ranges[i];
		    LP->BoundIndex[LP->NumberBounds] = sigcols;
		    LP->NumberBounds++;
		 }
	      sigcols++;
	      sigents++;
	   }
     }
  
  /* check that we added the right number of slacks */
  
  if (sigents != ents || sigcols != cols) 
     {
	printf("Got confused when adding slacks!:\n");
	printf("sigents = %d  ents = %d\n", sigents, ents);
	printf("sigcols = %d  cols = %d\n", sigcols, cols);
     }
  /* Find the transpose of the coefficient matrix */
  LP = MakeTranspose(LP);
  Free ((char *) RowMap);
  return (LP);
}


/************************************************************************/

/* given a solution to the LP, this routines inverts the shifts and sign
 * flips to produce a solution to the original MPS form.  SJW 12/21/94.  */

solution       *MPSsolution(LP, MPS, LPsolution, Changes, Inputs)
     LPtype         *LP;
     MPStype        *MPS;
     solution       *LPsolution;
     MPSchanges     *Changes;
     Parameters     *Inputs;
{
   solution       *Solution, *NewSolution();
   int             i, j, k, rows, cols, lprow, col, entry, row;
   double          cshift;
   
   /* Allocate space for new solution record 
      and copy across the common stuff */
   
   Solution = (solution *) Malloc(sizeof(solution), "Solution");

   rows = Solution->Rows = MPS->NumRows;
   cols = Solution->Columns = MPS->NumCols;
   
   Solution->x = NewDouble(cols, "Solution->x");
   Solution->pi = NewDouble(rows, "Solution->pi");
   Solution->DualUpper = NewDouble(cols, "Solution->DualUpper");
   Solution->DualLower = NewDouble(cols, "Solution->DualLower");
   
   Solution->ReadTime = LPsolution->ReadTime;
   Solution->PreprocessTime = LPsolution->PreprocessTime;
   Solution->SolutionTime = LPsolution->SolutionTime;
   Solution->FactorizationTime = LPsolution->FactorizationTime;
   Solution->SolveADATTime = LPsolution->SolveADATTime;
   Solution->InitTime = LPsolution->InitTime;
   Solution->LoopTime = LPsolution->LoopTime;
   Solution->PredictorTime = LPsolution->PredictorTime;
   Solution->FormADATtime = LPsolution->FormADATtime;
   Solution->CorrectorTime = LPsolution->CorrectorTime;
  
   Solution->PriorRows = LPsolution->PriorRows;
   Solution->PriorColumns = LPsolution->PriorColumns;
   Solution->ReducedRows = LPsolution->ReducedRows;
   Solution->ReducedColumns = LPsolution->ReducedColumns;
   Solution->Passes = LPsolution->Passes;
   
   Solution->Iterations = LPsolution->Iterations;
   Solution->RestoredIteration = LPsolution->RestoredIteration;
   
   Solution->Status = LPsolution->Status;
   Solution->IterationHistory = (IterationRecord *)
      Malloc((Solution->Iterations + 1) * sizeof(IterationRecord),
	     "Solution->IterationHistory");
   
   Solution->PrimalObjective = LPsolution->PrimalObjective + LP->cshift;
   Solution->DualObjective = LPsolution->DualObjective + LP->cshift;
   
   if (!Inputs->Minimize) 
      {
	 /* Change sign on objective value */
	 Solution->PrimalObjective = -Solution->PrimalObjective;
	 Solution->DualObjective = -Solution->DualObjective;
      }
   Solution->PrimalInfeasibility = LPsolution->PrimalInfeasibility;
   Solution->DualInfeasibility = LPsolution->DualInfeasibility;
   
   Solution->Complementarity = LPsolution->Complementarity;
   Solution->RelativeComplementarity = LPsolution->RelativeComplementarity;
   
   cshift = Changes->cshift;
   for (i = 0; i <= Solution->Iterations; i++) 
      {
	 Solution->IterationHistory[i].PrimalObjective =
	    LPsolution->IterationHistory[i].PrimalObjective + cshift;
	 Solution->IterationHistory[i].DualObjective =
	    LPsolution->IterationHistory[i].DualObjective + cshift;
	 
	 if (!Inputs->Minimize) 
	    {
	    Solution->IterationHistory[i].PrimalObjective =
	       -Solution->IterationHistory[i].PrimalObjective;
	    
	    Solution->IterationHistory[i].DualObjective =
	       -Solution->IterationHistory[i].DualObjective;
	    }
	 Solution->IterationHistory[i].PriInf =
	    LPsolution->IterationHistory[i].PriInf;
	 Solution->IterationHistory[i].DualInf =
	    LPsolution->IterationHistory[i].DualInf;
	 Solution->IterationHistory[i].logmu =
	    LPsolution->IterationHistory[i].logmu;
	 Solution->IterationHistory[i].NumCorrections =
	    LPsolution->IterationHistory[i].NumCorrections;
	 Solution->IterationHistory[i].phi =
	    LPsolution->IterationHistory[i].phi;
      }
   
   Solution->Factorizations = LPsolution->Factorizations;

   
   Solution->FactorizationCode = (char *) Malloc((strlen(LPsolution->FactorizationCode) + 1)*sizeof(char), "Solution->FactorizationCode");
   strcpy(Solution->FactorizationCode, LPsolution->FactorizationCode);

   Solution->FactorizationHistory = (FactorizationRecord *)
      Malloc(sizeof(FactorizationRecord), "Solution->FactorizationHistory");
   Solution->FactorizationHistory->Nonzeros =
      LPsolution->FactorizationHistory->Nonzeros;
   Solution->FactorizationHistory->Density =
      LPsolution->FactorizationHistory->Density;
   Solution->FactorizationHistory->Operations =
      LPsolution->FactorizationHistory->Operations;
   Solution->FactorizationHistory->NumDenseCols =
      LPsolution->FactorizationHistory->NumDenseCols;
   
   /* copy across the "significant" part of the primal solution, adjusting for
    * the shift where necessary */
   
   for (i = 0, lprow = 0; i < rows; i++)
      if (toupper(MPS->RowType[i]) != 'N')
	 Solution->pi[i] = LPsolution->pi[lprow++];
      else
	 Solution->pi[i] = 0.0;
   
   for (i = 0; i < MPS->NumCols; i++) 
      {
	 Solution->x[i] = LPsolution->x[i] + Changes->VarShifts[i];
	 Solution->DualUpper[i] = LPsolution->DualUpper[i];
	 Solution->DualLower[i] = LPsolution->DualLower[i];
	 if (Changes->SignChanges[i]) 
	    {
	       Solution->x[i] = -Solution->x[i];
	       Solution->DualUpper[i] = -Solution->DualUpper[i];
	       Solution->DualLower[i] = -Solution->DualLower[i];
	    }
      }
   
   /* Compute Row activity for MPS solution*/
   
   Solution->Activity = NewDouble(MPS->NumRows, "Solution->Activity");
   for (col = 0; col < MPS->NumCols; col++)
      for (entry = MPS->A.pBeginRow[col] - 1;
	   entry <= MPS->A.pEndRow[col] - 1; entry++) 
	 {
	    row = MPS->A.Row[entry] - 1;
	    Solution->Activity[row] += MPS->A.Value[entry] * Solution->x[col];
	 }
   
   /* Free old solution structure */
   FreeSolution(LPsolution);
   
   return (Solution);
}

/************************************************************************/
/* Fills in the Atranspose part of an LP data structure. Assumes that A is
 * already there, and that space has been allocated for Atranspose. */

LPtype         *MakeTranspose(LP)
     LPtype         *LP;
{
   int             i, j, k, index, rows, cols, *RowCount;
   
   rows = LP->Rows;
   cols = LP->Cols;
   
   /* Find the row lengths */
   RowCount = NewInt(rows, "RowCount");
   for (i = 0; i < cols; i++)
      for (j = LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++) 
	 {
	    k = (LP->A.Row[j - 1]) - 1;
	    RowCount[k]++;
	 }
   
   /* Adjust the begin/end pointers for the rows */
   LP->Atranspose.pBeginRow[0] = 1;
   LP->Atranspose.pEndRow[0] = RowCount[0];
   for (i = 1; i < rows; i++) 
      {
	 LP->Atranspose.pBeginRow[i] = LP->Atranspose.pEndRow[i - 1] + 1;
	 LP->Atranspose.pEndRow[i] = LP->Atranspose.pBeginRow[i] 
	    + RowCount[i] - 1;
      }
   
   /* Now put all the elements in place */
   for (j = 0; j < rows; j++)
      RowCount[j] = 0;
   for (i = 0; i < cols; i++)
      for (j = LP->A.pBeginRow[i] - 1; j <= LP->A.pEndRow[i] - 1; j++) 
	 {
	    k = LP->A.Row[j];
	    index = LP->Atranspose.pBeginRow[k - 1] + RowCount[k - 1];
	    LP->Atranspose.Row[index - 1] = i + 1;
	    LP->Atranspose.Value[index - 1] = LP->A.Value[j];
	    RowCount[k - 1]++;
	 }
   
   Free((char *) RowCount);
   return (LP);
}

/************************************************************************/

int DeleteLP(LP)
     LPtype         *LP;
{
   Free((char *) LP->A.pBeginRow);
   Free((char *) LP->A.pEndRow);
   Free((char *) LP->A.Row);
   Free((char *) LP->A.Value);
   
   Free((char *) LP->Atranspose.pBeginRow);
   Free((char *) LP->Atranspose.pEndRow);
   Free((char *) LP->Atranspose.Row);
   Free((char *) LP->Atranspose.Value);
   
   Free((char *) LP->b);
   Free((char *) LP->c);
   Free((char *) LP->VarType);
   Free((char *) LP->UpBound);
   
   Free((char *) LP->BoundIndex);
   Free((char *) LP->FreeIndex);
   
   Free((char *) LP->FreePlus);
   Free((char *) LP->FreeMinus);
   
   Free((char *) LP->ColScale);
   Free((char *) LP->RowScale);
   
   Free((char *) LP);
   return 0;
}
