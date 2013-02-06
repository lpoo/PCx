/* memory allocation
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
/* #include <malloc.h> */
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "memory.h"

void 
Free(c)
     char *c;
{
   /* TrFree(c,__LINE__,__FILE__); */
   free(c);
   return;
}

char  *
Malloc(size, message)
     int             size;
     char           *message;
{
   char           *ptr;
   
   if (size <= 0)
      size=1;

   ptr = (char *) malloc((unsigned) size);

  if (ptr == NULL) 
     {
	printf("Error in Malloc allocating %d for the variable %s.\n",
	       size, message);
	OutOfSpace();
     }
  return (ptr);
}

char   *
Calloc(num, size, message)
     int             num, size;
     char           *message;
{
   char           *ptr;
   
   if (size <= 0)
      size=1;
   
   ptr = (char *) calloc((unsigned) num, (unsigned) size);
   
#if (DEBUG > 0)
   {
      int i;
      for (i = 0; i < num*size; ++i)
	 ptr[i] = '\0';
   }
#endif
   
   if (ptr == NULL) 
      {
	 printf("Error in Calloc allocating %d for variable %s.\n", 
		num, message);
	 OutOfSpace();
      }
   return (ptr);
}

/* various allocation routines */

Iterate     *
NewIterate(NumRows, NumCols, NumBounds)
     int             NumRows, NumCols, NumBounds;
{
   Iterate        *ptr;
   
   ptr = (Iterate *) Malloc(sizeof(Iterate), "ptr in NewIterate");
   
   ptr->NumRows = NumRows;
   ptr->NumCols = NumCols;
   ptr->NumBounds = NumBounds;
   ptr->x = NewDouble(NumCols, "ptr->x");
   ptr->s = NewDouble(NumCols, "ptr->s");
   ptr->pi = NewDouble(NumRows, "ptr->pi");
   ptr->w = NewDouble(NumBounds, "ptr->w");
   ptr->r = NewDouble(NumBounds, "ptr->r");
   ptr->BoundIndex = NULL;
   return ptr;
}

solution   *
NewSolution(Rows, Cols, MaxIterations)
     int             Rows, Cols, MaxIterations;
{
   
   solution       *Solution;
   
   Solution = (solution *) Malloc(sizeof(solution), "Solution");
   Solution->Rows = Rows;
   Solution->Columns = Cols;
   Solution->x = NewDouble(Cols, "Solution->x");
   Solution->pi = NewDouble(Rows, "Solution->pi");
   Solution->DualUpper = NewDouble(Cols, "Solution->DualUpper");
   Solution->DualLower = NewDouble(Cols, "Solution->DualLower");
   Solution->PrimalObjective = 0.0;
   Solution->DualObjective = 0.0;
   Solution->Complementarity = 0.0;
   Solution->RelativeComplementarity = 0.0;
   Solution->PrimalInfeasibility = 0.0;
   Solution->DualInfeasibility = 0.0;
   Solution->Iterations = 0;
   Solution->IterationHistory = (IterationRecord *)
      Malloc((MaxIterations + 1)*sizeof(IterationRecord), "IterationHistory");
   Solution->Factorizations = 0;
   Solution->FactorizationHistory = (FactorizationRecord *)
      Malloc(sizeof(FactorizationRecord), "FactorizationHistory");
   Solution->FactorizationHistory->Operations = 0;
   Solution->FactorizationCode = (char *)Malloc(50*sizeof(char), 
						"Solution->FactorizationCode");
   strcpy(Solution->FactorizationCode, "not specified");
   Solution->RestoredIteration = -1;
   Solution->Activity = NULL;
   return Solution;
}

int
FreeSolution(Solution)
     solution       *Solution;
{
   Free((char *) Solution->x);
   Free((char *) Solution->pi);
   Free((char *) Solution->DualUpper);
   Free((char *) Solution->DualLower);
   Free((char *) Solution->IterationHistory);
   Free((char *) Solution->FactorizationHistory);
   Free((char *) Solution->Activity);
   Free((char *) Solution->FactorizationCode);
   Free((char *) Solution);
   return 0;
}

MMTtype    *
NewMMTtype(NumRows, NumCols, Nonzeros)
     int             NumRows, NumCols, Nonzeros;
{
   MMTtype        *MMT;
   
   MMT = (MMTtype *) Malloc(sizeof(MMTtype), "MMT");
   MMT->NumRows = NumRows;
   MMT->NumCols = NumCols;
   MMT->Nonzeros = Nonzeros;
   MMT->pBeginRow = NewInt(NumCols, "MMT->pBeginRow");
   MMT->pBeginRowT = NewInt(NumRows, "MMT->pBeginRow");
   MMT->pEndRow = NewInt(NumCols, "MMT->pBeginCol");
   MMT->pEndRowT = NewInt(NumRows, "MMT->pBeginCol");
   MMT->Row = NewInt(Nonzeros, "MMT->Row");
   MMT->RowT = NewInt(Nonzeros, "MMT->Row");
   MMT->Value = NewDouble(Nonzeros, "MMT->Value");
   MMT->ValueT = NULL;
   return MMT;
}

MPSchanges   *
NewChanges(Rows, Cols)
     int             Rows, Cols;
{
   MPSchanges     *Changes;
   
   Changes = (MPSchanges *) Malloc(sizeof(MPSchanges), "Changes");
   Changes->Rows = Rows;
   Changes->Cols = Cols;
   Changes->SignChanges = NewInt(Cols, "Changes->SignChanges");
   Changes->VarShifts = NewDouble(Cols, "Changes->VarShifts");
   Changes->cshift = 0.0;
   return Changes;
}

void		
DeleteChanges(Changes)
     MPSchanges	*Changes;
{
   Free((char*)Changes->VarShifts);
   Free((char*)Changes->SignChanges);
   Free((char*)Changes);
}

MPStype     *
NewMPS(rows, cols, ents)
     int             rows, cols, ents;
{
   MPStype        *MPS;
   
   MPS = (MPStype *) Malloc(sizeof(MPStype), "MPS");
   MPS->RowSize = rows;
   MPS->ColSize = cols;
   MPS->EntSize = ents;
   MPS->NumRows = 0;
   MPS->NumCols = 0;
   MPS->NumEnts = 0;
   
   MPS->A.pBeginRow = NewInt(cols, "MPS->A.pBeginRow");
   MPS->A.pEndRow = NewInt(cols, "MPS->A.pEndRow");
   MPS->A.Row = NewInt(ents, "MPS->A.Row");
   MPS->A.Value = NewDouble(ents, "MPS->A.Value");
   
   MPS->A.NumRows = 0;
   MPS->A.NumCols = 0; 
   MPS->A.Nonzeros = 0;
   
   MPS->b = NewDouble(rows, "MPS->b");
   MPS->c = NewDouble(cols, "MPS->c");
   MPS->cshift = 0.0;
   
   MPS->Ranges = NewDouble(rows, "MPS->Ranges");
   
   MPS->RowType = NewChar(rows, "MPS->RowType");
   
   MPS->BoundType = NewInt(cols, "MPS->BoundType");
   MPS->UpBound = NewDouble(cols, "MPS->UpBound");
   MPS->LowBound = NewDouble(cols, "MPS->LowBound");
   
   MPS->RowNames = NewCharPtr(rows, "MPS->RowNames");
   MPS->ColNames = NewCharPtr(cols, "MPS->ColNames");
   
   MPS->RowTable = NewHashTable(rows);
   MPS->ColTable = NewHashTable(cols);
   
   MPS->ProblemName   = NULL;
   MPS->ObjectiveName = NULL;
   MPS->RHSName       = NULL;
   MPS->RangeName     = NULL;
   MPS->BoundName     = NULL;
   
   return MPS;
}

LPtype      *
NewLP(Rows, Cols, Ents)
     int             Rows, Cols, Ents;
{
   LPtype         *LP;
   
   LP = (LPtype *) Malloc(sizeof(LPtype), "LP");

   LP->Atranspose.pBeginRow = NewInt(Rows, "LP->Atranspose.pBeginRow");
   LP->Atranspose.pEndRow = NewInt(Rows, "LP->Atranspose.pEndRow");
   LP->Atranspose.Row = NewInt(Ents, "LP->Atranspose.Row");
   LP->Atranspose.Value = NewDouble(Ents, "LP->Atranspose.Value");
   
   LP->Rows = Rows;
   LP->A.NumRows = Rows;
   LP->Atranspose.NumRows = Cols;
   LP->Cols = Cols;
   LP->A.NumCols = Cols;
   LP->Atranspose.NumCols = Rows;
   LP->Ents = Ents;
   LP->A.Nonzeros = Ents;
   LP->Atranspose.Nonzeros = Ents;
   LP->cshift = 0.0;
   
   LP->b = NewDouble(Rows, "LP->b");
   LP->c = NewDouble(Cols, "LP->c");

   LP->VarType = NewInt(Cols, "LP->BoundType");
   LP->UpBound = NewDouble(Cols, "LP->UpBound");
   LP->NumberBounds = 0;
   LP->BoundIndex = NewInt(Cols, "LP->BoundIndex");
   LP->NumberFree = 0;
   LP->FreeIndex = NewInt(Cols, "LP->FreeIndex");
   
   LP->A.pBeginRow = NewInt(Cols, "LP->A.pBeginRow");
   LP->A.pEndRow = NewInt(Cols, "LP->A.pEndRow");
   LP->A.Row = NewInt(Ents, "LP->A.Row");
   LP->A.Value = NewDouble(Ents, "LP->A.Value");
   
   LP->FreePlus  = NULL;
   LP->FreeMinus = NULL;
   LP->NumberSplit = 0;
   
   LP->ColScale  = NULL;
   LP->RowScale  = NULL;

   
   return LP;
}


char  *
StrDup(s1, message)
     char           *s1, *message;
{
   char           *ptr;
   
   ptr = strdup(s1);
   if (ptr == NULL) 
      {
	 printf("Error allocating space for %s\n", message);
	 OutOfSpace();
      }
   return (ptr);
   
}

double   *
NewDouble(size, message)
     int             size;
     char           *message;
{
   double         *ptr;
   
   if (size <= 0)
      size=1;
   ptr = (double *) Calloc(size, sizeof(double), message);
   return (ptr);
}

double        **
NewDouble2(size1, size2, message)
     int             size1, size2;
     char           *message;
{
   double        **ptr;
   int             i;
   
   if (size1 <= 0)
      size1=1;
   if (size2 <= 0)
      size2=1;
   ptr = (double **) malloc((unsigned) (size1 * sizeof(double *)));
   if (ptr == NULL) 
      {
	 printf("Error in NewDouble2 allocating the variable %s.\n", 
		message);
	 OutOfSpace();
      }
   
   ptr[0] = NewDouble(size1*size2, "ptr[]");
   if (ptr[0] == NULL) 
      {
	 printf("Error in NewDouble2 allocating the variable %s[].\n", 
		message);
	 OutOfSpace();
      }
   for (i = 1; i < size1; i++) 
      ptr[i] = ptr[i-1] + size2;
      
   return (ptr);
}

void   **
FreeDouble2(ptr)
     double** ptr;
{
   free((char*) ptr[0]);
   free((char*) ptr);
   return NULL;
}

int    *
NewInt(size, message)
     int             size;
     char           *message;
{
   int            *ptr;
   
   if (size <= 0)
      size=1;
   ptr = (int *) Calloc((unsigned) size, (unsigned) sizeof(int), message);
   return (ptr);
}

char  *
Realloc(oldptr, size, message)
     char           *oldptr;
     int             size;
     char           *message;
{
   char           *ptr;
   
   if (size <= 0)
      size=1;
   if (oldptr == NULL)
      ptr = Malloc(size, message);
   else
      ptr = (char *) realloc(oldptr, (unsigned) size);
   if (ptr == NULL) {
      printf("Error in Realloc allocating the variable %s.\n", message);
      OutOfSpace();
   }
   return (ptr);
}

double **
NewDoublePtr(size, message)
     int             size;
     char           *message;
     
{
   double        **ptr;
   int             i;
   
   if (size <= 0)
      size=1;
   ptr = (double **) Malloc(size * sizeof(double *), message);
   if (ptr == NULL) {
      printf("Error in NewDoublePtr allocating the variable %s.\n", message);
      OutOfSpace();
   }
   for (i = 0; i < size; i++)
      ptr[i] = NULL;
   return (ptr);
}

int   **
NewIntPtr(size, message)
     int             size;
     char           *message;
     
{
   int           **ptr, i;
   
   if (size <= 0)
      size=1;
   ptr = (int **) Malloc(size * sizeof(int *), message);
   if (ptr == NULL) 
      {
	 printf("Error in NewIntPtr allocating %d for %s.\n", size, message);
	 OutOfSpace();
      }
   for (i = 0; i < size; i++)
      ptr[i] = NULL;
   return (ptr);
}

char   *
NewChar(size, message)
     int             size;
     char           *message;
{
   char           *ptr;
   int             i;
   
   if (size <= 0)
      size=1;
   ptr = (char *) Malloc(size * sizeof(char), message);
   if (ptr == NULL) 
      {
	 printf("Error in NewChar allocating %d for %s.\n", size, message);
	 OutOfSpace();
      }
   for (i = 0; i < size; i++)
      ptr[i] = '\0';		/* each character is NULL */
  return (ptr);
}

char    **
NewCharPtr(size, message)
     int             size;
     char           *message;
{
   char          **ptr;
   int             i;
   
   if (size <= 0)
      size=1;
   ptr = (char **) Malloc(size * sizeof(char *), message);
   if (ptr == NULL) 
      {
	 printf("Error in NewCharPtr allocating %d for %s.\n", size, message);
	 OutOfSpace();
      }
   for (i = 0; i < size; i++)
      ptr[i] = NULL;
   return (ptr);
}

void 
OutOfSpace()
{
   printf("**Out of Memory**\n");
   exit(MEMORY_ERROR);
}


