/* readmps
 *
 * PCx 1.1 11/97
 *
 * Author: Joe Czyzyk, Steve Wright
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 *
 * last revised 2/27/2001
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "main.h"
#include "def.h"

/* Read the MPS file into an MPStype data structure, sorting the columns in
 * the sparse matrix structure so that the row indices appear in order in
 * each column. Also, keep track of the time needed to read and sort. */

MPStype        *ReadMPS(fp, filesize, Inputs, Time)
  FILE           *fp;
  int             filesize;
  Parameters     *Inputs;
  double         *Time;
{
  double          userTime0, systemTime0, userTime1, systemTime1;
  MPStype        *MPS, *NewMPS();
  int             ReadMPSFile();

  MPS = NewMPS(filesize / 150, filesize / 125, filesize / 30);

  GetTime(&userTime0, &systemTime0);

  /* set objective, rhs, range, and bound names */

  if (Inputs->ObjectiveName)
    MPS->ObjectiveName = StrDup(Inputs->ObjectiveName, "ObjectiveName");
  if (Inputs->RHSName)
    MPS->RHSName = StrDup(Inputs->RHSName, "RHSName");
  if (Inputs->RangeName)
    MPS->RangeName = StrDup(Inputs->RangeName, "RangeName");
  if (Inputs->BoundName)
    MPS->BoundName = StrDup(Inputs->BoundName, "BoundName");

  /* Read the file */
  if(ReadMPSFile(fp, MPS, Inputs) == BAD_INPUT) {
     printf("Error detected in MPS file\n");
     fclose(fp);
     return NULL;
  }

  /* Sort each column of the sparse matrix A so that the row indices appear
   * in order. (The preprocessor likes this.) */

  SortColumnRealSparseMatrix(MPS->A.Value, MPS->A.Row,
			     MPS->A.pBeginRow, MPS->A.pEndRow,
			     &(MPS->NumRows), &(MPS->NumCols));

  GetTime(&userTime1, &systemTime1);
  (*Time) = userTime1 - userTime0;

  fclose(fp);
  return MPS;
}

/* This routine reads in an MPS file.  Joe Czyzyk, May 1994. */

int ReadMPSFile(corefile, MPS, Inputs)
  FILE           *corefile;
  MPStype        *MPS;
  Parameters     *Inputs;
{

  int             col, row, oldrow, entry;
  char            line[200], name1[15], name2[15], name3[15], code[5],
    OldColumnName[15];
  int             return_getline;
  double          val1, val2;
  int             ColNumber, EntNumber, RowNumber;
  int             MPSQuickCheckBoundInfeasibility(),
    AddMPSRow(), AddMPSColumn(), AddMPSCoefficient(), 
    AddMPSRHS(), AddMPSBound(), AddMPSRange();

  /* read header line */

  if ((return_getline = GetLine(line, corefile, 200)) == HEADERLINE) {
    ParseHeaderLine(line, name1, name2);
    if (strncmp(name1, "NAME", 4) != 0) {
      printf("Error reading MPS file.\n");
      printf("Expected to find entry NAME, but found %s instead.\n", name1);
      return BAD_INPUT;
    }
    AddMPSProblemName(MPS, name2);
  } else {
    printf("Read error or possibly missing NAME line in MPS file.\n");
    return BAD_INPUT;
  }

  /*****************************************************************/

  /* read ROWS line */

  if ((return_getline = GetLine(line, corefile, 200)) == HEADERLINE) {
    ParseHeaderLine(line, name1, name2);
    if (strncmp(name1, "ROWS", 4) != 0) {
      printf("Missing ROWS line in MPS file.\n");
      return BAD_INPUT;
    } else if (DEBUG > 0)
      printf("Found ROWS line in MPS file.\n");
  } else {
    printf("Read error or possibly missing NAME line in MPS file.\n");
    return BAD_INPUT;
  }
  /*****************************************************************/

  while ((return_getline = GetLine(line, corefile, 200)) == DATALINE) {
    ParseDataLine(line, code, name1, name2, &val1, name3, &val2);

    RowNumber = AddMPSRow(MPS, name1, code);
    if(RowNumber == BAD_INPUT) {
      printf("Bad Row Number in MPS file\n"); 
      return BAD_INPUT; /* error return */
    }
  }				/* end while */

  if(return_getline == READERROR) {
    printf("Read error in MPS file.\n");
    return BAD_INPUT;
  }

  /* Have we read an objective row by this point? If not, error! */
  if(NoObjectiveRowName(MPS)) {
    printf("No objective row found in ROWS section - error!\n");
    return BAD_INPUT;
  }

  /*****************************************************************/

  /* Process COLUMNS section */

  ParseHeaderLine(line, name1, name2);

  if (strncmp(name1, "COLUMNS", 7) != 0) {
    printf("Missing COLUMNS row in MPS file.\n");
    return BAD_INPUT;
  }
  /*****************************************************************/

  strcpy(OldColumnName, "");

  while ((return_getline = GetLine(line, corefile, 200)) == DATALINE) {

    ParseDataLine(line, code, name1, name2, &val1, name3, &val2);

    if (strncmp(OldColumnName, name1, 8) != 0) {
      /* found new column */
      strcpy(OldColumnName, name1);
      ColNumber = AddMPSColumn(MPS, name1);
      if(ColNumber == BAD_INPUT) {
        printf("Bad Column Number in MPS file\n");
        return BAD_INPUT; 
      }
    }
    EntNumber = AddMPSCoefficient(MPS, name1, name2, val1, Inputs);
    if(EntNumber == BAD_INPUT) {
      printf("Bad Coefficient in MPS file\n");
      return BAD_INPUT; 
    }

    if (val2 != 0.0) {
      EntNumber = AddMPSCoefficient(MPS, name1, name3, val2, Inputs);
      if(EntNumber == BAD_INPUT) {
        printf("Bad Coefficient in MPS file\n");
        return BAD_INPUT; 
      }
    }
  }				/* end while */

  if(return_getline == READERROR) {
    printf("Read error in MPS file.\n");
    return BAD_INPUT;
  }

#if DEBUG > 1
  printf("Found %d columns in the constraint matrix.\n", ColNumber + 1);
  printf("Found %d entries in the constraint matrix.\n", EntNumber + 1);
#endif

  /*****************************************************************/

  /* Read in RHS section */

  ParseHeaderLine(line, name1, name2);

  if (strncmp(name1, "RHS", 3) != 0) {
    printf("Missing RHS statement in MPS file.\n");
    return BAD_INPUT;
  }
  while ((return_getline = GetLine(line, corefile, 200)) == DATALINE) {
    
    ParseDataLine(line, code, name1, name2, &val1, name3, &val2);

    if (NoRHSName(MPS))
      AddMPSRHSName(MPS, name1);

    if (strcmp(MPS->RHSName, name1) == 0) {
      /* only consider entries with the same RHS name */

      if(AddMPSRHS(MPS, name2, val1) == BAD_INPUT) {
        printf("Bad RHS element in MPS file\n");
        return BAD_INPUT; 
      }
      if (val2 != 0.0)
	if(AddMPSRHS(MPS, name3, val2) == BAD_INPUT) {
          printf("Bad RHS element in MPS file\n");
          return BAD_INPUT; 
	}
    }
  }				/* end while GetLine == DATALINE */

  if(return_getline == READERROR) {
    printf("Read error in MPS file.\n");
    return BAD_INPUT;
  }

  /*****************************************************************/

  ParseHeaderLine(line, name1, name2);

  if (strncmp(name1, "RANGES", 6) == 0) {
#if (DEBUG > 0)
    printf("Found RANGES: %s\n", name1);
#endif
    while ((return_getline = GetLine(line, corefile, 200)) == DATALINE) {

      ParseDataLine(line, code, name1, name2, &val1, name3, &val2);

      if (NoRangeName(MPS))
	AddMPSRangeName(MPS, name1);

      if (strcmp(MPS->RangeName, name1) == 0) {
	/* only consider ranges with matching range name */
	
	if(AddMPSRange(MPS, name2, val1) == BAD_INPUT) {
          printf("Bad Range in MPS file\n");
          return BAD_INPUT; 
	}
	if (val2 != 0.0)
	  if(AddMPSRange(MPS, name3, val2) == BAD_INPUT) {
            printf("Bad Range in MPS file\n");
            return BAD_INPUT; 
	  }
      }
    }

    if(return_getline == READERROR) {
      printf("Read error in MPS file.\n");
      return BAD_INPUT;
    }
    
  }
  /*****************************************************************/

  ParseHeaderLine(line, name1, name2);

  if (strncmp(name1, "BOUNDS", 6) == 0) {
#if (DEBUG > 0)
    printf("Found BOUNDS: %s\n", name1);
#endif

    while ((return_getline = GetLine(line, corefile, 200)) == DATALINE) {
      
      ParseDataLine(line, code, name1, name2, &val1, name3, &val2);

      if (NoBoundName(MPS))
	AddMPSBoundName(MPS, name1);

      if (strcmp(MPS->BoundName, name1) == 0) {
	/* only consider bounds with same Bound name */

	if(AddMPSBound(MPS, code, name2, val1) == BAD_INPUT) {
          printf("Bad bound in MPS file\n");
          return BAD_INPUT; 
	}
	if (val2 != 0.0)
	  if(AddMPSBound(MPS, code, name3, val2) == BAD_INPUT) {
            printf("Bad bound in MPS file\n");
            return BAD_INPUT; 
	  }
      } else
	printf("Ignoring bounds with name '%s'.\n", name1);
    }

    if(return_getline == READERROR) {
      printf("Read error in MPS file.\n");
      return BAD_INPUT;
    }

  }
  /*****************************************************************/

  ParseHeaderLine(line, name1, name2);

  if (strncmp(name1, "ENDATA", 6) == 0) {
    if (DEBUG > 0)
      printf("Found ENDATA\nSuccessful completion of MPS file.\n");
  } else {
    printf("Missing ENDATA statement in MPS file.\n");
    return BAD_INPUT;
  }

  /*****************************************************************/

  /* Resize MPS data structure to actual problem size */

  ResizeMPSMatrix(MPS, MPS->NumRows, MPS->NumCols, MPS->NumEnts);

  /*****************************************************************/

  /* Make a quick pass through the information to check for obvious
   * infeasibilities. */

  if(MPSQuickCheckBoundInfeasibility(MPS) == BAD_INPUT) return BAD_INPUT; 

  /*****************************************************************/

#if (DEBUG > 1)
  PrintHashTableStats(MPS->RowTable);
  PrintHashTableStats(MPS->ColTable);
#endif
#if (DEBUG > 10)
  printf("Print Row Table:\n");
  PrintHashTable(MPS->RowTable);

  printf("Print Col Table:\n");
  PrintHashTable(MPS->ColTable);
#endif

  /* Free hash table from MPS data structure */

  /* do not do this if you permit changes to the MPS data structure */

  /* DeleteHashTable(MPS->RowTable); DeleteHashTable(MPS->ColTable); */

  return 0;
}

/*****************************************************************/

int ResizeMPSMatrix(MPS, rows, cols, ents)
  MPStype        *MPS;
  int             rows, cols, ents;
{
  int             oldrows, oldcols, oldents, i;

  if (cols != MPS->ColSize) {
    MPS->A.pBeginRow = (int *) Realloc(MPS->A.pBeginRow, cols * sizeof(int),
				       "MPS->A.pBeginRow");
    MPS->A.pEndRow = (int *) Realloc(MPS->A.pEndRow, cols * sizeof(int),
				     "MPS->A.pEndRow");
    MPS->c = (double *) Realloc(MPS->c, cols * sizeof(double), "MPS->c");
    MPS->BoundType = (int *) Realloc(MPS->BoundType, cols * sizeof(int),
				     "MPS->BoundType");
    MPS->UpBound = (double *) Realloc(MPS->UpBound, cols * sizeof(double),
				      "MPS->UpBound");
    MPS->LowBound = (double *) Realloc(MPS->LowBound, cols * sizeof(double),
				       "MPS->LowBound");
    MPS->ColNames = (char **) Realloc(MPS->ColNames, cols * sizeof(char *),
				      "MPS->ColNames");
    oldcols = MPS->ColSize;
    for (i = oldcols; i < cols; i++) {
      MPS->A.pBeginRow[i] = 0;
      MPS->A.pEndRow[i] = 0;
      MPS->c[i] = 0.0;
      MPS->UpBound[i] = 0.0;
      MPS->LowBound[i] = 0.0;
      MPS->BoundType[i] = 0;
    }
    MPS->ColSize = cols;
  }

  if (ents != MPS->EntSize) {
    MPS->A.Row = (int *) Realloc(MPS->A.Row, ents * sizeof(int),
				 "MPS->A.Row");
    MPS->A.Value = (double *)
      Realloc(MPS->A.Value, ents * sizeof(double), "MPS->A.Value");

    oldents = MPS->EntSize;
    for (i = oldents; i < ents; i++) {
      MPS->A.Row[i] = 0;
      MPS->A.Value[i] = 0.0;
    }
    MPS->EntSize = ents;
  }

  if (rows != MPS->RowSize) {
    MPS->b = (double *) Realloc(MPS->b, rows * sizeof(double), "MPS->b");
    
    MPS->Ranges = (double *) Realloc(MPS->Ranges, rows * sizeof(double),
				     "MPS->Ranges");
    
    MPS->RowType = (char *) Realloc(MPS->RowType, rows * sizeof(char),
				    "MPS->RowType");
    
    MPS->RowNames = (char **) Realloc(MPS->RowNames, rows * sizeof(char *),
				      "MPS->RowNames");
    oldrows = MPS->RowSize;
    for (i = oldrows; i < rows; i++) {
      MPS->b[i] = 0.0;
      MPS->Ranges[i] = 0.0;
    }
    MPS->RowSize = rows;
  }
  return 0;
}

/*---------------------------------------------------*/

int DeleteMPS(MPS)
  MPStype        *MPS;
{
  int             i;

  Free((char *) MPS->A.pBeginRow);
  Free((char *) MPS->A.pEndRow);
  Free((char *) MPS->A.Row);
  Free((char *) MPS->A.Value);
  
  Free((char *) MPS->b);
  Free((char *) MPS->c);

  Free((char *) MPS->BoundType);
  Free((char *) MPS->UpBound);
  Free((char *) MPS->LowBound);

  Free((char *) MPS->Ranges);
  Free((char *) MPS->RowType);

  Free((char *) MPS->ProblemName);
  Free((char *) MPS->ObjectiveName);
  Free((char *) MPS->RangeName);
  Free((char *) MPS->BoundName);
  Free((char *) MPS->RHSName);

  for (i = 0; i < MPS->NumRows; i++)
    Free((char *) MPS->RowNames[i]);
  for (i = 0; i < MPS->NumCols; i++)
    Free((char *) MPS->ColNames[i]);

  Free((char *) MPS->RowNames);
  Free((char *) MPS->ColNames);

  DeleteHashTable(MPS->RowTable);
  DeleteHashTable(MPS->ColTable);
  Free((char *) MPS);

  return 0;

}

/*---------------------------------------------------*/

int AddMPSProblemName(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  MPS->ProblemName = StrDup(name, "MPS->ProblemName");
  return 0;
}

/*---------------------------------------------------*/

int AddMPSObjectiveRowName(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  MPS->ObjectiveName = StrDup(name, "MPS->ObjectiveName");
  return 0;
}

/*---------------------------------------------------*/

int             NoObjectiveRowName(MPS)
  MPStype        *MPS;
{
  if (MPS->ObjectiveName == NULL)
    return 1;			/* no objective name stored */
  else
    return 0;
}

/*---------------------------------------------------*/

int             AddMPSRow(MPS, name, type)
  MPStype        *MPS;
  char           *name, *type;
{
  int             oldrows, oldcols, oldents, newrows;
  int             RowNumber;
  char            RowType;

  oldrows = MPS->RowSize;
  oldcols = MPS->ColSize;
  oldents = MPS->EntSize;

  if (!isspace(type[0]))
    RowType = type[0];
  else if (!isspace(type[1]))
    RowType = type[1];
  else {
    printf("The row '%s' has a missing type in the ROWS section.\n", name);
    return BAD_INPUT;    
  }

  RowType = toupper(RowType);

  if (RowType == 'N')
    if (NoObjectiveRowName(MPS)) {
      AddMPSObjectiveRowName(MPS, name);
#if (PARSE_DEBUG > 0)
      printf("Found Objective Row: '%s'\n", name);
#endif
      return -1;
    }
  MPS->NumRows++;
  RowNumber = MPS->NumRows - 1;

  if (RowNumber >= oldrows) {
    newrows = FACTOR * oldrows;
    if (newrows < ROW_INCREMENT)
      newrows = ROW_INCREMENT;
    ResizeMPSMatrix(MPS, oldrows + newrows, oldcols, oldents);
  }
  if (Insert(MPS->RowTable, name, RowNumber) == 1) {
    printf("Duplicate row name '%s' in MPS file.\n", name);
    return BAD_INPUT;    
  }
  MPS->RowNames[RowNumber] = StrDup(name, "MPS->RowNames[]");
  MPS->RowType[RowNumber] = RowType;

  if (RowType != 'G' && RowType != 'L' && RowType != 'E' && RowType != 'N')
    printf("Unknown row type in MPS file: row = '%s' type = %c\n",
	   name, RowType);

  return (RowNumber);

}				/* end AddRow */

/*---------------------------------------------------*/

int             AddMPSColumn(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  int             oldrows, oldcols, oldents, newcols, ColNumber;

  oldrows = MPS->RowSize;
  oldcols = MPS->ColSize;
  oldents = MPS->EntSize;

  MPS->NumCols++;
  ColNumber = MPS->NumCols - 1;

  if (ColNumber >= oldcols) {
    newcols = FACTOR * oldcols;
    if (newcols < COL_INCREMENT)
      newcols = COL_INCREMENT;
    ResizeMPSMatrix(MPS, oldrows, oldcols + newcols, oldents);
  }
  if (Insert(MPS->ColTable, name, ColNumber) == 1) {
    printf("Duplicate column name '%s' in MPS file.\n", name);
    return BAD_INPUT; 
  }
  MPS->ColNames[ColNumber] = StrDup(name, "MPS->ColNames[]");

  MPS->A.pBeginRow[ColNumber] = MPS->NumEnts + 1;
  MPS->A.pEndRow[ColNumber]   = MPS->NumEnts;

  return (ColNumber);

}				/* end AddColumn */

/*---------------------------------------------------*/

int             AddMPSCoefficient(MPS, colname, rowname, value, Inputs)
  MPStype        *MPS;
  char           *colname, *rowname;
  double          value;
  Parameters     *Inputs;
{
  int             oldrows, oldcols, oldents, newents, RowNumber, ColNumber,
                  EntNumber;

  oldrows = MPS->RowSize;
  oldcols = MPS->ColSize;
  oldents = MPS->EntSize;

  ColNumber = GetColNumber(MPS, colname);
  if (!NoObjectiveRowName(MPS) && (strcmp(rowname, MPS->ObjectiveName) == 0)) {
    /* objective coefficient */
    if (Inputs->Minimize)
      MPS->c[ColNumber] = value;
    else
      MPS->c[ColNumber] = -value;

    return (MPS->NumEnts - 1);
  }
  RowNumber = GetRowNumber(MPS, rowname);
  if (RowNumber == -1) {
    printf("Unable to find row name '%s' in the row name array.\n", rowname);
    return BAD_INPUT; 
  }
  MPS->NumEnts++;
  EntNumber = MPS->NumEnts - 1;

  if (EntNumber >= oldents) {
    newents = FACTOR * oldents;
    if (newents < ENT_INCREMENT)
      newents = ENT_INCREMENT;
    ResizeMPSMatrix(MPS, oldrows, oldcols, oldents + newents);
  }
  MPS->A.Row[EntNumber] = RowNumber + 1;	/* +1 for Mehrotra's system */
  MPS->A.Value[EntNumber] = value;
  MPS->A.pEndRow[ColNumber] = EntNumber + 1;

  return (EntNumber);

}				/* end AddCoefficient */

/*---------------------------------------------------*/

int AddMPSRHSName(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  MPS->RHSName = StrDup(name, "MPS->RHSName");
  return 0;
}

/*---------------------------------------------------*/

int             NoRHSName(MPS)
  MPStype        *MPS;
{
  if (MPS->RHSName == NULL)
    return 1;			/* no objective name stored */
  else
    return 0;
}

/*---------------------------------------------------*/

int AddMPSRHS(MPS, rowname, value)
  MPStype        *MPS;
  char           *rowname;
  double          value;
{
  int             RowNumber;

  if (strcmp(rowname, MPS->ObjectiveName) == 0) {
    /* objective row in right-hand side: constant term in objective  */
    MPS->cshift = value;
    return 0;
  }
  RowNumber = GetRowNumber(MPS, rowname);
  /* RowNumber = GetNumberFromName(MPS->RowNames, rowname, MPS->RowSize); */

  if (RowNumber == -1) {
    printf("Unable to find row '%s' in row array in RHS section.\n", rowname);
    return BAD_INPUT; 
  }
  MPS->b[RowNumber] = value;
  return 0;
}				/* end AddRHS */

/*---------------------------------------------------*/

int AddMPSRangeName(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  MPS->RangeName = StrDup(name, "MPS->RangeName");
  return 0;
}

/*---------------------------------------------------*/

int             NoRangeName(MPS)
  MPStype        *MPS;
{
  if (MPS->RangeName == NULL)
    return 1;			/* no objective name stored */
  else
    return 0;
}

/*---------------------------------------------------*/

int AddMPSRange(MPS, rowname, value)
  MPStype        *MPS;
  char           *rowname;
  double          value;
{
  int             RowNumber;

  RowNumber = GetRowNumber(MPS, rowname);
  /* RowNumber = GetNumberFromName(MPS->RowNames, rowname, MPS->RowSize); */

  if (RowNumber == -1) {
    printf("Unable to find row '%s' in row array in RANGES section.\n",
	   rowname);
    return BAD_INPUT; 
  }
  MPS->Ranges[RowNumber] = value;

#if DEBUG > 1
  printf("Range[%d] = %f\n", RowNumber, value);
#endif
  return 0;
}				/* end AddMPSRange */

/*---------------------------------------------------*/

int AddMPSBoundName(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  MPS->BoundName = StrDup(name, "MPS->BoundName");
  return 0;
}

/*---------------------------------------------------*/

int             NoBoundName(MPS)
  MPStype        *MPS;
{
  if (MPS->BoundName == NULL)
    return 1;			/* no objective name stored */
  else
    return 0;
}

/*---------------------------------------------------*/

int AddMPSBound(MPS, code, colname, value)
  MPStype        *MPS;
  char           *code, *colname;
  double          value;
{
  int             ColNumber;

  ColNumber = GetColNumber(MPS, colname);
  if (ColNumber == -1) {
    printf("Unable to find col '%s' in col array in BOUNDS section.\n",
	   colname);
    return BAD_INPUT; 
  }
  if (strncmp(code, "LO", 2) == 0) {	/* lower bound */

    switch (MPS->BoundType[ColNumber]) {
    case UPPER:
      MPS->BoundType[ColNumber] = UPPERLOWER;
      MPS->LowBound[ColNumber] = value;
      break;
    case LOWER:
    case UPPERLOWER:
      printf("Found duplicate lower bound for variable '%s'.\n", colname);
      printf("Using more stringent value for bound.\n");
      if (MPS->LowBound[ColNumber] < value)
	MPS->LowBound[ColNumber] = value;
      break;
    case FREE:
      printf("The variable '%s' was previously declared to be FRee,\n",
	     colname);
      printf("but a lower bound has been indicated -- using lower bound.\n");
      MPS->BoundType[ColNumber] = LOWER;
      MPS->LowBound[ColNumber] = value;
      break;
    case PINFTY:
      MPS->BoundType[ColNumber] = LOWER;
      MPS->LowBound[ColNumber] = value;
      break;
    case MINFTY:
      printf("The variable '%s' was previously declared to be MI,\n", colname);
      printf("but a lower bound has been indicated -- using lower bound.\n");
      MPS->BoundType[ColNumber] = LOWER;
      MPS->LowBound[ColNumber] = value;
      break;
    case FIX:
      printf("The variable '%s' was previously declared to be FiXed,\n", colname);
      printf("but a lower bound has been indicated -- using lower bound.\n");
      MPS->BoundType[ColNumber] = LOWER;
      MPS->LowBound[ColNumber] = value;
    }
  }
  if (strncmp(code, "UP", 2) == 0) {	/* upper bound */
    switch (MPS->BoundType[ColNumber]) {
    case LOWER:
      MPS->BoundType[ColNumber] = UPPERLOWER;
      MPS->UpBound[ColNumber] = value;
      break;
    case UPPER:
    case UPPERLOWER:
      printf("Found duplicate upper bound for variable '%s'\n", colname);
      if (MPS->UpBound[ColNumber] > value)
	MPS->UpBound[ColNumber] = value;
      break;
    case FREE:
      printf("The variable '%s' was previously declared to be FRee,\n",
	     colname);
      printf("but an upper bound has been indicated -- using upper bound.\n");
      MPS->BoundType[ColNumber] = UPPER;
      MPS->UpBound[ColNumber] = value;
      break;
    case PINFTY:
      MPS->BoundType[ColNumber] = UPPER;
      MPS->UpBound[ColNumber] = value;
      break;
    case MINFTY:
      printf("The variable '%s' was previously declared to be MInfty,\n", colname);
      printf("but an upper bound has been indicated -- using this upper bound.\n");
      MPS->BoundType[ColNumber] = MINFTY;
      MPS->UpBound[ColNumber] = value;
      break;
    case FIX:
      printf("The variable '%s' was previously declared to be FiXed,\n", colname);
      printf("but an upper bound has been indicated -- using upper bound.\n");
      MPS->BoundType[ColNumber] = UPPER;
      MPS->UpBound[ColNumber] = value;
    }				/* end switch */
  }
  if (strncmp(code, "FX", 2) == 0) {	/* fixed bound */
    switch (MPS->BoundType[ColNumber]) {
    case LOWER:
      printf("The variable '%s' previously had a lower bound\n", colname);
      printf("Now it is fixed.\n");
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
      break;
    case UPPER:
    case UPPERLOWER:
      printf("The variable '%s' previously had an upper bound\n", colname);
      printf("Now it is fixed.\n");
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
      break;
    case FREE:
      printf("The variable '%s' was previously declared to be FRee,\n",
	     colname);
      printf("Now it is fixed.\n");
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
      break;
    case PINFTY:
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
      break;
    case MINFTY:
      printf("The variable '%s' was previously declared to be MI,\n", colname);
      printf("Now it is fixed.\n");
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
      break;
    case FIX:
      printf("The variable '%s' was previously declared to be FiXed.\n", colname);
      printf("It will now be fixed at a new value.\n");
      MPS->BoundType[ColNumber] = FIX;
      MPS->LowBound[ColNumber] = value;
      MPS->UpBound[ColNumber] = value;
    }				/* end switch */
  }
  if (strncmp(code, "FR", 2) == 0) {	/* free variable */
    switch (MPS->BoundType[ColNumber]) {
    case LOWER:
    case UPPER:
    case UPPERLOWER:
    case MINFTY:
    case FIX:
      printf("The variable '%s' previously had a bound\n", colname);
      printf("Now it is FRee.\n");
    case PINFTY:
      MPS->BoundType[ColNumber] = FREE;
      break;
    case FREE:
      printf("The variable '%s' has been declared FRee more than once.\n",
	     colname);
    }				/* end switch */
  }
  if (strncmp(code, "MI", 2) == 0) {	/* minus infinity */
    switch (MPS->BoundType[ColNumber]) {
    case LOWER:
    case UPPERLOWER:
      printf("The variable '%s' previously had a lower bound\n", colname);
      printf("Now it is MInus infinity.\n");
      MPS->BoundType[ColNumber] = MINFTY;
      break;
    case UPPER:
      MPS->BoundType[ColNumber] = MINFTY;
      break;
    case FREE:
      printf("The variable '%s' was previously declared to be FRee,\n",
	     colname);
      printf("Now it is MInus infinity.\n");
      MPS->BoundType[ColNumber] = MINFTY;
      MPS->UpBound[ColNumber] = 0.0;
      break;
    case MINFTY:
    case PINFTY:
      MPS->BoundType[ColNumber] = MINFTY;
      break;
    case FIX:
      printf("The variable '%s' was previously declared to be FiXed.\n", colname);
      printf("It is now MInus infinity.\n");
      MPS->BoundType[ColNumber] = MINFTY;
      MPS->LowBound[ColNumber] = 0.0;
      MPS->UpBound[ColNumber] = 0.0;
    }				/* end switch */
  }
  if (strncmp(code, "PL", 2) == 0) {	/* plus infinity */
    switch (MPS->BoundType[ColNumber]) {
    case MINFTY:
      printf("The variable '%s' used to be MInus infinity.\n", colname);
      printf("Now it is PLus infinity.\n");
      MPS->BoundType[ColNumber] = PINFTY;
      break;
    case UPPERLOWER:
      printf("The variable '%s' previously had upper and lower\n", colname);
      printf("bounds. Now it is lower bound (due to PL).\n");
      MPS->BoundType[ColNumber] = LOWER;
      break;
    case UPPER:
      printf("The variable '%s' previously had an upper bound\n", colname);
      printf("Now it is PLus infinity.\n");
      MPS->BoundType[ColNumber] = PINFTY;
      break;
    case FREE:
      printf("The variable '%s' was previously declared to be FRee,\n",
	     colname);
      printf("Now it is PLus infinity.\n");
      MPS->BoundType[ColNumber] = PINFTY;
      MPS->UpBound[ColNumber] = 0.0;
      MPS->LowBound[ColNumber] = 0.0;
      break;
    case FIX:
      printf("The variable '%s' was previously declared to be FiXed.\n", colname);
      printf("It is now LOWER (due to PL).\n");
      MPS->BoundType[ColNumber] = LOWER;
      MPS->UpBound[ColNumber] = 0.0;
      break;
    case PINFTY:
    case LOWER:
      break;
    }				/* end switch */
  }
  return 0;
}				/* end AddMPSBound */

/*---------------------------------------------------*/

int MPSQuickCheckBoundInfeasibility(MPS)
  MPStype        *MPS;
{
  int             col;

  for (col = 0; col < MPS->NumCols; col++)
    if (MPS->BoundType[col] != 0)
      switch (MPS->BoundType[col]) {
      case UPPER:
	if (MPS->UpBound[col] < 0.0) {
	  printf("Bound infeasibility:\n");
	  printf("The upper bound %f for column '%s' is less than zero.\n",
		 MPS->UpBound[col], MPS->ColNames[col]);
	  return BAD_INPUT; 
	}
	break;

      case UPPERLOWER:
	if (MPS->UpBound[col] < MPS->LowBound[col]) {
	  printf("Bound infeasibility:\n");
	  printf("The upper bound %f for column '%s' is less than the\n",
		 MPS->UpBound[col], MPS->ColNames[col]);
	  printf("lower bound %f.\n", MPS->LowBound[col]);
	  return BAD_INPUT; 
	}
	break;

      }				/* end switch */
      return 0;
}

/*---------------------------------------------------*/

int PrintMPS(MPS)
  MPStype        *MPS;
{
  int             row, col, entry;

  printf("-----------------------------------------------\n");
  printf("ProblemName   = %s\n", MPS->ProblemName);
  printf("ObjectiveName = %s\n", MPS->ObjectiveName);

  printf("NumRows = %d  NumCols = %d  NumEnts = %d\n",
	 MPS->NumRows, MPS->NumCols, MPS->NumEnts);

  printf("\nA:\n");
  for (col = 0; col < MPS->NumCols; col++) {

    printf(" Col %d:\n", col + 1);

    for (entry = MPS->A.pBeginRow[col] - 1;
	 entry <= MPS->A.pEndRow[col] - 1; entry++)
      printf("  Row %d  = %f\n", MPS->A.Row[entry], MPS->A.Value[entry]);
  }

  printf("\nRHS Name = '%s'\n", MPS->RHSName);

  printf("\nb:\n");
  for (row = 0; row < MPS->NumRows; row++)
    printf(" Row %d = %f\n", row + 1, MPS->b[row]);

  printf("\nObjective Name = '%s'\n", MPS->ObjectiveName);

  printf("\nc:\n");
  for (col = 0; col < MPS->NumCols; col++)
    printf(" Col %d = %f\n", col + 1, MPS->c[col]);

  printf("\nBounds:\n");

  for (col = 0; col < MPS->NumCols; col++)
    if (MPS->BoundType[col] != 0)
      printf("BoundType[%d] = %d, lower = %f, upper = %f\n",
       col + 1, MPS->BoundType[col], MPS->LowBound[col], MPS->UpBound[col]);

  if (MPS->RangeName != NULL) {
    printf("\nRanges:\n");
    printf(" RangeName: %s\n", MPS->RangeName);
    
    for (row = 0; row < MPS->NumRows; row++)
      if (MPS->Ranges[row] != 0.0)
	printf(" Range[%d] = %f\n", row+1, MPS->Ranges[row]);
  }

  printf("-----------------------------------------------\n\n");

  return 0;
}				/* end PrintMPS */

/******************************************************************/
int 
PrintMPSMatrixToFile(MPS, filename)
  MPStype        *MPS;
  char           *filename;
{
  FILE           *fp;
  int             row, col, entry;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Unable to open output file '%s'.\n", filename);
  }
  /* print dimensions */

  fprintf(fp, "%4d %4d %4d\n",
	  MPS->NumRows, MPS->NumCols, MPS->NumEnts);
  for (row = 0; row < MPS->NumRows; row++)
    fprintf(fp, "%s\n", MPS->RowNames[row]);

  for (col = 0; col < MPS->NumCols; col++)
    fprintf(fp, "%s\n", MPS->ColNames[col]);
  for (col = 0; col < MPS->NumCols; col++) {
    for (entry = MPS->A.pBeginRow[col] - 1;
	 entry <= MPS->A.pEndRow[col] - 1; entry++)
      fprintf(fp, "%4d %4d %f\n",
	      MPS->A.Row[entry], col + 1, MPS->A.Value[entry]);
  }
  fprintf(fp, "\nRHS Name = %s\n", MPS->RHSName);

  fprintf(fp, "\nb:\n");
  for (row = 0; row < MPS->NumRows; row++)
    fprintf(fp, " Row %d = %f\n", row + 1, MPS->b[row]);

  fprintf(fp, "\nObjective Name = '%s'\n", MPS->ObjectiveName);

  fprintf(fp, "\nc:\n");
  for (col = 0; col < MPS->NumCols; col++)
    fprintf(fp, " Col %d = %f\n", col + 1, MPS->c[col]);

  fprintf(fp, "Bounds:\n");

  for (col = 0; col < MPS->NumCols; col++)
    if (MPS->BoundType[col] != 0)
      fprintf(fp, "BoundType[%d] = %d, lower = %f, upper = %f\n",
       col + 1, MPS->BoundType[col], MPS->LowBound[col], MPS->UpBound[col]);
  fclose(fp);
  return 0;
}

/******************************************************************/

GetLine(line, file, length)
  char           *line;
  FILE           *file;
  int             length;
{
  char *return_line;
  int             i;

  do {

    for (i = 0; i < length; i++)
      line[i] = ' ';
    return_line = fgets(line, length, file);	/* get line from file */

  } while (return_line != NULL && line[0] == '*'); /* discard comment lines */

  if(return_line == NULL) return READERROR;

  if(line[0] == ' ') return DATALINE;
  else return HEADERLINE;
}

/******************************************************************/

int ParseHeaderLine(line, entry1, entry2)
  char           *entry1, *entry2, *line;
{
  string_copy(entry1, &line[0], 14);	/* characters  1 - 14 to entry1 */
  string_copy(entry2, &line[14], 10);	/* characters 15 - 24 to entry2 */

  return 0;
}

/******************************************************************/

int ParseDataLine(line, code, name1, name2, val1, name3, val2)
  char           *code, *line, *name1, *name2, *name3;
  double         *val1, *val2;
{
  int             i, j;
  char            valstr[15], NoBlank[15];

  *val1 = 0.0;
  *val2 = 0.0;

  string_copy(code, &line[1], 2);	/* characters  2 -  3 to code  */
  string_copy(name1, &line[4], 8);	/* characters  5 - 12 to name1 */
  string_copy(name2, &line[14], 8);	/* characters 15 - 22 to name2 */
  string_copy(valstr, &line[24], 12);

  j = 0;
  for (i = 0; i < 15 && valstr[i] != '\0'; i++)
    if (!isspace(valstr[i]))
      NoBlank[j++] = valstr[i];
  NoBlank[j] = '\0';

  *val1 = atof(NoBlank);	/* characters 25 - 36 to val1  */

  string_copy(name3, &line[39], 8);	/* characters 40 - 47 to name3 */
  string_copy(valstr, &line[49], 12);

  j = 0;
  for (i = 0; i < 15 && valstr[i] != '\0'; i++)
    if (!isspace(valstr[i]))
      NoBlank[j++] = valstr[i];
  NoBlank[j] = '\0';
  *val2 = atof(NoBlank);	/* characters 50 - 61 to val2  */

#if PARSE_DEBUG > 0
  printf("code = |%s|\nname1 = |%s|\nname2 = |%s|\nval1 = %f\n",
	 code, name1, name2, *val1);
  printf("name3 = |%s|\nval2 = %f\n", name3, *val2);
#endif

  return 0;
}

/******************************************************************/
int 
string_copy(dest, string, max)
  char           *dest, *string;
  int             max;
{
  /* Copy string to dest converting non-printable characters to spaces.
   * Null-terminate dest at max+1 */

  int             i;

  for (i = 0; i < max; i++)
    if (isprint(string[i]))
      dest[i] = string[i];
    else
      dest[i] = ' ';

  dest[max] = '\0';
  return 0;
}

/******************************************************************/

/* perform search through array for matching name, return index */

int             GetNumberFromName(NameList, name, size)
  char          **NameList, *name;
  int             size;
{
  int             Found, index, i;

  Found = 0;
  i = 0;

  if (name == NULL) {
    printf("No match for NULL name passed to GetNumberFromName().\n");
    return -1;
  }
  while (!Found && i < size) {
    if (NameList[i] != NULL &&
	strcmp(name, NameList[i]) == 0) {	/* a match */
      Found = 1;
      index = i;
    }
    i++;
  }

  if (Found)
    return index;
  else {
    printf("Unable to find name: '%s' in array.\n", name);
    return -1;
  }
}				/* end GetNumberFromName */

/******************************************************************/

char           *GetNameFromNumber(NameList, index)
  char          **NameList;
  int             index;
{
  return NameList[index];
}


/******************************************************************/

int             GetRowNumber(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  return GetIndex(MPS->RowTable, name);
}

/******************************************************************/

int             GetColNumber(MPS, name)
  MPStype        *MPS;
  char           *name;
{
  return GetIndex(MPS->ColTable, name);
}
