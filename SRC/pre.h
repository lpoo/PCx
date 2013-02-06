/* structures for the presolving routines in presolve.c
 *
 * PCx beta-2.0  10/31/96.
 *
 * Author: Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */


#ifndef PREPROCESS_INCLUDE_FILE
#define PREPROCESS_INCLUDE_FILE

#define NULL_RECORD		0
#define DUPLICATE_COLUMN	7
#define DUPLICATE_ROW   	8
#define ZERO_ROW		1
#define ZERO_COLUMN		2
#define FIXED_VARIABLE		3
#define SINGLETON_ROW		4
#define FREE_COLUMN_SINGLETON	5
#define FORCED_ROW              6

#define NORMAL			0
#define FREE			1

#define TRUE       1
#define FALSE      0
#define TINY       1.0e-10

#define STILL_ACTIVE		0
#define INFEASIBLE		-1
#define UNBOUNDED 		-2
#define PREPROCESSING_ERROR	-3

#define RANGELOWER              -1
#define RANGEUPPER              1

#define STACK_SIZE_FACTOR      0.1

typedef struct  {
  int size; /* size of the sparse vector storage structures in
		this record */
  int *Index;
  double *Value;
} SparseVector;

typedef struct  {
  int Split1, Split2; /* indices of split variable */
                      /* Column Split2 is the deleted one */
  int Type; /* Type is NORMAL if the two columns are positive
	       multiples of each other; FREE otherwise */
  double multiplier; /* value for which (Split2) = multiplier*(Split1) */
  double costdiff;   /* c(Split2) - multiplier * c(Split1) */
                     /* Any discrepancy affects the recovery strategy */
} DuplicateCol;

typedef struct  {
  int Split1, Split2; /* indices of duplicate row; Split1 < Split2 */
  /* Row Split2 is deleted */
} DuplicateRow;

typedef struct  {
  int RowIndex; /* index of deleted row */
} ZeroRow;

typedef struct  {
  int ColumnIndex; /* index of deleted column */
  double Value; /* value of fixed variable */
  double cshiftChange; /* change to the constant term in the primal objective */
  double cElement; /* the deleted member of the cost vector */
} ZeroColumn;

typedef struct  {
  int ColumnIndex; /* index of fixed column */
  double Value; /* value of fixed variable */
  double cshiftChange; /* change to the constant term in the primal objective */
  double cElement; /* the deleted member of the cost vector */
  SparseVector *DeletedColumn; /* the deleted column */
} FixedVariable;

typedef struct  {
  int ColumnIndex; /* column index */
  int RowIndex;    /* row index */
  double Value; /* value of fixed variable */
  double cshiftChange; /* change to the constant term in the primal objective */
  double cElement; /* the deleted member of the cost vector */
  double bElement; /* the deleted element from the rhs */
  double aElement; /* the deleted element from the matrix */
  SparseVector *DeletedColumn; /* the deleted column */
} SingletonRow;

typedef struct  {
  int ColumnIndex; /* column index */
  int RowIndex;    /* row index */
  double cshiftChange; /* change to the constant term in the primal objective */
  double bElement; /* the deleted element from the rhs */
  double aElement; /* the deleted element from the matrix */
  double cElement; /* the deleted element from the cost vector */
  SparseVector *DeletedRow; /* the deleted row */
} FreeColumnSingleton;

typedef struct {
  int RowIndex;    /* row index */
  int LowerUpper;  /* indicates whether row was forced to upper or lower
	end of its range RANGELOWER=-1, RANGEUPPER=+1 */
  double cshiftChange; /* change to the constant term in the primal objective */
  SparseVector *DeletedRow; /* the deleted row */
  SparseVector *Values;     /* values of the fixed variables */
  SparseVector *cElements;  /* costs for these fixed variables */
  SparseVector **DeletedColumns; /* the deleted columns */
  double bElement; /* the deleted element from the rhs */
} ForcedRow;

typedef struct  {
  int ChangeType; /* which preprocessing operation? */

/* only one of the following pointers is active within
                              each SingleChange record */
  DuplicateCol 		 *pDuplicateCol;
  DuplicateRow 		 *pDuplicateRow;
  ZeroRow                *pZeroRow;
  ZeroColumn             *pZeroColumn;
  FixedVariable          *pFixedVariable;
  SingletonRow           *pSingletonRow;
  FreeColumnSingleton    *pFreeColumnSingleton;
  ForcedRow              *pForcedRow;
} SingleChange;

typedef struct ChangeStack {
  int Rows; /* rows in the matrix */
  int Columns; /* cols in the matrix */
  double cshift;
  int Size; /* how many SingleChanges can we put into this stack? */
  int Top;  /* how many are in the stack right now? */

  SingleChange **StackOfChanges; /* stack of single-change records */

  int *RowMask;    /* mask vector indicating deleted status of rows.
		If nonzero, indicates the pass on which they were deleted. */
  int *ColumnMask; /* mask vector indicating deleted status of columns */

  int *VarType; /*  Variable Types; may change during preprocessing due to
		    combining of columns, etc */

  int Passes;      /* number of passes */
  int ReducedRows;  /* Rows remaining in the matrix after reduction */
  int ReducedColumns; /* Columns remaining in the matrix after reduction */
  int SizeUnit;    /* the amount of space added to the stack on each increment */
} ChangeStack;

extern ChangeStack *NewRecord();
extern LPtype   *Convert_MPS_LP();  
extern solution *MPSsolution();
#endif
