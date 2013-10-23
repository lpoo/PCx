typedef struct {
  int    NumSuperNodes;           /* NSUPER */
  int    *SuperPartitioning;      /* XSUPER */
  int    *mapColumnToSupernode;   /* SNODE  */
  int    *pSuperNodeCols;         /* XLINDX */
  int    *SuperNodeRows;          /* LINDX  */
  int    *pBeginRowL;             /* XLNZ   */
  double *L;                      /* LNZ    */
  int     NumCompressedCols;
} NgPeytonType;

int first = 0;
double *Fatorize_Tmp;
int * Fatorize_Work, *Fatorize_Split, Fatorize_WorkSize, Fatorize_TmpSize;
