
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
