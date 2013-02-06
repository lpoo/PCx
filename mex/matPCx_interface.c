
/************************************************************************

MATLAB - PCx interface

Michael Wagner
Old Dominion University

August 2000

 ***********************************************************************/

#include "PCx_mex.h"

/* generic service functions, adapted from Nathan Edwards' MATLAB-CPLEX
   interface */

double *get_double_vec(const mxArray* prhs[], int i)
{
  double *values;
  int buflen;
  int j;
  double *buf;

  if (mxIsEmpty(prhs[i])) 
    return(NULL);
  else
    {
      if (! ( mxIsDouble(prhs[i]) &&
	      ( mxGetM(prhs[i]) == 1 ||
		mxGetN(prhs[i]) == 1 ) ) ) 
	{
	  char msgbuf[1024];
	  sprintf(msgbuf,"Arg %d: Expected a vector of doubles",i+1);
	  AbortMsgTxt(msgbuf);
	}
      buflen = mxGetM(prhs[i]) * mxGetN(prhs[i]);
      buf    = (double*) mxCalloc(buflen,sizeof(double));
      values = mxGetPr(prhs[i]);
      memcpy(buf,values,buflen*sizeof(double));

      return buf;
    } 
}

void get_matrix(const mxArray*prhs[], int i,  int** beg, int** ind, double**val)
{
  int n;     /* number of columns */
  int nnz;   /* number of non-zeros */
  int j;
  int *matbeg, *matind;
  double *matval;

  if (! mxIsSparse(prhs[i])) 
    {
      char msgbuf[1024];
      sprintf(msgbuf,"Arg %d: Expected a sparse matrix",i+1);
      AbortMsgTxt(msgbuf);
    }
  
  n      = mxGetN(prhs[i]);
  nnz    = mxGetNzmax(prhs[i]);
  *beg   = (int *)mxCalloc(n+1,sizeof(int));
  matbeg = mxGetJc(prhs[i]);
  memcpy(*beg,matbeg,(n+1)*sizeof(int));
  *ind   = (int *)mxCalloc(nnz,sizeof(int));
  matind = mxGetIr(prhs[i]);
  memcpy(*ind,matind,nnz*sizeof(int));
  
  *val   = mxCalloc(nnz,sizeof(double));
  matval = mxGetPr(prhs[i]);
  memcpy(*val,matval,nnz*sizeof(double));
}

void return_double_vec(mxArray* plhs[], int i, double* value, int len)
{
/*   fprintf(stderr,"value: %X\n",value); */
/*   fprintf(stderr,"reallen: %d\n",reallen); */
/*   fprintf(stderr,"len: %d\n",len); */
  if (value != NULL) 
    {
      plhs[i] = mxCreateDoubleMatrix(len, 1, mxREAL);
      memcpy(mxGetPr(plhs[i]), value, len*sizeof(double));
    } 
  else 
    plhs[i] = mxCreateDoubleMatrix(0,1,mxREAL);
}

void signal_handler(int sig)
{
  switch (sig) {
  case SIGFPE:
    AbortMsgTxt("Signal SIGFPE received and handled...");
    break;
  case SIGBUS:
    AbortMsgTxt("Signal SIGBUS received and handled...");
    break;  
  case SIGSEGV:
    AbortMsgTxt("Signal SIGSEGV received and handled...");
    break;  
  default:
    /* do nothing */;
  }
}
