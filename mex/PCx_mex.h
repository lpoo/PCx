/************************************************************************

MATLAB - PCx interface

Michael Wagner
Old Dominion University

August 24, 2000

 ***********************************************************************/


#include <stdlib.h>
#include "mex.h"
#include "main.h"
#include "memory.h"
#include "pre.h"
#include "signal.h"

#define AbortMsgTxt(s) _AbortMsgTxt(s, __FILE__, __LINE__)
#define _AbortMsgTxt(s,f,l) fprintf(stderr,"Error at ==> %s:%d: %s\n",f,l,s);mexErrMsgTxt(s);

#define WarnMsgTxt(s) _WarnMsgTxt(s, __FILE__, __LINE__)
#define _WarnMsgTxt(s,f,l) fprintf(stderr,"Warning at ==> %s:%d: %s\n",f,l,s);

#define printf mexPrintf
#define malloc mxMalloc

double     *get_double_vec();
void        get_matrix();
int         PCx_main();
void        MexFunction();
Parameters *ReadOptions();
