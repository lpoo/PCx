/* global declaration of memory allocation routines in PCx()
 *
 * PCx beta-2.0  10/31/96.
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>

char    *Malloc();
char    *Calloc();
void     Free();
char    *StrDup();
double  *NewDouble();
double **NewDouble2();
int     *NewInt();
char    *NewChar();
char    *Realloc();


double  **NewDoublePtr();
int     **NewIntPtr();
char    **NewCharPtr();
void      OutOfSpace();
