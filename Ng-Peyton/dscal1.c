/* dscal1.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/*******     DSCAL1 .... SCALE A VECTOR                     ***************/
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE COMPUTES A <-- AX, WHERE A IS A */
/*               SCALAR AND X IS A VECTOR. */

/*     INPUT PARAMETERS - */
/*        N - LENGTH OF THE VECTOR X. */
/*        A - SCALAR MULIPLIER. */
/*        X - VECTOR TO BE SCALED. */

/*     OUTPUT PARAMETERS - */
/*        X - REPLACED BY THE SCALED VECTOR, AX. */

/* *********************************************************************** */

/* Subroutine */ int dscal1_(n, a, x)
integer *n;
doublereal *a, *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;


/* ***********************************************************************
 */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */

/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x[i] = *a * x[i];
/* L100: */
    }
    return 0;
} /* dscal1_ */

