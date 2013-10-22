/* smxpy1.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******     SMXPY1 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '1' SIGNIFIES NO LOOP UNROLLING, I.E., */
/*               LOOP-UNROLLING TO LEVEL 1. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int smxpy1_(m, n, y, apnt, a)
integer *m, *n;
doublereal *y;
integer *apnt;
doublereal *a;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal amult;
    static integer ii;


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
    --a;
    --apnt;
    --y;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ii = apnt[j + 1] - *m;
	amult = -a[ii];
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    y[i] += amult * a[ii];
	    ++ii;
/* L100: */
	}
/* L200: */
    }
    return 0;
} /* smxpy1_ */

