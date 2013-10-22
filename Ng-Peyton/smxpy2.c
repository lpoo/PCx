/* smxpy2.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******     SMXPY2 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '2' SIGNIFIES LEVEL 2 LOOP UNROLLING. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int smxpy2_(m, n, y, apnt, a)
integer *m, *n;
doublereal *y;
integer *apnt;
doublereal *a;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal a1, a2;
    static integer i1, i2, remain;


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
    remain = *n % 2;

    switch ((int)(remain + 1)) {
	case 1:  goto L2000;
	case 2:  goto L100;
    }

L100:
    i1 = apnt[2] - *m;
    a1 = -a[i1];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] += a1 * a[i1];
	++i1;
/* L150: */
    }
    goto L2000;

L2000:
    i__1 = *n;
    for (j = remain + 1; j <= i__1; j += 2) {
	i1 = apnt[j + 1] - *m;
	i2 = apnt[j + 2] - *m;
	a1 = -a[i1];
	a2 = -a[i2];
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    y[i] = y[i] + a1 * a[i1] + a2 * a[i2];
	    ++i1;
	    ++i2;
/* L3000: */
	}
/* L4000: */
    }

    return 0;
} /* smxpy2_ */

