/* mmpyi.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* *************     MMPYI  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       MATRIX X HAS ONLY 1 COLUMN. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN */
/*       RELIND(*)       -   RELATIVE INDICES. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int mmpyi_(m, q, xpnt, x, iy, y, relind)
integer *m, *q, *xpnt;
doublereal *x;
integer *iy;
doublereal *y;
integer *relind;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer isub;
    static doublereal a;
    static integer i, k, ylast, col;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --relind;
    --y;
    --iy;
    --x;
    --xpnt;

    /* Function Body */
    i__1 = *q;
    for (k = 1; k <= i__1; ++k) {
	col = xpnt[k];
	ylast = iy[col + 1] - 1;
	a = -x[k];
/* DIR$   IVDEP */
	i__2 = *m;
	for (i = k; i <= i__2; ++i) {
	    isub = xpnt[i];
	    isub = ylast - relind[isub];
	    y[isub] += a * x[i];
/* L100: */
	}
/* L200: */
    }
    return 0;

} /* mmpyi_ */

