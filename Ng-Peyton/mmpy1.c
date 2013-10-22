/* mmpy1.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* *************     MMPY1  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       LOOP UNROLLING: LEVEL 1 */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int mmpy1_(m, n, q, xpnt, x, y, ldy)
integer *m, *n, *q, *xpnt;
doublereal *x, *y;
integer *ldy;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer xcol, ycol, leny;
    static doublereal a1;
    static integer i1, mm, iy, iylast, iystop, iystrt;


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
    --y;
    --x;
    --xpnt;

    /* Function Body */
    mm = *m;
    iylast = 0;
    leny = *ldy;
/*       ------------------------------------ */
/*       TO COMPUTE EACH COLUMN YCOL OF Y ... */
/*       ------------------------------------ */
    i__1 = *q;
    for (ycol = 1; ycol <= i__1; ++ycol) {
	iystrt = iylast + 1;
	iystop = iystrt + mm - 1;
	iylast += leny;
/*           -------------------------------------------------- */
/*           ... PERFORM THE APPROPRATE MATRIX VECTOR MULTIPLY: */
/*               X * A(*,YCOL). */
/*           -------------------------------------------------- */
	i__2 = *n;
	for (xcol = 1; xcol <= i__2; ++xcol) {
	    i1 = xpnt[xcol + 1] - mm;
	    a1 = -x[i1];
	    i__3 = iystop;
	    for (iy = iystrt; iy <= i__3; ++iy) {
		y[iy] += a1 * x[i1];
		++i1;
/* L100: */
	    }
/* L200: */
	}
	--mm;
	--leny;
/* L300: */
    }

    return 0;
} /* mmpy1_ */

