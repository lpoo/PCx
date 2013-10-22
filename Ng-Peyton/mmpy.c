/* mmpy.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* **************     MMPY  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       SPLIT(*)        -   BLOCK PARTITIONING OF X. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY, */
/*                           WITH LEVEL N LOOP UNROLLING. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int mmpy_(m, n, q, split, xpnt, x, y, ldy, mmpyn)
integer *m, *n, *q, *split, *xpnt;
doublereal *x, *y;
integer *ldy;
/* Subroutine */ int (*mmpyn) ();
{
    static integer nn, fstcol, blk;


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
    --split;

    /* Function Body */
    blk = 1;
    fstcol = 1;
L100:
    if (fstcol <= *n) {
	nn = split[blk];
	(*mmpyn)(m, &nn, q, &xpnt[fstcol], &x[1], &y[1], ldy);
	fstcol += nn;
	++blk;
	goto L100;
    }
    return 0;

} /* mmpy_ */

