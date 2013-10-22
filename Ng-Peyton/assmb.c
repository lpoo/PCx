/* assmb.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ************     ASSMB .... INDEXED ASSEMBLY OPERATION     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS AN INDEXED ASSEMBLY (I.E., SCATTER-ADD) */
/*       OPERATION, ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE */
/*       CHOLESKY CODES. */

/*   INPUT PARAMETERS: */
/*       M               -   NUMBER OF ROWS IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN Y. */
/*       Y               -   BLOCK UPDATE TO BE INCORPORATED INTO FACTOR */
/*                           STORAGE. */
/*       RELIND          -   RELATIVE INDICES FOR MAPPING THE UPDATES */
/*                           ONTO THE TARGET COLUMNS. */
/*       XLNZ            -   POINTERS TO THE START OF EACH COLUMN IN THE */
/*                           TARGET MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   CONTAINS COLUMNS MODIFIED BY THE UPDATE */
/*                           MATRIX. */

/* *********************************************************************** */

/* Subroutine */ int assmb_(m, q, y, relind, xlnz, lnz, lda)
integer *m, *q;
doublereal *y;
integer *relind, *xlnz;
doublereal *lnz;
integer *lda;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer icol, ycol, lbot1, yoff1, ir, il1, iy1;


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
    --lnz;
    --xlnz;
    --relind;
    --y;

    /* Function Body */
    yoff1 = 0;
    i__1 = *q;
    for (icol = 1; icol <= i__1; ++icol) {
	ycol = *lda - relind[icol];
	lbot1 = xlnz[ycol + 1] - 1;
/* DIR$ IVDEP */
	i__2 = *m;
	for (ir = icol; ir <= i__2; ++ir) {
	    il1 = lbot1 - relind[ir];
	    iy1 = yoff1 + ir;
	    lnz[il1] += y[iy1];
	    y[iy1] = 0.;
/* L100: */
	}
	yoff1 = iy1 - icol;
/* L200: */
    }

    return 0;
} /* assmb_ */

