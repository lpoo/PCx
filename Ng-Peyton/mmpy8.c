/* mmpy8.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  May 26, 1995 */
/*   Authors:        Esmond G. Ng, Barry W. Peyton, and Guodong Zhang */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     MMPY8  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       LOOP UNROLLING: LEVEL 8 UPDATING TWO COLUMNS AT A TIME */

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

/* Subroutine */ int mmpy8_(m, n, q, xpnt, x, y, ldy)
integer *m, *n, *q, *xpnt;
doublereal *x, *y;
integer *ldy;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer leny, i, j, k, iybeg;
    static doublereal a1, a2, a3, a4, a5, a6, a7, a8;
    static integer i1, i2, i3, i4, i5, i6, i7, i8;
    static doublereal a9, b1, b2, b3, b4, b5, b6, b7, b8, y1, y2;
    static integer iybeg1, iybeg2;
    static doublereal a10, a11, a12, a13, a14, a15, a16;
    static integer mm, qq;
    extern /* Subroutine */ int smxpy8_();


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

/*       ---------------------------------------------------- */
/*       COMPUTE EACH DIAGONAL ENTRY OF THE ODD COLUMNS OF Y. */
/*       ---------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    --x;
    --xpnt;

    /* Function Body */
    mm = *m;
    qq = min(*m,*q);
    iybeg = 1;
    leny = *ldy - 1;
    i__1 = qq - 1;
    for (j = 1; j <= i__1; j += 2) {
/* DIR$   IVDEP */
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    i1 = xpnt[i + 1] - mm;
	    a1 = x[i1];
	    y[iybeg] -= a1 * a1;
/* L100: */
	}
	iybeg = iybeg + (leny << 1) + 1;
	leny += -2;
	mm += -2;
/* L200: */
    }

/*       ------------------------------------------------------- */
/*       UPDATE TWO COLUMNS OF Y AT A TIME,  EXCEPT THE DIAGONAL */
/*       ELEMENT. */
/*       NOTE: THE DIAGONAL ELEMENT OF THE ODD COLUMN HAS */
/*             BEEN COMPUTED, SO WE COMPUTE THE SAME NUMBER OF */
/*             ELEMENTS FOR THE TWO COLUMNS. */
/*       ------------------------------------------------------- */

    mm = *m;
    iybeg = 1;
    leny = *ldy - 1;

    i__1 = qq - 1;
    for (j = 1; j <= i__1; j += 2) {

	iybeg1 = iybeg;
	iybeg2 = iybeg + leny;

	i__2 = *n - 7;
	for (k = 1; k <= i__2; k += 8) {

/*               ----------------------------------- */
/*               EIGHT COLUMNS UPDATING TWO COLUMNS. */
/*               ----------------------------------- */

	    i1 = xpnt[k + 1] - mm;
	    i2 = xpnt[k + 2] - mm;
	    i3 = xpnt[k + 3] - mm;
	    i4 = xpnt[k + 4] - mm;
	    i5 = xpnt[k + 5] - mm;
	    i6 = xpnt[k + 6] - mm;
	    i7 = xpnt[k + 7] - mm;
	    i8 = xpnt[k + 8] - mm;
	    a1 = x[i1];
	    a2 = x[i2];
	    a3 = x[i3];
	    a4 = x[i4];
	    a5 = x[i5];
	    a6 = x[i6];
	    a7 = x[i7];
	    a8 = x[i8];
	    a9 = x[i1 + 1];
	    a10 = x[i2 + 1];
	    a11 = x[i3 + 1];
	    a12 = x[i4 + 1];
	    a13 = x[i5 + 1];
	    a14 = x[i6 + 1];
	    a15 = x[i7 + 1];
	    a16 = x[i8 + 1];

	    y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11 - 
		    a4 * a12 - a5 * a13 - a6 * a14 - a7 * a15 - a8 * a16;

	    y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11 - 
		    a12 * a12 - a13 * a13 - a14 * a14 - a15 * a15 - a16 * a16;


	    i__3 = mm - 1;
	    for (i = 2; i <= i__3; ++i) {
		y1 = y[iybeg1 + i];
		b1 = x[i1 + i];
		y1 -= b1 * a1;
		y2 = y[iybeg2 + i];
		b2 = x[i2 + i];
		y2 -= b1 * a9;
		y1 -= b2 * a2;
		b3 = x[i3 + i];
		y2 -= b2 * a10;
		y1 -= b3 * a3;
		b4 = x[i4 + i];
		y2 -= b3 * a11;
		y1 -= b4 * a4;
		b5 = x[i5 + i];
		y2 -= b4 * a12;
		y1 -= b5 * a5;
		b6 = x[i6 + i];
		y2 -= b5 * a13;
		y1 -= b6 * a6;
		b7 = x[i7 + i];
		y2 -= b6 * a14;
		y1 -= b7 * a7;
		b8 = x[i8 + i];
		y2 -= b7 * a15;
		y1 -= b8 * a8;
		y[iybeg1 + i] = y1;
		y2 -= b8 * a16;
		y[iybeg2 + i] = y2;
/* L300: */
	    }

/* L400: */
	}

/*           ----------------------------- */
/*           BOUNDARY CODE FOR THE K LOOP. */
/*           ----------------------------- */

	switch ((int)(*n - k + 2)) {
	    case 1:  goto L2000;
	    case 2:  goto L1700;
	    case 3:  goto L1500;
	    case 4:  goto L1300;
	    case 5:  goto L1100;
	    case 6:  goto L900;
	    case 7:  goto L700;
	    case 8:  goto L500;
	}

L500:

/*               ----------------------------------- */
/*               SEVEN COLUMNS UPDATING TWO COLUMNS. */
/*               ----------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	i3 = xpnt[k + 3] - mm;
	i4 = xpnt[k + 4] - mm;
	i5 = xpnt[k + 5] - mm;
	i6 = xpnt[k + 6] - mm;
	i7 = xpnt[k + 7] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a3 = x[i3];
	a4 = x[i4];
	a5 = x[i5];
	a6 = x[i6];
	a7 = x[i7];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];
	a11 = x[i3 + 1];
	a12 = x[i4 + 1];
	a13 = x[i5 + 1];
	a14 = x[i6 + 1];
	a15 = x[i7 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11 - a4 * 
		a12 - a5 * a13 - a6 * a14 - a7 * a15;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11 - a12 
		* a12 - a13 * a13 - a14 * a14 - a15 * a15;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    b3 = x[i3 + i];
	    y2 -= b2 * a10;
	    y1 -= b3 * a3;
	    b4 = x[i4 + i];
	    y2 -= b3 * a11;
	    y1 -= b4 * a4;
	    b5 = x[i5 + i];
	    y2 -= b4 * a12;
	    y1 -= b5 * a5;
	    b6 = x[i6 + i];
	    y2 -= b5 * a13;
	    y1 -= b6 * a6;
	    b7 = x[i7 + i];
	    y2 -= b6 * a14;
	    y1 -= b7 * a7;
	    y[iybeg1 + i] = y1;
	    y2 -= b7 * a15;
	    y[iybeg2 + i] = y2;
/* L600: */
	}

	goto L2000;

L700:

/*               --------------------------------- */
/*               SIX COLUMNS UPDATING TWO COLUMNS. */
/*               --------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	i3 = xpnt[k + 3] - mm;
	i4 = xpnt[k + 4] - mm;
	i5 = xpnt[k + 5] - mm;
	i6 = xpnt[k + 6] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a3 = x[i3];
	a4 = x[i4];
	a5 = x[i5];
	a6 = x[i6];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];
	a11 = x[i3 + 1];
	a12 = x[i4 + 1];
	a13 = x[i5 + 1];
	a14 = x[i6 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11 - a4 * 
		a12 - a5 * a13 - a6 * a14;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11 - a12 
		* a12 - a13 * a13 - a14 * a14;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    b3 = x[i3 + i];
	    y2 -= b2 * a10;
	    y1 -= b3 * a3;
	    b4 = x[i4 + i];
	    y2 -= b3 * a11;
	    y1 -= b4 * a4;
	    b5 = x[i5 + i];
	    y2 -= b4 * a12;
	    y1 -= b5 * a5;
	    b6 = x[i6 + i];
	    y2 -= b5 * a13;
	    y1 -= b6 * a6;
	    y[iybeg1 + i] = y1;
	    y2 -= b6 * a14;
	    y[iybeg2 + i] = y2;
/* L800: */
	}

	goto L2000;

L900:

/*               ---------------------------------- */
/*               FIVE COLUMNS UPDATING TWO COLUMNS. */
/*               ---------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	i3 = xpnt[k + 3] - mm;
	i4 = xpnt[k + 4] - mm;
	i5 = xpnt[k + 5] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a3 = x[i3];
	a4 = x[i4];
	a5 = x[i5];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];
	a11 = x[i3 + 1];
	a12 = x[i4 + 1];
	a13 = x[i5 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11 - a4 * 
		a12 - a5 * a13;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11 - a12 
		* a12 - a13 * a13;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    b3 = x[i3 + i];
	    y2 -= b2 * a10;
	    y1 -= b3 * a3;
	    b4 = x[i4 + i];
	    y2 -= b3 * a11;
	    y1 -= b4 * a4;
	    b5 = x[i5 + i];
	    y2 -= b4 * a12;
	    y1 -= b5 * a5;
	    y[iybeg1 + i] = y1;
	    y2 -= b5 * a13;
	    y[iybeg2 + i] = y2;
/* L1000: */
	}

	goto L2000;

L1100:

/*               ---------------------------------- */
/*               FOUR COLUMNS UPDATING TWO COLUMNS. */
/*               ---------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	i3 = xpnt[k + 3] - mm;
	i4 = xpnt[k + 4] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a3 = x[i3];
	a4 = x[i4];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];
	a11 = x[i3 + 1];
	a12 = x[i4 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11 - a4 * 
		a12;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11 - a12 
		* a12;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    b3 = x[i3 + i];
	    y2 -= b2 * a10;
	    y1 -= b3 * a3;
	    b4 = x[i4 + i];
	    y2 -= b3 * a11;
	    y1 -= b4 * a4;
	    y[iybeg1 + i] = y1;
	    y2 -= b4 * a12;
	    y[iybeg2 + i] = y2;
/* L1200: */
	}

	goto L2000;

L1300:

/*               ----------------------------------- */
/*               THREE COLUMNS UPDATING TWO COLUMNS. */
/*               ----------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	i3 = xpnt[k + 3] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a3 = x[i3];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];
	a11 = x[i3 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10 - a3 * a11;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10 - a11 * a11;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    b3 = x[i3 + i];
	    y2 -= b2 * a10;
	    y1 -= b3 * a3;
	    y[iybeg1 + i] = y1;
	    y2 -= b3 * a11;
	    y[iybeg2 + i] = y2;
/* L1400: */
	}

	goto L2000;

L1500:

/*               --------------------------------- */
/*               TWO COLUMNS UPDATING TWO COLUMNS. */
/*               --------------------------------- */

	i1 = xpnt[k + 1] - mm;
	i2 = xpnt[k + 2] - mm;
	a1 = x[i1];
	a2 = x[i2];
	a9 = x[i1 + 1];
	a10 = x[i2 + 1];

	y[iybeg1 + 1] = y[iybeg1 + 1] - a1 * a9 - a2 * a10;

	y[iybeg2 + 1] = y[iybeg2 + 1] - a9 * a9 - a10 * a10;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    b2 = x[i2 + i];
	    y2 -= b1 * a9;
	    y1 -= b2 * a2;
	    y[iybeg1 + i] = y1;
	    y2 -= b2 * a10;
	    y[iybeg2 + i] = y2;
/* L1600: */
	}

	goto L2000;

L1700:

/*               -------------------------------- */
/*               ONE COLUMN UPDATING TWO COLUMNS. */
/*               -------------------------------- */

	i1 = xpnt[k + 1] - mm;
	a1 = x[i1];
	a9 = x[i1 + 1];

	y[iybeg1 + 1] -= a1 * a9;

	y[iybeg2 + 1] -= a9 * a9;

	i__2 = mm - 1;
	for (i = 2; i <= i__2; ++i) {
	    y1 = y[iybeg1 + i];
	    b1 = x[i1 + i];
	    y1 -= b1 * a1;
	    y2 = y[iybeg2 + i];
	    y[iybeg1 + i] = y1;
	    y2 -= b1 * a9;
	    y[iybeg2 + i] = y2;
/* L1800: */
	}

	goto L2000;

/*           ----------------------------------------------- */
/*           PREPARE FOR NEXT PAIR OF COLUMNS TO BE UPDATED. */
/*           ----------------------------------------------- */

L2000:
	mm += -2;
	iybeg = iybeg2 + leny + 1;
	leny += -2;

/* L3000: */
    }

/*       ----------------------------------------------------- */
/*       BOUNDARY CODE FOR J LOOP:  EXECUTED WHENVER Q IS ODD. */
/*       ----------------------------------------------------- */

    if (j == qq) {
	smxpy8_(&mm, n, &y[iybeg], &xpnt[1], &x[1]);
    }

    return 0;
} /* mmpy8_ */

