/* smxpy8.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******     SMXPY8 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '8' SIGNIFIES LEVEL 8 LOOP UNROLLING. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  APNT(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int smxpy8_(m, n, y, apnt, a)
integer *m, *n;
doublereal *y;
integer *apnt;
doublereal *a;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal a1, a2, a3, a4, a5, a6, a7, a8;
    static integer i1, i2, i3, i4, i5, i6, i7, i8, remain;


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
    remain = *n % 8;

    switch ((int)(remain + 1)) {
	case 1:  goto L2000;
	case 2:  goto L100;
	case 3:  goto L200;
	case 4:  goto L300;
	case 5:  goto L400;
	case 6:  goto L500;
	case 7:  goto L600;
	case 8:  goto L700;
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

L200:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2];
	++i1;
	++i2;
/* L250: */
    }
    goto L2000;

L300:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    i3 = apnt[4] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    a3 = -a[i3];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3];
	++i1;
	++i2;
	++i3;
/* L350: */
    }
    goto L2000;

L400:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    i3 = apnt[4] - *m;
    i4 = apnt[5] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    a3 = -a[i3];
    a4 = -a[i4];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[i4];
	++i1;
	++i2;
	++i3;
	++i4;
/* L450: */
    }
    goto L2000;

L500:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    i3 = apnt[4] - *m;
    i4 = apnt[5] - *m;
    i5 = apnt[6] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    a3 = -a[i3];
    a4 = -a[i4];
    a5 = -a[i5];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[i4] + a5 *
		 a[i5];
	++i1;
	++i2;
	++i3;
	++i4;
	++i5;
/* L550: */
    }
    goto L2000;

L600:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    i3 = apnt[4] - *m;
    i4 = apnt[5] - *m;
    i5 = apnt[6] - *m;
    i6 = apnt[7] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    a3 = -a[i3];
    a4 = -a[i4];
    a5 = -a[i5];
    a6 = -a[i6];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[i4] + a5 *
		 a[i5] + a6 * a[i6];
	++i1;
	++i2;
	++i3;
	++i4;
	++i5;
	++i6;
/* L650: */
    }
    goto L2000;

L700:
    i1 = apnt[2] - *m;
    i2 = apnt[3] - *m;
    i3 = apnt[4] - *m;
    i4 = apnt[5] - *m;
    i5 = apnt[6] - *m;
    i6 = apnt[7] - *m;
    i7 = apnt[8] - *m;
    a1 = -a[i1];
    a2 = -a[i2];
    a3 = -a[i3];
    a4 = -a[i4];
    a5 = -a[i5];
    a6 = -a[i6];
    a7 = -a[i7];
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[i4] + a5 *
		 a[i5] + a6 * a[i6] + a7 * a[i7];
	++i1;
	++i2;
	++i3;
	++i4;
	++i5;
	++i6;
	++i7;
/* L750: */
    }
    goto L2000;

L2000:
    i__1 = *n;
    for (j = remain + 1; j <= i__1; j += 8) {
	i1 = apnt[j + 1] - *m;
	i2 = apnt[j + 2] - *m;
	i3 = apnt[j + 3] - *m;
	i4 = apnt[j + 4] - *m;
	i5 = apnt[j + 5] - *m;
	i6 = apnt[j + 6] - *m;
	i7 = apnt[j + 7] - *m;
	i8 = apnt[j + 8] - *m;
	a1 = -a[i1];
	a2 = -a[i2];
	a3 = -a[i3];
	a4 = -a[i4];
	a5 = -a[i5];
	a6 = -a[i6];
	a7 = -a[i7];
	a8 = -a[i8];
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    y[i] = y[i] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[i4] + 
		    a5 * a[i5] + a6 * a[i6] + a7 * a[i7] + a8 * a[i8];
	    ++i1;
	    ++i2;
	    ++i3;
	    ++i4;
	    ++i5;
	    ++i6;
	    ++i7;
	    ++i8;
/* L3000: */
	}
/* L4000: */
    }

    return 0;
} /* smxpy8_ */

