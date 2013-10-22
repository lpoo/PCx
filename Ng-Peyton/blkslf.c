/* blkslf.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Written:        October 6, 1996 by SJW. Based on routine BLKSLV of */
/*                   Esmond G. Ng and Barry W. Peyton. */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKSLF ... FORWARD TRIANGULAR SUBSTITUTION     ********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC */
/*       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE */
/*       FORWARD TRIANGULAR SUBSTITUTIOn.  IT USES OUTPUT FROM BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE. */
/*       (XLNZ,LNZ)      -   CHOLESKY FACTOR. */

/*   UPDATED PARAMETERS: */
/*       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON */
/*                           OUTPUT, CONTAINS THE SOLUTION. */

/* *********************************************************************** */

/* Subroutine */ int blkslf_(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
integer *nsuper, *xsuper, *xlindx, *lindx, *xlnz;
doublereal *lnz, *rhs;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer jcol, ipnt, jpnt, jsup, i;
    static doublereal t;
    static integer fjcol, ljcol, ix, ixstop, ixstrt;


/* ***********************************************************************
 */


/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --rhs;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --xsuper;

    /* Function Body */
    if (*nsuper <= 0) {
	return 0;
    }

/*       ------------------------ */
/*       FORWARD SUBSTITUTION ... */
/*       ------------------------ */
    fjcol = xsuper[1];
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {
	ljcol = xsuper[jsup + 1] - 1;
	ixstrt = xlnz[fjcol];
	jpnt = xlindx[jsup];
	i__2 = ljcol;
	for (jcol = fjcol; jcol <= i__2; ++jcol) {
	    ixstop = xlnz[jcol + 1] - 1;
	    t = rhs[jcol] / lnz[ixstrt];
	    rhs[jcol] = t;
	    ipnt = jpnt + 1;
/* DIR$           IVDEP */
	    i__3 = ixstop;
	    for (ix = ixstrt + 1; ix <= i__3; ++ix) {
		i = lindx[ipnt];
		rhs[i] -= t * lnz[ix];
		++ipnt;
/* L100: */
	    }
	    ixstrt = ixstop + 1;
	    ++jpnt;
/* L200: */
	}
	fjcol = ljcol + 1;
/* L300: */
    }

    return 0;
} /* blkslf_ */

