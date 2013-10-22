/* blkslb.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* *********     BLKSLB ... BACK TRIANGULAR SUBSTITUTION        ********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC */
/*       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE */
/*       BACKWARD TRIANGULAR SUBSTITUTION.  IT USES OUTPUT FROM BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE. */
/*       (XLNZ,LNZ)      -   CHOLESKY FACTOR. */

/*   UPDATED PARAMETERS: */
/*       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON */
/*                           OUTPUT, CONTAINS THE SOLUTION. */

/* *********************************************************************** */

/* Subroutine */ int blkslb_(nsuper, xsuper, xlindx, lindx, xlnz, lnz, rhs)
integer *nsuper, *xsuper, *xlindx, *lindx, *xlnz;
doublereal *lnz, *rhs;
{
    /* System generated locals */
    integer i__1, i__2;

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
/*       ------------------------- */
/*       BACKWARD SUBSTITUTION ... */
/*       ------------------------- */
    ljcol = xsuper[*nsuper + 1] - 1;
    for (jsup = *nsuper; jsup >= 1; --jsup) {
	fjcol = xsuper[jsup];
	ixstop = xlnz[ljcol + 1] - 1;
	jpnt = xlindx[jsup] + (ljcol - fjcol);
	i__1 = fjcol;
	for (jcol = ljcol; jcol >= i__1; --jcol) {
	    ixstrt = xlnz[jcol];
	    ipnt = jpnt + 1;
	    t = rhs[jcol];
/* DIR$           IVDEP */
	    i__2 = ixstop;
	    for (ix = ixstrt + 1; ix <= i__2; ++ix) {
		i = lindx[ipnt];
		t -= lnz[ix] * rhs[i];
		++ipnt;
/* L400: */
	    }
	    rhs[jcol] = t / lnz[ixstrt];
	    ixstop = ixstrt - 1;
	    --jpnt;
/* L500: */
	}
	ljcol = fjcol - 1;
/* L600: */
    }

    return 0;
} /* blkslb_ */

