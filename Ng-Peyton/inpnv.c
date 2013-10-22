/* inpnv.f -- translated by f2c (version of 25 March 1992  12:58:56).
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

/*     ------------------------------------------------------ */
/*     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ... */
/*     ------------------------------------------------------ */

/* Subroutine */ int inpnv_(neqns, xadjf, adjf, anzf, perm, invp, nsuper, 
	xsuper, xlindx, lindx, xlnz, lnz, offset)
integer *neqns, *xadjf, *adjf;
doublereal *anzf;
integer *perm, *invp, *nsuper, *xsuper, *xlindx, *lindx, *xlnz;
doublereal *lnz;
integer *offset;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer jlen, oldj, last, i, j, ii, jsuper;




    /* Parameter adjustments */
    --offset;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --xsuper;
    --invp;
    --perm;
    --anzf;
    --adjf;
    --xadjf;

    /* Function Body */
    i__1 = *nsuper;
    for (jsuper = 1; jsuper <= i__1; ++jsuper) {

/*           ---------------------------------------- */
/*           FOR EACH SUPERNODE, DO THE FOLLOWING ... */
/*           ---------------------------------------- */

/*           ----------------------------------------------- */
/*           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT. */
/*           ----------------------------------------------- */
	jlen = xlindx[jsuper + 1] - xlindx[jsuper];
	i__2 = xlindx[jsuper + 1] - 1;
	for (ii = xlindx[jsuper]; ii <= i__2; ++ii) {
	    i = lindx[ii];
	    --jlen;
	    offset[i] = jlen;
/* L100: */
	}

	i__2 = xsuper[jsuper + 1] - 1;
	for (j = xsuper[jsuper]; j <= i__2; ++j) {
/*               ----------------------------------------- */
/*               FOR EACH COLUMN IN THE CURRENT SUPERNODE, */
/*               FIRST INITIALIZE THE DATA STRUCTURE. */
/*               ----------------------------------------- */
	    i__3 = xlnz[j + 1] - 1;
	    for (ii = xlnz[j]; ii <= i__3; ++ii) {
		lnz[ii] = (float)0.;
/* L200: */
	    }

/*               ----------------------------------- */
/*               NEXT INPUT THE INDIVIDUAL NONZEROS. */
/*               ----------------------------------- */
	    oldj = perm[j];
	    last = xlnz[j + 1] - 1;
	    i__3 = xadjf[oldj + 1] - 1;
	    for (ii = xadjf[oldj]; ii <= i__3; ++ii) {
		i = invp[adjf[ii]];
		if (i >= j) {
		    lnz[last - offset[i]] = anzf[ii];
		}
/* L300: */
	    }
/* L400: */
	}

/* L500: */
    }
    return 0;
} /* inpnv_ */

