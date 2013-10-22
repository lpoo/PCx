/* getrhs.f -- translated by f2c (version of 25 March 1992  12:58:56).
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

/*     ------------------------------------ */
/*     CONSTRUCT RIGHT HAND SIDE VECTOR ... */
/*     ------------------------------------ */

/* Subroutine */ int getrhs_(n, xadjf, adjf, anzf, sol, rhs)
integer *n, *xadjf, *adjf;
doublereal *anzf, *sol, *rhs;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j;
    static doublereal t;
    static integer ii;




/*       --------------- */
/*       INITIALIZATION. */
/*       --------------- */
    /* Parameter adjustments */
    --rhs;
    --sol;
    --anzf;
    --adjf;
    --xadjf;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	rhs[j] = (float)0.;
/* L100: */
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*           ------------------- */
/*           FOR EACH COLUMN ... */
/*           ------------------- */
	t = sol[j];
	i__2 = xadjf[j + 1] - 1;
	for (ii = xadjf[j]; ii <= i__2; ++ii) {
	    rhs[adjf[ii]] += t * anzf[ii];
/* L200: */
	}
/* L300: */
    }

    return 0;
} /* getrhs_ */

