/* create.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/*     ------------------------------------------- */
/*     INSERT DIAGONAL NONZEROS INTO STRUCTURE AND */
/*     CREATE NUMERICAL VALUES. */
/*     ------------------------------------------- */

/* Subroutine */ int create_(n, xadj, adj, xadjf, adjf, anzf)
integer *n, *xadj, *adj, *xadjf, *adjf;
doublereal *anzf;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer next, j, nofnz, jstop, jstrt, ii;




    /* Parameter adjustments */
    --anzf;
    --adjf;
    --xadjf;
    --adj;
    --xadj;

    /* Function Body */
    next = 1;
    xadjf[1] = next;

    jstrt = xadj[1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jstop = xadj[j + 1] - 1;
	nofnz = jstop - jstrt + 1;

	adjf[next] = j;
	anzf[next] = (doublereal) (nofnz + 1);
	++next;

	i__2 = jstop;
	for (ii = jstrt; ii <= i__2; ++ii) {
	    adjf[next] = adj[ii];
	    anzf[next] = (float)-1.;
	    ++next;
/* L100: */
	}
	xadjf[j + 1] = next;
	jstrt = jstop + 1;
/* L200: */
    }

    return 0;
} /* create_ */

