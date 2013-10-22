/* ordnat.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ****     ORDNAT ..... NATURAL ORDERING                     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE RECORDS THE INITIAL ORDERING IN THE */
/*               ORDERING VECTORS PERM AND INVP. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE "NATURAL" ORDERING; I.E., THE INITIAL */
/*                 ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */

/* *********************************************************************** */

/* Subroutine */ int ordnat_(neqns, perm, invp)
integer *neqns, *perm, *invp;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;


/* ***********************************************************************
 */



/* ***********************************************************************
 */

    /* Parameter adjustments */
    --invp;
    --perm;

    /* Function Body */
    i__1 = *neqns;
    for (i = 1; i <= i__1; ++i) {
	perm[i] = i;
	invp[i] = i;
/* L700: */
    }
    return 0;

} /* ordnat_ */

