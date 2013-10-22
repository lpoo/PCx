/* invinv.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***********     INVINV ..... CONCATENATION OF TWO INVP     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO PERFORM THE MAPPING OF */
/*           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP */
/*       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION */
/*       VECTOR PERM IS ALSO COMPUTED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR. */

/*   UPDATED PARAMETERS: */
/*       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON */
/*                           OUTPUT, IT CONTAINS THE NEW INVERSE */
/*                           PERMUTATION. */

/*   OUTPUT PARAMETER: */
/*       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS */
/*                           INVP2). */

/* *********************************************************************** */

/* Subroutine */ int invinv_(neqns, invp, invp2, perm)
integer *neqns, *invp, *invp2, *perm;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, i, interm;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --perm;
    --invp2;
    --invp;

    /* Function Body */
    i__1 = *neqns;
    for (i = 1; i <= i__1; ++i) {
	interm = invp[i];
	invp[i] = invp2[interm];
/* L100: */
    }

    i__1 = *neqns;
    for (i = 1; i <= i__1; ++i) {
	node = invp[i];
	perm[node] = i;
/* L200: */
    }

    return 0;
} /* invinv_ */

