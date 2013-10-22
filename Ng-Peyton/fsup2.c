/* fsup2.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ****************    FSUP2  ..... FIND SUPERNODES #2   ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO */
/*       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE */
/*       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE */
/*       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS. */


/*   ASSUMPTIONS: */
/*       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT */
/*       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   OUTPUT PARAMETERS: */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           SUPERNODE PARTITIONING. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMEBER 22, 1994. */

/* *********************************************************************** */

/* Subroutine */ int fsup2_(neqns, nsuper, etpar, snode, xsuper)
integer *neqns, *nsuper, *etpar, *snode, *xsuper;
{
    static integer kcol, ksup, lstsup;


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

/*       ------------------------------------------------- */
/*       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*). */
/*       ------------------------------------------------- */
    /* Parameter adjustments */
    --xsuper;
    --snode;
    --etpar;

    /* Function Body */
    lstsup = *nsuper + 1;
    for (kcol = *neqns; kcol >= 1; --kcol) {
	ksup = snode[kcol];
	if (ksup != lstsup) {
	    xsuper[lstsup] = kcol + 1;
	}
	lstsup = ksup;
/* L100: */
    }
    xsuper[1] = 1;

    return 0;
} /* fsup2_ */

