/* fsup1.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ****************    FSUP1 ..... FIND SUPERNODES #1    ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE FIRST OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT RETURNS ONLY THE NUMBER OF */
/*       SUPERNODES NSUPER AND THE SUPERNODE MEMBERSHIP VECTOR SNODE(*), */
/*       WHICH IS OF LENGTH NEQNS.  THE VECTORS OF LENGTH NSUPER ARE */
/*       COMPUTED SUBSEQUENTLY BY THE COMPANION ROUTINE FSUP2. */

/*   METHOD AND ASSUMPTIONS: */
/*       THIS ROUTINE USES THE ELIMINATION TREE AND THE FACTOR COLUMN */
/*       COUNTS TO COMPUTE THE SUPERNODE PARTITION; IT ALSO ASSUMES A */
/*       POSTORDERING OF THE ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           FACTOR COLUMN COUNTS: I.E., THE NUMBER OF */
/*                           NONZERO ENTRIES IN EACH COLUMN OF L */
/*                           (INCLUDING THE DIAGONAL ENTRY). */

/*   OUTPUT PARAMETERS: */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMBER 11, 1994. */

/* *********************************************************************** */

/* Subroutine */ int fsup1_(neqns, etpar, colcnt, nofsub, nsuper, snode)
integer *neqns, *etpar, *colcnt, *nofsub, *nsuper, *snode;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kcol;


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

/*       -------------------------------------------- */
/*       COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION. */
/*       -------------------------------------------- */
    /* Parameter adjustments */
    --snode;
    --colcnt;
    --etpar;

    /* Function Body */
    *nsuper = 1;
    snode[1] = 1;
    *nofsub = colcnt[1];
    i__1 = *neqns;
    for (kcol = 2; kcol <= i__1; ++kcol) {
	if (etpar[kcol - 1] == kcol) {
	    if (colcnt[kcol - 1] == colcnt[kcol] + 1) {
		snode[kcol] = *nsuper;
		goto L300;
	    }
	}
	++(*nsuper);
	snode[kcol] = *nsuper;
	*nofsub += colcnt[kcol];
L300:
	;
    }

    return 0;
} /* fsup1_ */

