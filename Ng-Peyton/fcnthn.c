/* fcnthn.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************     FCNTHN  ..... FIND NONZERO COUNTS    *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN */
/*       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM. */

/*       TECHNIQUES: */
/*       1) SUPERNODE DETECTION. */
/*       2) PATH HALVING. */
/*       3) NO UNION BY RANK. */

/*   NOTES: */
/*       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH ROW OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                           THE DIAGONAL ENTRIES. */

/*   WORK PARAMETERS: */
/*       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE */
/*                           DISJOINT SETS (I.E., SUBTREES). */
/*       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS LEAF OF EACH ROW SUBTREE. */
/*       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL */
/*                           (DISTANCE FROM THE ROOT). */
/*       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS */
/*                           USED TO COMPUTE COLUMN COUNTS. */
/*       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT. */
/*       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           NUMBER OF CHILDREN. */
/*       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE. */

/*   FIRST CREATED ON    APRIL 12, 1990. */
/*   LAST UPDATED ON     JANUARY 12, 1995. */

/* *********************************************************************** */

/* Subroutine */ int fcnthn_(neqns, adjlen, xadj, adjncy, perm, invp, etpar, 
	rowcnt, colcnt, nlnz, set, prvlf, level, weight, fdesc, nchild, 
	prvnbr)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *etpar, *rowcnt, *
	colcnt, *nlnz, *set, *prvlf, *level, *weight, *fdesc, *nchild, *
	prvnbr;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer temp, xsup, last1, last2, j, k, lflag, pleaf, hinbr, jstop,
	     jstrt, ifdesc, oldnbr, parent, lownbr, lca;


/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

/*       -------------------------------------------------- */
/*       COMPUTE LEVEL(*), FDESC(*), NCHILD(*). */
/*       INITIALIZE XSUP, ROWCNT(*), COLCNT(*), */
/*                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*). */
/*       -------------------------------------------------- */
    /* Parameter adjustments */
    --prvnbr;
    --prvlf;
    --set;
    --colcnt;
    --rowcnt;
    --etpar;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    xsup = 1;
    level[0] = 0;
    for (k = *neqns; k >= 1; --k) {
	rowcnt[k] = 1;
	colcnt[k] = 0;
	set[k] = k;
	prvlf[k] = 0;
	level[k] = level[etpar[k]] + 1;
	weight[k] = 1;
	fdesc[k] = k;
	nchild[k] = 0;
	prvnbr[k] = 0;
/* L100: */
    }
    nchild[0] = 0;
    fdesc[0] = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	parent = etpar[k];
	weight[parent] = 0;
	++nchild[parent];
	ifdesc = fdesc[k];
	if (ifdesc < fdesc[parent]) {
	    fdesc[parent] = ifdesc;
	}
/* L200: */
    }
/*       ------------------------------------ */
/*       FOR EACH ``LOW NEIGHBOR'' LOWNBR ... */
/*       ------------------------------------ */
    i__1 = *neqns;
    for (lownbr = 1; lownbr <= i__1; ++lownbr) {
	lflag = 0;
	ifdesc = fdesc[lownbr];
	oldnbr = perm[lownbr];
	jstrt = xadj[oldnbr];
	jstop = xadj[oldnbr + 1] - 1;
/*           ----------------------------------------------- */
/*           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ... */
/*           ----------------------------------------------- */
	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    hinbr = invp[adjncy[j]];
	    if (hinbr > lownbr) {
		if (ifdesc > prvnbr[hinbr]) {
/*                       ------------------------- */
/*                       INCREMENT WEIGHT(LOWNBR). */
/*                       ------------------------- */
		    ++weight[lownbr];
		    pleaf = prvlf[hinbr];
/*                       --------------------------------
--------- */
/*                       IF HINBR HAS NO PREVIOUS ``LOW NE
IGHBOR'' */
/*                       THEN ... */
/*                       --------------------------------
--------- */
		    if (pleaf == 0) {
/*                           ------------------------
----------------- */
/*                           ... ACCUMULATE LOWNBR-->H
INBR PATH LENGTH */
/*                               IN ROWCNT(HINBR). */
/*                           ------------------------
----------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				hinbr];
		    } else {
/*                           ------------------------
----------------- */
/*                           ... OTHERWISE, LCA <-- FI
ND(PLEAF), WHICH */
/*                               IS THE LEAST COMMON A
NCESTOR OF PLEAF */
/*                               AND LOWNBR. */
/*                               (PATH HALVING.) */
/*                           ------------------------
----------------- */
			last1 = pleaf;
			last2 = set[last1];
			lca = set[last2];
L300:
			if (lca != last2) {
			    set[last1] = lca;
			    last1 = lca;
			    last2 = set[last1];
			    lca = set[last2];
			    goto L300;
			}
/*                           ------------------------
------------- */
/*                           ACCUMULATE PLEAF-->LCA PA
TH LENGTH IN */
/*                           ROWCNT(HINBR). */
/*                           DECREMENT WEIGHT(LCA). */

/*                           ------------------------
------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				lca];
			--weight[lca];
		    }
/*                       --------------------------------
-------------- */
/*                       LOWNBR NOW BECOMES ``PREVIOUS LEA
F'' OF HINBR. */
/*                       --------------------------------
-------------- */
		    prvlf[hinbr] = lownbr;
		    lflag = 1;
		}
/*                   ----------------------------------------
---------- */
/*                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' 
OF HINBR. */
/*                   ----------------------------------------
---------- */
		prvnbr[hinbr] = lownbr;
	    }
/* L500: */
	}
/*           ---------------------------------------------------- */
/*           DECREMENT WEIGHT ( PARENT(LOWNBR) ). */
/*           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP). */
/*           ---------------------------------------------------- */
	parent = etpar[lownbr];
	--weight[parent];
	if (lflag == 1 || nchild[lownbr] >= 2) {
	    xsup = lownbr;
	}
	set[xsup] = parent;
/* L600: */
    }
/*       --------------------------------------------------------- */
/*       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS. */
/*       --------------------------------------------------------- */
    *nlnz = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	temp = colcnt[k] + weight[k];
	colcnt[k] = temp;
	*nlnz += temp;
	parent = etpar[k];
	if (parent != 0) {
	    colcnt[parent] += temp;
	}
/* L700: */
    }

    return 0;
} /* fcnthn_ */

