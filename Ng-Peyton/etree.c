/* etree.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ****************     ETREE ..... ELIMINATION TREE     ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND */
/*       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */
/*       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */

/*   WORKING PARAMETERS: */
/*       ANCSTR          -   THE ANCESTOR VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int etree_(neqns, xadj, adjncy, perm, invp, parent, ancstr)
integer *neqns, *xadj, *adjncy, *perm, *invp, *parent, *ancstr;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer node, next, i, j, jstop, jstrt, nbr;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --ancstr;
    --parent;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (i = 1; i <= i__1; ++i) {
	parent[i] = 0;
	ancstr[i] = 0;
	node = perm[i];

	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	if (jstrt <= jstop) {
	    i__2 = jstop;
	    for (j = jstrt; j <= i__2; ++j) {
		nbr = adjncy[j];
		nbr = invp[nbr];
		if (nbr < i) {
/*                       --------------------------------
----------- */
/*                       FOR EACH NBR, FIND THE ROOT OF IT
S CURRENT */
/*                       ELIMINATION TREE.  PERFORM PATH C
OMPRESSION */
/*                       AS THE SUBTREE IS TRAVERSED. */
/*                       --------------------------------
----------- */
L100:
		    if (ancstr[nbr] == i) {
			goto L300;
		    }
		    if (ancstr[nbr] > 0) {
			next = ancstr[nbr];
			ancstr[nbr] = i;
			nbr = next;
			goto L100;
		    }
/*                       --------------------------------
------------ */
/*                       NOW, NBR IS THE ROOT OF THE SUBTR
EE.  MAKE I */
/*                       THE PARENT NODE OF THIS ROOT. */
/*                       --------------------------------
------------ */
		    parent[nbr] = i;
		    ancstr[nbr] = i;
		}
L300:
		;
	    }
	}
/* L400: */
    }

    return 0;
} /* etree_ */

