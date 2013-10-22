/* chordr.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* **********     CHORDR ..... CHILD REORDERING                *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       REARRANGE THE CHILDREN OF EACH VERTEX SO THAT THE LAST ONE */
/*       MAXIMIZES (AMONG THE CHILDREN) THE NUMBER OF NONZEROS IN THE */
/*       CORRESPONDING COLUMN OF L.  ALSO DETERMINE AN NEW POSTORDERING */
/*       BASED ON THE STRUCTURE OF THE MODIFIED ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE NEW */
/*                           POSTORDERING. */
/*       COLCNT          -   COLUMN COUNTS IN L UNDER INITIAL ORDERING; */
/*                           MODIFIED TO REFLECT THE NEW ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BTREE2, EPOST2, INVINV. */

/* *********************************************************************** */

/* Subroutine */ int chordr_(neqns, xadj, adjncy, perm, invp, colcnt, parent, 
	fson, brothr, invpos)
integer *neqns, *xadj, *adjncy, *perm, *invp, *colcnt, *parent, *fson, *
	brothr, *invpos;
{
    extern /* Subroutine */ int btree2_(), epost2_(), invinv_();


/* ***********************************************************************
 */



/* ***********************************************************************
 */

/*       ---------------------------------------------------------- */
/*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE, */
/*       SO THAT EACH "LAST CHILD" MAXIMIZES AMONG ITS SIBLINGS THE */
/*       NUMBER OF NONZEROS IN THE CORRESPONDING COLUMNS OF L. */
/*       ---------------------------------------------------------- */
    /* Parameter adjustments */
    --invpos;
    --brothr;
    --fson;
    --parent;
    --colcnt;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    btree2_(neqns, &parent[1], &colcnt[1], &fson[1], &brothr[1], &invpos[1]);

/*       ---------------------------------------------------- */
/*       POSTORDER THE ELIMINATION TREE (USING THE NEW BINARY */
/*       REPRESENTATION. */
/*       ---------------------------------------------------- */
    epost2_(neqns, &fson[1], &brothr[1], &invpos[1], &parent[1], &colcnt[1], &
	    perm[1]);

/*       -------------------------------------------------------- */
/*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
/*       -------------------------------------------------------- */
    invinv_(neqns, &invp[1], &invpos[1], &perm[1]);

    return 0;
} /* chordr_ */

