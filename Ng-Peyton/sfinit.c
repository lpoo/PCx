/* sfinit.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* **************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP */
/*       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION. */

/*   NOTE: */
/*       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E., */
/*       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES). */

/*   INPUT PARAMETERS: */
/*       NEQNS       -   NUMBER OF EQUATIONS. */
/*       NNZA        -   LENGTH OF ADJACENCY STRUCTURE. */
/*       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                       TO THE ADJACENCY STRUCTURE. */
/*       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                       THE ADJACENCY STRUCTURE. */
/*       PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       POSTORDERING. */
/*       INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       INVERSE OF THE POSTORDERING. */
/*       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                       OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                       INCLUDING THE DIAGONAL ENTRY. */
/*       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                       THE DIAGONAL ENTRIES. */
/*       NSUB        -   NUMBER OF SUBSCRIPTS. */
/*       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                       SUPERNODE MEMBERSHIP. */
/*       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                       SUPERNODE PARTITIONING. */
/*       IFLAG(*)    -   ERROR FLAG. */
/*                          0: SUCCESSFUL SF INITIALIZATION. */
/*                         -1: INSUFFICENT WORKING STORAGE */
/*                             [IWORK(*)]. */

/*   WORK PARAMETERS: */
/*       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3. */

/*   FIRST CREATED ON    NOVEMEBER 14, 1994. */
/*   LAST UPDATED ON     January 12, 1995. */

/* *********************************************************************** */

/* Subroutine */ int sfinit_(neqns, nnza, xadj, adjncy, perm, invp, colcnt, 
	nnzl, nsub, nsuper, snode, xsuper, iwsiz, iwork, iflag)
integer *neqns, *nnza, *xadj, *adjncy, *perm, *invp, *colcnt, *nnzl, *nsub, *
	nsuper, *snode, *xsuper, *iwsiz, *iwork, *iflag;
{
    extern /* Subroutine */ int fsup1_(), fsup2_(), fcnthn_(), chordr_(), 
	    etordr_();


/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* ***********************************************************************
 */

/*       -------------------------------------------------------- */
/*       RETURN IF THERE IS INSUFFICIENT INTEGER WORKING STORAGE. */
/*       -------------------------------------------------------- */
    /* Parameter adjustments */
    --iwork;
    --xsuper;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < *neqns * 7 + 3) {
	*iflag = -1;
	return 0;
    }

/*       ------------------------------------------ */
/*       COMPUTE ELIMINATION TREE AND POSTORDERING. */
/*       ------------------------------------------ */
    etordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], &
	    iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*neqns * 3 + 
	    1]);

/*       --------------------------------------------- */
/*       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS. */
/*       --------------------------------------------- */
    fcnthn_(neqns, nnza, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], 
	    &snode[1], &colcnt[1], nnzl, &iwork[*neqns + 1], &iwork[(*neqns <<
	     1) + 1], &xsuper[1], &iwork[*neqns * 3 + 1], &iwork[(*neqns << 2)
	     + 2], &iwork[*neqns * 5 + 3], &iwork[*neqns * 6 + 4]);

/*       --------------------------------------------------------- */
/*       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM */
/*       NUMBER OF NONZEROS IN ITS COLUMN OF L. */
/*       --------------------------------------------------------- */
    chordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[1], &
	    iwork[1], &iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*
	    neqns * 3 + 1]);

/*       ---------------- */
/*       FIND SUPERNODES. */
/*       ---------------- */
    fsup1_(neqns, &iwork[1], &colcnt[1], nsub, nsuper, &snode[1]);
    fsup2_(neqns, nsuper, &iwork[1], &snode[1], &xsuper[1]);

    return 0;
} /* sfinit_ */

