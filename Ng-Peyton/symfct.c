/* symfct.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFCT ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE CALLS SYMFC2 WHICH PERFORMS SUPERNODAL SYMBOLIC */
/*       FACTORIZATION ON A REORDERED LINEAR SYSTEM. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */
/*       (I) IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                              -1 - INSUFFICIENT INTEGER WORKING SPACE. */
/*                              -2 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) IWORK       -   WORKING ARRAY OF LENGTH NSUPER+2*NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int symfct_(neqns, adjlen, xadj, adjncy, perm, invp, colcnt, 
	nsuper, xsuper, snode, nofsub, xlindx, lindx, xlnz, iwsiz, iwork, 
	flag_)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *colcnt, *nsuper, *
	xsuper, *snode, *nofsub, *xlindx, *lindx, *xlnz, *iwsiz, *iwork, *
	flag_;
{
    extern /* Subroutine */ int symfc2_();


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    --iwork;
    --xlnz;
    --lindx;
    --xlindx;
    --snode;
    --xsuper;
    --colcnt;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    *flag_ = 0;
    if (*iwsiz < *nsuper + (*neqns << 1) + 1) {
	*flag_ = -1;
	return 0;
    }
    symfc2_(neqns, adjlen, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[
	    1], nsuper, &xsuper[1], &snode[1], nofsub, &xlindx[1], &lindx[1], 
	    &xlnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[*nsuper + *neqns 
	    + 2], flag_);
    return 0;
} /* symfct_ */

