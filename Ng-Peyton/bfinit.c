/* bfinit.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING */
/*       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */
/*       CACHSZ          -   CACHE SIZE (IN KBYTES). */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */

/* *********************************************************************** */

/* Subroutine */ int bfinit_(neqns, nsuper, xsuper, snode, xlindx, lindx, 
	cachsz, tmpsiz, split)
integer *neqns, *nsuper, *xsuper, *snode, *xlindx, *lindx, *cachsz, *tmpsiz, *
	split;
{
    extern /* Subroutine */ int fnsplt_(), fntsiz_();


/* ***********************************************************************
 */


/* ***********************************************************************
 */

/*       --------------------------------------------------- */
/*       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT. */
/*       --------------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --lindx;
    --xlindx;
    --snode;
    --xsuper;

    /* Function Body */
    fntsiz_(nsuper, &xsuper[1], &snode[1], &xlindx[1], &lindx[1], tmpsiz);

/*       ------------------------------- */
/*       PARTITION SUPERNODES FOR CACHE. */
/*       ------------------------------- */
    fnsplt_(neqns, nsuper, &xsuper[1], &xlindx[1], cachsz, &split[1]);

    return 0;
} /* bfinit_ */

