/* blkfct.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFCT .....  BLOCK GENERAL SPARSE CHOLESKY     ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE CALLS THE BLOCK GENERAL SPARSE CHOLESKY ROUTINE, */
/*       BLKFC2. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING */
/*                           IT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING */
/*                           THE DIAGONAL ELEMENTS). */
/*       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED. */
/*       IWSIZ           -   SIZE OF INTEGER WORKING STORAGE */
/*       TMPSIZ          -   SIZE OF FLOATING POINT WORKING STORAGE. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY. */
/*       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   ON OUTPUT, CONTAINS CHOLESKY FACTOR. */
/*       IFLAG           -   ERROR FLAG. */
/*                               0: SUCCESSFUL FACTORIZATION. */
/*                              -1: NONPOSITIVE DIAGONAL ENCOUNTERED, */
/*                                  MATRIX IS NOT POSITIVE DEFINITE. */
/*                              -2: INSUFFICIENT WORKING STORAGE */
/*                                  [TEMP(*)]. */
/*                              -3: INSUFFICIENT WORKING STORAGE */
/*                                  [IWORK(*)]. */

/*   WORKING PARAMETERS: */
/*       IWORK           -   INTEGER WORKING STORAGE OF LENGTH */
/*                           2*NEQNS + 2*NSUPER. */
/*       TMPVEC          -   DOUBLE PRECISION WORKING STORAGE OF LENGTH */
/*                           NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int blkfct_(neqns, nsuper, xsuper, snode, split, xlindx, 
	lindx, xlnz, lnz, iwsiz, iwork, tmpsiz, tmpvec, iflag, mmpyn, smxpy)
integer *neqns, *nsuper, *xsuper, *snode, *split, *xlindx, *lindx, *xlnz;
doublereal *lnz;
integer *iwsiz, *iwork, *tmpsiz;
doublereal *tmpvec;
integer *iflag;
/* Subroutine */ int (*mmpyn) (), (*smxpy) ();
{
    extern /* Subroutine */ int blkfc2_();


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/* ********************************************************************* 
*/

    /* Parameter adjustments */
    --tmpvec;
    --iwork;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --split;
    --snode;
    --xsuper;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < (*neqns << 1) + (*nsuper << 1)) {
	*iflag = -3;
	return 0;
    }
    blkfc2_(nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &lindx[1], &
	    xlnz[1], &lnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[(*nsuper 
	    << 1) + 1], &iwork[(*nsuper << 1) + *neqns + 1], tmpsiz, &tmpvec[
	    1], iflag, mmpyn, smxpy);
    return 0;
} /* blkfct_ */

