/* chlsup.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratoy */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     CHLSUP .... DENSE CHOLESKY WITHIN SUPERNODE   ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY */
/*               FACTORIZATION ON THE COLUMNS OF A SUPERNODE */
/*               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS */
/*               EXTERNAL TO THE SUPERNODE. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN). */
/*        N      - NUMBER OF COLUMNS IN THE SUPERNODE. */
/*        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END */
/*                 OF THE J-TH COLUMN OF THE SUPERNODE. */
/*        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO */
/*                 BE FACTORED. */
/*        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*     OUTPUT PARAMETERS - */
/*        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF */
/*                 THE SUPERNODE. */
/*        IFLAG  - UNCHANGED IF THERE IS NO ERROR. */
/*                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED. */

/* *********************************************************************** */

/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC */
/* Subroutine */ int chlsup_(m, n, split, xpnt, x, mxdiag, ntiny, iflag, 
	mmpyn, smxpy)
integer *m, *n, *split, *xpnt;
doublereal *x, *mxdiag;
integer *ntiny, *iflag;
/* Subroutine */ int (*mmpyn) (), (*smxpy) ();
{
    static integer jblk, jpnt, q;
    extern /* Subroutine */ int pchol_();
    static integer mm, nn, fstcol, nxtcol;


/* ***********************************************************************
 */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */




/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/

/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --x;
    --xpnt;
    --split;

    /* Function Body */
    jblk = 0;
    fstcol = 1;
    mm = *m;
    jpnt = xpnt[fstcol];

/*       ---------------------------------------- */
/*       FOR EACH BLOCK JBLK IN THE SUPERNODE ... */
/*       ---------------------------------------- */
L100:
    if (fstcol <= *n) {
	++jblk;
	nn = split[jblk];
/*           ------------------------------------------ */
/*           ... PERFORM PARTIAL CHOLESKY FACTORIZATION */
/*               ON THE BLOCK. */
/*           ------------------------------------------ */
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCx
PC */
	pchol_(&mm, &nn, &xpnt[fstcol], &x[1], mxdiag, ntiny, iflag, smxpy);
	if (*iflag == 1) {
	    return 0;
	}
/*           ---------------------------------------------- */
/*           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS */
/*               OF THE SUPERNODE REMAINING TO BE COMPUTED. */
/*           ---------------------------------------------- */
	nxtcol = fstcol + nn;
	q = *n - nxtcol + 1;
	mm -= nn;
	jpnt = xpnt[nxtcol];
	if (q > 0) {
	    (*mmpyn)(&mm, &nn, &q, &xpnt[fstcol], &x[1], &x[jpnt], &mm);
	}
	fstcol = nxtcol;
	goto L100;
    }

    return 0;
} /* chlsup_ */

