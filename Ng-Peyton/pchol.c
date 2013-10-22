/* pchol.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******     PCHOL .... DENSE PARTIAL CHOLESKY             ************** */
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
/* Subroutine */ int pchol_(m, n, xpnt, x, mxdiag, ntiny, iflag, smxpy)
integer *m, *n, *xpnt;
doublereal *x, *mxdiag;
integer *ntiny, *iflag;
/* Subroutine */ int (*smxpy) ();
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal diag;
    static integer jcol, jpnt;
    extern /* Subroutine */ int dscal1_();
    static integer mm;


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

/*       ------------------------------------------ */
/*       FOR EVERY COLUMN JCOL IN THE SUPERNODE ... */
/*       ------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --xpnt;

    /* Function Body */
    mm = *m;
    jpnt = xpnt[1];
    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {

/*           ---------------------------------- */
/*           UPDATE JCOL WITH PREVIOUS COLUMNS. */
/*           ---------------------------------- */
	if (jcol > 1) {
	    i__2 = jcol - 1;
	    (*smxpy)(&mm, &i__2, &x[jpnt], &xpnt[1], &x[1]);
	}

/*           --------------------------- */
/*           COMPUTE THE DIAGONAL ENTRY. */
/*           --------------------------- */
	diag = x[jpnt];
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCx
PC */
	if (diag <= *mxdiag * 1e-30) {
	    diag = 1e128;
	    ++(*ntiny);
	}
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCx
PC */
	diag = sqrt(diag);
	x[jpnt] = diag;
	diag = 1. / diag;

/*           ---------------------------------------------------- */
/*           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY. */
/*           ---------------------------------------------------- */
	--mm;
	++jpnt;
	dscal1_(&mm, &diag, &x[jpnt]);
	jpnt += mm;

/* L100: */
    }

    return 0;
} /* pchol_ */

