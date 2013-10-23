/* blkLVL.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;

/* Interface to the BLKFCT routine. */
/* PCx beta-0.1   5/1/96 */

/* Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright. */

/* (C) 1996 University of Chicago. See COPYRIGHT in main directory. */
/* It's called from the C routine with all the required */
/* parameters except for the matrix-matrix and matrix-vector */
/* multiply routines MMPYN and SMXPY. Instead, there's a */
/* flag "LEVEL" which indicates what level of loop unrolling */
/* is required. This routine plugs in the appropriate two */
/* routines and then calls BLKFCT */
/* Subroutine */ int blklvl_(neqns, nsuper, xsuper, snode, split, xlindx,
	lindx, xlnz, lnz, iwsiz, iwork, tmpsiz, tmpvec, iflag, level)
integer *neqns, *nsuper, *xsuper, *snode, *split, *xlindx, *lindx, *xlnz;
doublereal *lnz;
integer *iwsiz, *iwork, *tmpsiz;
doublereal *tmpvec;
integer *iflag, *level;
{
    /* Format strings */
    static char fmt_23[] = "(\002*** LOOP UNROLLING LEVEL = \002,i5,\002 **\
*\002/,\002*** SHOULD HAVE LEVEL = 1, 2, 4, OR 8 ***\002)";

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern /* Subroutine */ int mmpy1_(), mmpy2_(), mmpy4_(), mmpy8_(),
	    smxpy1_(), smxpy2_(), smxpy4_();
    extern /* Subroutine */ int blkfct_();
    extern /* Subroutine */ int smxpy8_();

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_23, 0 };



/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* ***********************************************************************
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
    if (*level != 1 && *level != 2 && *level != 4 && *level != 8) {
	s_wsfe(&io___1);
	do_fio(&c__1, (char *)&(*level), (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", 0L);
    }
    if (*level == 1) {
	blkfct_(neqns, nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &
		lindx[1], &xlnz[1], &lnz[1], iwsiz, &iwork[1], tmpsiz, &
		tmpvec[1], iflag, mmpy1_, smxpy1_);
    } else if (*level == 2) {
	blkfct_(neqns, nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &
		lindx[1], &xlnz[1], &lnz[1], iwsiz, &iwork[1], tmpsiz, &
		tmpvec[1], iflag, mmpy2_, smxpy2_);
    } else if (*level == 4) {
	blkfct_(neqns, nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &
		lindx[1], &xlnz[1], &lnz[1], iwsiz, &iwork[1], tmpsiz, &
		tmpvec[1], iflag, mmpy4_, smxpy4_);
    } else if (*level == 8) {
	blkfct_(neqns, nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &
		lindx[1], &xlnz[1], &lnz[1], iwsiz, &iwork[1], tmpsiz, &
		tmpvec[1], iflag, mmpy8_, smxpy8_);
    }
    return 0;
} /* blklvl_ */
