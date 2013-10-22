/* ordmmd.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ****     ORDMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE CALLS LIU'S MULTIPLE MINIMUM DEGREE */
/*               ROUTINE. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE. */
/*        IWSIZ  - SIZE OF INTEGER WORKING STORAGE. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE MINIMUM DEGREE ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO */
/*                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME. */
/*        IFLAG  - ERROR FLAG. */
/*                   0: SUCCESSFUL ORDERING */
/*                  -1: INSUFFICIENT WORKING STORAGE */
/*                      [IWORK(*)]. */

/*     WORKING PARAMETERS - */
/*        IWORK  - INTEGER WORKSPACE OF LENGTH 4*NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int ordmmd_(neqns, xadj, adjncy, invp, perm, iwsiz, iwork, 
	nofsub, iflag)
integer *neqns, *xadj, *adjncy, *invp, *perm, *iwsiz, *iwork, *nofsub, *iflag;
/*
 * iwszi = Worksize;
 * neqns = dimensions;
 *
 */
{
    static integer delta;
    extern /* Subroutine */ int genmmd_();
    static integer maxint;


/* ***********************************************************************
 */


/* ********************************************************************* 
*/

    /* Parameter adjustments */
    --iwork;
    --perm;
    --invp;
    --adjncy;
    --xadj;

    /* Function Body */
    *iflag = 0;
//    if (*iwsiz < *neqns << 2) {
//	*iflag = -1;
//	return 0;
//    }

/*       DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
/*       MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER */
/*                (ANY SMALLER ESTIMATE WILL DO) FOR MARKING */
/*                NODES. */

    delta = 0;
    maxint = 32767;
    genmmd_(neqns, &xadj[1], &adjncy[1], &invp[1], &perm[1], &delta, &iwork[1]
	    , &iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*neqns * 
	    3 + 1], &maxint, nofsub);
    return 0;

} /* ordmmd_ */

