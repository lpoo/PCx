/* igathr.f -- translated by f2c (version of 25 March 1992  12:58:56).
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
/* ******         IGATHR .... INTEGER GATHER OPERATION      ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A STANDARD INTEGER GATHER */
/*               OPERATION. */

/*     INPUT PARAMETERS - */
/*        KLEN   - LENGTH OF THE LIST OF GLOBAL INDICES. */
/*        LINDX  - LIST OF GLOBAL INDICES. */
/*        INDMAP - INDEXED BY GLOBAL INDICES, IT CONTAINS THE */
/*                 REQUIRED RELATIVE INDICES. */

/*     OUTPUT PARAMETERS - */
/*        RELIND - LIST RELATIVE INDICES. */

/* *********************************************************************** */

/* Subroutine */ int igathr_(klen, lindx, indmap, relind)
integer *klen, *lindx, *indmap, *relind;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;


/* ***********************************************************************
 */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */

/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */

/* ***********************************************************************
 */

/* DIR$ IVDEP */
    /* Parameter adjustments */
    --relind;
    --indmap;
    --lindx;

    /* Function Body */
    i__1 = *klen;
    for (i = 1; i <= i__1; ++i) {
	relind[i] = indmap[lindx[i]];
/* L100: */
    }
    return 0;
} /* igathr_ */

