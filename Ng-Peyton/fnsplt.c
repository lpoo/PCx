/* fnsplt.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  May 26, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     ***** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES */
/*       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE */
/*       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           POINTERS IN THE SUPERNODE INDICES. */
/*       CACHSZ          -   CACHE SIZE IN KILO BYTES. */
/*                           IF THERE IS NO CACHE, SET CACHSZ = 0. */

/*   OUTPUT PARAMETERS: */
/*       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE */
/*                           FINE PARTITIONING. */

/* *********************************************************************** */

/* Subroutine */ int fnsplt_(neqns, nsuper, xsuper, xlindx, cachsz, split)
integer *neqns, *nsuper, *xsuper, *xlindx, *cachsz, *split;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kcol, used, ksup, cache, ncols, width, height, curcol, 
	    fstcol, lstcol, nxtblk;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ******************************************************************* */

/*       -------------------------------------------- */
/*       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE. */
/*       -------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --xlindx;
    --xsuper;

    /* Function Body */
    if (*cachsz <= 0) {
	cache = 2000000000;
    } else {
	cache = (real) (*cachsz) * (float)1024. / (float)8. * (float).9;
    }

/*       --------------- */
/*       INITIALIZATION. */
/*       --------------- */
    i__1 = *neqns;
    for (kcol = 1; kcol <= i__1; ++kcol) {
	split[kcol] = 0;
/* L100: */
    }

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
/*           ----------------------- */
/*           ... GET SUPERNODE INFO. */
/*           ----------------------- */
	height = xlindx[ksup + 1] - xlindx[ksup];
	fstcol = xsuper[ksup];
	lstcol = xsuper[ksup + 1] - 1;
	width = lstcol - fstcol + 1;
	nxtblk = fstcol;
/*           -------------------------------------- */
/*           ... UNTIL ALL COLUMNS OF THE SUPERNODE */
/*               HAVE BEEN PROCESSED ... */
/*           -------------------------------------- */
	curcol = fstcol - 1;
L200:
/*               ------------------------------------------- */
/*               ... PLACE THE FIRST COLUMN(S) IN THE CACHE. */
/*               ------------------------------------------- */
	++curcol;
	if (curcol < lstcol) {
	    ++curcol;
	    ncols = 2;
	    used = (height << 2) - 1;
	    height += -2;
	} else {
	    ncols = 1;
	    used = height * 3;
	    --height;
	}

/*               -------------------------------------- */
/*               ... WHILE THE CACHE IS NOT FILLED AND */
/*                   THERE ARE COLUMNS OF THE SUPERNODE */
/*                   REMAINING TO BE PROCESSED ... */
/*               -------------------------------------- */
L300:
	if (used + height < cache && curcol < lstcol) {
/*                   -------------------------------- */
/*                   ... ADD ANOTHER COLUMN TO CACHE. */
/*                   -------------------------------- */
	    ++curcol;
	    ++ncols;
	    used += height;
	    --height;
	    goto L300;
	}
/*               ------------------------------------- */
/*               ... RECORD THE NUMBER OF COLUMNS THAT */
/*                   FILLED THE CACHE. */
/*               ------------------------------------- */
	split[nxtblk] = ncols;
	++nxtblk;
/*               -------------------------- */
/*               ... GO PROCESS NEXT BLOCK. */
/*               -------------------------- */
	if (curcol < lstcol) {
	    goto L200;
	}
/* L1000: */
    }

    return 0;
} /* fnsplt_ */

