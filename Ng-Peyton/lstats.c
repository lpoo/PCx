/* lstats.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */

/*     ----------------------------------------- */
/*     GATHER STATISTICS ABOUT FACTORIZATION ... */
/*     ----------------------------------------- */

/* Subroutine */ int lstats_(nsuper, xsuper, xlindx, lindx, xlnz, tmpsiz, 
	outunt)
integer *nsuper, *xsuper, *xlindx, *lindx, *xlnz, *tmpsiz, *outunt;
{
    /* Format strings */
    static char fmt_1[] = "(a40,i10)";
    static char fmt_2[] = "(a40,1pd20.10)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle(), s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer jlen, j, n, ncols, jsize, nofnz, nofsub;
    static doublereal fctops;
    static integer jsuper, maxsup;
    static doublereal slvops;
    static integer supsze;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, 0, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_2, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_2, 0 };





    /* Parameter adjustments */
    --xlnz;
    --lindx;
    --xlindx;
    --xsuper;

    /* Function Body */
    n = xsuper[*nsuper + 1] - 1;

    io___2.ciunit = *outunt;
    s_wsle(&io___2);
    do_lio(&c__9, &c__1, " ", 1L);
    e_wsle();
/*       ------------------------------------------------------- */
/*       DETERMINE THE NUMBER OF NONZEROS IN CHOLESKY FACTOR AND */
/*       THE NUMBER OF SUBSCRIPTS IN REPRESENTING THE SUPERNODAL */
/*       STRUCTURE. */
/*       ------------------------------------------------------- */
    nofnz = xlnz[n + 1] - 1;
    nofsub = xlindx[*nsuper + 1] - 1;
    io___5.ciunit = *outunt;
    s_wsfe(&io___5);
    do_fio(&c__1, "   NUMBER OF SUPERNODES               = ", 40L);
    do_fio(&c__1, (char *)&(*nsuper), (ftnlen)sizeof(integer));
    e_wsfe();
    io___6.ciunit = *outunt;
    s_wsfe(&io___6);
    do_fio(&c__1, "   NUMBER OF NONZEROS IN L            = ", 40L);
    do_fio(&c__1, (char *)&nofnz, (ftnlen)sizeof(integer));
    e_wsfe();
    io___7.ciunit = *outunt;
    s_wsfe(&io___7);
    do_fio(&c__1, "   NUMBER OF SUBSCRIPTS IN L          = ", 40L);
    do_fio(&c__1, (char *)&nofsub, (ftnlen)sizeof(integer));
    e_wsfe();

/*       ------------------------------------------------------- */
/*       DETERMINE THE LARGEST SUPERNODE IN THE CHOLESKY FACTOR. */
/*       ------------------------------------------------------- */
    maxsup = 0;
    supsze = 0;
    i__1 = *nsuper;
    for (jsuper = 1; jsuper <= i__1; ++jsuper) {
/*           --------------------------------------------------- */
/*           NCOLS IS THE NUMBER OF COLUMNS IN SUPERNODE JSUPER. */
/*           --------------------------------------------------- */
	ncols = xsuper[jsuper + 1] - xsuper[jsuper];
	if (ncols > maxsup) {
	    maxsup = ncols;
	}

/*           ---------------------------------------------------- */
/*           JSIZE IS THE NUMBER OF NONZEROS IN SUPERNDOE JSUPER. */
/*           ---------------------------------------------------- */
	jlen = xlindx[jsuper + 1] - xlindx[jsuper];
	jsize = ((jlen << 1) - ncols + 1) * ncols / 2;
	if (jsize > supsze) {
	    supsze = jsize;
	}
/* L100: */
    }
    io___14.ciunit = *outunt;
    s_wsfe(&io___14);
    do_fio(&c__1, "   LARGEST SUPERNODE BY COLUMNS       = ", 40L);
    do_fio(&c__1, (char *)&maxsup, (ftnlen)sizeof(integer));
    e_wsfe();
    io___15.ciunit = *outunt;
    s_wsfe(&io___15);
    do_fio(&c__1, "   LARGEST SUPERNODE BY NONZEROS      = ", 40L);
    do_fio(&c__1, (char *)&supsze, (ftnlen)sizeof(integer));
    e_wsfe();

    io___16.ciunit = *outunt;
    s_wsfe(&io___16);
    do_fio(&c__1, "   SIZE OF TEMPORARY WORK STORAGE     = ", 40L);
    do_fio(&c__1, (char *)&(*tmpsiz), (ftnlen)sizeof(integer));
    e_wsfe();

/*       --------------------------- */
/*       DETERMINE OPERATION COUNTS. */
/*       --------------------------- */
    slvops = (float)0.;
    fctops = (float)0.;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	jlen = xlnz[j + 1] - xlnz[j];
	slvops = slvops + (jlen << 1) - 1;
/* Computing 2nd power */
	i__2 = jlen;
	fctops = fctops + i__2 * i__2 - 1;
/* L400: */
    }
    slvops *= 2;
    io___20.ciunit = *outunt;
    s_wsfe(&io___20);
    do_fio(&c__1, "   FACTORIZATION OPERATION COUNT      = ", 40L);
    do_fio(&c__1, (char *)&fctops, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___21.ciunit = *outunt;
    s_wsfe(&io___21);
    do_fio(&c__1, "   TRIANGULAR SOLN OPERATION COUNT    = ", 40L);
    do_fio(&c__1, (char *)&slvops, (ftnlen)sizeof(doublereal));
    e_wsfe();

    return 0;
} /* lstats_ */

