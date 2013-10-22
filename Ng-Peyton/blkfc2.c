/* blkfc2.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratoy */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFC2 .....  BLOCK GENERAL SPARSE CHOLESKY     ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE FACTORS A SPARSE POSITIVE DEFINITE MATRIX. */
/*       THE COMPUTATION IS ORGANIZED AROUND KERNELS THAT PERFORM */
/*       SUPERNODE-TO-SUPERNODE UPDATES, I.E., BLOCK-TO-BLOCK UPDATES. */

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
/*       TMPSIZ          -   SIZE OF TEMPORARY WORKING STORAGE. */
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

/*   WORKING PARAMETERS: */
/*       LINK            -   LINKS TOGETHER THE SUPERNODES IN A SUPERNODE */
/*                           ROW. */
/*       LENGTH          -   LENGTH OF THE ACTIVE PORTION OF EACH */
/*                           SUPERNODE. */
/*       INDMAP          -   VECTOR OF SIZE NEQNS INTO WHICH THE GLOBAL */
/*                           INDICES ARE SCATTERED. */
/*       RELIND          -   MAPS LOCATIONS IN THE UPDATING COLUMNS TO */
/*                           THE CORRESPONDING LOCATIONS IN THE UPDATED */
/*                           COLUMNS.  (RELIND IS GATHERED FROM INDMAP). */
/*       TEMP            -   REAL VECTOR FOR ACCUMULATING UPDATES.  MUST */
/*                           ACCOMODATE ALL COLUMNS OF A SUPERNODE. */

/* *********************************************************************** */

/* Subroutine */ int blkfc2_(nsuper, xsuper, snode, split, xlindx, lindx, 
	xlnz, lnz, link, length, indmap, relind, tmpsiz, temp, iflag, mmpyn, 
	smxpy)
integer *nsuper, *xsuper, *snode, *split, *xlindx, *lindx, *xlnz;
doublereal *lnz;
integer *link, *length, *indmap, *relind, *tmpsiz;
doublereal *temp;
integer *iflag;
/* Subroutine */ int (*mmpyn) (), (*smxpy) ();
{
    /* Format strings */
    static char fmt_699[] = "(1x,\002 FOUND \002,i6,\002 TINY DIAGONALS; REP\
LACED WITH INF\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer ilen, jlen, klen, jsup, ksup;
    extern /* Subroutine */ int mmpy_();
    static integer i, fjcol, fkcol, ljcol;
    extern /* Subroutine */ int assmb_();
    static integer klast, ilpnt, jlpnt, klpnt, store;
    extern /* Subroutine */ int mmpyi_();
    static integer ntiny, jxpnt, kxpnt, inddif;
    static doublereal mxdiag;
    extern /* Subroutine */ int igathr_(), ldindx_();
    static integer njcols, nkcols;
    extern /* Subroutine */ int chlsup_();
    static integer ncolup, kfirst, nxtcol, nxksup, nxtsup;

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 6, 0, fmt_699, 0 };



/* ********************************************************************* 
*/

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/

/* ********************************************************************* 
*/

    /* Parameter adjustments */
    --temp;
    --relind;
    --indmap;
    --length;
    --link;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --split;
    --snode;
    --xsuper;

    /* Function Body */
    *iflag = 0;
    ntiny = 0;

/*       ----------------------------------------------------------- */
/*       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*). */
/*       ----------------------------------------------------------- */
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {
	link[jsup] = 0;
/* L100: */
    }
    i__1 = *tmpsiz;
    for (i = 1; i <= i__1; ++i) {
	temp[i] = 0.;
/* L200: */
    }
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/
/*       COMPUTE MAXIMUM DIAGONAL ELEMENT IN INPUT MATRIX */
    mxdiag = 0.;
    i__1 = xsuper[*nsuper + 1] - 1;
    for (i = 1; i <= i__1; ++i) {
	fjcol = xlnz[i];
/* Computing MAX */
	d__1 = mxdiag, d__2 = lnz[fjcol];
	mxdiag = max(d__1,d__2);
/* L201: */
    }
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/

/*       --------------------------- */
/*       FOR EACH SUPERNODE JSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {

/*           ------------------------------------------------ */
/*           FJCOL  ...  FIRST COLUMN OF SUPERNODE JSUP. */
/*           LJCOL  ...  LAST COLUMN OF SUPERNODE JSUP. */
/*           NJCOLS ...  NUMBER OF COLUMNS IN SUPERNODE JSUP. */
/*           JLEN   ...  LENGTH OF COLUMN FJCOL. */
/*           JXPNT  ...  POINTER TO INDEX OF FIRST */
/*                       NONZERO IN COLUMN FJCOL. */
/*           ------------------------------------------------ */
	fjcol = xsuper[jsup];
	njcols = xsuper[jsup + 1] - fjcol;
	ljcol = fjcol + njcols - 1;
	jlen = xlnz[fjcol + 1] - xlnz[fjcol];
	jxpnt = xlindx[jsup];
/*            print *, 'Super Node: ', JSUP, ' first: ', FJCOL, */
/*     .           ' last: ', LJCOL */


/*           ----------------------------------------------------- */
/*           SET UP INDMAP(*) TO MAP THE ENTRIES IN UPDATE COLUMNS */
/*           TO THEIR CORRESPONDING POSITIONS IN UPDATED COLUMNS, */
/*           RELATIVE THE THE BOTTOM OF EACH UPDATED COLUMN. */
/*           ----------------------------------------------------- */
	ldindx_(&jlen, &lindx[jxpnt], &indmap[1]);

/*           ----------------------------------------- */
/*           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ... */
/*           ----------------------------------------- */
	ksup = link[jsup];
L300:
	if (ksup > 0) {
	    nxksup = link[ksup];

/*               ------------------------------------------------
------- */
/*               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE. */

/*               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP. */
/*               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP. 
*/
/*               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FK
COL. */
/*               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN 
ACTIVE */
/*                           PORTION OF COLUMN FJCOL. */
/*               ------------------------------------------------
------- */
	    fkcol = xsuper[ksup];
	    nkcols = xsuper[ksup + 1] - fkcol;
	    klen = length[ksup];
	    kxpnt = xlindx[ksup + 1] - klen;

/*               ------------------------------------------- */
/*               PERFORM CMOD(JSUP,KSUP), WITH SPECIAL CASES */
/*               HANDLED DIFFERENTLY. */
/*               ------------------------------------------- */

	    if (klen != jlen) {

/*                   ----------------------------------------
--- */
/*                   SPARSE CMOD(JSUP,KSUP). */

/*                   NCOLUP ... NUMBER OF COLUMNS TO BE UPDATE
D. */
/*                   ----------------------------------------
--- */

		i__2 = klen - 1;
		for (i = 0; i <= i__2; ++i) {
		    nxtcol = lindx[kxpnt + i];
		    if (nxtcol > ljcol) {
			goto L500;
		    }
/* L400: */
		}
		i = klen;
L500:
		ncolup = i;

		if (nkcols == 1) {

/*                       --------------------------------
-------------- */
/*                       UPDATING TARGET SUPERNODE BY TRIV
IAL */
/*                       SUPERNODE (WITH ONE COLUMN). */

/*                       KLPNT  ...  POINTER TO FIRST NONZ
ERO IN ACTIVE */
/*                                   PORTION OF COLUMN FKC
OL. */
/*                       --------------------------------
-------------- */
		    klpnt = xlnz[fkcol + 1] - klen;
		    mmpyi_(&klen, &ncolup, &lindx[kxpnt], &lnz[klpnt], &xlnz[
			    1], &lnz[1], &indmap[1]);

		} else {

/*                       --------------------------------
------------ */
/*                       KFIRST ...  FIRST INDEX OF ACTIVE
 PORTION OF */
/*                                   SUPERNODE KSUP (FIRST
 COLUMN TO */
/*                                   BE UPDATED). */
/*                       KLAST  ...  LAST INDEX OF ACTIVE 
PORTION OF */
/*                                   SUPERNODE KSUP. */
/*                       --------------------------------
------------ */

		    kfirst = lindx[kxpnt];
		    klast = lindx[kxpnt + klen - 1];
		    inddif = indmap[kfirst] - indmap[klast];

		    if (inddif < klen) {

/*                           ------------------------
--------------- */
/*                           DENSE CMOD(JSUP,KSUP). */


/*                           ILPNT  ...  POINTER TO FI
RST NONZERO IN */
/*                                       COLUMN KFIRST
. */
/*                           ILEN   ...  LENGTH OF COL
UMN KFIRST. */
/*                           ------------------------
--------------- */
			ilpnt = xlnz[kfirst];
			ilen = xlnz[kfirst + 1] - ilpnt;
			mmpy_(&klen, &nkcols, &ncolup, &split[fkcol], &xlnz[
				fkcol], &lnz[1], &lnz[ilpnt], &ilen, mmpyn);

		    } else {

/*                           ------------------------
------- */
/*                           GENERAL SPARSE CMOD(JSUP,
KSUP). */
/*                           COMPUTE CMOD(JSUP,KSUP) U
PDATE */
/*                           IN WORK STORAGE. */
/*                           ------------------------
------- */
			store = klen * ncolup - ncolup * (ncolup - 1) / 2;
			if (store > *tmpsiz) {
			    *iflag = -2;
			    return 0;
			}
			mmpy_(&klen, &nkcols, &ncolup, &split[fkcol], &xlnz[
				fkcol], &lnz[1], &temp[1], &klen, mmpyn);
/*                           ------------------------
---------------- */
/*                           GATHER INDICES OF KSUP RE
LATIVE TO JSUP. */
/*                           ------------------------
---------------- */
			igathr_(&klen, &lindx[kxpnt], &indmap[1], &relind[1]);

/*                           ------------------------
-------------- */
/*                           INCORPORATE THE CMOD(JSUP
,KSUP) BLOCK */
/*                           UPDATE INTO THE TO APPROP
RIATE COLUMNS */
/*                           OF L. */
/*                           ------------------------
-------------- */
			assmb_(&klen, &ncolup, &temp[1], &relind[1], &xlnz[
				fjcol], &lnz[1], &jlen);

		    }

		}

	    } else {

/*                   ----------------------------------------
------ */
/*                   DENSE CMOD(JSUP,KSUP). */
/*                   JSUP AND KSUP HAVE IDENTICAL STRUCTURE. 
*/

/*                   JLPNT  ...  POINTER TO FIRST NONZERO IN C
OLUMN */
/*                               FJCOL. */
/*                   ----------------------------------------
------ */
		jlpnt = xlnz[fjcol];
		mmpy_(&klen, &nkcols, &njcols, &split[fkcol], &xlnz[fkcol], &
			lnz[1], &lnz[jlpnt], &jlen, mmpyn);
		ncolup = njcols;
		if (klen > njcols) {
		    nxtcol = lindx[jxpnt + njcols];
		}

	    }

/*               ------------------------------------------------ 
*/
/*               LINK KSUP INTO LINKED LIST OF THE NEXT SUPERNODE 
*/
/*               IT WILL UPDATE AND DECREMENT KSUP'S ACTIVE */
/*               LENGTH. */
/*               ------------------------------------------------ 
*/
	    if (klen > ncolup) {
		nxtsup = snode[nxtcol];
		link[ksup] = link[nxtsup];
		link[nxtsup] = ksup;
		length[ksup] = klen - ncolup;
	    } else {
		length[ksup] = 0;
	    }

/*               ------------------------------- */
/*               NEXT UPDATING SUPERNODE (KSUP). */
/*               ------------------------------- */
	    ksup = nxksup;
	    goto L300;

	}

/*           ---------------------------------------------- */
/*           APPLY PARTIAL CHOLESKY TO THE COLUMNS OF JSUP. */
/*           ---------------------------------------------- */
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCx
PC */
	chlsup_(&jlen, &njcols, &split[fjcol], &xlnz[fjcol], &lnz[1], &mxdiag,
		 &ntiny, iflag, mmpyn, smxpy);
	if (*iflag != 0) {
	    *iflag = -1;
	    return 0;
	}

/*           ----------------------------------------------- */
/*           INSERT JSUP INTO LINKED LIST OF FIRST SUPERNODE */
/*           IT WILL UPDATE. */
/*           ----------------------------------------------- */
	if (jlen > njcols) {
	    nxtcol = lindx[jxpnt + njcols];
	    nxtsup = snode[nxtcol];
	    link[jsup] = link[nxtsup];
	    link[nxtsup] = jsup;
	    length[jsup] = jlen - njcols;
	} else {
	    length[jsup] = 0;
	}

/* L600: */
    }

/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/
    if (ntiny != 0) {
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&ntiny, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/* SET IFLAG TO -1 TO INDICATE PRESENCE OF TINY DIAGONALS */

    if (ntiny != 0) {
	*iflag = -1;
    }
/* xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC 
*/
    return 0;
} /* blkfc2_ */

