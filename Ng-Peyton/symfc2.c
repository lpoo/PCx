/* symfc2.f -- translated by f2c (version of 25 March 1992  12:58:56).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include <f2c.h>

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.4 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFC2 ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A */
/*       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS */
/*       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE */
/*       ASSOCIATED WITH THE FACTOR MATRIX L. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                               1 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE */
/*                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST. */
/*       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                           CURRENT LINKED LIST OF MERGED INDICES (THE */
/*                           "REACH" SET). */
/*       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES */
/*                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S */
/*                           INDEX SET. */

/* *********************************************************************** */

/* Subroutine */ int symfc2_(neqns, adjlen, xadj, adjncy, perm, invp, colcnt, 
	nsuper, xsuper, snode, nofsub, xlindx, lindx, xlnz, mrglnk, rchlnk, 
	marker, flag_)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *colcnt, *nsuper, *
	xsuper, *snode, *nofsub, *xlindx, *lindx, *xlnz, *mrglnk, *rchlnk, *
	marker, *flag_;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer head, node, tail, pcol, newi, jptr, kptr, jsup, ksup, psup,
	     i, nzbeg, nzend, width, nexti, point, jnzbeg, knzbeg, length, 
	    jnzend, jwidth, fstcol, knzend, lstcol, knz;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    --marker;
    --mrglnk;
    --xlnz;
    --lindx;
    --xlindx;
    --snode;
    --xsuper;
    --colcnt;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    *flag_ = 0;
    if (*neqns <= 0) {
	return 0;
    }

/*       --------------------------------------------------- */
/*       INITIALIZATIONS ... */
/*           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX. */
/*           TAIL   : END OF LIST INDICATOR */
/*                    (IN RCHLNK(*), NOT MRGLNK(*)). */
/*           MRGLNK : CREATE EMPTY LISTS. */
/*           MARKER : "UNMARK" THE INDICES. */
/*       --------------------------------------------------- */
    nzend = 0;
    head = 0;
    tail = *neqns + 1;
    point = 1;
    i__1 = *neqns;
    for (i = 1; i <= i__1; ++i) {
	marker[i] = 0;
	xlnz[i] = point;
	point += colcnt[i];
/* L50: */
    }
    xlnz[*neqns + 1] = point;
    point = 1;
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
	mrglnk[ksup] = 0;
	fstcol = xsuper[ksup];
	xlindx[ksup] = point;
	point += colcnt[fstcol];
/* L100: */
    }
    xlindx[*nsuper + 1] = point;

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {

/*           ---------------------------------------------------------
 */
/*           INITIALIZATIONS ... */
/*               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP. */
/*               LSTCOL : LAST COLUMN OF SUPERNODE KSUP. */
/*               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL.
 */
/*               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL. */
/*           ---------------------------------------------------------
 */
	fstcol = xsuper[ksup];
	lstcol = xsuper[ksup + 1] - 1;
	width = lstcol - fstcol + 1;
	length = colcnt[fstcol];
	knz = 0;
	rchlnk[head] = tail;
	jsup = mrglnk[ksup];

/*           ------------------------------------------------- */
/*           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ... */
/*           ------------------------------------------------- */
	if (jsup > 0) {
/*               --------------------------------------------- */
/*               COPY THE INDICES OF THE FIRST CHILD JSUP INTO */
/*               THE LINKED LIST, AND MARK EACH WITH THE VALUE */
/*               KSUP. */
/*               --------------------------------------------- */
	    jwidth = xsuper[jsup + 1] - xsuper[jsup];
	    jnzbeg = xlindx[jsup] + jwidth;
	    jnzend = xlindx[jsup + 1] - 1;
	    i__2 = jnzbeg;
	    for (jptr = jnzend; jptr >= i__2; --jptr) {
		newi = lindx[jptr];
		++knz;
		marker[newi] = ksup;
		rchlnk[newi] = rchlnk[head];
		rchlnk[head] = newi;
/* L200: */
	    }
/*               ------------------------------------------ */
/*               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ... */
/*               ------------------------------------------ */
	    jsup = mrglnk[jsup];
L300:
	    if (jsup != 0 && knz < length) {
/*                   ---------------------------------------- 
*/
/*                   MERGE THE INDICES OF JSUP INTO THE LIST, 
*/
/*                   AND MARK NEW INDICES WITH VALUE KSUP. */
/*                   ---------------------------------------- 
*/
		jwidth = xsuper[jsup + 1] - xsuper[jsup];
		jnzbeg = xlindx[jsup] + jwidth;
		jnzend = xlindx[jsup + 1] - 1;
		nexti = head;
		i__2 = jnzend;
		for (jptr = jnzbeg; jptr <= i__2; ++jptr) {
		    newi = lindx[jptr];
L400:
		    i = nexti;
		    nexti = rchlnk[i];
		    if (newi > nexti) {
			goto L400;
		    }
		    if (newi < nexti) {
			++knz;
			rchlnk[i] = newi;
			rchlnk[newi] = nexti;
			marker[newi] = ksup;
			nexti = newi;
		    }
/* L500: */
		}
		jsup = mrglnk[jsup];
		goto L300;
	    }
	}
/*           --------------------------------------------------- */
/*           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET. */
/*           "SORT" ITS STRUCTURE INTO THE LINKED LIST, */
/*           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE */
/*           LIST. */
/*           --------------------------------------------------- */
	if (knz < length) {
	    node = perm[fstcol];
	    knzbeg = xadj[node];
	    knzend = xadj[node + 1] - 1;
	    i__2 = knzend;
	    for (kptr = knzbeg; kptr <= i__2; ++kptr) {
		newi = adjncy[kptr];
		newi = invp[newi];
		if (newi > fstcol && marker[newi] != ksup) {
/*                       -------------------------------- 
*/
/*                       POSITION AND INSERT NEWI IN LIST 
*/
/*                       AND MARK IT WITH KCOL. */
/*                       -------------------------------- 
*/
		    nexti = head;
L600:
		    i = nexti;
		    nexti = rchlnk[i];
		    if (newi > nexti) {
			goto L600;
		    }
		    ++knz;
		    rchlnk[i] = newi;
		    rchlnk[newi] = nexti;
		    marker[newi] = ksup;
		}
/* L700: */
	    }
	}
/*           --------------------------------------------------------
---- */
/*           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LI
ST. */
/*           --------------------------------------------------------
---- */
	if (rchlnk[head] != fstcol) {
	    rchlnk[fstcol] = rchlnk[head];
	    rchlnk[head] = fstcol;
	    ++knz;
	}

/*           -------------------------------------------- */
/*           COPY INDICES FROM LINKED LIST INTO LINDX(*). */
/*           -------------------------------------------- */
	nzbeg = nzend + 1;
	nzend += knz;
	if (nzend + 1 != xlindx[ksup + 1]) {
	    goto L8000;
	}
	i = head;
	i__2 = nzend;
	for (kptr = nzbeg; kptr <= i__2; ++kptr) {
	    i = rchlnk[i];
	    lindx[kptr] = i;
/* L800: */
	}

/*           --------------------------------------------------- */
/*           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S */
/*           "MERGE" LIST. */
/*           --------------------------------------------------- */
	if (length > width) {
	    pcol = lindx[xlindx[ksup] + width];
	    psup = snode[pcol];
	    mrglnk[ksup] = mrglnk[psup];
	    mrglnk[psup] = ksup;
	}

/* L1000: */
    }

    return 0;

/*       ----------------------------------------------- */
/*       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED. */
/*       ----------------------------------------------- */
L8000:
    *flag_ = -2;
    return 0;

} /* symfc2_ */

