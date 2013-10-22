/* gtimer.f -- translated by f2c (version of 25 March 1992  12:58:56).
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

doublereal gtimer_()
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern doublereal etime_();
    static real vec[2];

/*       -------------------------- */
/*       FOR IBM RS/6000 ... */
/*       INTEGER     MCLOCK */
/*       GTIMER = MCLOCK()/100.0 */
/*       -------------------------- */
/*       FOR MOST BERKELEY UNIX ... */
    ret_val = etime_(vec);
/*       -------------------------- */
/*       FOR CRAY ... */
/*       REAL        SECOND */
/*       GTIMER = SECOND() */
/*       -------------------------- */
    return ret_val;
} /* gtimer_ */

