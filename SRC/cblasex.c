#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]
#define DZ(I) dz[(I)-1]

void ddott(int *n, double *dx, double *dy, double *dz)
{
    /* Local variables */
    static int i, m;
    static int mp1;


    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DZ(i) = DX(i) * DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5)
    {
	DZ(i) = DX(i) * DY(i);
	DZ(i + 1) = DX(i + 1) * DY(i + 1);
	DZ(i + 2) = DX(i + 2) * DY(i + 2);
	DZ(i + 3) = DX(i + 3) * DY(i + 3);
	DZ(i + 4) = DX(i + 4) * DY(i + 4);
    }
    return;
} /* ddott */

void ddivt(int *n, double *dx, double *dy, double *dz)
{

    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DZ(i) = DX(i) / DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DZ(i) = DX(i) / DY(i);
	DZ(i + 1) = DX(i + 1) / DY(i + 1);
	DZ(i + 2) = DX(i + 2) / DY(i + 2);
	DZ(i + 3) = DX(i + 3) / DY(i + 3);
	DZ(i + 4) = DX(i + 4) / DY(i + 4);
    }
    return;
} /* ddivt */

void demy(int *n, double *dx, double *dy, double *dz)
{
    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DZ(i) = *dx - DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DZ(i) = *dx - DY(i);
	DZ(i + 1) = *dx - DY(i + 1);
	DZ(i + 2) = *dx - DY(i + 2);
	DZ(i + 3) = *dx - DY(i + 3);
	DZ(i + 4) = *dx - DY(i + 4);
    }
    return;
} /* demy */

void dety(int *n, double *dx, double *dy, double *dz)
{
    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DZ(i) = *dx * DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DZ(i) = *dx * DY(i);
	DZ(i + 1) = *dx * DY(i + 1);
	DZ(i + 2) = *dx * DY(i + 2);
	DZ(i + 3) = *dx * DY(i + 3);
	DZ(i + 4) = *dx * DY(i + 4);
    }
    return;
} /* dety */

void dxmy(int *n, double *dx, double *dy, double *dz)
{
    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DZ(i) = DX(i) - DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DZ(i) = DX(i) - DY(i);
	DZ(i + 1) = DX(i + 1) - DY(i + 1);
	DZ(i + 2) = DX(i + 2) - DY(i + 2);
	DZ(i + 3) = DX(i + 3) - DY(i + 3);
	DZ(i + 4) = DX(i + 4) - DY(i + 4);
    }
    return;
} /* dxmy */

void dxpy(int *n, double *dx, double *dy, double *dz)
{
    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;

    for (i = 1; i <= m; ++i)
	DZ(i) = DX(i) + DY(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DZ(i) = DX(i) + DY(i);
	DZ(i + 1) = DX(i + 1) + DY(i + 1);
	DZ(i + 2) = DX(i + 2) + DY(i + 2);
	DZ(i + 3) = DX(i + 3) + DY(i + 3);
	DZ(i + 4) = DX(i + 4) + DY(i + 4);
    }
    return;
} /* dxpy */

void dinv(int *n, double *dx, double *dy)
{
    /* Local variables */
    static int i, m;
    static int mp1;

    if (*n <= 0) return;

/*        code for both increments equal to 1
          clean-up loop */

    m = *n % 5;
    for (i = 1; i <= m; ++i)
	DY(i) = 1.0/DX(i);
    if (*n < 5) return;
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DY(i) = 1.0/DX(i);
	DY(i + 1) = 1.0/DX(i + 1);
	DY(i + 2) = 1.0/DX(i + 2);
	DY(i + 3) = 1.0/DX(i + 3);
	DY(i + 4) = 1.0/DX(i + 4);
    }
    return;
} /* dinv */

int idamax(int *n, double *dx)
{
    /* System generated locals */
    int ret_val;
    double d__1;

    /* Local variables */
    static double dmax__;
    static int i;

    if (*n < 1) return(-1);
    if (*n == 1) return(0);

    ret_val = 0;
    dmax__ = dx[0];
    for(i = 1; i < *n; ++i)
	if ((d__1 = dx[i]) > dmax__)
        {
	     ret_val = i;
	     dmax__ = d__1;
        }
    return(ret_val);
} /* idamax */

int idamin(int *n, double *dx)
{
    /* System generated locals */
    int ret_val;
    double d__1;

    /* Local variables */
    static double dmin__;
    static int i;

    if (*n < 1) return(-1);
    if (*n == 1) return(0);

    ret_val = 0;
    dmin__ = dx[0];
    for(i = 1; i < *n; ++i)
	if ((d__1 = dx[i]) < dmin__)
        {
	     ret_val = i;
	     dmin__ = d__1;
        }
    return(ret_val);
} /* idamin */

void daxpyx(int *n, double *da, double *dx, int *incx, double *dy, int *incy)
{
    /* Local variables */
    static int i, m, ix, iy, mp1;

    if (*n <= 0) return;
    if (*incx == 1 && *incy == 1) goto L20;

/*        code for unequal increments or equal increments
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) ix = (-(*n) + 1) * *incx + 1;
    if (*incy < 0) iy = (-(*n) + 1) * *incy + 1;
    for (i = 1; i <= *n; ++i) {
	DX(ix) = DY(iy) + *da * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1
          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) goto L40;
    for (i = 1; i <= m; ++i)
	DX(i) = DY(i) + (*da) * DX(i);
/* L30: */
    if (*n < 4) return;
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 4) {
	DX(i) = DY(i) + (*da) * DX(i);
	DX(i + 1) = DY(i + 1) + (*da) * DX(i + 1);
	DX(i + 2) = DY(i + 2) + (*da) * DX(i + 2);
	DX(i + 3) = DY(i + 3) + (*da) * DX(i + 3);
/* L50: */
    }
    return;
} /* daxpyx */

void daxpby(int *n, double *da, double *dx, int *incx, double *db, double *dy,
            int *incy)
{
    /* Local variables */
    static int i, m, ix, iy, mp1;

    if (*n <= 0) return;
    if (*incx == 1 && *incy == 1) goto L20;

/*        code for unequal increments or equal increments
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) ix = (-(*n) + 1) * *incx + 1;
    if (*incy < 0) iy = (-(*n) + 1) * *incy + 1;
    for (i = 1; i <= *n; ++i) {
	DX(ix) = (*db) * DY(iy) + (*da) * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1
          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) goto L40;
    for (i = 1; i <= m; ++i)
	DX(i) = (*db) * DY(i) + (*da) * DX(i);
/* L30: */
    if (*n < 4) return;
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 4) {
	DX(i) = (*db) * DY(i) + (*da) * DX(i);
	DX(i + 1) = (*db) * DY(i + 1) + (*da) * DX(i + 1);
	DX(i + 2) = (*db) * DY(i + 2) + (*da) * DX(i + 2);
	DX(i + 3) = (*db) * DY(i + 3) + (*da) * DX(i + 3);
/* L50: */
    }
    return;
} /* daxpby */
