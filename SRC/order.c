/*****************************************************************************

        Source      : order.c
	Author(s)   : January 12, 1991. Eva K. Lee

	Last Update : February 15, 1991. Eva


******************************************************************************/
#include <sys/types.h>
#include <sys/uio.h>
#include <errno.h>
#include <malloc.h>
#include <sys/time.h>
#include <sys/times.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))

#define  LPERR_ORDER              0


#define MALLOC(x,num,type) do {                                      \
        if ((x = (type *) malloc ((num) * sizeof (type))) == NULL) { \
                printf ("Out of memory?\n");        \
                return(LPERR_ORDER);                   \
        }} while (0)
/* A */

#define FREEL(x) do {                                                 \
        free ((char *) (x));                                         \
        x = NULL;                                                    \
} while (0)

#ifdef PROTOTYPE_MAX
int order (int rows, int *matbeg, int *matind, int *pivot, int *ipivot,
           int maxnzcnt, int *degree)
#else
int order (rows, matbeg, matind, pivot, ipivot, maxnzcnt, degree)
int        rows, *matbeg, *matind, *pivot, *ipivot, maxnzcnt, *degree;
#endif
{
#define  beg(i)       matbeg[i-1] + 1
#define  ind(i)       matind[i-1] + 1
#define  tpiv(i)      pivot[i-1]
#define  tipiv(i)     ipivot[i-1]

    int     min_degree, thresh_degree, u, v, v1, v2;
    int     vertex0, next0, vertex1, next1, vertex2, next2;
    int     nzcnt, tmpv;
    int     i, j, k, begin, end, deg1, kmax, *ip;
    int     flag = 0;
/*
    int     *degree = (int *) NULL;
*/
    int     *adj    = (int *) NULL;
/*
    MALLOC (degree, maxnzcnt, int);
    MALLOC (adj, maxnzcnt+1, int);
*/
    adj= degree + maxnzcnt;

    for (i = 0; i < rows; i++) {
        pivot[i]++;
        ipivot[i]++;
    }

    v = 1;
    for (i = rows; i <= maxnzcnt; i++) {
        adj[i] = i+1;
    }
    adj[maxnzcnt] = 0;
    nzcnt = 1;
    for (k = 1; k <=rows; k++) {
        v = tpiv(k);
        tipiv(v) = k;
        degree[k] = rows+1;
        adj[k] = k;
    }
    nzcnt = nzcnt + rows;

    for (v = 1; v <= rows; v++) {
        begin = beg(v);
        end = beg(v+1) - 1;
        if  (begin > end)  goto ERR1;
        for (j = begin; j <= end; j++) {
            u = ind(j);
            if  (u != v) {
                vertex0 =  v;
                do {
                    next0 = vertex0;
                    vertex0 = adj[next0];
                } while (degree[vertex0] < u);
                if  (degree[vertex0] > u) {
                    degree[v]++;
                    if  (nzcnt == 0)  goto ERR2;
                    vertex0 = nzcnt;
                    nzcnt = adj[nzcnt];
                    degree[vertex0] =  u;
                    adj[vertex0] = adj[next0];
                    adj[next0] = vertex0;
                    vertex1 =  u;
                    do {
                        next1 = vertex1;
                        vertex1 = adj[next1];
                    } while (degree[vertex1] < v);
                    if  (degree[vertex1] > v) {
                        degree[u]++;
                        if  (nzcnt == 0)  goto ERR2;
                        vertex1 = nzcnt;
                        nzcnt = adj[nzcnt];
                        degree[vertex1] = v;
                        adj[vertex1] = adj[next1];
                        adj[next1] = vertex1;
                    }
                }
            }
        }
    }
    i = 0;
    begin = 1;
    thresh_degree = 0;
    min_degree = rows+rows;

    do {
        do {
            begin = MAX (begin, i+1);
            for (j = begin; j <= rows; j++) {
                v1 = tpiv(j);
                if (degree[v1] <= thresh_degree)  goto MARK11;
                min_degree = MIN (min_degree, degree[v1]);
            }
            begin = 1;
            thresh_degree = min_degree;
            min_degree = rows+rows;
        } while (1);
MARK11: begin = j;
        i = i+1;
        u = tpiv(i);
        tpiv(j) =  u;
        tipiv(u) = j;
        tpiv(i) = v1;
        tipiv(v1) = i;
        deg1 = 1;
        next2 = v1;
        kmax = (degree[v1] - deg1) - rows;
        if  (kmax > 0)  {
            for (k = 1; k <= kmax; k++) {
                for (;;) {
                    vertex2 = next2;
                    next2 = adj[vertex2];
                    v = degree[next2];
                    if  (tipiv(v) <= i)  {
                        adj[vertex2] = adj[next2];
                        adj[next2] = nzcnt;
                        nzcnt = next2;
                        next2 = vertex2;
                    } else {
                        break;
                    }
                }
            }
        }

        next2 = v1;
        kmax = (degree[v1] - deg1) - rows;
        if  (kmax > 0) {
            for (k = 1; k <= kmax; k++) {
                vertex2 = next2;
                next2 = adj[vertex2];
                v = degree[next2];
                vertex0 = v;
                next1 = v1;
                end = (degree[v1] - deg1) - rows;
                if  (end > 0) {
                    for (j=1; j <= end; j++) {
                        next1 = adj[next1];
                        u = degree[next1];
                        if  ( u != v) {
                            for (;;) {
                                next0 = vertex0;
                                vertex0 = adj[next0];
                                v2 = degree[vertex0];
                                if  (u > v2) {
                                    if  (tipiv(v2) <= i) {
                                        adj[next0] = adj[vertex0];
                                        adj[vertex0] = nzcnt;
                                        nzcnt = vertex0;
                                        vertex0 = next0;
                                    }
                                } else {
                                    break;
                                }
                            }
                            if  ( u != v2) {
                                degree[v] = degree[v] + 1;
                                if  (nzcnt == 0)  goto ERR2;
                                    vertex0 = nzcnt;
                                    nzcnt = adj[nzcnt];
                                    degree[vertex0] =  u;
                                    adj[vertex0] = adj[next0];
                                    adj[next0] = vertex0;
                            }
                        }
                    }
                }
                if  (degree[v] <= degree[v1]) {
                    i = i+1;
                    j = tipiv(v);
                    u = tpiv(i);
                    tpiv(j) =  u;
                    tipiv(u) = j;
                    tpiv(i) = v;
                    tipiv(v) = i;
                    deg1 = deg1 + 1;
                    tmpv = adj[v];
                    adj[v] = nzcnt;
                    nzcnt = tmpv;
                    adj[vertex2] = adj[next2];
                    adj[next2] = nzcnt;
                    nzcnt = next2;
                    next2 = vertex2;
                }
            }
        }

        vertex2 = v1;
        kmax = (degree[v1] - deg1) - rows;
        if  (kmax > 0) {
            for (k=1; k <= kmax; k++) {
                vertex2 = adj[vertex2];
                v = degree[vertex2];
                degree[v] = degree[v] - deg1;
                if (degree[v] < min_degree) {
                    if (degree[v] <= thresh_degree) {
                        min_degree = thresh_degree;
                        thresh_degree = degree[v];
                        begin = MIN (begin, tipiv(v));
                    } else {
                        min_degree = degree[v];
                    }
                }
            }
        }
        tmpv = adj[v1];
        adj[v1] = nzcnt;
        nzcnt = tmpv;
    } while (i < rows);

    for (i = 0; i < rows; i++) {
        pivot[i]--;
        ipivot[i]--;
    }
/*
    FREEL (degree);
    FREEL (adj);
*/
    flag = 1;
    return (flag);

ERR1:   printf ("\nNull rows in A\n");
        error(1);
        return (LPERR_ORDER);
ERR2:   printf ("\nInsufficient Storage %d\n",v);
        //error(1);
        return (LPERR_ORDER);
} /* END ORDER */


/*****************************************************************************

        Source      : cheap_order.c
        Author(s)   : January 20, 1991. Eva K. Lee

        Last Update : February 27, 1991. Eva


******************************************************************************/

#ifdef PROTOTYPE_MAX
int cheap_order (int length, int *rmatbeg, int *rmatind, int *pivot, int *ipivot)
#else
int cheap_order (length, rmatbeg, rmatind, pivot, ipivot)
int     length, *rmatbeg, *rmatind, *pivot, *ipivot;
#endif
{
    int    *dghead;
    int    *linkfd;
    int    order = 0;
    int    i, rowindex, degree;

    MALLOC (dghead, length+1, int);
    MALLOC (linkfd, length, int);

    for (rowindex = 0; rowindex <= length; rowindex++) {
        dghead[rowindex] = -1;
    }

    for (rowindex = 0; rowindex < length; rowindex++) {
        degree = rmatbeg[rowindex+1] - rmatbeg[rowindex];
        linkfd[rowindex] = dghead[degree];
        dghead[degree] = rowindex;
    }

    for (degree = 0; degree <= length; degree++) {
        rowindex = dghead[degree];
        while (rowindex != -1) {
            ipivot[rowindex] = order;
            order++;
            rowindex = linkfd[rowindex];
        }
    }

    for (i = 0; i < length; i++) {
        pivot[ipivot[i]] = i;
    }

    FREEL (dghead);
    FREEL (linkfd);

    return (1);

} /* END CHEAP_ORDER */
