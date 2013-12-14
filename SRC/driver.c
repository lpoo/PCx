#define MAIN
#include "main.h"
#include "timer.h"

/* This file contains the main driver, the error handling
   functions and basic functions */

static TIMER xtime;

main(int argc, char *argv[])
{
/* Main driver for the interior point method.

   It solves

   (P) min c'x          (D) max b'y - u's

       Ax + w = b           A'y - s + z = c
        x + v = u

        x, w, v >=0               s, z >= 0

   s is the slack for primal inequalities constraints.
   z is the slack for dual inequalities constraints.

   Declare local variables */

   matrix *A, *AA, *schur, *L;
   char *rowt, *vart;
   double *b, *c, *u, *l, *r;
   double *C, *D, *E, *G, *H;
   double *x, *y, *z, *w, *s, *v;
   double *rx, *ry, *rz, *rw, *rs, *rv;
   double *dx, *dy, *dz, *dw, *ds, *dv;
   double *rhs, *rxy, *dwork;
   double cx, by, gap, gapp, cnt= 0.0;
   double phi, mu, passp, passd;
   double tol= 1.4901161e-8;
   int k, i, mn, nm, nz, np, box, ine, slack, nb, fil, maxit, *pos, stat, dual;
   int MD= TRUE, maxnz, *pivot, *ipivot, *diag, *iwork;
   float ops, cgops;

   options();

/* Allocate space to matrices' pointers */

   if((A= (matrix *) malloc(4*sizeof(matrix))) == NULL) error(4);
   AA=    A + 1;
   schur= A + 2;
   L=     A + 3;

   if(SAVE)
   {
      i= load(&A,&b,&c,&l,&u,&rowt,&vart,&cnt,&nm,&slack,&np,&phi,&AA,&schur,&L,
              &pos,&pivot,&k,&stat,&maxnz,&OPS,&x,&y,&z,&w,&v,&s,
              &xtime.UserTime,&xtime.SystemTime);
      if(i)
      {
         if(stat == 0) stat= -TRUE;
         mn= A->m + A->n;
         if((D = (double *) malloc((unsigned)
            (9*A->n + 8*A->m)*sizeof(double))) == NULL) error(4);
         E=  D  + A->n;

         dx= E  + A->m;
         dy= dx + A->n;
         dz= dy + A->m;
         dw= dz + A->n;
         dv= dw + A->m;
         ds= dv + A->n;

         ry= E  + A->m;
         rx= ry + A->n;
         rz= rx + A->m;
         rw= rz + A->n;
         rv= rw + A->m;
         rs= rv + A->n;

         rhs= ry;

         C=  ds + A->n;
         H=   C + A->n;
         rxy= H + A->n;
         dwork= rxy + mn;
         iwork= (int *) dwork;

         if(LU)
         {
            schur= NULL;
            if((pivot= (int *) malloc((3*L->m+1)*sizeof(int))) == NULL)
                error(4);
            diag= pivot + L->m;
            L->col= diag + L->m;
            if((L->row= (int *) malloc((maxnz)*sizeof(int))) == NULL) error(4);
            if((L->val= (double *) malloc(maxnz*sizeof(double))) == NULL)
                error(4);
         }
         RestartClock(&xtime);
         goto ITERATE;
      }
   }

/* Read problem using cplex callable library */

   input(&A,&b,&c,&u,&l,&r,&rowt,argc,argv);
/*
   ResetClock(&xtime);
   StartClock(&xtime);
*/
/* Determine variables type and shift lower bound to zero */

   if((vart= (char *) malloc((unsigned) (A->n+1))) == NULL) error(4);
   box= bounds(A,b,c,u,l,r,&cnt,rowt,vart);

/* preprocess problem */

   mn= A->n;
   nz= A->col[A->n];
   printf("%d X %d\n",A->m,A->n);
   fflush(stdout);
   preprocess(A,b,c,u,r,&cnt,rowt,vart, &box, (int *) l);
   ine= range(&A,b,&c,&u,r,mn,nz,&box,&rowt,&vart);

/* scale problem */
   if(r == NULL) r= (double *) malloc(A->m*sizeof(double));
   scale(A,b,c,u,l,r,vart);

   free((char *) r);

   if(SLACK)
   {
      slack= addslacks(&A,&c,&u,mn,nz,rowt,&vart);
      ine= 0;
   }
   else slack= 0;

   printf("%d X %d\n",A->m,A->n);
   printf("nzA= %d\n",A->col[A->n]);
   printf("bounded= %d\n",box);
   printf("inequalities= %d\n",ine);
   printf("slacks= %d\n",slack);
   fflush(stdout);

   if(A->m > A->n)
      if(ine == A->m && box == 0)
      {
         dual= TRUE;
         puts("Solving dual");
         change(A,b,c,rowt,vart);
         ine= A->m;
         error(9);
      }
      else error(9);
   else dual= FALSE;

   mn= A->m + A->n;

/* Allocate memory for vectors and work space*/
   maxit= MIN(200,A->m);

   if(AUG)
   {
      if(PCG)
      {
         if((D = (double *) malloc((unsigned)
         (16*A->n + 11*A->m)*sizeof(double))) == NULL) error(4);
      }
      else
      {
         if((D = (double *) malloc((unsigned)
           (12*A->m + 18*A->n + (maxit+A->n+A->m+8)*maxit)*sizeof(double)))
            == NULL) error(4);
      }
   }
   else
   {
      if(PCG)
      {
         if((D = (double *) malloc((unsigned)
            (14*A->n + 14*A->m)*sizeof(double))) == NULL) error(4);
      }
      else
      {
         if((D = (double *) malloc((unsigned)
            (8*A->m + 16*A->n + (maxit + A->m + 8)*maxit)*sizeof(double)))
             == NULL) error(4);
      }
   }

   E=  D  + A->n;
   x=  E  + A->m;
   y=  x  + A->n;
   z=  y  + A->m;
   w=  z  + A->n;
   v=  w  + A->m;
   s=  v  + A->n;

   dx= s  + A->n;
   dy= dx + A->n;
   dz= dy + A->m;
   dw= dz + A->n;
   dv= dw + A->m;
   ds= dv + A->n;

/* Atention: residual vectors (r) share space with directions vectors (d) */
   ry= s  + A->n;
   rx= ry + A->n;
   rz= rx + A->m;
   rw= rz + A->n;
   rv= rw + A->m;
   rs= rv + A->n;

/* rhs shares space with ry and rx */
   rhs= ry;

   C=  ds + A->n;
   H=   C + A->n;
   G=   H + A->n;
   rxy= G + A->n;
   dwork= rxy + mn;
   iwork= (int *) dwork;

/* Determine constraint type and change those from <= to >= */
   if(!SLACK) constr(A,b,rowt);

   L->m= L->n= A->m;
   if(LU)
   {
      schur= NULL;
      nz= A->col[A->n];
      if((pivot= (int *) malloc((3*A->m+A->n+1)*sizeof(int))) == NULL) error(4);
      diag= pivot + A->m;
      L->col= diag + A->m;
      pos= L->col + A->m + 1;
      maxnz= FILL*nz;
      if((L->row= (int *) malloc(maxnz*sizeof(int))) == NULL) error(4);
      if((L->val= (double *) malloc(maxnz*sizeof(double))) == NULL) error(4);
      for(i= A->col[k= 0]; k < A->n; k++)
      {
          pos[k]= k;
          for(G[k]= 0.0; i < A->col[k+1]; i++) G[k] -= fabs(A->val[i]);
      }
      colorder(A,pos,G,dwork,diag,L->row,(int *) L->val);
      stat= lufact(A,&L,&maxnz,pivot,pos,diag,iwork,H,&ops,FALSE,TRUE);
      if(stat)
      {
         printf("rank def, rank= %d\n",stat);
         nz= L->col[stat];
         for(k= stat-1; k >= 0; k--)
             if(iwork[k] != L->col[k+1])
             {
                np= L->col[k+1] - iwork[k];
                shrink(L,iwork[k],L->col[k+1],k,&nz);
                for(i= k+1; i < stat; i++) diag[i] -= np;
             }
         for(i= k= 0; k < A->m; k++)
             if(pivot[k] == A->m)
             {
                ldrow(A,k);
                iwork[i++]= k;
             }
         delrow(A,b,r,rowt,pivot,i,iwork,&iwork[i+1]);
         L->m= L->n= A->m;
         printf("%d X %d\n",A->m,A->n);
         nz= A->col[A->n];
         printf("nzA= %d\n",nz);
      }
      printf("nzLU= %d\n",L->col[L->m]);
      if(REFAC && L->col[L->m] > 4*L->m)
      {
         fflush(stdout);
         degree(A,pos,iwork,NULL,A->m);
         i= btb(A,pos,iwork,L->row,A->m,diag,maxnz);
         if(i)
            if(ORDER)
            {
               for(i= 0; i <A->m; i++) diag[i]= pivot[i]= i;
               i= order(A->m,iwork,L->row,diag,pivot,maxnz-1,(int *) L->val);
            }
            else i= cheap_order(A->m,iwork,L->row,diag,pivot);
         else puts("Insufficient Storage");
         if(i)
         {
            memcpy(pivot,pos,A->m*sizeof(int));
            for(i= 0; i < A->m; i++) pos[i]= pivot[diag[i]];
         }
         else MD= FALSE;
         lufact(A,&L,&maxnz,pivot,pos,diag,iwork,H,&ops,INC,FALSE);
         printf("nzLU= %d\n",L->col[L->m]);
      }
   }
   else
   {
      schur->m= schur->n= A->m;

/* Inicialize pointers for AA' */
      create_AA(A,&AA);

      if(dual) schur->col= (int *) malloc((unsigned) (A->m+1));
      else
      {
         schur->col= A->col;

      /* Change A from column oriented scheme into row oriented scheme */
         conv_row(&A);
      }

/* Determine nonzero structure of Schur complement A'A */

      create_schur(A,&schur,AA->n,iwork);
      nz= schur->col[schur->m];
      printf("nzL= %d\n",nz);

/* Allocate double space for the (in)complete Cholesky decomposition */

      if(INC)
      {
/*    This branch contains the code for incomplete Cholesky to be used with
      iterative methods */

         fil= 0;
         if((schur->val= (double *) malloc((unsigned) nz*sizeof(double)))
            == NULL) error(4);

/*    L and Schur share the same space */

         L->col= schur->row;
         L->row= schur->col;
         L->val= schur->val;
      }
      else
      {
/*    This branch contains the code for complete Cholesky to be used with
      direct methods methods. */

/* Inicialize pointers and allocate integer space for the
   Cholesky decomposition of the Schur complement  */

         if((L->row= (int *) malloc((unsigned)(nz + 3*L->m + 1)*sizeof(int))) ==
             NULL) error(4);
         pos= L->col= L->row + L->m + 1;

         pivot= pos + nz;
         ipivot= pivot + L->m;

         memcpy(L->col,schur->row,nz*sizeof(int));
         memcpy(L->row,schur->col,(L->m+1)*sizeof(int));

         if(ORDER)
            for(k= 0; k < L->m; k++) pivot[k]= ipivot[k]= k;

/*    reorder A's rows */

         reorder(L, pivot, ipivot,iwork);
         transf(A,b,pivot,rowt,iwork);
         memcpy(schur->row,L->col,nz*sizeof(int));
         memcpy(schur->col,L->row,(L->m+1)*sizeof(int));

/*    Compute the number of fill-ins */

         fil= fill_in(L,iwork);
         printf("fill-in= %d\n",fil);

         if((L->col= (int *) malloc((unsigned) (nz + fil)*sizeof(int))) == NULL)
             error(4);

/*    Determine L's structure */

         create_L(schur,L,iwork);
         statl(schur,L,pos);

         if((schur->val= (double *) malloc((unsigned)
            (2*nz + fil)*sizeof(double))) == NULL) error(4);
         L->val= schur->val + nz;

         ops= L->m;
         for(k= 0; k < L->m; k++)
             ops += (L->row[k+1]-L->row[k])*(L->row[k+1]-L->row[k]);
         printf("Flops for Cholesky= %.0f\n",ops);
      }

/* Determine the nonzero structure of each pair A_ik^' A_kj */

      compute_AA(A,&AA,iwork);
printf("Flops for Schur complement= %.0f\n",2.0*AA->row[AA->m]);
ops += 2*(2*(2*L->row[L->m] - L->m));
printf("Flops for all= %.0f\n",ops+2.0*AA->row[AA->m]);
   }
   np= A->n + box + ine;

/* Initialize diagonal matrices */

   for(k= 0; k < A->n; k++) D[k]= 1.0;
   for(k= 0; k < A->m; k++) E[k]= 0.0;

   if(!LU)
   {
/* Compute Schur S= ADA' + E and */
/* compute incomplete factorization for S */

      compute_schur(D,E,AA,schur,FALSE,rowt);
/*
ResetClock(&xtime);
StartClock(&xtime);
StopClock(&xtime);
error(0);
*/
      cholesky(schur,L,pos,dwork,(int *) (dwork + schur->m));
   }

   memset(v,0,A->n*sizeof(double));
   Ax(A,u,y);

   for(k= 0; k < A->m; k++)
       y[k]= 2.0*b[k] - y[k];
   memcpy(w,y,A->m*sizeof(double));
   nm= MIN(A->m, A->n - A->m);
   nm= A->m;

/* Solve positive definite system to get initial solution */

   if(AUG)
   {
      memset(x,0,A->n*sizeof(double));
      memset(z,0,A->n*sizeof(double));
      stat= (PCG) ? pcg(A,mn,D,D,x,z,L,pivot,diag,2*A->m,tol,dwork)
                  : minres(A,mn,D,x,z,L,mn,tol,dwork);
   }
   else
      stat= (PCG) ? pcg(A,A->m,D,D,y,w,L,pivot,diag,A->m,tol,dwork)
                  : minres(A,A->m,D,y,w,L,A->m,tol,dwork);

   if(stat < 0) error(3);
/*
cgops= 4*A->col[A->n] + A->n + 4*L->col[L->m] - L->m;
printf("pcgfl= %.0f\n",stat*(cgops+10.0*L->m+2) + cgops +4*L->m - 1);
if(stat*(cgops+10.0*L->m+2) + cgops +4*L->m - 1 < 5*ops) stat= FALSE;
fflush(stdout);
*/
   if(stat*20 < nm) stat= FALSE;

/* Compute initial solution */

   initial_sol(A,b,c,u,x,y,z,w,v,s,rowt,vart);
   memset(dw,0,A->m*sizeof(double));
   memset(dv,0,A->n*sizeof(double));
   memset(ds,0,A->n*sizeof(double));

/*   if(np < 1000) phi= 1.0/ ((double)(np * np));
   else*/ phi= 1.0/ ((double)(np * sqrt((double) np)));


   k= 0;
   if(SAVE)
   {
      StopClock(&xtime);
      savefix(A,b,c,l,u,rowt,vart,cnt,nm,slack,np,phi,AA,schur,L,pos,pivot);
      savevar(A->m,A->n,k,stat,maxnz,OPS,x,y,z,w,v,s,pos,
              xtime.UserTime*CLK_TCK,xtime.SystemTime*CLK_TCK);
      StartClock(&xtime);
   }

/* Iterate */
ITERATE:

   for( ; !converged(A,b,c,u,x,y,z,w,v,s,&cx,&by,&gap,rx,ry,rv,tol,rowt,vart,&i)
        ; k++)
   {
      printf("gap= %e\n",gap);
      if(k == 150) error(8);

      if(PC) mu= 0.0; else mu= compute_mu(gap,phi,i);

      memcpy(rxy,ry,mn*sizeof(double));

/*    Form linear system |-D A'| dx = rhs
                         | A E | dy        */
      form_ls(A,x,y,z,w,v,s,rx,ry,rz,rw,rv,rs,D,E,mu,dwork,rowt,vart);

      if(LU)
      {
         if(stat)
         {
            if(stat > 0)
            {
/*
ddivt(N,G,D,H);
               colorder(A,pos,H,dwork,diag,L->row,(int *) L->val);
*/
               colorder(A,pos,D,dwork,diag,L->row,(int *) L->val);
               lufact(A,&L,&maxnz,pivot,pos,diag,iwork,H,&ops,INC,TRUE);
            }
            else
            {
               degree(A,pos,iwork,NULL,A->m);
               lufact(A,&L,&maxnz,pivot,pos,diag,iwork,H,&ops,INC,FALSE);
            }
            printf("nzLU= %d\n",L->col[L->m]);
            fflush(stdout);
            if(REFAC && L->col[L->m] > 4*L->m)
            {
               degree(A,pos,iwork,NULL,A->m);
               if(MD) i= btb(A,pos,iwork,L->row,A->m,diag,maxnz); else i= FALSE;
               if(i)
                if(ORDER)
                {
                   for(i= 0; i <A->m; i++) diag[i]= pivot[i]= i;
                   i= order(A->m,iwork,L->row,diag,pivot,maxnz-1,(int *)L->val);
                }
                else i= cheap_order(A->m,iwork,L->row,diag,pivot);
               else puts("Insufficient Storage");
               if(i)
               {
                  memcpy(pivot,pos,A->m*sizeof(int));
                  for(i= 0; i < A->m; i++) pos[i]= pivot[diag[i]];
               }
               else MD= FALSE;
               lufact(A,&L,&maxnz,pivot,pos,diag,iwork,H,&ops,INC,FALSE);
               printf("nzLU= %d\n",L->col[L->m]);
               fflush(stdout);
             }
         }
         gather(A->m,pos,H,D);
      }

/*    Solve linear system */
      if(IMP)
         stat= solve_ls(A,D,E,H,AA,schur,L,pivot,diag,rhs,pos,tol,dwork,rowt);
      else
      {
/*
   compute_schur(D,E,AA,schur,TRUE,rowt);
   cholesky(schur,L,pos,dwork,(int *) (dwork + schur->m));
*/
         memcpy(dwork,rhs,mn*sizeof(double));
         stat= minres(A,mn,D,rhs,dwork,L,mn,tol,&dwork[mn]);
      }
      if(stat < 0) error(3);
      if(stat*20 < nm) stat= FALSE;
/*
cgops= 4*A->col[A->n] + A->n + 4*L->col[L->m] - L->m;
printf("pcgfl= %.0f\n",stat*(cgops+10.0*L->m+2) + cgops +4*L->m - 1);
if(stat*(cgops+10.0*L->m+2) + cgops +4*L->m - 1 < 5*ops) stat= FALSE;
*/

/*   Recover other directions */
     recover_dz(A,b,c,x,y,z,w,v,s,rz,rw,rv,rs,dx,dy,dz,dw,dv,ds,dwork,rowt,vart);

/*    Compute primal step lenght */
      passp= step_lenght(A,x,dx,w,dw,v,dv,rowt,vart);

/*    Compute dual step lenght */
      passd= step_lenght(A,z,dz,y,dy,s,ds,rowt,vart);

   if(PC)
   {
      gapp= compute_gp(A,x,y,z,w,v,s,dx,dy,dz,dw,dv,ds,rowt,vart,passp,passd);
      if(gap < 1.0) mu= gap*gap*phi;
      else mu= (gapp/gap)*(gapp/gap)*gapp/np;

      memcpy(ry,rxy,mn*sizeof(double));

/*    Form linear system |-D A'| dx = rhs
                         | A E | dy        */
      form_cls(A,u,x,y,z,w,v,s,rx,ry,rz,rw,rv,rs,mu,rowt,vart);

/*    Solve linear system */
      if(IMP)
         i= solve_ls(A,D,E,H,NULL,NULL,L,pivot,diag,rhs,pos,tol,dwork,rowt);
      else
      {
         memcpy(dwork,rhs,mn*sizeof(double));
         i= minres(A,mn,D,rhs,dwork,L,mn,tol,&dwork[mn]);
      }
      if(i < 0) error(3);
/*
      if(i*10 < nm) stat= FALSE;
*/
cgops= 4*A->col[A->n] + A->n + 4*L->col[L->m] - L->m;
if(i*(cgops+10.0*L->m+2) + cgops +4*L->m - 1 < 5*ops) stat= FALSE;

/*   Recover other directions */
     recover_dz(A,b,c,x,y,z,w,v,s,rz,rw,rv,rs,dx,dy,dz,dw,dv,ds,dwork,rowt,vart);

/*    Compute primal step lenght */
      passp= step_lenght(A,x,dx,w,dw,v,dv,rowt,vart);

/*    Compute dual step lenght */
      passd= step_lenght(A,z,dz,y,dy,s,ds,rowt,vart);
   }

/*    Compute next iterate */
      iterate(A,x,y,z,w,v,s,dx,dy,dz,dw,dv,ds,rowt,vart,passp,passd);
      fflush(stdout);
      if(SAVE)
      {
         StopClock(&xtime);
         savevar(A->m,A->n,k+1,stat,maxnz,OPS,x,y,z,w,v,s,pos,
                 xtime.UserTime*CLK_TCK,xtime.SystemTime*CLK_TCK);
         StartClock(&xtime);
      }
   }
/* Unscale(x,C); */
   ine= A->n - slack;
   if(SCALE) ddott(&ine,x,l,x);
/* ipivot */
   if(!INC && !LU) invpiv(A->m,y,w,pivot,dwork);

/* Recover solution for original problem */

   recover_sol();

   printf("%d iterations, cx= %e\n",k,cx+cnt);
   printf("rx= %e, ry= %e, rv= %e\n",dnrm2(M,rx,one)/(1.0+dnrm2(N,x,one)+dnrm2(M,w,one)),dnrm2(N,ry,one)/(1.0+dnrm2(M,y,one)+dnrm2(N,z,one)+dnrm2(M,w,one)),dnrm2(N,rv,one)/(1.0+dnrm2(N,v,one)));
   k++;
   OPS /= k;
   if(!LU) OPS += 2.0*AA->row[AA->m] + AA->row[AA->m]/k;
   printf("Flops for linear system= %f\n",OPS);
   error(-1);
}

/*-------------------------------------------------------------------*/
void error(int num)
/* Display error message and exit */
{
    switch(num)
    {
       case -1: num= 0;
                break;

       case  0: puts("Normal termination"); /* used for debugging */
                break;

       case  1: puts("Wrong Branch");
                break;

       case  2: puts("Input error");
                break;

       case  3: puts("Linear solver exceeds max iterations");
                break;

       case  4: puts("Out of memory");
                break;

       case  5: puts("Free variable");
                break;

       case  6: puts("Infinite lower bound");
                break;

       case  7: puts("Primal infeasible");
                break;

       case  8: puts("IP exceeds max iterations");
                break;

       case  9: puts("m > n");
                break;

       case 10: puts("Eig exceeds max iterations");
                break;

       case 11: puts("Primal unbounded");
                break;
    }
    StopClock(&xtime);
    printf("user= %fSec. system= %fSec.\n",xtime.UserTime,xtime.SystemTime);
/*
    remove("pdual.fix");
    remove("pdual.var");
*/
    exit(num);
}

/*-------------------------------------------------------------------*/
/* Display Math error message and exit */
/*
int matherr(struct exception *exc)
{
    switch(exc->type)
    {
       case DOMAIN:    puts("NaN");
                       break;

       case SING:      puts("Inf");
                       break;

       case OVERFLOW:  puts("Overflow");
                       break;

       case UNDERFLOW: puts("Underflow");
                       break;
     }
     error(1);
     return(0);
}
*/
int sign(int n)
{
   if(n == 0) return(0);
   if(n > 0)  return(1);
   return(-1);
}

int signd(double v)
{
   if(fabs(v) <= 1e-13) return(0);
   if(v > 0.0)  return(1);
   return(-1);
}

void options()
{

   PC= TRUE;
   LU= FALSE;
   PCG= TRUE;
   AUG= FALSE;
   RITZ= FALSE;
   IMP= TRUE;
   INC= FALSE;
   DIAG= TRUE;
   MR= TRUE;
   SCALE= TRUE;
   SLACK= TRUE;
   ORDER= TRUE;
   MARK= TRUE;
   REFAC= TRUE;
   SAVE= FALSE;
   if(!AUG || PCG) IMP= TRUE;
   if(AUG || !PCG) SAVE= FALSE;

   if(PC) printf("Predictor-Corrector, ");
   if(PCG) printf("PCG, "); else printf("MINRES, ");
   if(AUG) printf("Augmented system, "); else printf("Schur complement, ");
   if(INC) printf("Incomplete ");
   if(LU)  printf("LU"); else printf("Cholesky");
   if(!ORDER && !INC) printf(" cheap reordering");
   if(LU && MARK) printf(", Reorder columns");
   puts("");
}
