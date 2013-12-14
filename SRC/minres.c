#include "main.h"
//int pcg(int *flag, matrix *ADAt, matrix *LL, matrix *A, int n, double *D, double *H, double *x, double *b, decomp *S, int *pivot, int maxit, double tol, double *work, int piiter)
int minres(matrix *A, int n, double *D, double *x, double *b, matrix *L,
           int maxit, double tol, double *work)
{
   int i, j, k, o, p, p20, kend, kend2, info, MR;
   double nb, ap, theta, tau, c, s, tolb, nb1, tmp;
   double *u, *w, *v, *y, *alpha, *betha, *rho1, *rho2, *rho3, *q, *r;
   double *bt, *filter, *Q, *H;

   o= (n == A->m) ? 0 : A->n;
   if(!IMP) o= 0;

   maxit= MIN(200,A->m);
   k= p= maxit/2;
   maxit= k + p;
   u= work;
   w= work + MAX(n,A->n);
   bt= w + MAX(n,A->n);
   y= bt + n;
   q= y + maxit;
   alpha= q + maxit + 1;
   rho1= alpha + maxit;
   betha= rho1 + maxit;

   rho2= betha + maxit;
   rho3= rho2 + maxit;
   v= rho3 + maxit - 1;
   r= v + (maxit+1)*n;

   H= r + n;
   filter= H + A->n;
   Q= filter + maxit;

   if(!IMP)
   {
      memcpy(H,D,A->n*sizeof(double));
      dsort(H,A->n);
      c= H[0]*tol;
      s= H[A->n - 1]/tol;
      memcpy(H,D,A->n*sizeof(double));
         for(i= 0; i < A->n; i++)
                if(H[i] < s) H[i]= s; else ;
/*
             if(H[i] > c)    H[i]= c;
             else
      ap= s/c;
      if(ap > 1.0)
      else puts("erro");
*/
   }

   memcpy(bt,b,n*sizeof(double));
   left(A,L,&bt[o],H,D,w);
   memcpy(x,bt,n*sizeof(double));

   kend= maxit;
   k= 0; p= kend;
   kend2= kend*kend;

   nb1= nb= MAX(1.0,dnrm2(&n,bt,one));
   tolb= tol*nb;
   nb= MAX(1.0,dnrm2(&n,b,one));

   av(A,n,o,D,L,v,x,u,w,H);
   dxmy(&n,bt,v,v);
   theta = dnrm2(&n,v,one);

   if(theta < tolb)
   {
      memcpy(u,x,n*sizeof(double));
      right(A,L,&u[o],H,D,w);
      ADAtx(A,o,D,u,r,w);
      dxmy(&n,b,r,r);
      if(dnrm2(&n,r,one) < tol*nb)
      {
         memcpy(x,u,n*sizeof(double));
         return(0);
      }
   }
/*
      Mz(L,v,u);
*/

   ap = 1.0/theta;
   dscal(&n,&ap,v,one);

/* compute w = Av */

   av(A,n,o,D,L,v+n,v,u,w,H);
/*
      ADAtx(A,D,v,u,(v+n));
      Mz(L,(v+n),u);
*/
   ap= -ddot(&n,v,one,(v+n),one);

/* compute the residual in the second col of V */

   daxpy(&n,&ap,v,one,(v+n),one);
   alpha[0]= -ap;

/* one step of iterative refinement
     to correct the orthogonality */

   ap= -ddot(&n,v,one,(v+n),one);
   alpha[0] -= ap;
   daxpy(&n,&ap,v,one,(v+n),one);

   rho1[0]= alpha[0];
   betha[0]= rho2[0]= dnrm2(&n,(v+n),one);

   info= FALSE;

   for(j= 0; ;)
   {
      givens(rho1[0],betha[0],&c,&s);
      rho1[0]= s*betha[0] - c*rho1[0];
      q[0]= -c;
      q[1]= s;
      if(q[1] == 0.0) info= TRUE;

      for(i= 1; i < kend; i++)
      {
         if(i >= k && !info)
         {
            if(j + i - k > 8*n) return(-1);
            lanczos(i,n,o,A,D,L,H,v,alpha,betha,tolb,u,w,&info);
            betha[i]= rho2[i]= dnrm2(&n,(v+n*(i+1)),one);
            rho1[i]= alpha[i];
         }
         if(info > 0) i--;
         else
         {
            tmp= qr(i,betha,rho1,rho2,rho3,&c,&s);
            q[i+1]= s*q[i];
            if(!MR)
            {
               ap= q[i];
               tau= fabs(betha[i]*ap/tmp*theta);
            }
            else tau= fabs(theta*q[i+1]);
            q[i] *= -c;
            if(tau == 0.0) info= TRUE;
            if(fabs(betha[i]) < 1e-16*fabs(alpha[i])) info= TRUE;
         }
         if(tau < tolb || info)
         {
            if(MR) triaR(i,y,rho1,rho2,rho3,q);
            else
            {
               q[i]= ap;
               SWAP(rho1[i],tmp,ap);
               triaR(i,y,rho1,rho2,rho3,q);
               SWAP(rho1[i],tmp,ap);
               q[i] *= -c;
            }
            i++;
            dscal(&i,&theta,y,one);

            Vy(n,i,v,y,u);
            dxpy(&n,x,u,u);
            right(A,L,&u[o],H,D,w);
            ADAtx(A,o,D,u,r,w);

            dxmy(&n,b,r,w);
            ap= dnrm2(&n,w,one);

            if(AUG)
            printf("rx= %e, ry= %e\n",dnrm2(N,w,one)/nb,dnrm2(M,&w[A->n],one)/nb);
/*
            left(A,L,&r[o],H,D,w);
            dxmy(&n,bt,r,r);
            printf("diff= %e\n",(dnrm2(&n,r,one) - tau)/MAX(tau,1.0));
*/

            if(tau < tolb*10e-2)
               if(ap < tol*10*nb) info= TRUE;
            if(tau < tolb*10e-4) info= TRUE;
/*
*/
            if(ap < tol*nb || info)
            {
               memcpy(x,u,n*sizeof(double));
             if(AUG)
                if(dnrm2(N,w,one) > dnrm2(M,&w[A->n],one))
                {
                   Aty(A,&x[A->n],w);
                   dxmy(N,w,b,x);
                   ddott(N,x,D,x);
                }
               printf("minres= %d\n",i + j - k);
               return(i);
            }
            else i--;
         }
      }

      printf("taup= %e\n",tau/nb1);

      memset(Q,0,kend2*sizeof(double));
      if(RITZ)
      {
         memcpy(filter,alpha,kend*sizeof(double));
         memcpy(rho2,betha,kend*sizeof(double));
         eig(kend,filter,rho2,tol);
      }
      else harmonic(kend,alpha,betha,filter,rho1,rho2,rho3,Q,work,tol);
      dsort(filter,kend);

      memset(y,0,kend*sizeof(double));
      memset(Q,0,kend2*sizeof(double));
      for(i= 0; i < kend2; i += kend + 1) Q[i]= 1.0;
      memcpy(rho2,betha,kend*sizeof(double));

/*
      for(i= 0, ap= 1.0; i < kend; i++)
          ap *= filter[i];
      ap= pow(fabs(ap),1.0/((double) kend));
      for(i= 0, ap= 0.0; i < kend; i++)
          ap += fabs(filter[i]);
      ap /= ((double) kend);
      for(p= 0; p < kend-1 && fabs(filter[p]) > ap; p++)
          shift(kend,alpha,betha,rho2,rho3,filter[p],y,Q,&theta,u,w);
      p20= p + 10;
      p20= MIN(p20,kend-1);
*/

      j += p;
      p20= kend/2;
      for(p= 0; p < p20; p++)
      {
          i= shift(kend,alpha,betha,rho2,rho3,filter[p],y,Q,&theta,u,w);
          if(i < kend - p) { p= kend-i-1; break;}
      }

      if(p != p20) printf("p= %d\n",p);
      k= kend - p;

      Vy(n,kend,v,y,u);
      dxpy(&n,x,u,x);
      VQ(n,kend,k,v,Q,w);
      daxpby(&n,&betha[k-1],(v+n*k),one,(Q+kend*k-1),(v+n*kend),one);

      if(k < kend/2) printf("zero= %e= %e\n",betha[k-1],Q[kend*k-1]);
      memcpy(rho1,alpha,k*sizeof(double));
      memcpy(rho2,betha,k*sizeof(double));
   }
}

void lanczos(int i, int n, int o, matrix *A, double *D, matrix *L, double *H,
             double *v, double *alpha, double *betha, double tol, double *u,
             double *w, int *info)
{
   int j, jm1, jp1;
   double beta, s1, s2, normv, normw, k= 0.5;

   j= n*i;
   jm1 = j - n;
   jp1 = j + n;

   beta= betha[i-1];
   if(fabs(beta) < 1e-8)
   {
      puts("beta");
      error(1);
      *info= TRUE;
      v[i+j]= 1.0;
      orthog(n,(v+jm1-n),(v+jm1),(v+j),&s2,&s1);
      orthog(n,(v+jm1-n),(v+jm1),(v+j),&s2,&s1);
      beta = dnrm2(&n,(v+j),one);
      beta = 1.0/beta;
   }
   else beta = 1.0/beta;
   dscal(&n,&beta,(v+j),one);

/* compute w = Av and store w in j+1 -st col of V */

   av(A,n,o,D,L,v+jp1,v+j,u,w,H);
/*
   ADAtx(A,D,(v+j),u,(v+jp1));
   Mz(L,(v+jp1),u);
*/
   normw= dnrm2(&n,(v+jp1),one);

/* compute the next (j-th)  column of H and
        compute the residual in the j+1 -st col of V */

   orthog(n,(v+jm1),(v+j),(v+jp1),&alpha[i],&betha[i]);
   normv= dnrm2(&n,(v+jp1),one);
   if(normv < k*normw)
   {

/* one step of iterative refinement
        to correct the orthogonality         */

      orthog(n,(v+jm1),(v+j),(v+jp1),&s2,&s1);
      if(dnrm2(&n,(v+jp1),one) > k*normv)
      {
         betha[i] += s1;
         alpha[i] += s2;
      }
      else
      {
         puts("zero");
         *info= -TRUE;
         memset((v+jp1),0,n*sizeof(double));
      }

   }
}

void av(matrix *A, int n, int o, double *D, matrix *L, double *v, double *v1,
        double *u, double *w, double *H)
{
   memcpy(u,v1,n*sizeof(double));
   right(A,L,&u[o],H,D,w);
   ADAtx(A,o,D,u,v,w);
   left(A,L,&v[o],H,D,w);
}

void orthog(int n, double *v1, double *v2, double *w, double *alpha,
            double *betha)
{
   *alpha= -ddot(&n,v2,one,w,one);
   *betha= -ddot(&n,v1,one,w,one);
   daxpy(&n,alpha,v2,one,w,one);
   daxpy(&n,betha,v1,one,w,one);
   *alpha= -(*alpha);
   *betha= -(*betha);
}

void Vy(int n, int m, double *v, double *y, double *u)
{
   int j;

   for(j= 0; j < n; j++) u[j]= v[j]*y[0];
   for(j= 1; j < m; j++)
       daxpy(&n,&y[j],(v+j*n),one,u,one);
}

void VQ(int n, int m, int k, double *v, double *Q, double *w)
{
   int i ,j, k1= k+1;

   for(i=0; i < n; i++)
   {
       memset(w,0,k1*sizeof(double));
       for(j= 0; j < m; j++)
           daxpy(&k1,&v[i+j*n],(Q+j),&m,w,one);
       for(j= 0; j <= k; j++)
           v[i+j*n]= w[j];
   }
}

int shift(int m, double *rho1, double *betha, double *rho2,
          double *rho3, double filter, double *q0, double *Q, double *theta,
          double *c, double *s)
{
    int i, j, i1;
    int mi, mi1;
    double rt;

    rt= *theta/filter;
    daxpy(&m,&rt,Q,one,q0,one);
    for(i= 0; i < m; i++)
        rho1[i] -= filter;

    givens(rho1[0],betha[0],c,s);
    rho1[0]= s[0]*betha[0] - c[0]*rho1[0];
    *theta = -rho1[0]*rt;
    for(i= 1; i < m - 1; i++)
    {
        i1= i - 1;
        rot(rho2[i1],rho1[i],c[i1],s[i1]);
        rho3[i1]= s[i1]*rho2[i];
        rho2[i] *= c[i1];
        if(fabs(betha[i]) < 1e-8*(fabs(rho1[i])+fabs(rho1[i+1]))) betha[i]= 0.0;
        givens(rho1[i],betha[i],&c[i],&s[i]);
        rho1[i]= s[i]*betha[i] - c[i]*rho1[i];
    }
    i1= i - 1;
    rot(rho2[i1],rho1[i],c[i1],s[i1]);

    memset(betha,0,(m-1)*sizeof(double));
    rot(betha[0],rho1[1],c[0],s[0]);
    rot(rho1[0],rho2[0],c[0],s[0]);
    for(j= 0; j < m; j++)
    {
        rot(Q[j],Q[m+j],c[0],s[0]);
    }
    for(i= 1; i < m-1; i++)
    {
        i1= i - 1;
        mi= m*i;
        mi1= mi + m;
        rot(betha[i],rho1[i+1],c[i],s[i]);
        rot(rho1[i],rho2[i],c[i],s[i]);
/*
        rot(rho2[i1],rho3[i1],c[i],s[i]);
*/
        rho2[i1]= betha[i1];
        for(j= 0; j < m; j++)
        {
            rot(Q[mi+j],Q[mi1+j],c[i],s[i]);
        }
        rho1[i1] += filter;
    }
    rho1[i1= i-1] += filter;
    rho1[i] += filter;
    rho2[i1]= betha[i1];
    for(i= 0; i < m; i++)
        if(betha[i] == 0.0) break;
    return(i);
}

void triaR(int m, double *y, double *rho1, double *rho2, double *rho3,double *q)
{
   int j;

   for(j= m; j >= 0; j--)
   {
       y[j]= q[j];
       if(j+2 <= m) y[j] -= (y[j+1]*rho2[j] + y[j+2]*rho3[j]);
       else
          if(j+1 <= m) y[j] -= y[j+1]*rho2[j];
       y[j] /= rho1[j];
   }
}

void triaL(int m, double *y, double *rho1, double *rho2, double *rho3,double *q)
{
   int j;

   for(j= 0; j <= m; j++)
   {
       y[j]= q[j];
       if(j > 1)
          y[j] -= (y[j-1]*rho2[j-1] + y[j-2]*rho3[j-2]);
       else
          if(j > 0) y[j] -= y[j-1]*rho2[j-1];
       y[j] /= rho1[j];
   }
}

void multL(int m, double *y, double *rho1, double *rho2, double *rho3,double *q)
{
   int j;

   for(j= 0; j <= m; j++)
   {
       y[j]= q[j]*rho1[j];
       if(j > 1)
          y[j] += q[j-1]*rho2[j-1] + q[j-2]*rho3[j-2];
       else
          if(j > 0) y[j] += q[j-1]*rho2[j-1];
   }
}

void harmonic(int m, double *alpha, double *betha, double *filter, double *rho1,
              double *rho2, double *rho3, double *Q, double *w, double tol)
{
   int i, mi, m2, m1= m - 1;

   rho1[m1]= sqrt(rho1[m1]*rho1[m1] + betha[m1]*betha[m1]);
   memset(w,0,m*sizeof(double));
   w[0]= alpha[0];
   w[1]= betha[0];
   triaL(m1,Q,rho1,rho2,rho3,w);
   for(i= 1; i < m; i++)
   {
       mi= i-1;
       w[mi]= betha[mi];
       w[i]= alpha[i];
       w[i+1]= betha[i];
       triaL(m-i,Q+i*m+mi,&rho1[mi],&rho2[mi],&rho3[mi],&w[mi]);
   }

   memset(w,0,m*sizeof(double));
   w[0]= Q[0];
   w[1]= Q[m];
   triaL(m1,Q,rho1,rho2,rho3,w);
   for(m2= m, i= 1; i < m; i++, m2 += m)
   {
       mi= m2 + i;
       w[i]= Q[mi] - rho2[i-1]*Q[mi-m];
       if(i > 1) w[i] -= rho3[i-2]*Q[mi-2*m];
       w[i+1]= Q[mi+m] - rho3[i-1]*Q[mi-m];
       triaL(m1-i,Q+mi,&rho1[i],&rho2[i],&rho3[i],&w[i]);
   }

   tridiag(m,Q,filter,w,w+m);
   eig(m,filter,w,tol);
   dinv(&m,filter,filter);
}
