#include "main.h"

double house(int n, double *v)
{
   double betha, mu;

   mu= (v[0] > 0.0) ? dnrm2(&n,v,one) : -dnrm2(&n,v,one);
   betha= 1.0/(v[0] + mu);
   n--;
   dscal(&n,&betha,&v[1],one);
   v[0]= 1.0;
   return(-mu);
}

void tridiag(int m, double *Q, double *alpha, double *betha, double *p)
{
   int n, i, j, k;
   double dd, s, *v, *q;

   for(i= 0; i < m - 2; i++)
   {
       n= m - i - 1;
       v= &Q[i+1+i*m];
       alpha[i]= v[-1];
       betha[i]= house(n,v);
       dd= ddot(&n,v,one,v,one);
       memset(p,0,n*sizeof(double));
       for(j= 0, q= v+m; j < n; j++, q += m)
       {
           daxpy(&j,&v[j],&v[j+m],&m,p,one);
           k= n - j;
           daxpy(&k,&v[j],&q[j],one,&p[j],one);
       }
       s= 2/dd;
       dscal(&n,&s,p,one);
       s= -ddot(&n,p,one,v,one)/dd;
       daxpy(&n,&s,v,one,p,one);
       for(j= 0, q= v+m; j < n; j++, q += m)
           for(k= j; k < n; k++)
               q[k] -= (v[k]*p[j] + v[j]*p[k]);
   }
   v= &Q[(m+1)*(m-2)];
   alpha[m-2]= v[0];
   betha[m-2]= v[1];
   alpha[m-1]= Q[m*m-1];
}

void eig(int m, double *alpha, double *betha, double tol)
{
   int i, i1, j, k, l, l1;
   double d, mu, z, c, s, beta2, nt;

   for(i= m-1; i > 0; i--)
   {
       i1= i - 1;

       l1= 0; l= 1;
       for(j= 0; fabs(betha[i1]) > tol*(fabs(alpha[i]) + fabs(alpha[i1])); j++)
       {
          for(l1--; l < i; l++)
              if(fabs(betha[l-1]) < tol*(fabs(alpha[l]) + fabs(alpha[l-1])))
                 betha[l1= l - 1]= 0.0;
          l1++;
          l= l1 + 1;
          if(l > i) break;
          if(j == 100) error(10);

          if(l == i)
          {
             d= (alpha[i1] - alpha[i])*0.5;
             beta2= betha[i1]*betha[i1];
             z= sqrt(d*d + beta2);
             alpha[i1]= alpha[i] + d - z;
             alpha[i] += d + z;
             betha[i1]= 0.0;
             i--;
             break;
          }
/*
          if(d > 0.0)
             mu= alpha[i] - beta2/(d + sqrt(d*d + beta2));
          else
             mu= alpha[i] - beta2/(d - sqrt(d*d + beta2));
          if(d > 0.0)
             mu= alpha[i] + d - sqrt(d*d + beta2);
          else
             mu= alpha[i] + d + sqrt(d*d + beta2);
*/
          mu= alpha[i];

          z= betha[l1];
          givens(alpha[l1] - mu,z,&c,&s);
          rot(alpha[l1],betha[l1],c,s);
          rot(z,alpha[l],c,s);

          rot(betha[l1],alpha[l],c,s);
          alpha[l1]= s*z - c*alpha[l1];

          for(k= l; k < i; k++)
          {
             z= s*betha[k];
             betha[k] *= c;
             d= betha[k];

             givens(betha[k-1],z,&c,&s);
             rot(alpha[k],d,c,s);
             rot(betha[k],alpha[k+1],c,s);

             rot(alpha[k],betha[k],c,s);
             alpha[k+1]= s*d + c*alpha[k+1];
             betha[k-1]= s*z - c*betha[k-1];
           }

       }
   }
}

double qr(int m, double *betha, double *rho1, double *rho2, double *rho3,
               double *c, double *s)
{
   int m1= m-1;
   double aux;

   rot(rho2[m1],rho1[m],*c,*s);
   aux= rho1[m];
   rho3[m1]= (*s)*rho2[m];
   rho2[m] *= *c;
   givens(rho1[m],betha[m],c,s);
   rho1[m]= -(*c)*rho1[m] + (*s)*betha[m];
   return(aux);
}

void givens(double a, double b, double *c, double *s)
{
   double tau;

   if(b == 0.0)
   {
      *c= 1.0;
      *s= 0.0;
   }
   else
      if(fabs(b) > fabs(a))
      {
         tau= -a/b;
         *s= 1.0/(sqrt(1.0 + tau*tau));
         *c= (*s)*tau;
      }
      else
      {
         tau= -b/a;
         *c= 1.0/(sqrt(1.0 + tau*tau));
         *s= (*c)*tau;
      }
}
/*
void rot(double *r1, double *r2, double c, double s)
{
   double tau1;

   tau1= *r1;
   *r1= s*(*r2) - c*tau1;
   *r2= c*(*r2) + s*tau1;
}
*/
