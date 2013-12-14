#include "main.h"

double norm1(matrix *A)
{
   int i, j, nz= A->row[A->m];
   double *temp, tmp, max= 0.0;

   if(LU)
   {
      for(i= A->col[j= 0]; j < A->n; j++)
      {
          for(tmp= 0.0; i < A->col[j+1]; i++)
              tmp += fabs(A->val[i]);
          max= MAX(max,tmp);
      }
   }
   else
   {
      if((temp= (double *) malloc(A->n*sizeof(double))) == NULL) error(4);
      memset(temp,0,A->n*sizeof(double));

      for(j= 0; j < nz; j++)
          temp[A->col[j]] += fabs(A->val[j]);

      max= temp[idamax(N,temp)];
      free((char *)temp);
   }
   return(max);
}	

void Ax(matrix *A, double *x, double *v)
{
   int i, j, k;

   if(LU)
   {
      memset(v,0,A->m*sizeof(double));
      for(i= A->col[j= 0]; j < A->n; j++)
      {
          k= A->col[j+1];
          for( ; i < k; i++)
              v[A->row[i]] += A->val[i]*x[j];
      }
   }
   else
      for(i= A->row[j= 0]; j < A->m; j++)
      {
          k= A->row[j+1];
          for(v[j]= 0.0 ; i < k; i++)
              v[j] += A->val[i]*x[A->col[i]];
      }
}

void Aty(matrix *A, double *y, double *w)
{
   int i, j, k;

   if(LU)
      for(i= A->col[j= 0]; j < A->n; j++)
      {
          k= A->col[j+1];
          for(w[j]= 0.0; i < k; i++)
              w[j] += A->val[i]*y[A->row[i]];
      }
   else
   {
      memset(w,0,A->n*sizeof(double));
      for(i= A->row[j= 0]; j < A->m; j++)
      {
          k= A->row[j+1];
          for( ; i < k; i++)
              w[A->col[i]] += A->val[i]*y[j];
      }
   }
}

void both(matrix *A, double *x, double *v, double *y, double *w)
{
   int i, j, k;

   if(LU)
   {
      memset(v,0,A->m*sizeof(double));
      for(i= A->col[j= 0]; j < A->n; j++)
      {
          k= A->col[j+1];
          for(w[j]= 0.0; i < k; i++)
          {
              w[j] += A->val[i]*y[A->row[i]];
              v[A->row[i]] += A->val[i]*x[j];
          }
      }
   }
   else
   {
      memset(w,0,A->n*sizeof(double));
      for(i= A->row[j= 0]; j < A->m; j++)
      {
          k= A->row[j+1];
          for(v[j]= 0.0; i < k; i++)
          {
              w[A->col[i]] += A->val[i]*y[j];
              v[j] += A->val[i]*x[A->col[i]];
          }
      }
   }
}

void compute_schur(double *D, double *E, matrix *AA, matrix *s, int d1,
                   char *type)
{
   int i, j, nz= s->col[s->m];

   memset(s->val,0,nz*sizeof(double));
   if(d1)
   {
      for(i= 0; i < s->m; i++)
          if(type[i] != 'E')
             s->val[s->col[i]]= E[i];
      for(j= AA->row[i= 0]; i < nz; i++)
          for( ; j < AA->row[i+1]; j++)
              s->val[i] += AA->val[j]*D[AA->col[j]];
   }
   else
   {
/*
      for(i= 0; i < s->m; i++)
          if(type[i] != 'E')
             s->val[s->col[i]]= 1.0;
*/
      for(j= AA->row[i= 0]; i < nz; i++)
          for( ; j < AA->row[i+1]; j++)
              s->val[i] += AA->val[j];
   }
}

void cholesky(matrix *s, matrix *L, int *pos, double *diag, int *row)
{
   int i, j, k, l, p, ri, rr= 0;
   int nz= s->col[s->m];

   if(DIAG) memset(diag,0,L->m*sizeof(double));
   else for(i= 0; i < s->m; i++) diag[i]= s->val[s->col[i]];

   if(!INC)
   {
      memset(L->val,0,L->row[L->m]*sizeof(double));
      for(k= 0; k < nz; k++) L->val[pos[k]]= s->val[k];
   }

   for(k= 0; k < L->m; k++)
   {
       p= L->row[k];
       OPS += 2;
/*
       if(L->val[p] > 1e-8*diag[k]) L->val[p]= 1.0/sqrt(L->val[p]);
*/
       if(L->val[p] > 1e-16*diag[k]) L->val[p]= 1.0/L->val[p];
       else
          {
             rr++;
             if(DIAG) L->val[p]= 1.0/diag[k];
             else
             {
                L->val[p]= 1.0;
                for(j= p+1; j < L->row[k+1]; j++) L->val[j]= 0.0;
                for(j= 0; j < p; j++)
                    if(L->col[j] == k) L->val[j]= 0.0;
             }
          }

       memset(row,0,L->m*sizeof(int));
for(j= p+1; j < L->row[k+1]; j++) row[L->col[j]]= j;
/*
*/
       for(j= p+1; j < L->row[k+1]; j++)
       {
           row[L->col[j]]= j;
           L->val[j] *= L->val[p];
           OPS += 1;
           if(DIAG) diag[L->col[j]] += fabs(L->val[j]);
/*
       }
       for(j= p+1; j < L->row[k+1]; j++)
*/
           for(i= L->row[l= L->col[j]]; i < L->row[l+1]; i++)
               if((ri= row[L->col[i]]) > 0)
               {
                  OPS += 2;
                  L->val[i] -= L->val[ri]*L->val[j];
               }
/*
               else
                L->val[L->row[L->col[j]]] -= L->val[ri]*L->val[j];
*/
       }
   }
   if(rr) printf("%d smal pivots\n",rr);
/*
for(j= 0; j < L->m; j++)
 for(i= L->row[j+1]-1; i >= L->row[j]; i--)
  L->val[i] *= L->val[L->row[j]];
*/
}

void Mz(matrix *L, double *H, int *pivot, int *diag, int *p, double *z,
        double *b, double *u)
{
   int i;

   if(LU)
   {
      for(i= 0; i < L->n; i++) u[pivot[i]]= b[i];
      lub(L,diag,u);
      ddivt(&L->m,u,H,u);
      lutb(L,diag,u);
      gather(L->n,pivot,z,u);
   }
   else
   {
      if(*p) memcpy(z,b,(L->m+*p)*sizeof(double));
      else  memcpy(z,b, L->m   *sizeof(double));
      lower(L,&z[*p]);
      upper(L,&z[*p]);
   }
}

void left(matrix *A, matrix *L, double *z, double *H, double *D, double *w)
{
   int i;
   double eps= 1e-4;

   if(IMP) lower(L,z);
   else
   {
/*
      ddott(N,z,H,z);
      ddott(N,z,H,w);
      Ax(A,w,&w[A->n]);
      dxpy(M,&z[A->n],&w[A->n],&z[A->n]);
      lower(L,&z[A->n]);

      EQUIV
      FIRST
*/
/*
      ddott(N,z,H,w);
      lower(L,&z[A->n]);
      upper(L,&z[A->n]);
      Aty(A,&z[A->n],z);
      ddivt(N,z,H,z);
      dxpy(N,z,w,z);
      for(i= 0; i < A->m; i++) z[A->n+i]= w[ind[i+(A->n-A->m)/2]];

      Second
      Ax(A,w,&z[A->n]);
      lower(L,&z[A->n]);

*/
      for(i= 0; i < A->n; i++) w[i]= z[i]*sqrt(D[i]);
      Ax(A,w,&w[A->n]);
/*
      lower(L,&z[A->n]);
      upper(L,&z[A->n]);
*/

      Aty(A,&z[A->n],w);
      memcpy(&z[A->n],&w[A->n],A->m*sizeof(double));
      lower(L,&z[A->n]);
      for(i= 0; i < A->n; i++) z[i]= z[i]*sqrt(H[i]);
      dxpy(N,z,w,z);
   }
}

void right(matrix *A, matrix *L, double *z, double *H, double *D, double *w)
{
   int i;
   double eps= 1e-4;

   if(IMP) upper(L,z);
   else
   {
/*
      upper(L,&z[A->n]);
      Aty(A,&z[A->n],w);
      ddott(N,w,H,w);
      dxpy(N,z,w,z);
      ddott(N,z,H,z);
      EQUIV

      FIRST
*/

/*
      memcpy(&w[A->n],&z[A->n],A->m*sizeof(double));

      ddivt(N,z,H,w);
      Ax(A,w,&z[A->n]);
      lower(L,&z[A->n]);
      upper(L,&z[A->n]);

      upper(L,&w[A->n]);
      Aty(A,&w[A->n],w);
      dxpy(N,z,w,z);
      Second ^

      for(i= 0; i < A->m; i++) z[ind[i+(A->n-A->m)/2]] += z[A->n+i];
      ddott(N,z,H,z);
*/

      upper(L,&z[A->n]);
      Aty(A,&z[A->n],w);
      Ax(A,z,&z[A->n]);

/*
      lower(L,&z[A->n]);
      upper(L,&z[A->n]);
*/

      for(i= 0; i < A->n; i++) z[i]= z[i]*sqrt(H[i]);
      for(i= 0; i < A->n; i++) w[i]= w[i]*sqrt(D[i]);
      dxpy(N,z,w,z);
   }
}

void lower(matrix *L, double *z)
{
     int i, j, ind;

     for(i= 0; i < L->m; i++)
     {
        ind= L->row[i+1];
        for(j= L->row[i] + 1; j < ind; j++)
            z[L->col[j]] -= z[i]*L->val[j];
        z[i] *= L->val[L->row[i]];
     }
}
/*
void lower(matrix *L, double *z)
{
     int i, j, ind;

     for(i= 0; i < L->m; i++)
     {
        z[i] *= L->val[L->row[i]];
        ind= L->row[i+1];
        for(j= L->row[i] + 1; j < ind; j++)
            z[L->col[j]] -= z[i]*L->val[j];
     }
}
*/
void upper(matrix *L, double *z)
{
     int i, j, ind;

     for(i= L->m - 1; i >= 0; i--)
     {
        ind= L->row[i+1];
        for(j= L->row[i] + 1; j < ind; j++)
            z[i] -= z[L->col[j]]*L->val[j];
/*
        z[i] *= L->val[L->row[i]];
*/
     }
}
/*
void ADAtx(matrix *A, double *D, double *x, double *y, double *w)
{
   int i,j;

   memset(y,0,A->m*sizeof(double));
   for(i= A->col[j= 0]; j < A->m; j++)
       for( ; i < A->col[j+1]; i++)
       {
           y[A->row[i]] += A->val[i]*x[j];
           if(A->row[i] > j) y[j] += A->val[i]*x[A->row[i]];
       }
}
*/
void ADAtx(matrix *A, int p, double *D, double *x, double *y, double *w)
{
   int i;
   double *E;

   if(IMP)
   {
      E= D + A->n;
      Aty(A,&x[p],w);
      ddott(N,w,D,w);
      Ax(A,w,&y[p]);
      ddott(M,&x[p],E,w);
      dxpy(M,&y[p],w,&y[p]);
      if(p) memcpy(y,x,A->n*sizeof(double));
   }
   else
      Mx(A,D,x,y);
}

void Mx(matrix *A, double *D, double *x, double *y)
{
   int i;

   both(A,x,&y[A->n],&x[A->n],y);
   for(i= 0; i < A->n; i++)
       y[i] -= x[i]/D[i];
   for( ; i < A->m+A->n; i++)
       y[i] += x[i]*D[i];
}

int solve_ls(matrix *A, double *D, double *E, double *H, matrix *AA, matrix *s,
             matrix *L, int *pivot, int * diag, double *rhs, int *pos,
             double tol, double *work, char *type)
{
   int i, stat, mn= A->m + A->n;
double *sol, norm, er1;
if((sol= (double *) malloc(mn*sizeof(double))) == NULL) error(4);

   if(s != NULL)
   {
      compute_schur(D,E,AA,s,TRUE,type);
      cholesky(s,L,pos,work,(int *) (work + s->m));
   }

norm= dnrm2(&mn,rhs,one);
norm= MAX(norm,1.0);
memcpy(sol,rhs,mn*sizeof(double));

   ddott(N,rhs,D,work);
/*
   Ax(A,work,&work[A->n]);
   dxpy(M,&rhs[A->n],&work[A->n],work);
*/
   Ax(A,work,&work[mn]);
   dxpy(M,&rhs[A->n],&work[mn],work);
   memcpy(&rhs[A->n],work,A->m*sizeof(double));

   if(AUG)
   {
      for(i= 0; i < A->n; i++)
          rhs[i] *= -sqrt(D[i]);

      memcpy(work,rhs,mn*sizeof(double));

      if(PCG)
      {
         stat= pcg(A,mn,D,H,rhs,work,L,pivot,diag,mn*4,tol,&work[mn]);
      }
      else
         stat= minres(A,mn,D,rhs,work,L,mn,tol,&work[mn]);

      for(i= 0; i < A->n; i++)
          rhs[i] *= sqrt(D[i]);
   }
   else
      if(PCG)
      {
      stat= pcg(A,A->m,D,H,&rhs[A->n],work,L,pivot,diag,A->m*2,tol,&work[A->m]);
      }
      else
        stat= minres(A,A->m,D,&rhs[A->n],work,L,A->m,tol,&work[A->m]);
   if(stat < 0) return(stat);

   Aty(A,&rhs[A->n],work);

   if(AUG)
   {
      ddott(N,work,D,work);
      dxpy(N,work,rhs,rhs);
   }
   else
   {
      dxmy(N,work,rhs,rhs);
      ddott(N,rhs,D,rhs);
Mx(A,D,rhs,work);
dxmy(&mn,sol,work,work);
er1= dnrm2(&mn,work,one);
/*
printf("error1= %e\n",er1/norm);
*/
   }
memcpy(work,sol,mn*sizeof(double));
free(sol);
   return(stat);
}
