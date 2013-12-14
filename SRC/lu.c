#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "splitting.h"
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SMALL 1e-20
#define u 0.1

int luret(matrix *A, matrix **L, int si, int *maxnz, int *perm_r, int *perm_c,
           int *diag, int *inv, double *dense, float *ops, int inc, int SKIP,
           int RECO)
{
   matrix *LU;
   int i, j, k, l, r, ind, e, skip;
   int *nz, *row, *w;
   double piv, thresh;
   char *temp;
   int c, unm, ii, jj, ll, rec= FALSE;
   int *pr, *below, *next, *perm;
   unsigned int max;

   if(RECO) max= *maxnz; else max= *maxnz/2 + A->m;
   if(SKIP) rec= TRUE;
   // SKIP= FALSE; /*no key
   LU= *L;
   *ops= 0.0;
   nz= inv + A->m;
   row= nz + A->m;
   w= (SKIP) ? diag : perm_c;
   if(!SKIP)
   {
     memset(row,0,A->m*sizeof(int));
     for(k= 0; k < A->m; k++)
       for(i= A->col[w[k]]; i < A->col[w[k]+1]; i++) row[A->row[i]]++;
   }
   w= row + A->m;

   if(MATCH && SKIP)
   {
      pr= w + A->n - A->m;
      below= pr + A->m;
      next= below + A->m;
      perm= next + A->m;
      for(i= 0; i < A->m; i++) perm[i]= nz[i]= -1;
   }

   for(k= 0; k < LU->m; k++) perm_r[k]= LU->m;
   LU->col[0]= 0;

   skip= e= unm= 0;
   for(k= c= r= 0; k < LU->m; )
   {
     j= r + k;
     if(j == A->n)
     {
       puts("Acabaram as colunas");
       break;
     }
     if(SKIP)
       if(j < diag[k])
       {
         i= diag[k] - j;
         skip += i;
         memcpy(&w[r],&perm_c[j],i*sizeof(int));
         r= diag[k] - k;
       }

     if(MATCH && SKIP) if(r+k != diag[k])
     {
       perm_c[k]= perm_c[k+r];
       for(jj= r+c-1; c <= k; c++)
       {
         pr[j= c]= -1;
         below[j]= A->col[perm_c[j]];
         jj++;
         do
         {
           for(l= below[j]; l < A->col[perm_c[j]+1]; l++)
             if(perm[i= A->row[l]] == -1) break;
             if(l < A->col[perm_c[j]+1])
             {
               if(c != k)
               {
                 perm[i]= j;
                 below[j]= l + 1;
                 for(j= pr[j]; j != -1; j= pr[j]) perm[A->row[next[j]-1]]= j;
               }
               else
               {
                 jj= j;
                 ii= i;
                 ll= l;
                 j= TRUE;
               }
               break;
             }
             else
             {
               below[j]= l;
               next[j]= A->col[perm_c[j]];
               do
               {
                 ind= FALSE;
                 for(l= next[j]; l < A->col[perm_c[j]+1]; l++)
                   if(nz[i= A->row[l]] != jj)
                   {
                     nz[i]= jj;
                     next[j]= l + 1;
                     ind= j;
                     j= perm[i];
                     pr[j]= ind;
                     ind= TRUE;
                     break;
                   }
                 if(ind) break;
                 j= pr[j];
               } while(j != -1);
             }
         } while(j != -1);
       }

       if(j == -1)
       {
         c--;
         unm++;
         w[r++]= perm_c[k];
         continue;
       }
     }

     perm_c[k]= perm_c[k+r];
     i= A->col[perm_c[k]];
     j= A->col[perm_c[k]+1] - i;
     l= LU->col[k];
     LU->col[k+1]= l + j;
     scatter(j,LU->m,&A->row[i],nz,&A->val[i],dense);

     //if((!SKIP) || (k + r == diag[k]))


     if(!SKIP)
       for(j= A->col[perm_c[k]+1]; i < j; i++) row[A->row[i]]--;

     ind= MIN(k,si);
     for(i= 0; i < ind; i++)
        if(nz[inv[i]])
          for(l= LU->col[i]+1; l < LU->col[i+1]; l++)
          {
            j= LU->row[l];
            if(nz[j]) dense[j] -= LU->val[l]*dense[inv[i]];
            else
            {
              nz[j]= TRUE;
              dense[j]= -LU->val[l]*dense[inv[i]];
            }
            *ops += 2;
          }

     for(ind= TRUE; i < k; i++)
       if(nz[inv[i]])
       {
         for(l= LU->col[i]+1; l < LU->col[i+1]; l++)
         {
           j= LU->row[l];
           if(nz[j]) dense[j] -= LU->val[l]*dense[inv[i]];
           else
           {
             nz[j]= TRUE;
             dense[j]= -LU->val[l]*dense[inv[i]];
           }
           if(ind)
             if(perm_r[j] > i)
               if(fabs(dense[j]) > 5e-8) ind= FALSE;
           *ops += 2;
         }
         if(ind)
         {
           for(j= 0; j < LU->m; j++)
             if(perm_r[j] > i)
               if(nz[j])
                 if(fabs(dense[j]) > 5e-8) break;
           if(j == LU->m) break;
         }
         else ind= TRUE;
       }

       ind= -1;
       if(i == k)
       {
         piv= -1.0;
         for(i= l= 0; i < LU->m; i++)
           if(nz[i])
             if(perm_r[i] == LU->m)
             {
               l++;
               if((thresh= fabs(dense[i])) > piv)
               {
                 piv= thresh;
                 ind= i;
               }
             }
       }

       if(ind == -1 || piv < 1e-6)
       {
         //if(ind != -1) printf("Pivo pequeno %e\n",piv);
         if(MATCH && SKIP) if(r+k != diag[k]) c--;
         w[r++]= perm_c[k];
         continue;
       }

       if(MATCH && SKIP) if(r+k != diag[k])
       {
         perm[ii]= jj;
         below[jj]= ll + 1;
         for(jj= pr[jj]; jj != -1; jj= pr[jj]) perm[A->row[next[jj]-1]]= jj;
       }

       if(!SKIP)
       {
         thresh= piv*u;
         for(i= 0; i < LU->m; i++)
             if(nz[i])
               if(perm_r[i] == LU->m)
                 if(row[i] < row[ind])
                if(fabs(dense[i]) >= thresh) ind= i;
       }

       //if(piv < 1e-5) printf("piv= %e\n",piv);
       if(inc) thresh= piv*SMALL;
       perm_r[ind]= k;
       inv[k]= ind;

       while(LU->col[k] + l >= max)
       {
         if(rec && k <= si)
         {
           max= max - A->m;
           rec= FALSE;
          // puts("dobra");
           //printf("Factor flops= %.0f, %d columns",*ops,k+r-skip-unm);
          // if(r) printf(", %d skipped, %d unmatched, %d elim\n",skip,unm,r); else puts("");
           continue;
         }

         if(RECO)
         {
           // puts("Recomeca");
            j= k+r;
            OPS += *ops;
            //printf("Factor flops= %.0f, %d columns",*ops,j-skip-unm);
            //if(r) printf(", %d skipped, %d unmatched, %d elim\n",skip,unm,r);
            //else puts("");
            fflush(stdout);
            memcpy(&perm_c[k+1],&perm_c[j+1],(A->n-j-1)*sizeof(int));
            memcpy(&perm_c[A->n-r],w,r*sizeof(int));
            //printf("Valor de k + f = %d \n", k+r);
            return(-k);
         }
         puts("Expand");
         //Teste feito no dia 04/12/2006
         RECO= TRUE;

         //error(1);

         fflush(stdout);
         i= *maxnz;
         *maxnz *= 2;
         max= *maxnz;
         if((temp= Malloc((*maxnz)*sizeof(double))) == NULL) error(4);
         memcpy(temp,LU->val,i*sizeof(double));

         LU->val= (*L)->val= (double *) temp;
         if((temp= Malloc((*maxnz)*sizeof(int))) == NULL) error(4);
         memcpy(temp,LU->row,i*sizeof(int));

         LU->row= (*L)->row= (int *) temp;
       }

       *ops += l;
       j= LU->col[k] + 1;
       LU->row[LU->col[k]]= k;
       //if(fabs(dense[ind]) < 1e-16) printf("pivret= %e\n",dense[ind]);
       LU->val[LU->col[k]]= piv= 1.0/dense[ind];

       for(i= 0; i < LU->m; i++)
         if(nz[i]) /* if(dense[i] != 0.0) */
           if(perm_r[i] == LU->m)
             if(fabs(thresh= dense[i]*piv) > 1e-16)
             {
               LU->row[j]= i;
               LU->val[j++]= thresh;
             }
             else
             {
               if(inc) if( i != ind)
               {
                 if(fabs(dense[i]*dense[i]*LU->val[LU->col[perm_r[i]]]) < thresh)
                 { e++; continue; }
               }
             }
       LU->col[k+1]= j;

       if(inc)
       {
         thresh /= sqrt(SMALL);
         for(l= i= LU->col[k]+1; i < LU->col[k+1]; i++)
         if(fabs(LU->val[i]) < thresh) e++;
         else
         {
           LU->row[l]= LU->row[i];
           LU->val[l++]= LU->val[i];
         }
         LU->col[k+1]= l;
       }

       k++;
     }

     if(k < LU->m)
       if(MIX)
         for(i= 0; i < k; i++)
         {
           for(j= l= LU->col[i]+1; l < LU->col[i+1]; l++)
             if(perm_r[LU->row[l]] < LU->m)
             {
                LU->val[j]= LU->val[l];
                LU->row[j++]= perm_r[LU->row[l]];
             }
             inv[i]= j;
         }
       else ;
     else
       if(!REFAC)
         for(i= 0; i < k; i++)
           for(l= LU->col[i]+1; l < LU->col[i+1]; l++)
             LU->row[l]= perm_r[LU->row[l]];

     //if(inc) printf("%d entries elim.\n",e);

     OPS += *ops;
     //printf("Factor flops= %.0f, %d columns",*ops,k+r-skip-unm);
     //if(SKIP && r) printf(", %d skipped, %d unmatched\n",skip,unm); else puts("");
     memcpy(&perm_c[k],w,r*sizeof(int));
     //printf("Valor de k + f = %d \n", k+r);
     return(k);
}



/*
int desfaz(matrix *A, matrix *LU, int k, int ii, int *perm_r, int *perm_c, int *diag, int *inv, int *nz, double *dense)
{
int i, j, l, li;

li= FALSE;
      i= A->col[perm_c[k]];
      j= A->col[perm_c[k]+1] - i;
      scatter(j,LU->m,&A->row[i],nz,&A->val[i],dense);
      for(i= 0; i < k; i++)
          if(nz[inv[i]])
             for(l= diag[i]+1; l < LU->col[i+1]; l++)
             {
                j= LU->row[l];
if(i != ii)                dense[j] -= LU->val[l]*dense[inv[i]];
                nz[j]= TRUE;
             }
      for(i= 0; i < LU->m; i++)
          if(nz[i])
             if((perm_r[i] == LU->m) || (perm_r[i] == ii))
                if(fabs(dense[i]) >= 5e-8) { li= TRUE; break; }
return(li);
}
*/
int cld(matrix *LU, int *pos, int k, int ind, int *perm_r, int *diag, int *inv, int *nz, double *dense, int *w)
{
   int i, j, l, m, *p;
   p= pos + LU->m;

   memcpy(w,pos,k*sizeof(int));
   for(i= 0; i < k; i++) p[i]= i;
   sorti(w,p,k);
   for(i= k-1; i >= 0; i--) if(pos[ind] > w[i]) break;
   m= i;

   for(j= k-1; j > m; j--)
   {
     i= p[j];
     if(nz[inv[i]])
       for(l= diag[i]+1; l < LU->col[i+1]; l++)
         if((perm_r[LU->row[l]] == LU->m))
           if(fabs(dense[LU->row[l]] + LU->val[l]*dense[inv[i]]) >= 5e-8)  return(i);
   }
   return(ind);
}

int lufact(matrix *A, matrix **L, int *maxnz, int *perm_r, int *perm_c, int cc,
           int cr, int *diag, int *nz, double *dense, float *ops/*, int cont*/,
           int *compr)
{

   matrix *LU;
   int i, j, k, l, r, ind, first;
   int *inv, *row, *w;
   double piv, thresh;
   char *temp;
   // int *pos;

   LU= *L;
   *ops= 0.0;
   inv= nz + A->m;
   row= inv + A->m;
   w= row + A->m;

   memset(row,0,A->m*sizeof(int));
   // if(cont < 0) k= 0; else k= cont;

   for(k= cr; k < A->m; k++)
     for(i= A->col[perm_c[k]]; i < A->col[perm_c[k]+1]; i++)
       if(compr[diag[k]] == compr[A->row[i]]) row[A->row[i]]++;


   for(k= 0; k < LU->m; k++) perm_r[k]= LU->m;
   LU->col[0]= 0;

   ind= first= 0;
   for(k= 0; k < cr; k++)
   {
      ind= inv[k]= diag[k];
      perm_r[ind]= k;
      i= A->col[perm_c[k]];
      j= A->col[perm_c[k]+1];
      LU->col[k+1]= LU->col[k] + j - i;
      if(k < cc)
      {
         diag[k]= LU->col[k];
         l= diag[k] + 1;
      }
      else
      {
         diag[k]= LU->col[k+1]-1;
         l= LU->col[k];
      }
      for( ; i < j; i++)
      {
        if((r= A->row[i]) != ind)
        {
          LU->row[l]= (k < cc) ? r : perm_r[r];
          LU->val[l++]= A->val[i];
        }
        else
        {
          LU->row[diag[k]]= k;
          //if(fabs(A->val[i]) < 1e-5) printf("piv= %e\n",fabs(A->val[i]));
          LU->val[diag[k]]= 1.0/A->val[i];
        }
      }
   }

   for(r= 0; k + r < LU->m; )
   {
      perm_c[k]= perm_c[k+r];
      i= A->col[perm_c[k]];
      j= A->col[perm_c[k]+1] - i;
      l= LU->col[k];
      LU->col[k+1]= l + j;

      scatter(j,LU->m,&A->row[i],nz,&A->val[i],dense);

      for(j= A->col[perm_c[k]+1]; i < j; i++)
          if(compr[diag[k]] == compr[A->row[i]]) row[A->row[i]]--;

      if(compr[diag[k]] != compr[ind]) w[compr[ind]]= first= k;
      else
        for(i= first; i < k; i++)
          if(nz[inv[i]])
            for(l= diag[i]+1; l < LU->col[i+1]; l++)
            {
              if(compr[diag[k]] != compr[j= LU->row[l]]) continue;
              if(nz[j]) dense[j] -= LU->val[l]*dense[inv[i]];
              else
              {
                 nz[j]= TRUE;
                 dense[j]= -LU->val[l]*dense[inv[i]];
              }
              *ops += 2;
            }
        piv= -1.0;
        ind= -1;
        for(i= j= l= 0; i < LU->m; i++)
          if(nz[i])
            if(perm_r[i] == LU->m)
            {
              l++;
              if(compr[diag[k]] != compr[i]) continue;
                if((thresh= fabs(dense[i])) > piv)
                {
                   piv= thresh;
                   ind= i;
                }
             }
             else j++;
        if(ind == -1/* || piv < 5e-6*/)
        {
          w[r++]= perm_c[k];
          printf("funcao lufact\n");
          error(1);
          continue;
        }

        thresh= piv*u;
        for(i= 0; i < LU->m; i++)
          if(nz[i]) if(perm_r[i] == LU->m) if(compr[diag[k]] == compr[i])
            if(row[i] < row[ind]) if(fabs(dense[i]) >= thresh) ind= i;

        //if(piv < 1e-5) printf("piv= %e\n",piv);
        perm_r[ind]= k;
        inv[k]= ind;

        while(LU->col[k] + j + l >= *maxnz)
        {
          //fflush(stdout);
          i= *maxnz;
          *maxnz *= 2;
          if((temp= Malloc((*maxnz)*sizeof(double))) == NULL) error(4);
          memcpy(temp,LU->val,i*sizeof(double));
          Free(LU->val);
          LU->val= (*L)->val= (double *) temp;
          if((temp= Malloc((*maxnz)*sizeof(int))) == NULL) error(4);
          memcpy(temp,LU->row,i*sizeof(int));
          Free(LU->row);
          LU->row= (*L)->row= (int *) temp;
        }

        l= LU->col[k];
        diag[k]= l + j;
        j= diag[k] + 1;
        LU->row[diag[k]]= k;
        if(dense[ind] == 0.0) LU->val[diag[k]]= piv= 1.0;
        else LU->val[diag[k]]= piv= 1.0/dense[ind];

        for(i= 0; i < LU->m; i++)
          if(nz[i]) /* if(dense[i] != 0.0) */
             if(perm_r[i] == LU->m)
                if(compr[ind] != compr[i])
                {
                   LU->row[j]= i;
                   LU->val[j++]= dense[i];
                }
                else
                   if(fabs(thresh= dense[i]*piv) > 1e-18)
                   {
                      *ops += 1;
                      LU->row[j]= i;
                      LU->val[j++]= thresh;
                   }
                   else ;
             else
                if(i != ind)
                {
                   LU->row[l]= perm_r[i];
                   LU->val[l++]= dense[i];
                }
        LU->col[++k]= j;
   }

   // if(cont < 0)

   for(i= 0; i < k; i++)
     for(l= diag[i]+1; l < LU->col[i+1]; l++)
        LU->row[l]= perm_r[LU->row[l]];

   OPS += *ops;
   // if(cont < 0)

   //printf("Factor flops= %.0f, %d columns\n",*ops,k+r);

   /*if(cont < 0)

   else
     if(k != A->m)
       memcpy(perm_c,pos,k*sizeof(int));
       memcpy(&perm_c[k],w,r*sizeof(int));
   */
   i= compr[ind]*sizeof(int);
   memcpy(compr,w,i);
   return(k);
}
#undef u

void lub(decomp *S, double *z)
{
   int i, j, k, ind;
   int *diag, *comp;
   matrix *LU;

   LU= S->L;
   diag= S->diag;
   comp= S->comp;

   for(i= 0; i < S->cc; i++)
   {
       z[i] *= LU->val[diag[i]];
       ind= LU->col[i+1];
       for(j= diag[i]+1; j < ind; j++)
           z[LU->row[j]] -= z[i]*LU->val[j];
   }

   for(k= S->scc; k > 0; k--)
   {
      for(i= comp[k-1]; i < comp[k]; i++)
      {
        ind= LU->col[i+1];
        for(j= diag[i]+1; j < ind; j++)
           z[LU->row[j]] -= z[i]*LU->val[j];
      }

      for(i= comp[k] - 1; i >= comp[k-1]; i--)
      {
        z[i] *= LU->val[ind= diag[i]];
        for(j= LU->col[i]; j < ind; j++)
           z[LU->row[j]] -= z[i]*LU->val[j];
      }
   }

   for(i= comp[0] - 1; i >= S->cc; i--)
   {
      z[i] *= LU->val[ind= diag[i]];
      for(j= LU->col[i]; j < ind; j++)
         z[LU->row[j]] -= z[i]*LU->val[j];
   }
}

void lutb(decomp *S, double *z)
{
   int i, j, k, ind;
   int *diag, *comp;
   matrix *LU;

   LU= S->L;
   diag= S->diag;
   comp= S->comp;

   for(i= S->cc; i < comp[0]; i++)
   {
      ind= diag[i];
      for(j= LU->col[i]; j < ind; j++)
         z[i] -= z[LU->row[j]]*LU->val[j];
      z[i] *= LU->val[j];
   }

   for(k= 1; k <= S->scc; k++)
   {
     for(i= comp[k-1]; i < comp[k]; i++)
     {
       ind= diag[i];
       for(j= LU->col[i]; j < ind; j++)
          z[i] -= z[LU->row[j]]*LU->val[j];
       z[i] *= LU->val[j];
     }

     for(i= comp[k] - 1; i >= comp[k-1]; i--)
     {
       ind= LU->col[i+1];
       for(j= diag[i]+1; j < ind; j++)
           z[i] -= z[LU->row[j]]*LU->val[j];
     }
   }

   for(i= S->cc - 1; i >= 0; i--)
   {
      ind= LU->col[i+1];
      for(j= diag[i]+1; j < ind; j++)
         z[i] -= z[LU->row[j]]*LU->val[j];
      z[i] *= LU->val[diag[i]];
   }
}

int forma(matrix *A, int *pos, int *perm_c, int stat, int *key, int *d, int *w)
{
   int i, j, jj, li, *r;

   r= d + A->n;
   for(j= 0, i= stat; i < A->m; i++)
   {
      if(pos[perm_c[i]] == -1) continue;
      w[j++]= pos[perm_c[i]];
      pos[perm_c[i]]= -1;
   }
   jj= A->n - j - 1;
   li= stat;
   for(i= 0; i < stat; i++) if(pos[perm_c[i]] == -1) { li= i; break;}
   for(j= i= 0; i < stat; i++)
     if(pos[perm_c[i]] != -1) perm_c[j++]= pos[perm_c[i]];
   stat= j;
   if(li > stat) { printf("funçao forma li > stat"); error(1);};
   for(j= i= 0; i < jj; i++)
   {
      while(pos[i+j] == -1) j++;
      pos[i]= pos[i+j];
   }
   memcpy(&pos[i],w,j*sizeof(int));
   for(i= 0; i < A->n; i++)
   {
      d[i]= A->col[pos[i]];
      r[i]= A->col[pos[i]+1];
   }
   i= match(A->m,A->n,d,A->row,r,key,w);
   memset(d,0,A->n*sizeof(int));
   for(i= 0; i < stat; i++) d[perm_c[i]]= TRUE;
   for(j= stat, i= 0; i < A->m; i++)
      if(!d[pos[key[i]]]) perm_c[j++]= pos[key[i]];
   if(j != A->m) { printf("funcao forma: j diferente de m"); error(1); };
   isort(&perm_c[stat],A->m-stat);
   for(i= 0; i < A->n; i++) d[pos[i]]= i;
   for(i= 0; i < A->m; i++) key[i]= d[perm_c[i]];
   return(li);
}

int colorder(matrix *A, int *pos, double *H, double *w, int *key, int *r, int *ant,
              int *o)
{
   int i, j, k, l, *p, first, second, *c, *col, *d= (int *) w;
   int m, s, ii, jj, mm, *e;

   if(H != NULL)
   {
     /*i= r[0];
     j= r[1];
     k= r[2];
     l= A->n - i - j - k;*/
     gather(A->n,pos,w,H);
     sortd(w,pos,A->n);
     /*
     sortd(w,pos,i+j+k);
     sortd(&w[i],&pos[i],j+k+l);
     sortd(&w[i+j],&pos[i+j],k);
     sortd(&w[i+j+k],&pos[i+j+k],l);
     */
     memset(r,0,A->n*sizeof(int));
     for(i= 0; i < A->n; i++) r[pos[i]]= TRUE;
     for(i= 0; i < A->n; i++) if(!r[i]) error(0);
   }

   //return; /*no key columns
   // Foi eliminado para testar as colunas de B não LI?//
   for(i= 0; i < A->n; i++)
   {
      d[i]= A->col[pos[i]];
      r[i]= A->col[pos[i]+1];
   }
   i= match(A->m,A->n,d,A->row,r,key,&r[A->n]);
   //if(i != A->m) printf("Symbolically singular, rank= %d\n",i);
   /*
   if(H != NULL) return; /*no independent
   */
   c= r + A->m + 1;
   col= c + A->n;
   p= col + A->n;
   e= d + A->m;
   for(i= 0; i < A->m; i++)
   {
      d[i]= A->col[pos[key[i]]];
      r[i]= A->col[pos[key[i]]+1];
   }
   m= scc(A->m,A->row,d,r,c,col,o);
   memcpy(c,&o[A->m],A->m*sizeof(int));
   for(i= 0; i < A->m; i++) r[i]= i;
   sortii(key,r,c,A->m);
   for(i= 0; i < m; i++) o[i]= --col[i];
   for(j= i= 0; i < A->m; i++)
   {
      d[i]= j;
      if(i == col[j]) j++;
   }
   mm= 0;
   for(s= i= 0; i < A->m; i++)
     if(col[j= d[c[i]]] > 0)
     {
        col[j]= -key[i];
        /*
        s++;
        */
     }
     else
     {
       if((l= -col[j]) == A->m) break;
       ii= pos[key[i]];
       for(k= A->col[ii]; k < A->col[ii+1]; k++)
          if(A->row[k] >= r[l]) break;
       if(k < A->col[ii+1])
         if(A->row[k] > r[l])
         {
            col[j]= -A->m;
            s++;
            continue;
         }
       ii= r[key[i]];
       l= pos[l];
       for(k= A->col[l]; k < A->col[l+1]; k++)
         if(A->row[k] >= ii) break;
       if(k < A->col[l+1])
         if(A->row[k] > ii)
         {
            col[j]= -A->m;
            s++;
            continue;
         }
         break;
     }
     for(i= 0; i < A->m; i++)
       if(o[j= d[c[A->m-i-1]]] > 0) o[j]= -key[i];
       else
       {
          if((l= -o[j]) == A->m) break;
          ii= pos[key[i]];
          for(k= A->col[ii]; k < A->col[ii+1]; k++)
              if(A->row[k] >= r[l]) break;
          if(k < A->col[ii+1])
             if(A->row[k] > r[l])
             {
                o[j]= -A->m;
                mm++;
                continue;
             }
          ii= r[key[i]];
          l= pos[l];
          for(k= A->col[l]; k < A->col[l+1]; k++)
              if(A->row[k] >= ii) break;
          if(k < A->col[l+1])
             if(A->row[k] > ii)
             {
                o[j]= -A->m;
                mm++;
                continue;
             }
          break;
       }
     /*
     isort(key,A->m);
     */
     for(i= 0; i < A->m; i++) c[i]= pos[key[i]];
     block(&i,&j,A,c,col,o,FALSE);
     memset(o,0,A->n*sizeof(int));
     for(i= 0; i < j; i++) o[c[i]]= TRUE;

     memset(r,0,A->m*sizeof(int));
     memset(col,0,A->n*sizeof(int));
     key[-1]= -1;

     for(i= 0; i < A->m; i++) if(!ant[pos[i]]) break;
     //printf("numero de scc= %d, s= %d, mm= %d, repetidos= %d\n",m,s,mm,i);
     if(s < mm) s= mm;
     if(s < i) s= i;

     for(mm= k= l= m= 0; m < A->m; m++)
     {
       j= key[m] - key[m-1] - 1;
       if(j)
       {
          memcpy(&c[l],&pos[key[m-1]+1],j*sizeof(int));
          l += j;
       }

       j= pos[key[m]];
       if(o[j]) first= FALSE;
       else
       {
         first= TRUE;
         second= FALSE;
         for(i= A->col[j]; i < A->col[j+1]; i++)
         {
           if(r[A->row[i]] == 0)
           {
              first= FALSE;
              d[A->row[i]]= k;
              col[k]++;
           }
           else
              if(r[A->row[i]] == 1 && col[d[A->row[i]]] > 1) p[second++]= i;
           r[A->row[i]]++;
       }

       if(first && second)
          for(i= 0; i < second; i++)
          {
              jj= d[A->row[p[i]]];
              for(ii= A->col[pos[jj]]; ii < A->col[pos[jj]+1]; ii++)
                  if(r[A->row[ii]] == 1)
                  {
                     col[jj]--;
                     col[k]= 1;
                     d[A->row[p[i]]]= k;
                     first= FALSE;
                     break;
                  }
              if(!first) break;
          }
      }

      if(m < s)
         if(first)
         {
            first= FALSE;
            mm++;
         }

       if(first)
       {
          e[m-k]= l;
          c[l++]= j;
       }
       else pos[k++]= j;
   }

   //printf("Symbolically independent= %d, (mais %d)\n",k,mm);
   memcpy(&pos[k],c,l*sizeof(int));

   if(MARK)
   {
      for(i= 0; i < k; i++) r[i]= i;
      degree(A,pos,c,r,k);
/*
      btb(A,pos,col,p,k,key);
      if(ORDER)
      {
         for(i= 0; i < k; i++) key[i]= key[i+k]= i;
         i= order(k,col,p,key,&key[k],10*A->col[A->n],o);
      }
      else i= cheap_order(k,col,p,key,&key[k]);
      if(i) error(-1);
      memcpy(p,pos,k*sizeof(int));
      memcpy(col,r,k*sizeof(int));
      for(i= 0; i < k; i++)
      {
          pos[i]= p[key[i]];
          r[i]= col[key[i]];
      }
/* com emerge deixe
*/
/*
      c[k]= A->n;
      for(i= k; i < A->n; i++) c[i+1]= A->col[pos[i]+1] - A->col[pos[i]];
      for(i= k; i < A->m; i++) key[i]= e[i-k] + k;
      for(i= 0; i < k; i++)
      {
         e[r[i]]= i;
         r[i]= pos[i];
      }
      for(i= 0; i < A->m; i++) d[i]= e[d[i]];
      emerge(r,k,pos,A->n,key,c,d,A);
   }
   else
   {
retire para emerge */
      for(i= 0; i < k; i++) key[i]= i;
      for(    ; i < A->m; i++) key[i]= e[i-k] + k;
   }
   return(k);
}

void block(int *m, int *n, matrix *A, int *key, int *perm_r, int *row, int MAIS)
{
   int i, j, k, l, p, mais, *col, *w;

   col= row + A->m;
   w= col + A->m;
   *n= 0;
   do
   {
      memset(row,0,A->m*sizeof(int));
      //for (i = 0; i < A->m; i++) row[i] = -1;
      for(i= 0; i < A->m; i++)
      {
          k= key[i];
          if(k != -1)
             for(j= A->col[k]; j < A->col[k+1]; j++)
             {
                 row[A->row[j]]++;
                 col[A->row[j]]= i;
             }
      }
      j= *n;
      for(i= 0; i < A->m; i++)
          if(row[i] == 1) //&& (key[col[i]]!=-1))
          {
             //if (key[col[i]]==-1) puts("Cuidado");
             perm_r[*n]= i;
             w[(*n)++]= key[col[i]];
             key[col[i]]= -1;
          }
   }  while(j != *n);
   *m= *n;
   //Begin Teste 1
   /*for (i = 0; i < A->m; i++) col[i] = FALSE;
   for (i = 0; i < *n; i++)
     col[perm_r[i]] = TRUE;
   for (i=0,j=0; i < A->m; i++)
     if (col[i]) j++;
   if (j!=*n)
   {
     puts("Error antes do mais diag não tem todos");
     printf("valor do j %d e valor de *n %d\n", j, *n);
     error(1);
   }*/
   //End do teste 1
   //Begin Teste 2
//   j=0;
////   for (i=0; i < A->m; i++)
//     if (row[i]==0) //(!col[i]) &&
//       j++;
  // if (j > 0)
//   {
  //   printf("Foram encontradas %d linhas de zeros\n",j);
     //error(1);
//   }
   //End Teste 2
   mais= MAIS;
   while(mais)
   {
      mais= FALSE;
      for(i= 0; i < A->m; i++)
      {
          k= key[i];
          if(k != -1)
          {
             l= 0;
             for(j= A->col[k]; j < A->col[k+1]; j++)
                 if(row[A->row[j]] > 0)
                 {
                    l++;
                    if(l > 1) break;
                    p= A->row[j];
                 }
             if(l == 1)
             {
                perm_r[*n]= p;
                row[p]= -1;
                w[(*n)++]= k;
                mais= TRUE;
                key[i]= -1;
             }
          }
      }
   }
   //Begin Teste 3
//   for (i = 0; i < A->m; i++)
//     col[i] = FALSE;
//   for (i = 0; i < *n; i++)
//     col[perm_r[i]] = TRUE;
//   for (i=0,j=0; i < A->m; i++)
////     if (col[i]) j++;
//   if (j!=*n)
//   {
//     puts("Error depois do mais diag não tem todos");
//     error(1);
//   }
   //End do teste 3
//   printf("Valor do *n %d e valor do *m %d\n",*n,*m);
   //Begin Teste 2
 //  for (i=0; i < A->m; i++)
 //    if ((!col[i]) && (row[i]==0))
  //   {
 //      puts("Tem uma linha de zeros no restante");
 //      error(1);
//     }
//   //Begin Teste 2

   if (*n == 0) return;
   /*j = 0;
   for(i= A->m - 1; i >= *n; i--)
   {
       while(key[i-j] == -1) j++;
       key[i]= key[i-j];
   }*/
   for (i=0, j= *n; i < A->m; i++)
     if (key[i] != -1) w[j++] = key[i];
   memcpy(key,w,(A->m)*sizeof(int));
   if (j!= A->m)
   {
     printf("j e diferente de m no block");
     //error(0);
   }
   //printf(" mais %d colunas\n",*n-*m);
}
/*
void block(int *m, int *n, matrix *A, int *key, int *org, int *pos, int *perm_r,
           int *row)
{
   int i, j, k, l, p, mais, *col, *w;

   col= row + A->m;
   w= col + A->m;
   *n= 0;
   do
   {
      memset(row,0,A->m*sizeof(int));
      for(i= 0; i < A->m; i++)
      {
          k= key[i];
          if(k != -1)
             for(j= A->col[k]; j < A->col[k+1]; j++)
             {
                 row[A->row[j]]++;
                 col[A->row[j]]= i;
             }
      }
      j= *n;
      for(i= 0; i < A->m; i++)
          if(row[i] == 1)
          {
             perm_r[*n]= i;
             w[(*n)++]= org[col[i]];
             key[col[i]]= -1;
          }
   }  while(j != *n);
   printf("%d colunas",*n);
   *m= *n;

do
{
   mais= FALSE;
   for(i= 0; i < A->m; i++)
   {
       k= key[i];
       if(k != -1)
       {
          l= 0;
          for(j= A->col[k]; j < A->col[k+1]; j++)
              if(row[A->row[j]])
              {
                 l++;
                 if(l > 1) break;
                 p= A->row[j];
              }
          if(l == 1)
          {
             perm_r[*n]= p;
             row[p]= 0;
             w[(*n)++]= org[i];
             mais= TRUE;
             key[i]= -1;
          }
       }
   }
}  while(mais);

   j= 0;
   for(i= A->m - 1; i >= *n; i--)
   {
       while(key[i-j] == -1) j++;
       key[i]= key[i-j];
org[i]= org[i-j];
   }
   memcpy(org,w,(*n)*sizeof(int));
if(pos != NULL) for( ; i >= 0; i--) key[i]= pos[org[i]];
   printf(" mais %d colunas\n",*n-*m);
}
*/
int blocktri(int n, matrix *A, matrix *L, int *pos, int *pivot, int *diag,
             int *key, int *w, int maxnz)
{
   int m, i, j, k;
   int *r, *c;

   c= key + A->m;
   r= c + A->m;
   m= A->m - n;
  /* printf("Valor do m no blocktri = %d\n",m);
   printf("valor do n no blocktri = %d\n",n);
   printf("Valor da A->m no blocktri = %d\n",A->m);*/
   memset(key,0,A->m*sizeof(int));
   for(i= 0; i < n; i++) key[diag[i]]= A->n;
   for(j= i= 0; i < A->m; i++)
       if(!key[i])
       {
          key[i]= j;
          c[j++]= i;
       }
   if (j!=m)
   {
     printf("No blocktri j saiu diferente de m\n");
     error(1);
   }
   for(j= 0, i= 0; i < m; i++)
   {
       L->col[i]= j;
       for(k= A->col[pos[i]]; k < A->col[pos[i]+1]; k++)
           if(key[A->row[k]] < A->n) L->row[j++]= key[A->row[k]];
   }
   L->col[m]= j;
   for(i= 0; i < m; i++)
   //Begin teste
//   if(L->col[i]== L->col[i+1])
////   {
 //      printf("Estou no blocktri saido do novo teste");
//       error(1);
 //  }
   // End teste
   /*printf("Valor do m no blocktri antes do match = %d\n",m);
   printf("valor do n no blocktri antes do match = %d\n",n);
   printf("Valor da A->m no blocktri antes do match = %d\n",A->m);*/
   i= match(m,m,L->col,L->row,&L->col[1],key,r);
   //printf("Valor do i depois do match no Blocktri = %d\n");
   /*printf("Valor do m no blocktri = %d\n",m);
   printf("valor do n no blocktri = %d\n",n);
   printf("Valor da A->m no blocktri = %d\n",A->m);*/
   //if(i != m) {printf("funcao blocktri: i diferente de m\n"); error(1);};
   for(i= 0; i < m; i++)
   {
      w[i]= L->col[key[i]];
      r[i]= L->col[key[i]+1];
   }
   k= scc(m,L->row,w,r,&w[m],pivot,&r[m]);
   if(k == A->m) {printf("funcao blocktri: k diferente de m\n");error(1);} else pivot[k]= m+1;
   r= w + m;
   for(i= 0; i < m; i++)
   {
       diag[i+n]= c[r[i]];
       w[i]= pos[r[i]= key[r[i]]];
   }
   memcpy(pos,w,m*sizeof(int));

   for(j= n-1, i= 0; i < k; i++)
   {
       if(pivot[i+1]-pivot[i] > 3)
       {
          degree(L,r,c,&pos[pivot[i]-1],pivot[i+1]-pivot[i]);
          if(MD) MD= ordena(A,&pos[pivot[i]-1],&w[A->m],L->row,pivot[i]+j,
                 pivot[i+1]+j,diag,L->col,w,maxnz,key);
       }
       pivot[i] += j;
   }
   pivot[k]= A->m;
   return(k);
}

/*
void degree(matrix *A, int *pos, int *r)
{
   int i, j;

   memset(r,0,A->m*sizeof(int));
   for(j= 0; j < A->m; j++)
       for(i= A->col[pos[j]]; i < A->col[pos[j]+1]; i++)
           r[A->row[i]]++;
}
*/
void degree(matrix *A, int *pos, int *w, int *r, int n)
{
   int j;

   for(j= 0; j < n; j++)
       w[j]= A->col[pos[j]+1] - A->col[pos[j]];
   if(r == NULL) sorti(w,pos,n); else sortii(w,pos,r,n);
}

void gather(int m, int *perm_c, double *v, double *w)
{
     int i;

     for(i= 0; i < m; i++) v[i]= w[perm_c[i]];
}

void scatter(int m, int n, int *perm, int *d, double *v, double *w)
{
     int i;

     memset(d,0,n*sizeof(int));
/*
     memset(w,0,n*sizeof(double));
*/
     for(i= 0; i < m; i++)
     {
        d[perm[i]]= TRUE;
        w[perm[i]]= v[i];
/*
        r[perm[i]]--;
*/
     }
}
