#include <stdio.h>
#include <math.h>
#include "main.h"
#include "pre.h"
#include "memory.h"
#include "string.h"

#include "splitting.h"
#define MIN(a,b) (((a)<(b)) ? (a) : (b))

extern        char             P_name[100];

/*-------------------------------------------------------------------*/
int merge(int *p1, int n1, int *p2, int n2, int *p3)
/* Merge sorted int vectors p1 and p2 into p3 eliminating
   repeteated entries 

   n1= size of p1
   n2= size of p2

   returns size of p3 */
{
    int i, j, k, t, *pt;

    i= j= k= 0;
    do
    {
       while(p1[i] < p2[j]) p3[k++]= p1[i++];
       if(p1[i] == p2[j])
       {
          p3[k++]= p1[i++];
          j++;
       }
       else
       {
          SWAP(p1,p2,pt);
          SWAP(n1,n2,t);
          SWAP(i,j,t);
       }
    }  while((i < n1) || (j < n2));
    return(k);
}

/*-------------------------------------------------------------------*/
void merge2(int *pi1, double * pd1, int n1, int *pi2, double *pd2, int n2, 
int *pi3, double *pd3)
/* Merge sorted int vectors pi1 and pi2 into pi3
             double         pd1     pd2      pd3

   n1= size of p1
   n2= size of p2 */

{
    int i, j, k;

    i= j= k= 0;
    do
    {
       while(pi1[i] < pi2[j])
       {
          pd3[k]= pd1[i];
          pi3[k++]= pi1[i++];
       }
       if(j < n2)
       {
          pd3[k]= pd2[j];
          pi3[k++]= pi2[j++];
       }
    }  while(j < n2);
    while(i < n1)
    {
       pd3[k]= pd1[i];
       pi3[k++]= pi1[i++];
    }
}

/*-------------------------------------------------------------------*/
void emerge(int *p1, int n1, int *pt, int n2, int *key, int *d, int *f,
     matrix *A)
/* Merge int vectors p1 and pt

   n1= size of p1
   n2= size of pt */

{
    int i, j, k, l, m, n, t, *p2, *d2;

    p2= pt + n1;
    d2= d + n1 + 1;
    i= j= k= m= 0;
    n= n1;
    do
    {
       while(d[i] <= d2[j]) { key[m++]= k; pt[k++]= p1[i++]; }
       if(i == n1) return;
       t= TRUE;
       for(l= A->col[p2[j]]; l < A->col[p2[j]+1]; l++) 
           if(f[A->row[l]] >= i)
           {
              t= FALSE;
              break;
           }
       if(t)
       {
          if(key[n] == j + n1)
          {
             key[m++]= k;
             n++;
          }
          pt[k++]= p2[j++];
       }
       else { key[m++]= k; pt[k++]= p1[i++]; }
    }  while(i < n1);
}

/*------------------------------------------------------------------*/
void iiddownheap(int *v, int *v2, double *w, int n, int k)
{
   int j, temp, temp2;
   double vtemp;

   temp= v[k];
   temp2= v2[k];
   vtemp= w[k];
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(v[j] >= v[j+1]) 
            if(v[j] == v[j+1])
               if(v2[j] < v2[j+1]) j++; else ;
            else ;
         else j++;
      if(temp > v[j]) break;
      if(temp == v[j]) if(temp2 > v2[j]) break;
      v[k]= v[j];
      v2[k]= v2[j];
      w[k]= w[j];
      k= j;
   }
   v[k]= temp;
   v2[k]= temp2;
   w[k]= vtemp;
}

void sort(int *v, int *v2, double *va, int n)
/* Sort int vector v, int v2 and double va by v's values 
   Uses heapsort method
   method eventually)

   n= vectors size */

{
   int k, temp;
   double vtemp;

   for(k= n/2; k >= 1; k--) iiddownheap(v-1,v2-1,va-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(v[0],v[n],temp);
      SWAP(v2[0],v2[n],temp);
      SWAP(va[0],va[n],vtemp);
      iiddownheap(v-1,v2-1,va-1,n,1);
   }
}

/*------------------------------------------------------------------*/
void ddownheap(double *v, int n, int k)
{
   int j;
   double temp, abtemp;

   temp= v[k];
   abtemp= fabs(temp);
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(fabs(v[j]) < fabs(v[j+1])) j++;
      if(abtemp >= fabs(v[j])) break;
      v[k]= v[j];
      k= j;
   }
   v[k]= temp;
}

void dsort(double *v, int n)
/* Sort double vector v by descending absolute value
   Uses heapsort method

   n= vector size */

{
   int k;
   double temp;

   for(k= n/2; k >= 1; k--) ddownheap(v-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(v[0],v[n],temp);
      ddownheap(v-1,n,1);
   }
}

/*------------------------------------------------------------------*/
void idownheap(int *v, int n, int k)
{
   int j;
   int temp;

   temp= v[k];
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(v[j] < v[j+1]) j++;
      if(temp >= v[j]) break;
      v[k]= v[j];
      k= j;
   }
   v[k]= temp;
}

void isort(int *v, int n)
/* Sort int vector v by ascending value
   Uses heapsort method

   n= vector size */

{
   int k, temp;

   for(k= n/2; k >= 1; k--) idownheap(v-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(v[0],v[n],temp);
      idownheap(v-1,n,1);
   }
}

/*------------------------------------------------------------------*/
void diownheap(double *v, int *w, int n, int k)
{
   int j, temp;
   double vtemp;

   vtemp= v[k];
   temp= w[k];
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(v[j] > v[j+1]) j++;
      if(vtemp <= v[j]) break;
      v[k]= v[j];
      w[k]= w[j];
      k= j;
   }
   v[k]= vtemp;
   w[k]= temp;
}

void sortd(double *va, int *v, int n)
/* Sort double va by va's descending values 
   Return original position on v
   Uses heapsort method

   n= vectors size */

{
   int k, temp;
   double vtemp;

   for(k= n/2; k >= 1; k--) diownheap(va-1,v-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(va[0],va[n],vtemp);
      SWAP(v[0],v[n],temp);
      diownheap(va-1,v-1,n,1);
   }
}
/*------------------------------------------------------------------*/
void iidownheap(int *v, int *w, int n, int k)
{
   int j;
   int temp, tempw;

   temp= v[k];
   tempw= w[k];
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(v[j] < v[j+1]) j++;
      if(temp >= v[j]) break;
      v[k]= v[j];
      w[k]= w[j];
      k= j;
   }
   v[k]= temp;
   w[k]= tempw;
}

void sorti(int *v, int *w, int n)
/* Sort int v by ascending value 
   Uses heapsort method

   n= vectors size */

{
   int k, temp;

   for(k= n/2; k >= 1; k--) iidownheap(v-1,w-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(v[0],v[n],temp);
      SWAP(w[0],w[n],temp);
      iidownheap(v-1,w-1,n,1);
   }
}

/*------------------------------------------------------------------*/
void iiidownheap(int *v, int *w, int *u, int n, int k)
{
   int j;
   int temp, tempw, tempu;

   temp= v[k];
   tempw= w[k];
   tempu= u[k];
   while(k <= n/2)
   {
      j= k + k;
      if(j < n)
         if(v[j] < v[j+1]) j++;
      if(temp >= v[j]) break;
      v[k]= v[j];
      w[k]= w[j];
      u[k]= u[j];
      k= j;
   }
   v[k]= temp;
   w[k]= tempw;
   u[k]= tempu;
}

void sortii(int *v, int *w, int *u, int n)
/* Sort int v by ascending value 
   Uses heapsort method

   n= vectors size */

{
   int k, temp;

   for(k= n/2; k >= 1; k--) iiidownheap(v-1,w-1,u-1,n,k);
   while(n > 1)
   {
      n--;
      SWAP(v[0],v[n],temp);
      SWAP(w[0],w[n],temp);
      SWAP(u[0],u[n],temp);
      iiidownheap(v-1,w-1,u-1,n,1);
   }
}

int match(int m, int n, int *col, int *row, int *end, int *perm, int *work)
{
    int i, j, k, c, aux, mn, rank= 0; 
    int *pr, *below, *vc, *next;

    mn= MIN(m,n);
    pr= work;
    below= pr + n;
    next= below + n;
    vc= next + n;
    for(i= 0; i < m; i++) perm[i]= vc[i]= -1;

    for(c= 0; c < n; c++)
    {
        pr[j= c]= -1;
        below[j]= col[j];
        do
        {
           for(k= below[j]; k < end[j]; k++)
               if(perm[i= row[k]] == -1) break;
           if(k < end[j])
           {
              perm[i]= j;
              below[j]= k + 1;
              for(j= pr[j]; j != -1; j= pr[j]) perm[row[next[j]-1]]= j;
              if(++rank == mn) return(rank);
              break; /* do */
           }
           else
           {
              below[j]= k;
              next[j]= col[j];
              do
              {
                 aux= FALSE;
                 for(k= next[j]; k < end[j]; k++)
                     if(vc[i= row[k]] != c)
                     {
                        vc[i]= c;
                        next[j]= k + 1;
                        aux= j;
                        j= perm[i];
                        pr[j]= aux;
                        aux= TRUE;
                        break;
                     }
                 if(aux) break;
                 j= pr[j];
              }  while(j != -1);
           }
        }  while(j != -1);
    }
    return(rank);
}

int scc(int n, int *icn, int *ip, int *ip1, int *arp, int *ib, int *work)
{
   int i, ii, icnt, num, nnm1, isn, ist, iv, iw, dummy;
   int *lowl, *numb, *prev;

   lowl= work;
   numb= lowl + n;
   prev= numb + n;

   num= icnt= 0;
   nnm1= n + n - 1;
   for(i= 0; i < n; i++)
   {
       numb[i]= -1;
       arp[i]= ip[i];
   }

   for(isn= 0; isn < n; isn++)
   {
      if(numb[isn] >= 0) continue;
      iv= isn;
      ist= 1;
      lowl[iv]= numb[iv]= 0;
      ib[n-1]= iv;

      for(dummy= 0; dummy < nnm1; dummy++)
      {
         if(arp[iv] >= ip[iv])
         {
            ii= ip1[iv];
            for(i= arp[iv]; i < ii; i++)
            {
               if(numb[iw= icn[i]] == -1)
               {
                  arp[iv]= i + 1;
                  prev[iw]= iv;
                  iv= iw;
                  lowl[iv]= numb[iv]= ist++;
                  ib[n-ist]= iv;
                  break;
               }
               if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
            }
            if(i < ii) continue;
            arp[iv]= -1;
         }
         if(lowl[iv] < numb[iv])
         {
            iw= iv;
            iv= prev[iv];
            if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
            continue;
         }
         ii= icnt + 1;

         for(i= n - ist; i < n; i++)
         {
            iw= ib[i];
            lowl[iw]= n;
            numb[iw]= icnt++;
            if(iw == iv) break;
         }
         ib[num++]= ii;
         if((ist= n - i - 1) <= 0) break;
         iw= iv;
         iv= prev[iv];
         if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
      }
      if(icnt >= n) break;
   }

   for(i= 0; i < n; i++) arp[numb[i]]= i;
   return(num);
}
void fscanf_vetor(char name[200],int size,double *vetor )
{
  FILE *file_vetor;
  int i;
  int *vetor_aux;

  
  file_vetor = fopen(name,"r");
  if (file_vetor == NULL)
  {
     printf("Nao foi possivel abrir o arquivo %s \n",name);
     error(0);
  }
  for(i = 0; i < size; i++) 
    fscanf(file_vetor, "%lf", &(vetor[i]));
  fclose(file_vetor);
}

void fprintf_vetor(double *vetor, int size, char name_file[200] )
{
  FILE      *file_vetor;
  int       i;
  char      vetor_file[200];
  // Create a file for saving a vetor, name= problem
  strcpy(vetor_file,"./OUT/");
  strcat(vetor_file,P_name);
  strcat(vetor_file,name_file);
  strcat(vetor_file,".txt");
  file_vetor=fopen(vetor_file,"w");
  if (file_vetor == NULL)
  {
    printf("Nao foi possivel abrir o arquivo %s \n",vetor_file);
    error(0);
  }
  //Save in file file_vetor the vetor 
  fprintf(file_vetor, "\n");
       for (i = 0; i < size; i++)
         fprintf(file_vetor,"%3.15e \n ",vetor[i]);
  fclose(file_vetor);
}

void fprintf_vetor_int(int *vetor, int size, char name_file[200] )
{
  FILE      *file_vetor;
  int       i;
  char      vetor_file[200];
  // Create a file for saving a vetor, name= problem
  strcpy(vetor_file,"./OUT/");
  strcat(vetor_file,P_name);
  strcat(vetor_file,name_file);
  strcat(vetor_file,".txt");
  file_vetor=fopen(vetor_file,"w");
  if (file_vetor == NULL)
  {
    printf("Nao foi possivel abrir o arquivo %s \n",vetor_file);
    error(0);
  }
  //Save in file file_vetor the vetor 
  fprintf(file_vetor, "\n");
       for (i = 0; i < size; i++)
         fprintf(file_vetor,"%i \n ",vetor[i]);
  fclose(file_vetor);
}

void fprintf_col(int *vetor, int size, int nonzeros, char name_file[200] )
{
  FILE      *file_vetor;
  int       i, j;
  char      vetor_file[200];
  // Create a file for saving a vetor, name= problem
  strcpy(vetor_file,"./OUT/");
  strcat(vetor_file,P_name);
  strcat(vetor_file,name_file);
  strcat(vetor_file,".txt");
  file_vetor=fopen(vetor_file,"w");
  if (file_vetor == NULL)
  {
    printf("Nao foi possivel abrir o arquivo %s \n",vetor_file);
    error(0);
  }
  
  //Save in file file_vetor the vetor col
  fprintf(file_vetor, "\n");
  for (i = 0; i < size; i++)
  {
    if (i == size-1)
      for (j = vetor[i]; j < nonzeros; j++)
        fprintf(file_vetor,"%i \n ",i);
    else 
      for (j = vetor[i]; j < vetor[i+1]; j++)
        fprintf(file_vetor,"%i \n ",i);
  }
  fclose(file_vetor);
}

void fprintf_A_sparse(matrix *A)
{
  fprintf_col(A->col, A->n, A->nnulos,"_col");
  fprintf_vetor_int(A->row,A->nnulos,"_row");
  fprintf_vetor(A->val,A->nnulos,"_val");  
}
void fprintf_A(matrix *A)
    
{
  int i, j, h, col, nul;
  char            A_file[200]; // Saving A
  FILE            *Pre_A;

  strcpy(A_file,"./OUT/");
  strcat(A_file,P_name);
  strcat(A_file,"_A.txt");
  Pre_A=fopen(A_file,"w");
  if (Pre_A == NULL)
   {
      printf("Nao foi possivel abrir o arquivo %s \n",A_file);
      error(1);
   }
  
  for (i = 0; i < A -> n; i++)
  {
    nul = 0;
    j = A->col[i];
    col = A->row[j];
    do
    {
      if (nul == col)
      {
        fprintf(Pre_A, "%7.4f ", A->val[j]);
        nul = nul + 1;
      }
      else
      {
        if (nul < col)
        {
          for (h = nul; h < col; h++) fprintf(Pre_A,"%7.4f ", 0);
          fprintf(Pre_A, "%7.4f ", A->val[j]);
          nul = h+1;
        }
        else
        {
          for (h = nul; h < A->m; h++) fprintf(Pre_A,"%7.4f ", 0);
          nul = h+1;
        }
      }
      if (j < A->nnulos - 1)
      {  
        if (i == A->n-1)
        {
          j = j + 1;
          col = A->row[j];
        }
        else
          if (j+1 < A->col[i+1])
          {
            j = j + 1;
            col = A->row[j];
          }
      } 
    } 
    while (nul < A->m);
    fprintf(Pre_A," \n");
  }
  fclose(Pre_A);
}
/*
int scc(int n, int *icn, int *ip, int * ip1, int *arp, int *ib, int *work)
{
   int i, ii, icnt, num, nnm1, isn, ist, iv, iw, dummy;
   int *lowl, *numb, *prev;

   lowl= work;
   numb= lowl + n + 1;
   prev= numb + n + 1;

   num= icnt= 0;
   nnm1= n + n - 1;
   for(i= 1; i <= n; i++)
   {
       numb[i]= 0;
       arp[i]= ip1[i] - ip[i] - 1;
   }

   for(isn= 1; isn <= n; isn++)
   {
      if(numb[isn] > 0) continue;
      iv= isn;
      ist= 1;
      lowl[iv]= numb[iv]= 1;
      ib[n]= iv;

      for(dummy= 1; dummy <= nnm1; dummy++)
      {
         if(arp[iv] >= 0)
         {
            ii= ip1[iv] - 1;

            for(i= ii - arp[iv]; i <= ii; i++)
            {
               if(numb[iw= icn[i]+1] == 0)
               {
                  arp[iv]= ii - i - 1;
                  prev[iw]= iv;
                  iv= iw;
                  lowl[iv]= numb[iv]= ++ist;
                  ib[n+1-ist]= iv;
                  goto I70;
               }
               if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
            }
            arp[iv]= -1;
         }
         if(lowl[iv] < numb[iv])
         {
            iw= iv;
            iv= prev[iv];
            if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
            continue;
         }
         ii= icnt + 1;

         for(i= n + 1 - ist; i <= n; i++)
         {
            iw= ib[i];
            lowl[iw]= n+1;
            numb[iw]= ++icnt;
            if(iw == iv) break;
         }
         ib[++num]= ii;
         if((ist= n - i) <= 0) break;
         iw= iv;
         iv= prev[iv];
         if(lowl[iw] < lowl[iv]) lowl[iv]= lowl[iw];
I70: ;
      }
      if(icnt >= n) break;
   }

   for(i= 1; i <= n; i++) arp[numb[i]]= i-1;
   return(num);
}
*/
