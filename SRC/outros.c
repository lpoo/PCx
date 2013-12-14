#include "splitting.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define ccf ccf_
#define dscale dscale_
//#define ritzva ritzva_



extern FILE *eigout;
extern int *trocou,Eta;


/* Display error message and exit */

void error(num)
int num;
{
    switch(num)
    {
       case -1: num= 0;
                break;

       case  0: puts("Normal termination");
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
    exit(num);
}

int signd(double v)
{
   if(fabs(v) <= 1e-13) return(0);
   if(v > 0.0)  return(1);
   return(-1);
}

void options()
{
   TROCA= FALSE;
   LU= TRUE;
   MARK= TRUE;
   REFAC= TRUE;
   MATCH= FALSE;
   MIX= TRUE;
   TOL= TRUE;
   INE= TRUE;
   MD= TRUE;

   LLt= FALSE;

   INC= FALSE;
   ORDER= TRUE;// Deve ser atualizado a TRUE  teste feito no dia 31/01/07

   if(!LU) MIX= FALSE;
}

int ordena(matrix *A, int *pos, int *cols, int *rows, int i1, int i2,
           int *perm, int *pivot, int *iwork, int maxnz, int *work)
{
    int i, nm= i2-i1;

    if(nm == 0) return(TRUE);
    i= btb(A,pos,cols,rows,i1,i2,perm,pivot,maxnz);
    if(!i) return(FALSE);
    if(ORDER)
    {
       for(i= 0; i < nm; i++) iwork[i]= pivot[i]= i;
       i= order(nm,cols,rows,iwork,pivot,maxnz-1,work);
       if(!i) ORDER= FALSE;
    }
    if(!ORDER) i= cheap_order(nm,cols,rows,iwork,pivot);
    if(!i) return(FALSE);
    memcpy(pivot,pos,(nm)*sizeof(int));
    for(i= 0; i < nm; i++) pos[i]= pivot[iwork[i]];
    return(TRUE);
}

int btb(matrix *A, int *pos, int *cols, int *rows, int i1, int i2,
        int *perm, int *adr, int maxnz)
{
   int i, j, k, cc, kk, mm, ind= 0;

   mm= i2 - i1;
   memset(adr,0,A->m*sizeof(int));
   for(i= 0; i < i1; i++) adr[perm[i]]= -1;
   for(i= i2; i < A->m; i++) adr[perm[i]]= -1;
   for(i= 0; i < mm; i++)
   {
      cols[i]= ind;
      cc= A->col[pos[i]+1];
      for(j= A->col[pos[i]]; j < cc; j++)
          if(adr[A->row[j]] != -1) adr[A->row[j]]= TRUE;
      for(k= 0; k < mm; k++)
      {
          kk= A->col[pos[k]+1];
          for(j= A->col[pos[k]]; j < kk; j++)
              if(adr[A->row[j]] == TRUE)
              {
                 rows[ind++]= k;
                 if(ind < maxnz) break;
                 puts("btb Insufficient Storage");
                 return(FALSE);
              }
      }
      for(j= A->col[pos[i]]; j < cc; j++)
          if(adr[A->row[j]] == TRUE) adr[A->row[j]]= 0;
   }
   cols[i]= ind;
   return(TRUE);
}

void Mz(decomp *S, double *H, int *pivot, int *p, double *z, double *b,
        double *u)
{
   int i;
   matrix *L;

   L= S->L;
   if(LU)
   {
      if(INC) for(i= 0; i < L->n; i++) u[pivot[i]]= b[i];
      else    for(i= 0; i < L->m; i++) u[pivot[i]]= b[i];
      lub(S,u);
      ddivt(&L->m,u,H,u);
      lutb(S,u);
      if(INC) gather(L->n,pivot,z,u); else gather(L->m,pivot,z,u);
   }
   else error(1);
}

void Axv(matrix *A, double *x, double *v)
{
   int i, j, k;

      memset(v,0,A->m*sizeof(double));
      for(i= A->col[j= 0]; j < A->n; j++)
      {
          k= A->col[j+1];
          for( ; i < k; i++)
              v[A->row[i]] += A->val[i]*x[j];
      }
}

void Aty(matrix *A, double *y, double *w)
{
   int i, j, k;

      for(i= A->col[j= 0]; j < A->n; j++)
      {
          k= A->col[j+1];
          for(w[j]= 0.0; i < k; i++)
              w[j] += A->val[i]*y[A->row[i]];
      }
}

void ADAtx(matrix *A, int p, double *D, double *x, double *y, double *w)
{
   int i;

   Aty(A,&x[p],w);
   ddott(&A->n,w,D,w);
   Axv(A,w,&y[p]);
   if(p) memcpy(y,x,A->n*sizeof(double));
}

int reduz(matrix *A, int *row, int *col, int *rowm, int *colm)
{
   int i, j, k, m, tot;

   for(m= i= 0; i < A->m; i++)
       if(rowm[i] == STILL_ACTIVE) { m++; row[i]= -1;} else row[i]= A->m;
   for(i= 0; i < A->n; i++)
       if(colm[i] == STILL_ACTIVE) col[i]= -1; else col[i]= A->n;
   j= 0;
   do
   {
     tot= j;
     for(i= 0; i < A->n; i++)
     {
        if(col[i] >= 0) continue;
        if(row[A->row[k= A->col[i]]] < 0) col[i]= A->row[k];
        for(k++; k < A->col[i+1]; k++)
            if(row[A->row[k]] < 0)
               if(col[i] >= 0)
               {
                  col[i]= -1;
                  break;
               }
               else col[i]= A->row[k];
        if(col[i] != -1)
        {
           row[col[i]]= 0;
           j++;
        }
     }
   } while(tot != j);
   //printf("Independentes= %d\n",tot);
   return(tot == m);
}

/******************************************************************************
/*************************************************************/
unsigned long int
calculannulosAAt(A,At)
     matrix *A,*At;
{
    int i,kbeg,kend,IR,k,j,ii;
    int *marca;
    int status,ancol,aancol;
    int *Aendcol,*Atendcol;
    unsigned long int nzaa;

    /* numera linhas e colunas a partir de 1 */
    for (i=0; i<= A->n; i++) A->col[i]=A->col[i]+1;

    for (i=0; i< A->nnulos; i++) A->row[i]=A->row[i]+1;

    /*....................................................*/

    /*aloca memoria para vetores temporarios*/
    if(!(Aendcol=(int*)malloc(A->n*sizeof(int))))
      {printf("sem memoria para A->col \n");
	exit(1);
      }
    if(!(Atendcol=(int*)calloc(A->m,sizeof(int))))
      {printf("sem memoria para A->row  \n");
	exit(1);
      }

    /*..................................................*/



    for (i=0; i< A->n; i++) Aendcol[i]=A->col[i+1]-1;

    /*calula transposta da A*/

    //printf(" Calcula transposta \n");
    status= TransposeSparseRealMatrix(A->val, A->col, Aendcol, A->row,
                    At->val, At->col, Atendcol, At->row, &A->m, &A->n);
    if(status)
       printf("Error: TransposeSparseRealMatrix = %d\n", status);

    /*armazena a posicao m+1 de Atcol*/
    At->col[A->m]=Atendcol[A->m-1]+1;


    /*libera memoria de vetores temporarios */
    free(Aendcol);
    free(Atendcol);

    /* aloca memoria para vetores locais */
    if(!(marca=(int*)calloc(A->n, sizeof(int))))
    {
      printf("sem memoria para A->row  \n");
      exit(1);
    }

    nzaa=0;



    for (i=0; i< A->m; i++)
    {
       nzaa+=1;
       kbeg=At->col[i];
       kend=At->col[i+1];

       /* marca quais colunas da linha i que possuem elementos nao nulos */
       for (k=kbeg; k < kend; k++)
       {
         IR=At->row[k-1];
         marca[IR-1]=1;
       }

       /* Procura linhas adjacentes a linha i*/
       for (j=i+1; j< A->m ; j++)
       {
          kbeg=At->col[j];
          kend=At->col[j+1];
          for (k= kbeg ; k< kend; k++)
	  {
             IR=At->row[k-1];
             if (marca[IR-1]==1)
             {
                nzaa+=1;
                break;
	     }
          }
       }


      /*reinicia vetor para marcar as linhas que contem elementos em cada coluna */
      for (ii=0; ii< A->n ; ii++) marca[ii]=0;

    }

    free(marca);

    /*retorna estrutura de dados comecando de zero */
    for (i=0; i<= A->n; i++) A->col[i]=A->col[i]-1;

    for (i=0; i< A->nnulos; i++) A->row[i]=A->row[i]-1;

    return nzaa;
}





/*************************************************************/
int
calculaestruturaAAt(A,At,AAt)
     matrix *A,*At,*AAt;
{
    int i,kbeg,kend,IR,k,j,sil,ii;
    int *marca;
    int status,ancol,aancol;
    int *Aendcol,*Atendcol;
    unsigned long int nzaa;

    /* numera linhas e colunas a partir de 1 */

    for (i=0; i<= A->n; i++) A->col[i]=A->col[i]+1;

    for (i=0; i< A->nnulos; i++) A->row[i]=A->row[i]+1;

    /* aloca memoria para vetores locais */
    if(!(marca=(int*)calloc(A->n, sizeof(int))))
    {
      printf("sem memoria para A->row  \n");
      exit(1);
    }

    nzaa=0;
    AAt->col[0]=1;

    for (i=0; i< A->m; i++)
    {
      AAt->row[nzaa]=i+1;
      kbeg=At->col[i];
      kend=At->col[i+1];

      /* marca quais colunas da linha i que possuem elementos nao nulos */
      for (k=kbeg; k <kend; k++)
      {
         IR=At->row[k-1];
         marca[IR-1]=1;
      }


      /* Procura linhas adjacentes a linha i*/
      for (j=i+1; j< A->m ; j++)
      {
         kbeg=At->col[j];
         kend=At->col[j+1];
         for (k= kbeg ; k< kend; k++)
	 {
           IR=At->row[k-1];
           if (marca[IR-1]==1)
	   {
              nzaa+=1;
              AAt->row[nzaa]=j+1;
              break;
	   }
         }
      }

      nzaa+=1;
      //printf("nzaa: %4lu \n",nzaa);

      /*armazena no do 1o elemento que aparece na coluna  */
      AAt->col[i+1]= nzaa+1;

      /*reinicia vetor para marcar as linhas que contem elementos em cada coluna */
      for (ii=0; ii< A->n ; ii++) marca[ii]=0;


    }


    free(marca);

    AAt->nnulos =nzaa;
    AAt->m=A->m;

    /*retorna estrutura de dados comecando de zero */
    for (i=0; i<= A->n; i++) A->col[i]=A->col[i]-1;

    for (i=0; i< A->nnulos; i++) A->row[i]=A->row[i]-1;

    return 0;
}



/******************************************************************/
/*Calula matriz AAT */
int
calculaAAt(A,At,d,AAt,iwork,dwork)
     matrix *A,*At,*AAt;
     int *iwork;
     double *d, *dwork;
{
   int cont,j,i,kbeg,kend,k,IR,kk,IRR,kkbeg,kkend,lbeg,lend,l,lrow,ii;
   int ncol,nrow;
   double temp, *temp2;


   ncol=A->n;
   nrow=A->m;

   temp2=dwork;
   memset(temp2,0,ncol*sizeof(double));


   cont=-1;
   for(j=0; j < nrow; j++)
   {
     kbeg=At->col[j]-1;
     kend=At->col[j+1]-1;
     /*calcula coluna j da matriz AAt */
     for(k=kbeg;k<kend;k++)
     {
       IR=At->row[k]-1;
       temp=At->val[k]*d[IR];
       kkbeg=A->col[IR]-1 +1;
       kkend=A->col[IR+1]-1 +1;
       for(kk=kkbeg;kk<kkend;kk++)
       {
         IRR=A->row[kk]-1 +1;
	 if(IRR>=j) temp2[IRR]+=temp*A->val[kk];
       }
       temp=0;
     }

     /*armazena coluna j na estrutura AAt */
     lbeg=AAt->col[j]-1;
     lend=AAt->col[j+1]-1;

     for(l=lbeg;l<lend;l++)
     {
       cont+=1;
       lrow=AAt->row[l]-1;
       AAt->val[cont]=temp2[lrow];
       temp2[lrow]=0;
     }

   }
   return 0;

}

/*-----------------------------------------------------------------------*/

/******************************************************************/
/*  Calcula precondicionador                                      */
/*  Faz a chamada para a subrotina FCC em fortran                 */


int
precond(AAt,L,Eta,dwork)
     matrix *AAt,*L;
     double *dwork;
     int Eta;
{
   int       Info, DimL ,*LinRow, *LogRow, *RowCol, ncol1, i;
   double    Eps8, *T, *W;

   L->m=AAt->m;
   DimL=AAt->nnulos+Eta*AAt->m;
   Eps8 = 1e-8;
   ncol1=AAt->m+1;
   Info=0;

   if(DimL<AAt->m) DimL=AAt->m;

   /*aloca memoria para vetores temporarios usados na CCF*/
   if(!(LinRow=(int*)malloc(DimL*sizeof(int))))
   {
     printf("sem memoria para Linrow \n");
     exit(1);
   }


   T=dwork;
   W=T+AAt->m;
   LogRow=(int*)(W+AAt->m);
   RowCol=(int*)(LogRow + AAt->m);


   /*------------------------------------------------------------*/

   /*chamada p/ funcao em fortran */
   /* faz diagonal scaling na matriz AAt */
   dscale( AAt->val, AAt->col, AAt->row, &(AAt->m),&(AAt->nnulos), &Info);

   /*calcula matriz precondicionadora L */
   /* Fatoracao Controlada de Cholesky */


   ccf(&(AAt->m), &(AAt->nnulos), AAt->val, AAt->col, AAt->row,
       &(L->nnulos), L->val, L->col, L->row,
       &Eta, &DimL, &Eps8, LinRow, LogRow, RowCol, T, W, &Info);

   //printf("Precondicionador foi calculado \n ");
   //printf("Info: %4i \n",Info);
   fprintf(eigout," ETA: %5i   NZFCC: %8i,   Info: %3i \n",Eta, L->nnulos , Info);


	
   /*libera memoria de vetores temporarios*/
   free(LinRow);
   /*
   free(LogRow);
   free(RowCol);
   free(T);
   free(W);
   */
   return(Info);

}

/*......................................................*/



/******************************************************************/
/*  Dado precondicionador L, com D armazenado na diagonal         */
/*  Resolve o sistema LDLtx=b                                      */


int
solsis(AAt,L,b,x)
     matrix *L, *AAt;
     double *b, *x;
{

   double *v;
   int    i,j,k,kbeg,kend,I,n;

   n=L->m;

   if(!(v=(double *)malloc(n*sizeof(double))))
   {
     printf("sem memoria para v de LDLtx=b  \n");
     exit(1);
   }


   /*Resolve o sistema D^(-1/2)(L*d*L')D^(1/2)x = D^(-1/2)b
   D= diagonal scaling - esta armazenado na diagonal de AAt */


   /* Solve the equation  L * v = b */

   /* Subsuc(L,b,v) */
   /* Multiplica vetor b por D^(-1/2) */
   for(i=0;i < n; i++) v[i]=b[i]*(AAt->val[AAt->col[i]-1]);

   /* Solve the equation  L * v = b */
   /* Subsuc(L,b,v) */

   for(j=0; j<n; j++)
   {
     kbeg=L->col[j]-1;
     kend=L->col[j+1]-1;
     for(k=kbeg+1; k<kend; k++)
     {
       I=L->row[k]-1;
       v[I]=v[I]- L->val[k]* v[j];
     }
   }

   /* Solve the equation  D * y = v */
   /* D esta armazenada na diagonal de L */
   /* armazena y em v para economizar memoria */

   for(i=0;i<n;i++) v[i]=v[i]/(L->val[L->col[i]-1]);

   /* Solve the equation  Lt * x  = y */
   /* y esta armazenado em v */
   /* Subret(L,v,x) */

   for(j=n-1; j>=0 ; j--)
   {
     x[j]=v[j];
     kbeg=L->col[j]-1;
     kend=L->col[j+1]-1;
     for(k=kbeg+1;k<kend;k++)
     {
       I=L->row[k]-1;
       x[j]=x[j]-x[I]*L->val[k];
     }
   }


   /* faz a operacao inversa do diagonal scaling p/ recuperar a solucao */
   for(i=0;i < n; i++)
     x[i]=x[i]*(AAt->val[AAt->col[i]-1]); // x=D^(-1/2)b D=diagonal scaling

   free(v);
   return(0);
}

/*------------------------------------------------------------------------*/
//Rotina para calcular valores de Ritz
/*
int VaRitz(valpha, vbeta, iter, piiter)
     double *valpha, *vbeta;
     int iter,piiter;

{
   double  rmax,rmin,raio;
   int     info,i;
   int     *C;

   C=(int*)NewInt(21,"C:rotina Varitz");
   //C=NewInt(21,"C:rotina Varitz");

   if (iter==-1)
   {
     //salva no arquivo numero de iteracoes do CG e valores de ritz
     fprintf(eigout," IT:%4i  Eta:%3i   GC:%6i  \n\n", piiter, Eta, 0);
     return(0);
   }

   // Rotina em fortran para calculo dos valores de Ritz
   ritzva(valpha, vbeta, &iter, &info, &rmax, &rmin, &raio,C);

   //libera memoria alocada em pcg
   Free((char *) valpha);
   Free((char *) vbeta);

   if(info == -5)
   printf("Problemas no calculo dos valores de Ritz: dsterf falhou \n");

   //imprime  numero de iteracoes do CG e valores de ritz

   if(trocou==0)
   {
      //printf("IT: %4i,   , GC:%6i,  RMAX:%12.5e,  RMIN:%12.5e, RAIO:%12.5e \n", piiter, iter, rmax,rmin,raio);

     //salva no arquivo numero de iteracoes do CG e valores de ritz
     fprintf(eigout,"IT %4i,  GC:%6i, RMAX:%12.5e, RMIN:%12.5e, RAIO:%12.5e \n\n", piiter,  iter, rmax,rmin,raio);


     fprintf(eigout, "  -inf   1e-8   1e-7   1e-6   1e-5   1e-4   1e-3   1e-2   1e-1   .316   1.00 \n  %6i %6i %6i %6i   %6i %6i %6i %6i %6i %6i \n", C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]);

    fprintf(eigout, "  1.00   3.16   1e01   1e02   1e03   1e04   1e05   1e06   1e07   1e08   +inf  \n  %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i \n---------------------------------------------------------------------------------------------------\n", C[10],  C[11], C[12], C[13], C[14], C[15], C[16], C[17], C[18], C[19]);
    }

  else
    {
     //printf("IT: %4i, %4s ,  GC:%6i,  RMAX:%12.5e,  RMIN:%12.5e, RAIO:%12.5e \n", piiter,"      ",iter,rmax,rmin,raio);

     fprintf(eigout,"IT: %4i,  GC:%6i,  RMAX:%12.5e,  RMIN:%12.5e, RAIO:%12.5e \n",	    piiter,iter,rmax,rmin,raio);


     fprintf(eigout, "  -inf   1e-8   1e-7   1e-6   1e-5   1e-4   1e-3   1e-2   1e-1   .316   1.00 \n  %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i \n", C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]);

     fprintf(eigout, "  1.00   3.16   1e01   1e02   1e03   1e04   1e05   1e06   1e07   1e08   +inf  \n  %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i \n---------------------------------------------------------------------------------------------------\n", C[10],  C[11], C[12], C[13], C[14], C[15], C[16], C[17], C[18], C[19]);
   }


   Free((char *) C);
   return(0);
}*/

