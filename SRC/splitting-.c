/* splitting routines
 *
 * Author: Aurelio R. L. Oliveira
 * 
 * 2000 Unicamp.
 */  

#define FILL 16


#include <stdio.h>
#include <math.h>
#include "main.h"
#include "pre.h"
#include "memory.h"
#include "string.h"
#define MAIN
#include "splitting.h" 

#define transforma transforma_  

extern int      Eta, Etamax, aumenta,marca,marcatroca,marcares,GC_Split;
extern FILE     *eigout;//*OUT_POS;

/*****************************************************************/
/* Allocation and deallocation routines for the FactorType data  */
/* structure                                                     */
/* These are specific to the solver since it might require the   */
/* use of the space in FactorType->ptr                           */
/*****************************************************************/

FactorType * 
SplitFactorType(A, NumRows, NumCols)
     MMTtype *A;
     int      NumRows, NumCols;
{
  int             N;
  FactorType     *FactorSpace;
  double          EndUserTime, EndSysTime, StartUserTime, StartSysTime;
  SplittingType   *Splitting, *NewSplittingType();

  FactorSpace = (FactorType *) Malloc(sizeof(FactorType), "FactorSpace");
  
  
  if(LLt)
  {
     FactorSpace->AAT = (MMTtype *) Malloc(sizeof(MMTtype), "AAT");
     ComputeStructureAAT(A, FactorSpace->AAT);
     N = FactorSpace->AAT->NumCols;
  }

  else N= NumRows;

  FactorSpace->Ndense = 0;
  FactorSpace->N = N;
  FactorSpace->Perm = NewInt(N, "Perm");
  FactorSpace->InvPerm = NewInt(N, "InvPerm");
  FactorSpace->SmallDiagonals = NO;

  if(LLt) Splitting= NULL; else Splitting = NewSplittingType(N,NumCols);

  FactorSpace->ptr = Splitting;
  
  FactorSpace->FactorizationCode = (char *) Malloc(50*sizeof(char), 
						  "FactorizationCode");
  strcpy(FactorSpace->FactorizationCode,"Splitting library");

  return FactorSpace;
}

void
FreeSplitType(Factor)
     FactorType *Factor;
{
   void FreeSplittingType();
   
   Free((char *) Factor->Perm);
   Free((char *) Factor->InvPerm);
   if(LLt)
   {
      Free((char *) Factor->maskDense);
      FreeDouble2(Factor->W);
      FreeDouble2(Factor->Ldense);
      Free((char *) Factor->AAT->pBeginRow);
      Free((char *) Factor->AAT->pEndRow);
      Free((char *) Factor->AAT->Row);
      Free((char *) Factor->AAT->Value);
      Free((char *) Factor->AAT);
   }

   Free((char *) Factor->FactorizationCode);
   
   /* THIS IS THE ONLY SOLVER SPECIFIC LINE */
   if(!LLt) FreeSplittingType(Factor->ptr);
   
   Free((char *) Factor);
}

/*****************************************************************/
/* Allocation and deallocation routines for Splitting data       */
/*****************************************************************/

SplittingType *
NewSplittingType(N,NumCols)
     int  N, NumCols;
{
   SplittingType *SplittingSpace;
  

   SplittingSpace = (SplittingType *) Malloc(sizeof(SplittingType),"Splitting");
   SplittingSpace->A = (matrix *) Malloc(sizeof(matrix),"A");
  
    /*acrescentei */
    SplittingSpace->ptr = (fcc *) Malloc(sizeof(fcc),"FCC");  
    SplittingSpace->ptr->At = (matrix *) Malloc(sizeof(matrix),"At");
    SplittingSpace->ptr->ADAt = (matrix *) Malloc(sizeof(matrix),"ADAt"); 
    SplittingSpace->ptr->L = (matrix *) Malloc(sizeof(matrix),"L"); 
   /**/

   SplittingSpace->L = (matrix *) Malloc(sizeof(matrix),"L");
   SplittingSpace->S = (decomp *) Malloc(sizeof(decomp),"S");
   SplittingSpace->S->L= SplittingSpace->L;
   SplittingSpace->A->m = N;
   SplittingSpace->A->n = NumCols;
   SplittingSpace->L->m = SplittingSpace->L->n = N;
   SplittingSpace->L->col= NewInt(N + 1,"L->col");
   SplittingSpace->pos= NewInt(NumCols,"pos");
   SplittingSpace->ant= NewInt(NumCols,"ant");
   SplittingSpace->diag= NewInt(N + 1,"diag");
   SplittingSpace->S->diag = SplittingSpace->diag + 1;
   SplittingSpace->G= NewDouble(NumCols,"G");
   SplittingSpace->H= NewDouble(NumCols,"H");
   SplittingSpace->dwork= NewDouble(NumCols + 4*N,"dwork");
   SplittingSpace->compr= NewInt(N,"compr");
   SplittingSpace->S->comp= SplittingSpace->compr;


   /* space for L will be allocated after the ordering phase */
   /* espaco para ADAt sera alocado depois de montar a estrutura */
   return SplittingSpace;
}

void 
FreeSplittingType(Splitting)
     SplittingType *Splitting;
{
   Free((char *) Splitting->A);
   Free((char *) Splitting->L->col);
   Free((char *) Splitting->L);
   Free((char *) Splitting->S);
   Free((char *) Splitting->pos);
   Free((char *) Splitting->ant);
   Free((char *) Splitting->diag);
   Free((char *) Splitting->G);
   Free((char *) Splitting->H);
   Free((char *) Splitting->dwork);
   Free((char *) Splitting->compr);

   /*acrescentei */
  
   /*
    Free((char *) Splitting->ptr->At->row);
    Free((char *) Splitting->ptr->At->col);
    Free((char *) Splitting->ptr->At->val);
    Free((char *) Splitting->ptr->ADAt->col); 
    Free((char *) Splitting->ptr->ADAt->row); 
    Free((char *) Splitting->ptr->ADAt->val); 
    Free((char *) Splitting->ptr->L->col); 
    Free((char *) Splitting->ptr->L->row); 
    Free((char *) Splitting->ptr->L->val); 
   */

    Free((char *) Splitting->ptr->At); 
    Free((char *) Splitting->ptr->ADAt); 
    Free((char *) Splitting->ptr->L); 




   Free((char *) Splitting->ptr); 
   Free((char *) Splitting);

}

void parametros(k,gap, Factor)
    int k;
    double gap;
    FactorType *Factor;
{
   double *gapa, *tol, *tola;
   SplittingType   *Splitting;
   Splitting= (SplittingType *) Factor->ptr;
   gapa= &Splitting->gapa;
   tol= &Splitting->tol;
   tola= &Splitting->tola;

   Splitting->gap= gap;

   if(k == 0) Splitting->gap0= gap*1e-4;
   //if(k == 0) Splitting->gap0= gap*1e-6;
   else
      if(TOL || INE)
        if(*gapa < gap*2.0)
        {
           if(INE)
           {
              //puts("sqrt(eps)");
              *tola= *tol= 1.4901161e-8;
              INE= FALSE;
              // MIX= FALSE; Splitting->stat= LU= TROCA= TRUE;
           }
           if(TOL == TRUE)
           {
              //puts("1e-4");
              *tola= 1e-4;
              TOL= -TRUE;
           }
           else
           {
              //puts("1e-8");
              *tola= 1.4901161e-8;
              TOL= FALSE;
           }
        }
   *gapa= gap;
}
/*****************************************************************/
/* Space Allocation for matrix LU                                */
/*****************************************************************/

void create_LU(M, Factor)
     MMTtype *M;
     FactorType *Factor;
{
   int N, i, j, k, maxnz, *pivot, *pos, *ant, *diag, *iwork, DimmaxL;
   double *G, *H, *dwork, Mel;
   unsigned long int nnaat;
   matrix *L, *A, *At, *ADAt, *LL;
   decomp *S;
   float ops;
   

   SplittingType   *Splitting;
   Splitting= (SplittingType *) Factor->ptr;

   printf("Create LU \n "); 

   L= Splitting->L;
   A= Splitting->A;
   A->row= NewInt(M->Nonzeros,"A->row");
   pos= Splitting->pos;
   ant= Splitting->ant;
   diag= Splitting->diag + 1;
   G= Splitting->G;
   H= Splitting->H;
   dwork= Splitting->dwork;
   iwork= (int *) dwork;
   
   pivot= Factor->Perm;

   memset(Splitting->ant,0,M->NumCols*sizeof(int));
  
   A->col= M->pEndRow - 1;
   A->col[0]= 0;
   for(maxnz= i= k= 0; k < M->NumCols; k++)
   {
       j= A->col[k+1] - A->col[k];
       maxnz += j*(j + 1);
       pos[k]= k;
       for(G[k]= 0.0; i < M->pEndRow[k]; i++)
       {
           //G[k] -= fabs(M->Value[i]); /* Norma 1 */
           G[k] += M->Value[i] * M->Value[i]; /* Norma 2 */
           //if (G[k] < fabs(M->Value[i])) G[k] = fabs(M->Value[i]); //*Norma infinito/ 
           A->row[i] = M->Row[i] - 1;
       }
       G[k] = -sqrt(G[k]); /* Norma 2*/
       //G[k] = -G[k]; //Norma Infinito
   }
   printf("\n Norma 2 \n");
   fprintf(eigout," Ordenacao pela Norma 2 \n");
   
   maxnz /= 2;
   A->val= M->Value;
   //acresecentei
   maxnz *=10;
 
   /*
   maxnz= M->Nonzeros;
   maxnz *= FILL;
   */
   Splitting->maxnz= maxnz;

   //printf("MAXNZ: %4i \n",Splitting->maxnz);
   /*
   L->row= NewInt(maxnz,"L->row");
   L->val= NewDouble(maxnz,"L->val");
   */

   /* verificar estrutura de dados
    printf("A->m: %4i ,A->n: %4i \n ", A->m,A->n);
    for (i=0;i<10;i++)
      printf("A->row[i]: %4i \n",A->row[i]);
    for (i=0;i<10;i++)
      printf("A->col[i]: %4i \n",A->col[i]);
    for (i=0;i<10;i++)
      printf("A->val[i]: %4f \n",A->val[i]);
   */

   
   
   if(MIX)
   {
     LU= FALSE;
     /*acrescentei*/
     N=Splitting->A->m;
     //printf("N: %4i\n",N);

     At=Splitting->ptr->At;
     ADAt=Splitting->ptr->ADAt;
     LL=Splitting->ptr->L; 
   
     //printf("Aloca \n ");

     /*aloca espaco para matriz At e ADAt */
     At->col = NewInt((N + 1),"At->col");
     At->row = NewInt(M->Nonzeros,"At->row");
     At->val = NewDouble(M->Nonzeros,"At->val");

     /* calcula numero de elementos nnulos em ADAt */
     //printf("Calcula \n ");
 
     A->nnulos=M->Nonzeros;
   
     nnaat=0;
     nnaat=calculannulosAAt(A,At);
         
     /*aloca espaco para matriz ADAt */
     ADAt->col= NewInt(N+1,"ADAt->col"); 
     ADAt->row= NewInt(nnaat,"ADAt->row");         
   
     /*calcula estrutura da matriz ADAt */
     calculaestruturaAAt(A,At,ADAt);  

     /*aloca espaco para armazenar valores numericos para ADAt */     
     ADAt->val= NewDouble(ADAt->nnulos,"ADAt->val");
     
     /*calcula valor inicial para Eta*/
     Mel=(nnaat/N);
     printf("Mel: %7.4lf \n",Mel);
     //Eta= 20 - floor(Mel);
     //if (N < 1000) Eta=N;

    
     //printf("AAt->m: %4i,   Eta: %4i \n", ADAt->m,Eta);
     //if (Eta < -50) Eta=-50;
     //printf("AAt->m: %4i,   Eta: %4i \n", ADAt->m,Eta);

     
     //if (Eta < -250) Eta= -ADAt->m; 

     fprintf(eigout,"NZAAT:%4i   N:%4i   Mel:%7.4f \n",nnaat, N, Mel);
     fprintf(eigout,"---------------------------------------------------------------------------------- \n"); 

//Eta=-N; //  modifiquei aqui jair

     if (Mel<10) Eta=Mel;  //mudei aqui jair, troquei Mel por -N
     else  Eta=-Mel;
     //Eta=10;
     //Eta+=100;
     //Eta+=5;
	     //Eta += 70; //nug15
     //Eta=5;
     //Eta -= 10;
     printf("Eta inicial: %7.4lf \n",Eta); 
   
     Etamax=(100000000 - ADAt->nnulos)/ADAt->m;
     if (Etamax > ADAt->m) Etamax=ADAt->m;
     if (Etamax < Eta) Etamax=Eta; 

     //printf("Eta:%4i,    ETAMAX: %4i \n",Eta, Etamax); 
     fprintf(eigout,"Eta:%4i,    ETAMAX: %4i \n",Eta, Etamax);
     DimmaxL=ADAt->nnulos + Etamax*ADAt->m;
   
     if (DimmaxL <=0) DimmaxL=ADAt->nnulos;
     //DimmaxL *= 2; // Multiplica o espa�o alocado por 2
     LL->col=L->col;
     LL->row=NewInt(DimmaxL,"L->row");
     LL->val=NewDouble(DimmaxL,"L->val");      
         
   }
  
   Splitting->tol= (INE) ? 1e-4 : 1.4901161e-8;
   Splitting->tola= (TOL) ? 1e-2 : Splitting->tol;
   // Create a file for saving A matrix, name= problem_A.txt
   //fprintf_A_sparse(Splitting->A);
   
}

/*****************************************************************/
/* Driver for solving the linear system                          */
/*****************************************************************/

void splitting(D, M, Factor, dx, dz, x, z,Iteration)
     double *D;
     MMTtype *M;
     FactorType *Factor;
     double *dx, *dz, *x, *z;
     int Iteration;
{
   matrix *A, *L, *At, *ADAt,*LL;
   decomp *S;
   int i, j, stat, maxnz, *pos, *ant, *diag, *iwork, *pivot;
   int *temp_pos; // Variavel para verificar se pos esta correta
   int rec;
   int k,info,Mel;
   int *compr;
   double *G, *H, *dwork, *D_Sqr;
   float ops, luops;
   SplittingType *Splitting;
   Splitting= (SplittingType *) Factor->ptr;
   


   if(Eta > Etamax)
   {
     Eta=Etamax;
     printf(" \n Aumente o valor de Etamax \n");  
   }    


  /*acrescentei*/
   At=Splitting->ptr->At;
   ADAt=Splitting->ptr->ADAt;
   LL=Splitting->ptr->L;
   
   Splitting->ptr->iter=Iteration;  
   maxnz= Splitting->maxnz;
   A= Splitting->A;
   L= Splitting->L;
   S= Splitting->S;
   pos= Splitting->pos;
   ant= Splitting->ant;
   diag= Splitting->diag + 1;
   G= Splitting->G;
   H= Splitting->H;
   compr= Splitting->compr;
   dwork= Splitting->dwork;
   iwork= (int *) dwork;

   pivot= Factor->Perm;
   
   /*
    for (i=0;i<A->n;i++)
       printf("D[i]: %12.8f \n",D[i]);

      for (i=0; i<=A->n ; i++)
        printf("A->col[i]: %4i \n",A->col[i]);


     for (i=0; i < A->nnulos ; i++)
        printf("A->row[i]: %4i \n",A->row[i]);&& (aumenta>=3))


     for (i=0; i < A->nnulos ; i++)
        printf("A->val[i]: %4f \n",A->val[i]);


      getchar();
  */
   
   Mel=floor(ADAt->nnulos / ADAt->m);
   
   if(MIX)
   { 
      /* calcular ADAt chamada na rutina outros.c*/
      calculaAAt(A,At,D,ADAt,iwork,dwork);
      /* CALCULAR PRECOND chamada na rutina outros.c*/  
      info=precond(ADAt,LL,Eta,dwork);
      if ( info < 0) info=precond(ADAt,LL,-Mel,dwork);
      Splitting->ptr->info=info;
      Splitting->corretor= FALSE;
      return;
   }
   if (!Splitting->stat) return; // Se o precondicionar esta ok, n�o precisa organizar e criar uma nova organizacao do B
   
   fprintf(eigout," Troca  %8i \n",Iteration);
   printf("\n Reordenamento das colunas\n");
   //D_Sqr = (double*) malloc((A->n)*sizeof(double));
   //for (i = 0; i < A->n; i++) 
     // D_Sqr[i] = D[i]*sqrt(D[i]);
      //D_Sqr[i] = D[i]*D[i];
   //ddivt(&A->n,G,D_Sqr,H);
   //free(D_Sqr);
   ddivt(&A->n,G,D,H); //Norma 2, Norma Infinito
   /*if(dx == NULL) { L->row[0]= A->n; 2L->row[1]= L->row[2]= 0;}
   else
    {
     for(i= j= k= stat= rec= 0; i < A->n; i++)
      if(signd(dx[pos[i]]) != signd(dz[pos[i]])) // testes para troca de fases
         if(z[pos[i]] <= x[pos[i]]) iwork[j++]= pos[i];    /* B1 
         else iwork[A->n+rec++]= pos[i];                   /* N1 
      else
         if(z[pos[i]] > x[pos[i]]) L->row[stat++]= pos[i]; /* B2 
         else L->row[A->n+k++]= pos[i];                    /* N2 
     memcpy(pos,iwork,j*sizeof(int));
     memcpy(&pos[j],&L->row[A->n],k*sizeof(int));
     memcpy(&pos[j+k],L->row,stat*sizeof(int));
     memcpy(&pos[j+stat+k],&iwork[A->n],rec*sizeof(int));
     L->row[0]= j; L->row[1]= k; L->row[2]= stat;
     printf("B1: %d, B2: %d, N2: %d, N1: %d\n",j,stat,k,rec);

    }*/
   stat= colorder(A,pos,H,dwork,diag,L->row,ant,(int *) L->val); /*Ordena as colunas (lu.c)*/
   if(stat < A->m)
   {
      i= maxnz/2 + A->m;
      stat= luret(A,&L,stat,&i,pivot,pos,diag,iwork,H,&luops,FALSE,TRUE,TRUE);
      rec= TRUE;
      while(stat < 0)
      {
         stat= 1 - stat;
         degree(A,pos,iwork,NULL,stat);
         i= stat;
         //puts("Estamos no while, socorro!");
         stat= luret(A,&L,stat,&maxnz,pivot,pos,diag,iwork,H,&luops,INC,FALSE,
                     rec);
         if(i >= stat) rec= FALSE; else rec= TRUE;
         ops += luops;
      }
      //printf("nzL= %d\n",L->col[L->m]);
      fflush(stdout);
   }
   //printf("Colunas LI encontradas %d\n", stat);
   /*//Begin teste
   for(j= 0, i= 0; i < A->m; i++)
   {
       L->col[i]= j;
       for(k= A->col[pos[i]]; k < A->col[pos[i]+1]; k++)
           L->row[j++]= A->row[k];
   }
   L->col[A->m]= j;
   temp_pos = (int*) malloc((A->n)*sizeof(int));
   i= match(A->m,A->m,L->col,L->row,&L->col[1],temp_pos,iwork);
   if(i != A->m) {printf("funcao luret: i= %d diferente de m\n",i); error(1);};
   // End Teste*/
   //degree(A,pos,iwork,NULL,A->m);
   
   // Teste Begin 01/02/07
   // Teste para verificar que pos tem todos 01/02/07
   /*temp_pos = (int*) malloc((A->n)*sizeof(int));
   for (i = 0; i < A->n; i++) temp_pos[i] = FALSE;
   for (i = 0; i < A->n; i++) temp_pos[pos[i]] = TRUE; 
   for (i = 0; i < A->n; i++) 
     if (!temp_pos[i]) printf("Nao tem a coluna %d depois de block\n",i);
  
   for (i = 0; i < A->n; i++)
         fprintf(OUT_POS,"%d  \n",pos[i]);
   fclose(OUT_POS);
   printf("Estou fechando OUT_POS");
   error(1);
   // Teste End 01/02/07*/
   //printf("Estou entrando no block...\n");    
   block(&(S->cc),&stat,A,pos,diag,(int *)L->val,TRUE);
   //printf("Estou saido do block...\n");
   /*for (i = 0; i< A->m; i++)
     if (!temp_pos[pos[i]]) printf("Mudou Key \n");
   printf("O valor de *n %d depois de block\n", stat);
   for (i = 0; i < A->n; i++) temp_pos[i] = FALSE;
   for (i = 0; i < A->n; i++) temp_pos[pos[i]] = TRUE; 
   for (i = 0; i < A->n; i++) 
     if (!temp_pos[i]) printf("Nao tem a coluna %d depois de block\n",i);
   //free(temp_pos);*/

   S->scc= (stat < A->m) ? blocktri(stat,A,L,&pos[stat],pivot,diag,(int *)L->val,
                   iwork,maxnz) : 0;
   //S->scc= 0; stat = 0;
   //printf("%d scc\n",S->scc);
   for(j= 0; j < stat; j++) compr[diag[j]]= 0;
   for(rec= 1; j < A->m; j++)
   {
      if(j >= pivot[rec]) rec++;
      compr[diag[j]]= rec;
   }
   //j = maxnz; Teste 30/01/07
   lufact(A,&L,&maxnz,pivot,pos,S->cc,stat,diag,iwork,H,&ops,compr);
   
   /*/Teste Begin 30/01/07
   if (maxnz < 0) 
   {
     maxnz= j;
     stat=luret(A,&L,0,&maxnz,pivot,pos,diag,iwork,H,&luops,FALSE,TRUE,FALSE);
     rec= TRUE;
      while(stat < 0)
      {
         stat= 1 - stat;
         printf("Colunas LI %d",stat);
         degree(A,pos,iwork,NULL,stat);
         i= stat;
         puts("Estamos no while, socorro!");
         stat= luret(A,&L,stat,&maxnz,pivot,pos,diag,iwork,H,&luops,FALSE,FALSE,
                     rec);
         if(i >= stat) rec= FALSE; else rec= TRUE;
         ops += luops;
      }
     printf("Colunas LI %d",stat);
     error(1);
   }
   //Teste End 30/01/07 */
   compr[0]= S->cc;
   compr[S->scc]= A->m;
   memcpy(iwork,ant,A->n*sizeof(int));
   memset(ant,0,A->n*sizeof(int));
   for(j= 0; j < A->m; j++) ant[pos[j]]= TRUE;
   for(i= j= 0; j < A->m; j++) { ant[pos[j]]= TRUE;
       if(iwork[pos[j]]) i++; }
   /*printf("Repetidos= %d\n",i);
   printf("nzLU= %d\n",L->col[L->m]);*/
   Splitting->maxnz= maxnz;
   gather(A->m,pos,H,D);
   fprintf(eigout," NZSPLIT:  %8i ",L->col[L->m]);
   
}
/************************************************************************/
/*
/* Fun��o Solve que chama ao gradiente conjugado, resolve o sistema ADAt 
/* A chamada � feita por SolveADAT (chamada por SolveAugmented, chamada 
/* por ComputeGondzioCorrections, ComputeCorrector e ComputePredictor)   
/*
/************************************************************************/

void solve(Factor, rhs, Solution, D, corretor,teste)
      FactorType *Factor;
      double *rhs, *Solution, *D;
      int corretor, teste;
{
   matrix *A, *L, *LL,*ADAt;
   decomp *S;
   int stat, *pivot, i,iter, maxit,info,maxnz,Mel,stat1, flag;
   double tola, *H, *dwork;
 

   SplittingType *Splitting;
   Splitting= (SplittingType *) Factor->ptr;

   A= Splitting->A;
   L= Splitting->L;
   S= Splitting->S;
   tola= Splitting->tola;
   H= Splitting->H;
   dwork= Splitting->dwork;
   maxnz=Splitting->maxnz;

   /*acrescentei */
   LL=Splitting->ptr->L;
   ADAt=Splitting->ptr->ADAt;
   iter=Splitting->ptr->iter;   
   info=Splitting->ptr->info;   
   maxit=A->m;
 
   Mel=floor(ADAt->nnulos / ADAt->m);

   /* Tolerancia para o gradiente conjugado. Se a dire��o � corretora a tolerancia � menor*/
   if (MIX)
      tola=1e-4;
   else
      tola=1e-4;
   if (marca==1 && teste==1) //FriInf < 10e-5 e direcao corretora
	tola=1e-8;
        //tola=1e-7; //nug15
   
   pivot= Factor->Perm; 

   //printf("tolerancia  %7.2e : \n",tola);
   
   stat= pcg(&(flag),ADAt,LL,A,A->m,D,H,Solution,rhs,S,pivot,maxit,tola,dwork,iter);
   //printf("STAT %4i , FLAG: %4i \n",stat, flag);  
   fprintf(eigout,"\n GC:%6i, \n",  stat);
   fprintf(eigout,"__________________________________________________________________________________________________________\n");

   // Refina a solucao na direcao corretora calculada pelo split se tol>1e-7.
   if (flag==-2 && teste==1) 
   {
      printf("Refinamento da solucao  \n");
      stat1= refina(ADAt,LL,A,A->m,D,H,Solution,rhs,S,pivot,(maxit/2),tola,dwork,iter);
      //printf("STAT1 refina  %4i \n",stat1);
      fprintf(eigout,"STAT1 refina  %4i \n",stat1);
      marcares=0;
   }

  
   /*
   if ((teste==1) && (iter >=1))
     transforma(&iter, &(A->m), &(ADAt->nnulos), ADAt->row, ADAt->col, ADAt->val, rhs, Solution); 
   */ 
 
   if(stat < 0) error(3);
   Splitting->stat= (stat*8 < A->m) ? FALSE : TRUE;
   if(corretor) corretor= Splitting->corretor;
   if(corretor && MIX) 
   {
     
     //Troca da Marta
     // testes para troca de fases
     // ITCG >= A->m/4 ou PriInf < 10^-5
 //Eta=12 ;//mudei aqui jair
     

     if(stat*5>= A->m) //|| marca==1) 
       if (Eta <10) //aumenta < 3 && (stat*6 >= A->m &&) // mudei aqui jair troquei < por >
       {
         Eta=Eta+10;             
         aumenta++;
         printf("Aumenta = %d \n",aumenta);
       }
       else
         
       {
          
          fprintf(eigout,"Parametros da troca IGC*5>A->m e Eta < 10 \n");			
          MIX= FALSE; 
          Splitting->stat= LU= TROCA= TRUE;
          //libera memoria usado no prec fcc
          Free((char *) Splitting->ptr->At->row);
          Free((char *) Splitting->ptr->At->col);
          Free((char *) Splitting->ptr->At->val);
          Free((char *) ADAt->col); 
          Free((char *) ADAt->row); 
          Free((char *) ADAt->val); 
          Free((char *) LL->row); 
          Free((char *) LL->val); 
          // aloca memoria p/ split 
          L->row=NewInt(maxnz,"L->row");
          L->val=NewDouble(maxnz,"L->val");    
       }
       /*/
       // testes para troca de fases da Silvana
       // ITCG >= A->m/4 ou PriInf < 10^-5 //////////////O que � PriInf?
       if(stat*4>= A->m || marca==1) 
       {

	  if((( stat*2>= A->m) && (aumenta>=3))  ||
	     (stat*4>=A->m && A->n < 16000 && info >=10)  ||
	     (marca==1 && stat*4 >=A->m)) 
          {			
              MIX= FALSE; 
              Splitting->stat= LU= TROCA= TRUE;
              //libera memoria usado no prec fcc
              Free((char *) Splitting->ptr->At->row);
              Free((char *) Splitting->ptr->At->col);
              Free((char *) Splitting->ptr->At->val);
              Free((char *) ADAt->col); 
              Free((char *) ADAt->row); 
              Free((char *) ADAt->val); 
              Free((char *) LL->row); 
              Free((char *) LL->val);
              // aloca memoria p/ split 
              L->row=NewInt(maxnz,"L->row");
              L->val=NewDouble(maxnz,"L->val");    
          }
           
	  if (stat*4 >= A->m && (Eta <=100) )
          {
	     aumenta=aumenta+1;
             printf("Aumentando o valor de Eta: \n");
             if (info >= 12 &&  Eta>=0)	
	        Eta=0;
	     else
	        Eta=Eta+10;      
	     //printf("%4i \n ",Eta);             
	  }
      	}*/
     else
       aumenta = 0;     
   }

     Splitting->corretor= MIX;
}


#define INDEPENDENTE 9
 
void rank(LP, Record, top, Pass)
     LPtype *LP;
     ChangeStack *Record;
     int *top, *Pass;
{
   int            i, j, stat, maxnz, *pos, *ant, *pivot, *diag, *iwork;
   int            cols, rows, inx, ents, *RowMap;
   int            *RowMask, *ColumnMask;
   double         *H, *dwork;
   matrix         *A, *L;
   float          ops;
   DuplicateRow   *pDuplicateRow;
   //char           POS_file[100]; // Variave teste
   //FILE           *Pos_in; // Variavel teste
   //int *temp_pos; // Variavel para verificar se pos esta correta; teste

   RowMask= Record->RowMask + 1;
   ColumnMask= Record->ColumnMask + 1;

   options();

   A= (matrix *) Malloc(2*sizeof(matrix),"A");
   L= A + 1;

   ents= LP->Ents;
   rows= LP->Rows;
   cols= LP->Cols;
   L->m= L->n= A->m= rows;
   A->n= cols;
   A->row= NewInt(A->n+A->m+1+ents,"A->row");
   A->col = LP->A.pEndRow - 1;
   A->col[0]= 0;
   for(i= 0; i < ents; i++) A->row[i]= LP->A.Row[i] - 1;
   A->val= LP->A.Value;
   diag= A->row + ents + 1;
   pos= diag + A->m;
  
   if(reduz(A,diag,pos,RowMask,ColumnMask))
   {
      Free((char *) A->row);
      Free((char *) A);
      return;
   }
   
   RowMap = NewInt(LP->Rows, "RowMap");
   
   for(i= rows= 0; i < LP->Rows; i++)
   {
       if(diag[i] > -1)
          if(RowMask[i] == STILL_ACTIVE) RowMask[i]= INDEPENDENTE; else ;
       else
          if(RowMask[i] == STILL_ACTIVE) RowMap[i]= rows++;
   }

   ents = 0;
   for(i= 0, cols= 0; i < LP->Cols; i++)
   {
       if(pos[i] > -1)
          if(ColumnMask[i] == STILL_ACTIVE) ColumnMask[i]= INDEPENDENTE; else ;
       else if(ColumnMask[i] == STILL_ACTIVE)
            {
               cols++;
               for(j= LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++)
                   if(RowMask[LP->A.Row[j-1]-1] == STILL_ACTIVE) ents++;
            }
   }

   L->m= L->n= A->m= rows;
   A->n= cols;
   for(maxnz= i= 0; i < LP->Cols; i++)
   {
       j= LP->A.pEndRow[i] - LP->A.pBeginRow[i] + 1;
       maxnz += j*(j + 1);
   }
   if(maxnz > 8*A->n) maxnz /= 2;
   // Multiplicado por 4 o valor de maxnz dia 05/12/2006
   maxnz*=4;
   A->col= NewInt(A->n+1,"A->col");
   A->val= NewDouble(ents,"A->val");
   
   for(i= cols= ents= 0; i < LP->Cols; i++)
       if(ColumnMask[i] == STILL_ACTIVE)
       {
          A->col[cols]= ents;
          // add this column to the reduced LP 
          for(j= LP->A.pBeginRow[i]; j <= LP->A.pEndRow[i]; j++)
          {
              inx= LP->A.Row[j-1]-1;
              if(RowMask[inx] == STILL_ACTIVE)
              {
                 // add this element to the reduced matrix 
                 A->row[ents]= RowMap[inx];
                 A->val[ents++]= LP->A.Value[j-1];
              }
          }
          cols++;
       }
       else if(ColumnMask[i] == INDEPENDENTE) ColumnMask[i]= STILL_ACTIVE;
   A->col[cols]= ents;
   
   for(maxnz= i= 0; i < LP->Cols; i++)
   {
       j= LP->A.pEndRow[i] - LP->A.pBeginRow[i] + 1;
       maxnz += j*(j + 1);
   }
   if(maxnz > 8*A->n) maxnz /= 2;
   // Multiplicado por 4 o valor de maxnz dia 05/12/2006
   maxnz*=4;
   L->col= NewInt(A->n+2*A->m+1+maxnz,"L->col");
  
   L->row= L->col + A->m + 1;
   ant= L->row + maxnz;
   pivot= ant + A->n;
  
   L->val= NewDouble(3*A->n+8*A->m+maxnz,"L->val");
   H= L->val + maxnz;
   dwork= H + A->n;
   iwork= (int *) dwork;
   stat= 0;
   for(i= 0; i < A->n; i++)
   {
       pos[i]= i;
       H[i]= A->col[i] - A->col[i+1];
       ant[i]= 0;
   }
   /*// Teste Begin 01/02/07: Leer o pos do arquivo e fazer em H a ordena�ao
   strcpy(POS_file,"./OUT/Pos_file");
   Pos_in = fopen(POS_file,"r");
   for(i= 0; i < LP->Cols; i++)
   {  
      fscanf(Pos_in, "%d", &(pos[i]));
      //printf("pos %d %d\n",i,pos[i]);
      H[pos[i]] = -i;
       
   }
   fclose(Pos_in);
   printf("\n Valor de LP->Cols %d", LP->Cols);
   printf("\n Valor de A-n %d", A->n);
   printf("\n Lei o pos");
   /*temp_pos = (int*) malloc((A->n)*sizeof(int));
   for (i = 0; i < A->n; i++) temp_pos[i] = FALSE;
   for (i = 0; i < A->n; i++) temp_pos[pos[i]] = TRUE; 
   for (i = 0; i < A->n; i++) 
     if (!temp_pos[i]) printf("Nao tem a coluna %d depois da leitura\n",i);
   // Teste End 01/02/07*/
   stat= colorder(A,pos,H,dwork,diag,L->row,ant,(int *) L->val);
   if(stat != A->m)
   {
      i= maxnz/2 + A->m;
      printf("\n no comeco 1\n\n");
      stat= luret(A,&L,stat,&i,pivot,pos,diag,iwork,H,&ops,FALSE,TRUE,FALSE);
      if(i > maxnz) maxnz= i;
      while(stat < 0)
      {
         stat= 1 - stat;
         degree(A,pos,iwork,NULL,stat);
         printf("\n no comeco 2\n\n");
         stat= luret(A,&L,stat,&maxnz,pivot,pos,diag,iwork,H,&ops,FALSE,FALSE,
                     FALSE);
      }
      if(stat != A->m)
      {
         //printf("rank def, rank= %d --- nzL= %d\n",stat,L->col[stat]);
         for(i= j= 0; i < LP->Rows; i++)
         {
             if(RowMask[i] == STILL_ACTIVE) RowMap[j++]= i;
             if(RowMask[i] == INDEPENDENTE) RowMask[i]= STILL_ACTIVE;
         }
         for(j= 0; j < A->m; j++)
         {
             if(pivot[j] == A->m)
             {
                i= RowMap[j];
                Record->StackOfChanges[*top] = (SingleChange *)
                Malloc(sizeof(SingleChange), "Record->StackOfChanges[]");
                Record->StackOfChanges[*top]->ChangeType = DUPLICATE_ROW;
                pDuplicateRow = (DuplicateRow *)
                Malloc(sizeof(DuplicateRow), "pDuplicateRow");
                pDuplicateRow->Split1= 1;
                pDuplicateRow->Split2= i+1;
                Record->StackOfChanges[*top]->pDuplicateRow = pDuplicateRow;
                // increment top-of-stack pointer, resizing the stack if necessary 
                (*top)++;
                Record->Top = *top;
                if (Record->Top >= Record->Size) ResizeRecord(Record);
                RowMask[i] = *Pass;
             }
         }
      }
   }
   if(stat == A->m)
       for(i= 0; i < LP->Rows; i++)
          if(RowMask[i] == INDEPENDENTE) RowMask[i]= STILL_ACTIVE;
   
   Free((char *) RowMap);
   
   Free((char *) A->row);
   Free((char *) A->val);
   Free((char *) L->col);
   Free((char *) L->val);
   if(stat != A->m) (*Pass)++; 
   Free((char *) A);
}
