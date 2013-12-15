#include "splitting.h"

/* #define ccf ccf_ */

/*************************************************************/
int
calculaestruturaAAt(A,At,AAt)
     matrix *A,*At,*AAt;
{
    int i,kbeg,kend,IR,k,j,sil,ii; 
    int *marca; 
    int  nzaa,status,ancol,aancol;
    int *Aendcol,*Atendcol;


    ancol=A->m+1;
    aancol= (A->m*A->m + A->m)/2;
    /*
     printf("Entrei na cal estru \n");   
     printf("A->n:  %4i \n" , A->n); 
     printf("A->m:  %4i \n" , A->m);
     printf("A->nnulos:  %4i \n" , A->nnulos);
    */
     
    /* numera linhas e colunas a partir de 1 */
    for (i=0; i<= A->n; i++)
      {
       A->col[i]=A->col[i]+1;
       /*printf(" %4i" , A->col[i]); */
      }

     for (i=0; i< A->nnulos; i++)
      {
       A->row[i]=A->row[i]+1;
       /* printf(" %4i" , A->row[i]); */
      }

   



   /*aloca memoria para os ponteiros de colunas e linhas e  valores  de At */
    if(!(At->col=calloc(ancol,sizeof(int))))
      {printf("sem memoria para At->col \n");
	exit(1);
      } 
    if(!(At->row=calloc(A->nnulos,sizeof(int))))
      {printf("sem memoria para At->row  \n");
	exit(1);
      } 
    if(!(At->val=calloc(A->nnulos,sizeof(double))))
      {printf("sem memoria para At->val  \n");
	exit(1);
      } 

    /*....................................................*/

    /*aloca memoria para vetores temporarios*/
    if(!(Aendcol=malloc(A->n*sizeof(int))))
      {printf("sem memoria para A->col \n");
	exit(1);
      } 
    if(!(Atendcol=calloc(A->m,sizeof(int))))
      {printf("sem memoria para A->row  \n");
	exit(1);
      } 
   
    /*..................................................*/
    
    printf("alocada na cal estru \n");   
   



     for (i=0; i< A->n; i++)
      {
       Aendcol[i]=A->col[i+1]-1;
       /*printf(" %4i" , Aendcol[i]); */
      }

   
    /*calula transposta da A*/
    status= TransposeSparseRealMatrix(A->val, A->col, Aendcol, A->row,
            At->val, At->col, Atendcol, At->row, &A->m, &A->n);

   


    /*armazena a prosicao m+1 de Atcol*/ 
    At->col[A->m]=Atendcol[A->m-1]+1;

  
    
    /*libera memoria de vetores temporarios */
    free(Aendcol);
    free(Atendcol);



   /*aloca memoria para os ponteiros de linha e colunas de AAt*/
   if(!(AAt->col=malloc(ancol*sizeof(int))))
      {
      printf("sem memoria para AAt->col \n");
	exit(1);
      } 
   if(!(AAt->row=malloc(aancol*sizeof(int))))
      {
      printf("sem memoria para AA->row  \n");
	exit(1);
      } 


       

   /* aloca memoria para vetores locais */
   if(!(marca=calloc(A->n, sizeof(int))))
      {
      printf("sem memoria para A->row  \n");
	exit(1);
      } 
    
   nzaa=-1;
   
   AAt->col[0]=1;

   for (i=0; i< A->m; i++)
     {
       nzaa+=1;
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
           /*printf(" j: %4i ",j); */
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

    
         /*armazena no do 1o elemento que aparece na coluna  */
         AAt->col[i+1]= nzaa+2;

         /*reinicia vetor para marcar as linhas que contem elementos 
         em cada coluna */        
         for (ii=0; ii< A->m ; ii++)
             marca[ii]=0;
             
       }
 

   free(marca);
 
   AAt->nnulos =nzaa+1;

    printf("AAt-col");
   for (ii=0; ii< A->m+1 ; ii++)
   printf(" %4i", AAt->col[ii]); 
   printf("\n");
  
   printf("AAt-row");
   for (ii=0; ii< nzaa+1 ; ii++)
   printf(" %4i", AAt->row[ii]); 
   printf("\n");
   
   
 
  


   /*retorna estrutura de dados comecando de zero */
     for (i=0; i<= A->n; i++)
      {
       A->col[i]=A->col[i]-1;
       /*printf(" %4i" , A->col[i]); */
      }

     for (i=0; i< A->nnulos; i++)
      {
       A->row[i]=A->row[i]-1;
       /*  printf(" %4i" , A->row[i]);*/ 
      }



   return 0;
} 



/******************************************************************/
/*Calula matriz AAT */
int
calculaAAt(A,At,d,AAt)
     matrix *A,*At,*AAt;
     double *d;
{
  int cont,j,i,kbeg,kend,k,IR,kk,IRR,kkbeg,kkend,lbeg,lend,l,lrow,sil,ii;
  int ncol,nrow;
  double temp,*temp2; 


  ncol=A->n;
  nrow=A->m;  

 /*aloca memoria para valores numericos de AAt */
   if(!(AAt->val=malloc(AAt->nnulos*sizeof(double))))
      {
      printf("sem memoria para AAt->col \n");
	exit(1);
      } 
      

  /*aloca memoria para vetor temporario */
    if(!(temp2=calloc(ncol,sizeof(double))))
      {
      printf("sem memoria para temp  \n");
	exit(1);
      } 
  
   
  
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
	     if(IRR>=j)
	       temp2[IRR]+=temp*A->val[kk];
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

  free(temp2);

  printf("\n ADAt->val \n");
  for (i=0; i<= cont; i++)
     {
      printf(" %4f" , AAt->val[i]);
     }
  printf("\n"); 

  return 0;

}

/*-----------------------------------------------------------------------*/

/******************************************************************/
/*  Calcula precondicionador                                      */   
/*  Faz a chamada para a subrotina FCC em fortran                 */ 


int
precond(AAt,L,Eta)
     matrix *AAt,*L;
{
  int Eta, Info, DimL ,*LinRow, *LogRow, *RowCol,ncol1,i;
  double Eps8, *T, *W;


 
  DimL=AAt->nnulos+Eta*AAt->m;
  Eps8 = 0.0000000000001;
  ncol1=AAt->m+1;
 
  /*-----------------------------------------------------------------*/
  /*aloca memoria para os ponteiros de linhas, colunas e valores de L*/
   if(!(L->col=malloc(AAt->m*sizeof(int))))
      {
      printf("sem memoria para L->col \n");
      exit(1);
      } 
 
   if(!(L->row=malloc(DimL*sizeof(int))))
      {
      printf("sem memoria para L>row  \n");
	exit(1);
      } 
   if(!(L->val=malloc(DimL*sizeof(double))))
      {
      printf("sem memoria para L->value  \n");
	exit(1);
      } 


/*------------------------------------------------------------*/
/*aloca memoria para vetores temporarios usados na CCF*/
   if(!(LinRow=malloc(DimL*sizeof(int))))
      {
      printf("sem memoria para Linrow \n");
	exit(1);
      } 
   if(!(LogRow=malloc(AAt->m*sizeof(int))))
      {
      printf("sem memoria para LogRow  \n");
	exit(1);
      } 
   if(!(RowCol=malloc(3*AAt->m*sizeof(int))))
      {
      printf("sem memoria para RowCol  \n");
	exit(1);
      } 

  if(!(T=malloc(AAt->m*sizeof(double))))
      {
      printf("sem memoria para T  \n");
	exit(1);
      } 
   if(!(W=malloc(AAt->m*sizeof(double))))
      {
      printf("sem memoria para W  \n");
	exit(1);
      } 
   /*------------------------------------------------------------*/
   /*calcula matriz precondicionadora L */
   /* Fatoracao Controlada de Cholesky */

   
   printf("CCF \n");   
   /*
   ccf(&(AAt->m), &(AAt->nnulos), AAt->val, AAt->col, AAt->row,
       &(L->nnulos), L->val, L->col, L->row, 
       &Eta, &DimL, &Eps8, LinRow, LogRow, RowCol, T, W, &Info);

   */
	    
     /*libera memoria de vetores temporarios*/
     free(LinRow);
     free(LogRow);
     free(RowCol); 
     free(T); 
     free(W); 

     return(0);
 
}


/*......................................................*/
