#include <stdlib.h>
#include <string.h>
#include "memory.h"
#include "splitting.h"
#define MAX(a,b) (((a)>(b)) ? (a) : (b))


int refina(matrix *ADAt, matrix *LL, matrix *A, int n, double *D, double *H, double *x, double *b, decomp *S,
        int *pivot, int maxit, double tol, double *work, int piiter)

{
   int     i, o, info;
   double  nb, nr, rho, rho1, ak, bk;
   double  *r, *p, *w, *u;
   matrix  *L;
   float   nz;


   o= (n == A->m) ? 0 : A->n;
   L= S->L;


   // variaveis auxiliares
   r= work;
   p= r + n;
   w= p + n;
   u= w + n;

   //considera vetor inicial igual a x
   // calcula residuo
   ADAtx(A,o,D,x,r,u);  // produto ADAT por x - armazena em r
   dxmy(&n,b,r,r); //diferenca entre b e r - armazena em r
   nr= dnrm2(&n,r,one); //norma do residuo ||b-ADATx||
   //printf("residuo refina  %7.2e : \n",nr);

   //testa se o residuo inicial satisfaz
   if(nr < tol)
   {
     printf("(O residuo inicial do GG < tol  )\n");
     return(0);
   }

   if(MIX == TRUE)
       solsis(ADAt,LL,r,w);    // x=inv(LDLT)*r
   else Mz(S,H,pivot,&o,w,r,u);

   if(!INC && !LU && MIX == FALSE) return(0);

   rho= ddot(&n,r,one,w,one);// faz produto interno r.w
   memcpy(p,w,n*sizeof(double)); // p=w

   //laco principal

   for(i= 0; i < maxit; i++)
   {
     ADAtx(A,o,D,p,w,u); // w= ADAT*p
     ak= rho/ddot(&n,p,one,w,one);    // r.w / p.w
     daxpy(&n,&ak,p,one,x,one); // x= x + ak*p
     ak= -ak;
     daxpy(&n,&ak,w,one,r,one); // r = r-ak*w
     nr= dnrm2(&n,r,one);
     //printf(" norma res:  %8.4e \n", nr);
     // verifica se ja'  e' satisfeita a condicao de parada
     if(nr < tol || i + 1 == maxit )
     {
       //calcula a norma do residuo verdadeiro
       ADAtx(A,o,D,x,r,u);   // r=ADAT*x
       dxmy(&n,b,r,r);   // r=b-r
       nr= dnrm2(&n,r,one); //||r||
       //printf(" res final refina :  %8.4e \n", nr);
       return(i);
     }
     rho1= rho;
     if(MIX)
       solsis(ADAt,LL,r,w);    //w=inv(LDL)r
     else
       Mz(S,H,pivot,&o,w,r,u);

     rho= ddot(&n,r,one,w,one); // rho=r*w
     bk= rho/rho1;
     daxpyx(&n,&bk,p,one,w,one);  // p=w+bk*p

   }

   //printf("nr= %e\n",nr/nb);
   return(-1);

}
