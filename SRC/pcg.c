#include <stdlib.h>
#include <string.h>
#include "memory.h"
#include "splitting.h"
#define MAX(a,b) (((a)>(b)) ? (a) : (b))


//extern int marcares; //marcatroca,

int pcg(int *flag, matrix *ADAt, matrix *LL, matrix *A, int n, double *D, double *H, double *x, double *b, decomp *S, int *pivot, int maxit, double tol, double *work, int piiter)
{
   int i, o, info;
   double nb, nr, rho, rho1, ak, bk;
   double *r, *p, *w, *u;
   matrix *L;
   float nz;
   double tol1,*r1,nr1;


   *flag=0;

   o= (n == A->m) ? 0 : A->n;
   L= S->L;

   // variaveis auxiliares
   r= work;
   p= r + n;
   w= p + n;
   u= w + n;


   tol1=tol;


   //considera vetor inicial igual a zero
   // dessa forma r=b


   if(MIX == TRUE)
     solsis(ADAt,LL,b,x); // x=inv(LDLT)*b
   else Mz(S,H,pivot,&o,x,b,u);

   if(!INC && !LU && MIX == FALSE) return(0);



   nb= dnrm2(&n,b,one);
   nb= MAX(1.0, nb);
   tol *= nb;

   ADAtx(A,o,D,x,r,u);  // produto ADAT por x - armazena em r


   if(MIX == TRUE) nz= 4*A->col[A->n] + A->n + L->m;
   else if(LU) nz= 4*A->col[A->n] + A->n + 4*L->col[L->m] - L->m;
        else nz= 4*A->col[A->n] + A->n + 4*L->row[L->m] - 3*L->m;



   //residuo
   dxmy(&n,b,r,r); //diferenca entre b e r - armazena em r


   nr= dnrm2(&n,r,one); //norma do residuo ||b-ADATx||

   if(nr < tol && tol<=1e-8)
   {
     OPS += nz + 6*n;
     /*printf("(O residuo inicial do GG < tol  )\n");
     VaRitz(valpha, vbeta, -1 ,piiter);*/
     return(0);
   }


   if(MIX == TRUE)
    solsis(ADAt,LL,r,w); // w=inv(LDL)r
   else Mz(S,H,pivot,&o,w,r,u);


   rho= ddot(&n,r,one,w,one);// faz produto interno r.w
   memcpy(p,w,n*sizeof(double)); // p=w



   //laco principal

   for(i= 0; i < maxit; i++)

   {
     ADAtx(A,o,D,p,w,u); // w= ADAT*p
     ak= rho/ddot(&n,p,one,w,one);       // rw/pw
     daxpy(&n,&ak,p,one,x,one); // x= x + ak*p
     ak= -ak;
     daxpy(&n,&ak,w,one,r,one); // r = r-ak*w
     nr= dnrm2(&n,r,one);

     // verifica se ja e' satisfeita a condicao de parada

     if(nr < tol || i + 1 == maxit )
     {
       if(nr < tol*1e-4 /*|| (rho1 < tol*1e-4 && nr >= tol)*/) info= TRUE;
       else info= FALSE;
       //calcula a norma do residuo verdadeiro
       ADAtx(A,o,D,x,r,u);   // r=ADAT*x
       dxmy(&n,b,r,r);   // r=b-r
       nr= dnrm2(&n,r,one); //||r||
       printf("Iter do gradiente conjugado: %4i \n", i);
       //acrescentei teste para refinar a solucao no precondicionador separador
       if ((MIX == FALSE) &&  (nr > 1.0e-7) ) *flag=-2;
       if(info || nr < tol)
       {
         OPS += i*(nz+10*n+2) + nz +4*n - 1;
         return(i);
       }
       OPS += i*(nz+10*n+2) + nz +4*n - 1;

       return(i);
     }

     //termina teste da condicao de parada


     rho1= rho;
     if(MIX)
       /* ddivt(&n,r,H,w); */
       solsis(ADAt,LL,r,w);  //w=inv(LDL)r
     else Mz(S,H,pivot,&o,w,r,u);

     rho= ddot(&n,r,one,w,one); // rho=r*w
     bk= rho/rho1;
     daxpyx(&n,&bk,p,one,w,one);  // p=w+bk*p
   }
   return(-1);

}
