#include <stdio.h>
#include <math.h>
#include "main.h"
#include "pre.h"
#include "memory.h"
#include "string.h"
#include <time.h>

double dnrm2(int *n, double *x, int *incx);
double ddot(int *n, double *dx, int *incx, double *dy, int *incy);
void   sortd(double *va, int *v, int n);
void   both(MMTtype *M, int *Acol, double *v1, double *v2, double *x, double *v, double *y, double *w);
void   cholesky(int aux6, int pp2, int p, double s1, double v, double valgaux[], double valg[], double valg1[], double wwaux[], double l1[], double rr1[], double rr2[], double rr4[], double *rr3, double u[], double dd12[]);
void   resolucao_sistema(int aux6, int p, double s1, double sti[], double ssf[], double ssf2[], double dd12[], double valg1[], double rr4[]);

extern int param_p;
extern double param_t;

void jair1(MMTtype *M, FactorType *Factor, double b[], double c[], double *px, double *pi, double *ps, double *pw, double *pr, double upbound[], int NumBounds, int BoundIndex[])
{
	int    i, aux, auxx, aux1, aux2, aux3, aux4, aux5, aux6, aux7,ccc;
	int    j, k, zz, cc, p, mxzero, mxzero1, xzero, l, cci, pp2, p_i, pj, pk, onem6, ONE, tam, wwa;
	int    *dl1posi, *duposi, *r1r2, *ww, *Acol;
	double *x, *s, *w, *r, *vv1;
	double c1, c11, s1, somat,  v, alpu, alpl1, auxl1, ep1, ep2, resiff, resif1, dif, rr3, normabc,d__1,absxi,scale, auxvn,  resiff2, l1, resiff0;
	double *B1, *y1, *y2, *B2, *B,*somalin,*xgrande,*wb, *Paux1, *Paux2, *rr1, *rr2, *rr4, *u, alpp;
	double *val1,*val2,*resi1,*resi2,*resiv,*resiu,*w2,*uu1,*uu2, *valg, *valg1, *wwaux, *dd12,*valgaux, *resif, *normau;
	double *sti, *sts, *ssf, *du, *dl1, *N1, *sti2, *sts2, *ssf2, *dl1mzero, *dumzero, *fornang;
	clock_t start, end;

        start  = clock();
//param_p=2;
//param_t=10;

/*	if((M->NumRows) < 1000)
	{
		param_p = 2;
		ep1     = 0.0000000000000001;
		ep2     = 0.000000000000000001;
	}
	else
	{
		if((M->NumRows) < 10000)
		{
			param_p = 4;
			ep1     = 0.000000000000001;
			ep2     = 0.00000000000000001;
		}
		else
		{	   
			if((M->NumRows) < 20000)
			{
				param_p = 8;
				ep1     = 0.0000000000000001;
				ep2     = 0.000000000000000001;
			}
			else
			{
				if((M->NumRows) < 100000)	
				{
					param_p = 20;
					ep1     = 0.000000000001;
					ep2     = 0.00000000000001;
				}
				else
				{
					if((M->NumRows) < 150000)	
					{
						param_p = 40;
						ep1     = 0.00000000000001;
						ep2     = 0.0000000000000001;
					}
					else
					{
						param_p = 80;
						ep1     = 0.00000000001;
//						          0.0000000000001;
						ep2     = 0.0000000000001;
					}
				}
			}
		}
	}
*/

//	param_p = 20;
	p       = param_p;
//	p       = 2;
//	ccc     = param_ccc;
//	ccc     = 11;
//	param_p = 10;
//	ep1     = 0.000000000000001;
//	ep2     = 0.0000000000000001;
	ep1     = 0.0000000001;
	ep2     = 0.00000000001;
	resiff  = 0.000001;
	aux     = M->NumCols + NumBounds;
	auxx    = M->NumRows + NumBounds;
	aux1    = aux + M->NumRows;
	aux2    = aux1 + M->NumRows;
	aux3    = aux2 + NumBounds;
	aux4    = aux3 + M->NumCols;
	aux5    = p / 2;
	aux6    = p + 1;
	aux7    = aux1 + 1;
	pp2     = (aux5 + 1) * aux6;
	onem6   = M->NumRows + 2;

//	printf(" p = %d \n", p);
//	printf(" ccc = %d \n", ccc);

	ww         = NewInt(aux1 ,"ww");
	//wwa        = NewInt(aux1,"wwa");
	Acol       = NewInt(M->NumCols + 1, "Acol");
	r1r2       = NewInt(p,"r1r2");
	duposi     = NewInt(aux6,"duposi");      
	dl1posi    = NewInt(aux6,"dl1posi");     
	normau     = NewDouble(NumBounds,  "normau");
	fornang    = NewDouble(aux1,"fornang1");
	somalin    = NewDouble(M->NumRows,"somalin");
	B          = NewDouble(M->NumCols,"B");
	B2         = NewDouble(M->NumRows,"B2");
	B1         = NewDouble(onem6,"B1");
	xgrande    = NewDouble(aux4+1,"xgrande");
	uu1        = NewDouble(M->NumRows,"uu1");
	uu2        = NewDouble(M->NumCols,"uu2");
	val1       = NewDouble(M->Nonzeros,"val1");
	val2       = NewDouble(M->Nonzeros,"val2");
	w2         = NewDouble(M->NumRows,"w2");
	valgaux    = NewDouble(pp2,"valgaux");
	valg       = NewDouble(pp2,"valg");
	valg1      = NewDouble(pp2,"valg1");
	dumzero    = NewDouble(aux6,"dumzero");  
	dl1mzero   = NewDouble(aux6,"dl1mzero"); 
	N1         = NewDouble(aux6,"N1");
	du         = NewDouble(aux6,"du");
	dl1        = NewDouble(aux6,"dl1");
	ssf2       = NewDouble(aux6,"ssf2");
	sti        = NewDouble(aux6,"sti");
	ssf        = NewDouble(aux6,"ssf");
	wwaux      = NewDouble(aux6,"wwaux");
	dd12       = NewDouble(aux6,"dd12");
	//l1         = NewDouble(aux6,"l1");
	u          = NewDouble(aux6,"u");
	rr1        = NewDouble(aux6,"rr1");
	rr2        = NewDouble(aux6,"rr2");
	rr4        = NewDouble(aux6,"rr4");
	wb         = NewDouble(aux7,"wb");
	Paux1      = NewDouble(aux7,"Paux1");
	Paux2      = NewDouble(aux7,"Paux2");
	resif      = NewDouble(aux7,"resif");

	resi1      = resif;
	resiu      = resif + M->NumRows;
	resi2      = resif + auxx;
	resiv      = resif + aux1;	
	ONE        = 1;

	//vduration             = 0.0;

/*	for(i = 0; i < M->NumCols; i++)
	{
		xgrande[i] = px[i];
	}

	for(i = M->NumCols; i < aux; i++)
	{
		xgrande[i] = pw[i - M->NumCols];
	}

	for(i = aux2; i < aux3; i++)
	{
		xgrande[i] = pr[i - aux2];
	}

	for(i = aux3; i < aux4; i++)
	{
		xgrande[i] = ps[i - aux3];
	}
*/

	for(i = 0; i < (aux4+1); i++)
	{
		xgrande[i] = (double)1/(aux4+1);
	}


	x  = xgrande;
	w  = xgrande + M->NumCols;
	y1 = xgrande + aux;
	y2 = xgrande + aux1;
	r  = xgrande + aux2;
	s  = xgrande + aux3;
        vv1 = xgrande + aux4;
	for(i = k = 0; k < M->NumCols; k++)
	{
		if(M->pEndRow[k] > i)
			i += M->pEndRow[k] - i;

		Acol[k + 1] = i;
	}

	/*****************************************************************************/
	/*                     Fazendo a primeira transformação                      */
	/*****************************************************************************/

//	vv1 = 0;

/*	for(i = 0; i < M->NumCols; i++)
	{
		vv1 += x[i] + s[i];  
	}

	for(i = 0; i < M->NumRows; i++)
	{
		if(pi[i] >= 0)
		{
			y1[i] = pi[i];
		}
		else
		{
			y2[i] = -pi[i];
		}

		vv1 += y1[i] + y2[i]; 
	}

	for (i = 0; i < NumBounds; i++)
	{
		vv1 += w[i] + r[i];
	}

	vv1 = 1 / (vv1 + 1);	// constante da primeira transformação

	for(i = 0; i < NumBounds; i++)
	{
		w[i] *= vv1;	// fazendo a primeira transformação em w
		r[i] *= vv1;	// fazendo a primeira transformação em r
	}					// final da primeira transformação

*/	
	c1 = 0;
   	k  = 0;
	
	for(j = 0; j < M->NumCols; j++)         // achando a norma das colunas de  A | c
	{
		tam = 0;		
//		x[j] *= vv1;			// fazendo a primeira transformacao em x
//		s[j] *= vv1;			// fazendo a primeira transformacao em s = z

		c1 += c[j] * c[j];

		B1[tam++] = c[j];

		for(i = Acol[j]; i < Acol[j + 1]; i++)
		{
			B1[tam++] = M->Value[i];
		}

		if(NumBounds != 0)
		{
			for(i = 0; i < NumBounds; i++)
			{
				if(j == (BoundIndex[i] - 1))
				{
					B1[tam++] = 1;
				}
			}
		}

		B[j] = dnrm2(&tam, B1, &ONE);
	}

	
	for (i=0; i< M->NumRows; i++) /*inicializando o vetor somalin*/
	{
		somalin[i]=0;
	}

	for(i = 0; i < M->Nonzeros; i++)	// somando as linhas de A
	{
		somalin[M->Row[i] - 1] += M->Value[i] * M->Value[i];
	}

	c11 = 0;

	for(i = 0; i < M->NumRows; i++)	// a norma das colunas de A'|b'
	{
	//	y1[i] *= vv1;
	//	y2[i] *= vv1;						// fazendo a primeira transformação em y2=y-, vtil=vv1
		c11    = b[i] * b[i];
		B2[i]  = sqrt(somalin[i] + c11);	// B2 é um vetor com as normas das colunas de A'|b'
		c1    += c11;
	}

	for(i = 0; i < NumBounds; i++)
	{
		normau[i]  = upbound[i] * upbound[i];
		c1        += normau[i];
		normau[i]  = sqrt(normau[i] + 1);
	//	r[i]      *= normau[i];					
	}

	normabc  = sqrt(c1);	// norma de bc
	//vv1     *= normabc;		// fazendo transformação na última variavel v
	//somat    = vv1;

//	for(i = 0; i < M->NumCols; i++)	// fazendo a segunda transformação em x
//	{
//		x[i]  *= B[i];
//		somat += x[i] + s[i];	// somando os vetores x e s=z
//	}

//	for(i = 0; i < M->NumRows; i++)	// fazendo a segunda transformação em y1
//	{
//		y1[i] *= B2[i];
//		y2[i] *= B2[i];
//		somat += y1[i] + y2[i];
//	}

//	for(i = 0; i < NumBounds; i++)	// somando os vetores w e r
//	{
//		somat += w[i] + r[i];
//
//	}

//	vv1 /= somat;

	for(j = 0; j < M->NumCols; j++)	// criando um novo vetor de valor para Atil
	{
//		x[j] /= somat;
//		s[j] /= somat;

		for(i = Acol[j]; i < Acol[j+1]; i++)	
		{
			val1[i] = M->Value[i] / B[j];
		}
	}

	for(i = 0; i < M->Nonzeros; i++)	// criando um novo vetor de valor para Atil1
	{
		val2[i] = M->Value[i] / B2[M->Row[i] - 1];
		//val2[i] = M->Value[i] / B2[M->Row[i]];
	}

	/*****************************************************************************/
	/* CALCULANDO O RESIDUO                                                      */
	/*****************************************************************************/
	//vv1 = (double)1/(aux4+1);
	c11    = -*vv1/normabc;
	*resiv = 0;

	for(i = 0; i < M->NumRows; i++)
	{
//		y1[i]  /= somat;
//		y2[i]  /= somat;
		pi[i]   = y1[i] - y2[i];
		c1      = b[i]/B2[i];
		*resiv += (y2[i] - y1[i]) * c1;	// calculando o produto b'*y+=b'*y1
	}

	both(M, Acol, val1, val2, x, uu1, pi, uu2);

	for(i = 0; i < M->NumRows; i++)
	{
		resi1[i]  = uu1[i] + c11 * b[i];	// calculando o primeiro residuo
	}

	for(i = 0; i < M->NumCols; i++)	// calculando o produto vv1*c
	{
		resi2[i] = uu2[i] + s[i] + c11 * c[i];	// calculando o residou 2
		*resiv   += (c[i] / B[i]) * x[i];		// calculando o produto c'*x
	}

	for(i = 0; i < NumBounds; i++)
	{
//		w[i]                     /= somat;
//		r[i]                     /= somat;
		resiu[i]                  = c11 * upbound[i] + w[i] + (x[BoundIndex[i] - 1] / B[BoundIndex[i] - 1]);
		resi2[BoundIndex[i] - 1] -= r[i] / normau[i];
		*resiv                    += (upbound[i] / normau[i]) * r[i];	// ultimo residuo, residuo da variavel v
	}

	/*****************************************************************************/
	/* RESOLVENDO O PROBLEMA                                                     */
	/*****************************************************************************/

	cc     = 0;
	zz     = 1;
	resif1 = 5;
	resiff = dnrm2(&aux1, resif, &ONE);
        resiff0=resiff; 
	//printf("resi0  %30.30lf \n ", resiff);

	while(zz == 1)
	{
		cc++;
                 //  dif = fabs((resif1 - resiff)/resif1);
	

	//	printf( "resiff:      %2.8f \n", resiff);


		/*****************************************************************************/
		/* Criterio de parada                                                        */
		/*****************************************************************************/



		//if(((double)(clock() - start)/CLOCKS_PER_SEC) > param_t)
		 //if(dif<0.005)	
		 if(cc > p)
		{
			end = clock();

			printf( "tempo cc resiffim  %30.30g %d %11.10lf\n",(double)(end - start) / CLOCKS_PER_SEC,cc-1,resiff);
			//printf("resiff fim  %11.10lf \n ", resiff);

			break;
		}


			
/********************************************************************************************
			somat = vv1 / normabc;			

			for(i = 0; i < M->NumCols; i++)
			{
				somat += x[i] / B[i] + s[i];
			}

			for(i = 0; i < M->NumRows; i++)
			{
				somat += (y1[i] + y2[i])/ B2[i];
			}

			for(i = 0; i < NumBounds; i++)
			{
				somat += w[i] + r[i] / normau[i];
			}

			vv1 /= normabc * somat;
			c11  = somat * vv1;

			for(i = 0; i < M->NumCols; i++)
			{
				x[i] /= B[i] * c11;
				s[i] /= c11;

				if(x[i] == 0)
				{
					x[i] = 0.001;
				}

				if(s[i] == 0)
				{
					s[i] = 0.001;
				}


			}

			for(i = 0; i < NumBounds; i++)
			{
				w[i] /= c11;
				r[i] /= normau[i] * c11;

				if(w[i] == 0)
				{
					w[i] = 0.001;
				}

				if(r[i] == 0)
				{
					r[i] = 0.001;	
				}
			}

			for(i = 0; i < M->NumRows; i++)
			{
				c1     = B2[i] * c11;
				y1[i] /= c1;
				y2[i] /= c1;
				pi[i]  = y1[i] - y2[i];
                        }


			both(M, Acol, val1, val2, x, uu1, pi, uu2);

   			for( i=0; i<M->NumRows; i++)
			{
			  uu1[i]=uu1[i]-b[i];

			}

   			for( i=0; i<M->NumCols; i++)
			{
			  uu2[i]=uu2[i]+s[i]-b[i];

			}



			break;
		}
*/
		//xgrande[aux4] = vv1;
		both(M, Acol, val2, val1, resi2, w2, resi1, uu2);

		for(i = 0; i < aux1; i++)
		{
			ww[i]  = i; 
			//wwa[i] = 0;
		}

		if((cc % 2) == 0)
		{
			for(i = 0; i < aux1; i++)
			{
				fornang[i] = 0; 				
			}
			
			for(i = 0; i < M->NumCols; i++)
			{
				fornang[i] = uu2[i] + (c[i] / B[i]) * *resiv;	// vetor que fornece angulos P'b=alg
			}

			for(i = 0; i < NumBounds; i++)
			{
				fornang[BoundIndex[i] - 1] += resiu[i];
			}

			for(i = M->NumCols; i < aux; i++)
			{
				fornang[i] = resiu[i - M->NumCols];
			}

			for(i = aux; i < aux1; i++)
			{
				fornang[i] = w2[i - aux] - (b[i - aux] / B2[i - aux]) * *resiv;
			}
			
			wwa=0;
			sortd(fornang, ww, aux1);	// ordena em ordem descendente, guarda a posição ori em ww		
   
                        wwa   = ww[aux1-1];
			auxvn = fornang[aux1-1];
			//printf("jairpar wwa cc auxvn %d %d %g %g \n ", wwa, cc,auxvn,fornang[aux1-2]);
   
			r1r2[0] = wwa;

		}
		else
		{
			for(i = 0; i < aux1; i++)
			{
				fornang[i] = 0; 				
			}

			for(i = 0; i < M->NumRows; i++)
			{
				fornang[i] = -w2[i] + (b[i] / B2[i]) * *resiv;
			}

			// referente as variaveis canalizadas

			for(i = M->NumRows; i < auxx; i++)
			{
				fornang[i] = resi2[BoundIndex[i - M->NumRows] - 1] / (normau[i - M->NumRows]) + (*resiv * upbound[i - M->NumRows]) / (normau[i - M->NumRows]);
			}

			// referente a parte da Id

			for(i = auxx; i < aux1; i++)
			{
				fornang[i] = resi2[i - auxx];
			}

			wwa=0;	
			sortd(fornang, ww, aux1); 			// ordena em ordem descendente, guarda a posição ori em ww
			xzero = 0;

                        wwa   = ww[aux1-1];
			auxvn = fornang[aux1-1];
			//	printf("jairpar wwa cc auxvn %d %d %g %g \n ", wwa, cc,auxvn,fornang[aux1-2]);
   
			r1r2[0] = wwa+aux1;
	
		}

		
	//	resif1 = resiff;
	//	resiff = dnrm2(&aux1, resif, &ONE);

		//printf(" resiff  %30.30lf \n", resiff);
		/*****************************************************************************/
		/* SUBPROBLEMA PARA ACHAR A MENOR NORMA DE b E ASSIM ATUALIZAR X             */
		/*****************************************************************************/


                resiff2 = 0;
                resiff2 = resiff * resiff;
                l1      = 0 ;
                l1      = (1 - auxvn) / ((resiff2 - 2*auxvn) + 1);
                auxvn   = 0;

		for(i = 0; i < aux4+1; i++)
		{
			xgrande[i] *= l1;
			auxvn      += xgrande[i];
		}

                auxvn = auxvn + (1 - l1);

		for(i = 0; i < 1; i++)
		{
			xgrande[r1r2[i]] = xgrande[r1r2[i]]+(1-l1);
		}

	        //auxvn = 0;

		//for(i = 0; i < aux4+1; i++)
		//{
		//	auxvn += xgrande[i];
		//}






			/*****************************************************************************/
			/*           	CALCULANDO O RESIDUO EM CADA ITERAÇÃO                        */
			/*****************************************************************************/
			
			//vv1 = xgrande[aux4];
			c11 = -*vv1/normabc;

			for(i = 0; i < M->NumRows; i++)
			{
				pi[i]    = y1[i] - y2[i];
			}

			both(M, Acol, val1, val2, x, uu1, pi, uu2);

			*resiv = 0;

			for(i = 0; i < M->NumCols; i++)
			{
				resi2[i] = uu2[i] + s[i] + c11 * c[i];			/* calculando o residou 2*/
				*resiv   += (c[i] / B[i]) * x[i];			/*calculando o produto c'*x */
			}

			for(i = 0; i < M->NumRows; i++)
			{
				c1     = b[i] / B2[i];
				*resiv += (y2[i] - y1[i]) * c1;				/*calculando o produto b'*y+=b'*y1 */
				resi1[i] = uu1[i] + c11 * b[i]; 				/*calculando o primeiro residuo*/
			}

			for(i = 0; i < NumBounds; i++)
			{
				resiu[i]                  = c11 * upbound[i] + w[i] + (x[BoundIndex[i] - 1] / B[BoundIndex[i] - 1]);
				resi2[BoundIndex[i] - 1] -= r[i] / normau[i];
				*resiv                    += (upbound[i] / normau[i]) * r[i];
			}
		resif1 = resiff;
		resiff = dnrm2(&aux1, resif, &ONE);
		//printf(" resiff  %30.30lf \n", resiff);


		}//fim do loço while

	
	//printf( "cc =  %d \n", cc-1);
	
	Free((char *)resif);
	Free((char *)Paux2);
	Free((char *)Paux1);

	Free((char *)wb);
	Free((char *)rr4);
	Free((char *)rr2);
	Free((char *)rr1);
	Free((char *)u);
	//Free((char *)l1);
	Free((char *)dd12);
	Free((char *)wwaux);
	Free((char *)ssf);
	Free((char *)sti);
	Free((char *)ssf2);
	Free((char *)dl1);
	Free((char *)du);
	Free((char *)N1);
	Free((char *)dl1mzero);
	Free((char *)dumzero);
	Free((char *)valg1);

	Free((char *)valg);
	Free((char *)valgaux);
	Free((char *)w2);
	Free((char *)val2);
	Free((char *)val1);
	Free((char *)uu2);
	Free((char *)uu1);
	Free((char *)xgrande);
	Free((char *)B1);
	Free((char *)B2);
	Free((char *)B);
	Free((char *)somalin);
	Free((char *)fornang);
	Free((char *)normau);
	Free((char *)dl1posi);
	Free((char *)duposi);
	Free((char *)r1r2);
	Free((char *)Acol);
	//Free((char *)wwa);
	Free((char *)ww);
}

void both(MMTtype *M, int *Acol, double *v1, double *v2, double *x, double *v, double *y, double *w)
{
  int i, j, k;

     memset(v, 0, M->NumRows*sizeof(double));

     for(i = Acol[j= 0]; j < M->NumCols; j++)
     {
         k = Acol[j+1];

         for(w[j] = 0.0; i < k; i++)
         {
             w[j]         += v2[i]*y[M->Row[i] - 1];
             v[M->Row[i] - 1] += v1[i]*x[j];
         }
     }
}


void cholesky(int aux6, int pp2, int p, double s1, double v, double valgaux[], double valg[], double valg1[], double wwaux[], double l1[], double rr1[], double rr2[], double rr4[], double *rr3, double u[], double dd12[])
{
    int    i, j, k, p_i, pj, pk, pw;
    double prodl1u, prodl1, gaux;

    for(i = 0; i < pp2; i++)
    {
        valgaux[i] = valg[i];
    }

    p_i = 0;
    i  = 0;

    while(i < aux6)
    {
        wwaux[i]  = 0;
        wwaux[0] += valg[p_i] * l1[i];

        i++;
        p_i += i;
    }

    pj = 1;
    j  = 1;

    while (j < aux6)
    {
        pk = pj;

        for(k = 0; k < j; k++)
        {
            wwaux[j] += valg[pk] * l1[k];
            pk++;
        }

        p_i = (j * (j + 1)) / 2 + j;
        i  = j;

        while(i < aux6)
        {
            wwaux[j] += valg[p_i] * l1[i];

            i++;
            p_i += i;
        }

        j++;
	pj += j;
    }

	rr1[0] = u[0] - s1 * v - wwaux[0]; // calculo do primeiro residuo

	prodl1u = l1[0] * u[0];

    	for(i = 1; i < aux6; i++)
    	{
        	prodl1u += l1[i] * u[i];
        	rr1[i]   = u[i] - v - wwaux[i]; // calculo do primeiro residuo
	}

	prodl1u /= (aux6 * aux6);
	prodl1   = 0;

	rr2[0] = -l1[0] * u[0] + prodl1u; // calculo do segundo residuo

	for(i = 1; i < aux6; i++)
	{
		rr2[i]  = -l1[i] * u[i] + prodl1u;    // calculo do segundo residuo
        	prodl1 += l1[i];
    }

    *rr3 = 1 - l1[0] * s1 - prodl1;

    for(i = 0; i < p + 1; i++)
    {
        rr4[i] = rr1[i] + (1 / l1[i]) * rr2[i];
    }

    // Somando U*L1^(-1) na diagonal                             

    p_i = 0;
    i  = 0;

    while (i < aux6)
    {
        valgaux[p_i + i] += u[i] / l1[i];
        dd12[i] = 1 / sqrt(valgaux[p_i + i]);

        i++;
        p_i += i;
    }
   
    p_i = 0;
    i  = 0;
   
    while (i < aux6)
    {
        pj = p_i;
        k  = i + 1;
       
        for(j = 0; j < k; j++)
        {
            valgaux[pj] *= dd12[i] * dd12[j];
            pj++;
        }
       
        i++;
        p_i += i;
    }

    valg1[0] = sqrt(valgaux[0]);

    p_i = 1;
    i  = 1;

    while (i < aux6)
    {
        valg1[p_i] = valgaux[p_i] / valg1[0];

        i++;
        p_i += i;
    }

    for(j = 1; j < aux6; j++)
    {
        gaux = 0;
        pj   = j * (j + 1) / 2;

        for(k = 0; k < j; k++)
        {
            gaux += valg1[pj] * valg1[pj];
            pj++;
        }

        pj = j * (j + 1) / 2 + j;

        valg1[pj] = sqrt(valgaux[pj] - gaux);

        for(i = j + 1; i < aux6; i++)
        {
        	gaux = 0;
            p_i   = i * (i + 1) / 2;
            pj   = j * (j + 1) / 2;
            pk   = p_i;
            pw   = pj;

            for(k = 0; k < j; k++)
            {
		gaux += valg1[pk++] * valg1[pw++];
            }

            p_i += j;
            valg1[p_i] = (valgaux[p_i] - gaux) / valg1[pj + j];
        }
    }
}

void resolucao_sistema(int aux6, int p, double s1, double sti[], double ssf[], double ssf2[], double dd12[], double valg1[], double rr4[])
{
	int    i, k, p_i, pk, pw;
	double saux;

	sti[0] = (s1 * dd12[0]) / valg1[0];

	p_i = 1;
	i  = 1;

	while (i < aux6)
	{
		saux = 0;
		pk   = p_i;

		for(k = 0; k < i; k++)
		{
			saux += valg1[pk] * sti[k];
			pk++;
		}

		sti[i] = (dd12[i] - saux) / valg1[p_i + i];

		i++;
		p_i += i;
	}

	ssf[p] = sti[p] / valg1[(p * aux6) / 2 + p];
	
	p_i = p;

	for(i = 1; i < aux6; i++) // resolucao do primeiro s. linear
	{
		saux = 0;
		pk   = p;
		pw   = pk * (pk + 1) / 2 + --p_i;

		for(k = 0; k < i; k++)
		{
			saux += valg1[pw] * ssf[pk];
			
			pw -= pk; 
			pk--;
		}

		ssf[p_i] = (sti[p_i] - saux) / valg1[(p_i * (p_i + 1)) / 2 + p_i];
	}

	for(i = 0; i < aux6; i++) // solucao final do primeiro sistema
	{
		ssf[i] = dd12[i] * ssf[i]; 				
	}

	sti[0] = (rr4[0] * dd12[0]) / valg1[0];

	p_i = 1;
	i  = 1;

	while (i < aux6)
	{
		saux = 0;
		pk   = p_i;

		for(k = 0; k < i; k++)
		{
			saux += valg1[pk] * sti[k];
			pk++;			
		}

		sti[i] = (rr4[i] * dd12[i] - saux) / valg1[p_i + i];

		i++;
		p_i += i;
	}

	ssf2[p] = sti[p] / valg1[(p * aux6) / 2 + p];

	p_i = p;

	for(i = 1; i < aux6; i++) // resolucao do segundo s. linear
	{
		saux = 0;
		pk   = p;
		pw   = pk * (pk + 1) / 2 + --p_i;

		for(k = 0; k < i; k++)
		{
			saux += valg1[pw] * ssf2[pk];

			pw -= pk;
			pk--;
		}

		ssf2[p_i] = (sti[p_i] - saux) / valg1[(p_i * (p_i + 1)) / 2 + p_i];
	}

	for(i = 0; i < aux6; i++)
	{
		ssf2[i] = dd12[i] * ssf2[i]; // solucao final do segundo sistema
	}
}

