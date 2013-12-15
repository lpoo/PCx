
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
	int    i, aux, auxx, auxx02, aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, ccc;
	int    j, k, zz, cc, p, mxzero, mxzero1, xzero, l, cci, pp2, p_i, pj, pk, onem6, ONE, tam, wwa, wwa02;
	int    *dl1posi, *duposi, *r1r2, *ww, *ww02, *Acol;
	double *x, *s, *w, *r;
	double vv1, c1, c11, s1, somat,  v, alpu, alpl1, auxl1, resiff, resif1, dif, rr3, normabc, d__1,absxi, scale, resiff2, auxvn,l1;
	double *B1, *y1, *B2, *B,*somalin,*xgrande,*wb, *Paux1, *Paux2, *rr1, *rr2, *rr4, *u, *alpp;
	double *val1,*val2,*resi1,*resi2,*resiv,*resiu,*w2,*uu1,*uu2, *valg, *valg1, *wwaux, *dd12,*valgaux, *resif, *normau;
	double *sti, *sts, *ssf, *du, *dl1, *N1, *sti2, *sts2, *ssf2, *dl1mzero, *dumzero, *fornang, *fornang02;
	double	vduration;
	clock_t start, end;

	start  = clock();
	p      = param_p;
	resiff = 0.000001;
        auxx02 = M->NumCols+NumBounds;
	aux    = M->NumCols + NumBounds;
	auxx   = M->NumRows + NumBounds;
	aux1   = aux + M->NumRows;
	aux2   = aux1 + M->NumRows;
	aux3   = aux2 + NumBounds;
	aux4   = aux3 + M->NumCols;
	aux5   = p / 2;
	aux6   = p + 1;
	aux7   = aux1 + 1;
        aux9   = aux1+ NumBounds;
	aux8   = aux9+M->NumCols+1;
	pp2    = (aux5 + 1) * aux6;
	onem6  = M->NumRows + 2;

	ww         = NewInt(auxx02 ,"ww");
        ww02       = NewInt(aux1 ,"ww02");
	Acol       = NewInt(M->NumCols + 1, "Acol");
	r1r2       = NewInt(p,"r1r2");
	duposi     = NewInt(aux6,"duposi");      
	dl1posi    = NewInt(aux6,"dl1posi");     
	normau     = NewDouble(NumBounds,  "normau");
	fornang    = NewDouble(auxx02,"fornang");
        fornang02  = NewDouble(aux1,"fornang02");
	somalin    = NewDouble(M->NumRows,"somalin");
	B          = NewDouble(M->NumCols,"B");
	B2         = NewDouble(M->NumRows,"B2");
	B1         = NewDouble(onem6,"B1");
	xgrande    = NewDouble(aux8,"xgrande");
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

	vduration  = 0.0;
        auxvn      = 0;
 
	for(i = 0; i < aux8; i++)
	{
	   xgrande[i] = 10 * (0.1 / aux8);
           auxvn      = auxvn + xgrande[i];
 	}

	x  = xgrande;
	w  = xgrande + M->NumCols;
	y1 = xgrande + aux;
	r  = xgrande + aux1;
	s  = xgrande + aux9;

	for(i = k = 0; k < M->NumCols; k++)
	{
		if(M->pEndRow[k] > i)
			i += M->pEndRow[k] - i;

		Acol[k + 1] = i;
	}

	/*****************************************************************************/
	/*                     Fazendo a primeira transformação                      */
	/*****************************************************************************/

	vv1 = 0;
	c1  = 0;
   	k   = 0;

	for(j = 0; j < M->NumCols; j++)         // achando a norma das colunas de  A | c
	{
		tam = 0;		
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

	for(i = 0; i < M->Nonzeros; i++)	// somando as linhas de A
	{
		somalin[M->Row[i] - 1] += M->Value[i] * M->Value[i];
	}

	c11 = 0;


	for(i = 0; i < M->NumRows; i++)	// a norma das colunas de A'|b'
	{
		c11    = b[i] * b[i];
        	B2[i]  = sqrt(somalin[i] + c11);	// B2 é um vetor com as normas das colunas de A'|b'
		c1    += c11;
	}

	for(i = 0; i < NumBounds; i++)
	{
		normau[i]  = upbound[i] * upbound[i];
		c1        += normau[i];
		normau[i]  = sqrt(normau[i] + 1);
	}

	normabc  = sqrt(c1);	// norma de bc

	for(j = 0; j < M->NumCols; j++)	// criando um novo vetor de valor para Atil
	{
		for(i = Acol[j]; i < Acol[j+1]; i++)	
		{
			val1[i] = M->Value[i] / B[j];
		}
	}

	for(i = 0; i < M->Nonzeros; i++)	// criando um novo vetor de valor para Atil1
	{
		val2[i] = M->Value[i] / B2[M->Row[i] - 1];
	}

	/*****************************************************************************/
	/* CALCULANDO O RESIDUO                                                      */
	/*****************************************************************************/

        vv1 = 10*(0.1/aux8);
        c11 = -vv1/normabc;

	*resiv = 0;

	for(i = 0; i < M->NumRows; i++)
	{
                  pi[i]  = y1[i];
		  c1     = b[i]/B2[i];
                 *resiv += (- y1[i]) * c1;
	}

	both(M, Acol, val1, val2, x, uu1, pi, uu2);

	for(i = 0; i < M->NumRows; i++)
	{
		resi1[i]  = uu1[i] + c11 * b[i];	// calculando o primeiro residuo
	}

	for(i = 0; i < M->NumCols; i++)			// calculando o produto vv1*c
	{
		resi2[i] = uu2[i] + s[i] + c11 * c[i];	// calculando o residou 2
		*resiv   += (c[i] / B[i]) * x[i];	// calculando o produto c'*x
	}

	for(i = 0; i < NumBounds; i++)
	{
		resiu[i]                  = c11 * upbound[i] + w[i] + (x[BoundIndex[i] - 1] / B[BoundIndex[i] - 1]);
		resi2[BoundIndex[i] - 1] -= r[i] / normau[i];
		*resiv                    += (upbound[i] / normau[i]) * r[i];	// ultimo residuo, residuo da variavel v
	}

	//printf("B \n");
	//for (i = 0; i < M->NumCols; i++)	
	//{
	//	printf("%.4lf ", B[i]);
	//}
	//printf("\n");
	//getchar();

	/*****************************************************************************/
	/* RESOLVENDO O PROBLEMA                                                     */
	/*****************************************************************************/

	cc     = 0;
	zz     = 1;
	resif1 = 1;

	while(zz == 1)
	{
		cc++;

		dif = fabs((resif1 - resiff)/resif1);

		/*****************************************************************************/
		/* Criterio de parada                                                        */
		/*****************************************************************************/

		// if(cc > 10)
		if(((double)(clock() - start)/CLOCKS_PER_SEC) > param_t)
               	{
			end = clock();
			//printf( "tempo: %6.4g seconds\n",(double)(end - start) / CLOCKS_PER_SEC);
			break;
		}

		xgrande[aux8-1] = vv1;


                 both(M, Acol, val2, val1, resi2, w2, resi1, uu2);
		//both(M, Acol, val1, val2, resi2, w2, resi1, uu2);(troquei aki, porque temos (val1)^t e val2 não )

		if((cc % 2) == 0)
                {
			auxvn=0;

			for(i = 0; i < auxx02; i++)
			{
				ww[i]  = i;
			}


			for(i = 0; i < auxx02; i++)
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

			wwa   = 0;
			sortd(fornang, ww, auxx02);	// ordena em ordem descendente, guarda a posição ori em ww
                        wwa   = ww[auxx02-1];
			auxvn = fornang[auxx02-1];
			//	printf("jairpar wwa cc auxvn %d %d %g %g \n ", wwa[0], cc,auxvn,fornang[auxx02-2]);
   
			r1r2[0] = wwa;

		}
		else
		{

		for(i = 0; i < aux1 ; i++)
		{
			ww02[i]  = i;
		}


		for(i = 0; i < aux1; i++)
		{
			fornang02[i] =0;
		}


		for(i = 0; i < M->NumRows; i++)
		{
			fornang02[i] = w2[i] - (b[i] / B2[i]) * *resiv;
		}

		// referente as variaveis canalizadas

		for(i =  M->NumRows; i < auxx; i++)
		{
			fornang02[i] = resi2[BoundIndex[i- M->NumRows] - 1] / (normau[i- M->NumRows]) + (*resiv * upbound[i- M->NumRows]) / (normau[i- M->NumRows]);
		}

		// referente a parte da Id

		for(i =auxx ; i < aux1; i++)
		{
			fornang02[i] = resi2[i -auxx];
		}

		wwa02 = 0;
		auxvn = 0;
		sortd(fornang02, ww02, aux1); 			// ordena em ordem descendente, guarda a posição ori em ww
		xzero = 0;
                wwa02 = ww02[aux1-1]+aux;
		auxvn = fornang02[aux1-1];
		//	printf("jairimpar wwa02 cc auxvn %d %d %g \n ", wwa02[0], cc,auxvn);

		r1r2[0] = wwa02;
		}

		resif1  = resiff;
		resiff  = dnrm2(&aux7, resif, &ONE);
                 printf("jair resif %g %d\n ", resiff, cc);

		/*****************************************************************************/
		/* SUBPROBLEMA PARA ACHAR A MENOR NORMA DE b E ASSIM ATUALIZAR X             */
		/*****************************************************************************/

                resiff2 = 0;
                resiff2 = resiff * resiff;
                l1      = 0 ;
                l1      = (1 - auxvn) / ((resiff2 - 2*auxvn) + 1);
                auxvn   = 0;

		for(i = 0; i < aux8; i++)
		{
			xgrande[i] *= l1;
			auxvn      += xgrande[i];
		}

                auxvn = auxvn + (1 - l1);

		for(i = 0; i < 1; i++)
		{
			xgrande[r1r2[i]] = xgrande[r1r2[i]]+(1-l1);
		}

	        auxvn = 0;

		for(i = 0; i < aux8; i++)
		{
			auxvn += xgrande[i];
		}

		/*****************************************************************************/
		/*           	CALCULANDO O RESIDUO EM CADA ITERAÇÃO                        */
		/*****************************************************************************/

		vv1 = xgrande[aux8-1];
		c11 = -vv1/normabc;

		for(i = 0; i < M->NumRows; i++)
		{
			pi[i] = y1[i];
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
			c1       = b[i] / B2[i];
                        *resiv  += (- y1[i]) * c1;
			resi1[i] = uu1[i] + c11 * b[i]; 				/*calculando o primeiro residuo*/
		}

		for(i = 0; i < NumBounds; i++)
		{
			resiu[i]                   = c11 * upbound[i] + w[i] + (x[BoundIndex[i] - 1] / B[BoundIndex[i] - 1]);
			resi2[BoundIndex[i] - 1]  -= r[i] / normau[i];
			*resiv                    += (upbound[i] / normau[i]) * r[i];
		}

	} //fim primeiro laço while

	for(i = 0; i < M->NumCols; i++)
	{
		px[i] = xgrande[i];
	}

	for(i = M->NumCols; i < aux; i++)
	{
		pw[i - M->NumCols] = xgrande[i];
	}

	for(i = aux2; i < aux3; i++)
	{
		pr[i - aux2] = xgrande[i];
	}

	for(i = aux3; i < aux4; i++)
	{
		ps[i - aux3] = xgrande[i];
	}

	printf( "cc =  %d \n", cc-1);

	Free((char *)resif);
	Free((char *)Paux2);
	Free((char *)Paux1);
	Free((char *)wb);
	Free((char *)rr4);
	Free((char *)rr2);
	Free((char *)rr1);
	Free((char *)u);
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
        Free((char *)fornang02);
	Free((char *)normau);
	Free((char *)dl1posi);
	Free((char *)duposi);
	Free((char *)r1r2);
	Free((char *)Acol);
	Free((char *)ww02);
	Free((char *)ww);
} //programa acaba aki

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

