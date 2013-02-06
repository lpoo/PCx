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

extern int param_p, param_ccc;

void jair1(MMTtype *M, FactorType *Factor, double b[], double c[], double *px, double *pi, double *ps, double *pw, double *pr, double upbound[], int NumBounds, int BoundIndex[])
{
	int    i, aux, auxx, aux1, aux2, aux3, aux4, aux5, aux6, aux7,ccc;
	int    j, k, zz, cc, p, mxzero, mxzero1, xzero, l, cci, pp2, p_i, pj, pk, onem6, ONE, tam;
	int    *wwa, *dl1posi, *duposi, *r1r2, *ww, *Acol, param_p;
	double *x, *s, *w, *r;
	double vv1, c1, c11, s1, somat,  v, alpu, alpl1, auxl1, ep1, ep2, resiff, resif1, dif, rr3, normabc,d__1,absxi,scale;
	double *B1, *y1, *y2, *B2, *B,*somalin,*xgrande,*wb, *Paux1, *Paux2, *rr1, *rr2, *rr4, *l1, *u, alpp;
	double *val1,*val2,*resi1,*resi2,*resiv,*resiu,*w2,*uu1,*uu2, *valg, *valg1, *wwaux, *dd12,*valgaux, *resif, *normau;
	double *sti, *sts, *ssf, *du, *dl1, *N1, *sti2, *sts2, *ssf2, *dl1mzero, *dumzero, *fornang;
	
	if((M->NumRows) < 1000)
	{
		param_p = 2;
//		ep1     = 0.0000000000000001;
//		ep2     = 0.000000000000000001;
	}
	else
	{
		if((M->NumRows) < 10000)
		{
			param_p = 4;
//			ep1     = 0.000000000000001;
//			ep2     = 0.00000000000000001;
		}
		else
		{	   
			if((M->NumRows) < 20000)
			{
				param_p = 8;
//				ep1     = 0.0000000000000001;
//				ep2     = 0.000000000000000001;
			}
			else
			{
				if((M->NumRows) < 100000)	
				{
					param_p = 20;
//					ep1     = 0.000000000001;
//					ep2     = 0.00000000000001;
				}
				else
				{
					if((M->NumRows) < 150000)	
					{
						param_p = 40;
//						ep1     = 0.00000000000001;
//						ep2     = 0.0000000000000001;
					}
					else
					{
						param_p = 80;
//						ep1     = 0.00000000001;
//						          0.0000000000001;
//						ep2     = 0.0000000000001;
					}
				}
			}
		}
	}


//	param_p = 20;
	p       = param_p;
//	p       = 2;
	ccc     = param_ccc;
//	ccc     = 11;
//	param_p = 10;
//	ep1     = 0.000000000000001;
//	ep2     = 0.0000000000000001;
	ep1     = 0.00000000001;
	ep2     = 0.000000000001;
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

	printf(" p = %d \n", p);
	//printf(" ccc = %d \n", ccc);

	ww         = NewInt(aux1 ,"ww");
	wwa        = NewInt(aux1,"wwa");
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
	l1         = NewDouble(aux6,"l1");
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

	for(i = 0; i < M->NumCols; i++)
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

	x  = xgrande;
	w  = xgrande + M->NumCols;
	y1 = xgrande + aux;
	y2 = xgrande + aux1;
	r  = xgrande + aux2;
	s  = xgrande + aux3;

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

	for(i = 0; i < M->NumCols; i++)
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

	
	c1 = 0;
   	k  = 0;
	
	for(j = 0; j < M->NumCols; j++)         // achando a norma das colunas de  A | c
	{
		tam = 0;		
		x[j] *= vv1;			// fazendo a primeira transformacao em x
		s[j] *= vv1;			// fazendo a primeira transformacao em s = z

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
		y1[i] *= vv1;
		y2[i] *= vv1;						// fazendo a primeira transformação em y2=y-, vtil=vv1
		c11    = b[i] * b[i];
		B2[i]  = sqrt(somalin[i] + c11);	// B2 é um vetor com as normas das colunas de A'|b'
		c1    += c11;
	}

	for(i = 0; i < NumBounds; i++)
	{
		normau[i]  = upbound[i] * upbound[i];
		c1        += normau[i];
		normau[i]  = sqrt(normau[i] + 1);
		r[i]      *= normau[i];					
	}

	normabc  = sqrt(c1);	// norma de bc
	vv1     *= normabc;		// fazendo transformação na última variavel v
	somat    = vv1;

	for(i = 0; i < M->NumCols; i++)	// fazendo a segunda transformação em x
	{
		x[i]  *= B[i];
		somat += x[i] + s[i];	// somando os vetores x e s=z
	}

	for(i = 0; i < M->NumRows; i++)	// fazendo a segunda transformação em y1
	{
		y1[i] *= B2[i];
		y2[i] *= B2[i];
		somat += y1[i] + y2[i];
	}

	for(i = 0; i < NumBounds; i++)	// somando os vetores w e r
	{
		somat += w[i] + r[i];

	}

	vv1 /= somat;

	for(j = 0; j < M->NumCols; j++)	// criando um novo vetor de valor para Atil
	{
		x[j] /= somat;
		s[j] /= somat;

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

	c11    = -vv1/normabc;
	*resiv = 0;

	for(i = 0; i < M->NumRows; i++)
	{
		y1[i]  /= somat;
		y2[i]  /= somat;
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
		w[i]                     /= somat;
		r[i]                     /= somat;
		resiu[i]                  = c11 * upbound[i] + w[i] + (x[BoundIndex[i] - 1] / B[BoundIndex[i] - 1]);
		resi2[BoundIndex[i] - 1] -= r[i] / normau[i];
		*resiv                    += (upbound[i] / normau[i]) * r[i];	// ultimo residuo, residuo da variavel v
	}

	/*****************************************************************************/
	/* RESOLVENDO O PROBLEMA                                                     */
	/*****************************************************************************/

	cc     = 0;
	zz     = 1;

	resiff = dnrm2(&aux1, resif, &ONE);
	resif1 = resiff;
//	printf("resi0  %30.30lf \n", resiff);
	dif = 1;
	while(zz == 1)
	{
		cc++;

	//	printf( "resif1:      %2.8f \n", resif1);
	//	printf( "resiff:      %2.8f \n", resiff);
	//	printf( "dif:      %20.20lf \n", dif);

		/*****************************************************************************/
		/* Criterio de parada                                                        */
		/*****************************************************************************/

		//if((cc > 50) || cc  == param_p+1)
		if((cc > 50) || ((cc > 2) &&  (dif < 0.001)))
//		if((cc > 50) || ((cc > 2) && ((-0.001 dif)&&(dif < 0.001))))
	  // 	if((cc > 200) || ((cc > 2) && (cc == ccc)))
		//if(cc == ccc)
		{

			
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

			break;
		}

		xgrande[aux4] = vv1;
		both(M, Acol, val2, val1, resi2, w2, resi1, uu2);

		for(i = 0; i < aux1; i++)
		{
			ww[i]  = i; 
			wwa[i] = 0;
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

			sortd(fornang, ww, aux1);	// ordena em ordem descendente, guarda a posição ori em ww		
   
			xzero = 0;

			for(i = 0; i < aux1; i++)
			{
				if(xgrande[ww[i]] == 0)
				{
					xzero++;
				}
				else
				{
					wwa[i - xzero] = ww[i];
				}
			}

			for(i = 0; i < aux5; i++)
			{
				r1r2[i]      = wwa[aux1 - xzero - i - 1];	// ind dos (p/2) vetores que formam o maior com residuo
				r1r2[i+aux5] = wwa[i];						// indices dos (p/2) vetores que formam o menor ang com resi
			}

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

			sortd(fornang, ww, aux1); 			// ordena em ordem descendente, guarda a posição ori em ww
			xzero = 0;

			for(i = 0; i < aux1; i++)
			{
				if(xgrande[ww[i] + aux1] == 0.0)
				{
					xzero++;
				}
				else
				{
					wwa[i - xzero] = ww[i];
				}
			}
		
			for(i = 0; i < aux5; i++)
			{
				r1r2[i]      = wwa[aux1 - xzero - i -1] + aux1;	// ind dos (p/2) vetores que formam o maior com residuo
				r1r2[i+aux5] = wwa[i] + aux1;						// indices dos (p/2) vetores que formam o menor ang com resi
			}
	
		}

		
		resif1 = resiff;
		resiff = dnrm2(&aux1, resif, &ONE);
		dif = fabs((resif1 - resiff)/resiff);
		//printf(" resiff  %30.30lf \n", resiff);

		/*****************************************************************************/
		/* SUBPROBLEMA PARA ACHAR A MENOR NORMA DE b E ASSIM ATUALIZAR X             */
		/*****************************************************************************/

		if(zz == 1)
		{
			s1 = 1;

			for(i = 0; i < aux5; i++)
			{
				s1 -= (xgrande[r1r2[i]] + xgrande[r1r2[i+aux5]]);
			}

			/*****************************************************************************/
 			/* Construindo o vetor w = b - sum(xiPi)                                     */
			/*****************************************************************************/

			for(i = 0; i < aux7; i++)
			{
				wb[i] = 0;
			}

			for(i = 0; i < p; i++)
			{
				if(r1r2[i] < M->NumCols)
				{
					wb[aux1] += xgrande[r1r2[i]] * (c[r1r2[i]] / B[r1r2[i]]);

					for(j = Acol[r1r2[i]]; j < Acol[r1r2[i] + 1]; j++)
					{
						wb[M->Row[j]-1] += val1[j] * xgrande[r1r2[i]];
					}

					for(k = 0; k < NumBounds; k++)
					{
						if(r1r2[i] == (BoundIndex[k] - 1))
						{
							wb[k + M->NumRows] += (xgrande[r1r2[i]] / (B[r1r2[i]]));
						}
					}
				}
				else
					if(r1r2[i] < aux)
					{
						wb[r1r2[i] - M->NumCols + M->NumRows] += xgrande[r1r2[i]];
					}
					else
						if(r1r2[i] < aux1)
						{
							wb[aux1] += xgrande[r1r2[i]] * (-b[r1r2[i] - aux] / B2[r1r2[i] - aux]);

							for(j = Acol[k = 0]; k < M->NumCols; k++)
							{
								for( ; j < Acol[k + 1]; j++)
								{
									if((r1r2[i] - aux) == M->Row[j] - 1)
									{
										wb[k + auxx] += val2[j] * xgrande[r1r2[i]];
									}
								}
							}
						}
						else
							if(r1r2[i] < aux2)
							{
								wb[aux1] += xgrande[r1r2[i]] * (b[r1r2[i] - aux1] / B2[r1r2[i] - aux1]);

								for(j = Acol[k = 0]; k < M->NumCols; k++)
								{
									for( ; j < Acol[k + 1]; j++)
									{
										if((r1r2[i] - aux1) == M->Row[j] - 1)
										{
											wb[k + auxx] += val2[j] * (-xgrande[r1r2[i]]);
										}
									}
								}
							}
							else
								if(r1r2[i] < aux3)
								{
									wb[aux1] += (xgrande[r1r2[i]] * upbound[r1r2[i] - aux2]) / (normau[r1r2[i] - aux2]);
									wb[auxx + r1r2[i] - aux2] += xgrande[r1r2[i]] / (normau[r1r2[i] - aux2]);
								}
								else
									if(r1r2[i] < aux4)
									{
										wb[auxx + r1r2[i] - aux3] += xgrande[r1r2[i]];
									}
									else
									{
										for(j = 0; j < M->NumRows; j++)
										{
											wb[j] += xgrande[r1r2[i]] * (-b[j] / normabc);
										}

										for(j = M->NumRows; j < auxx; j++)
										{
											wb[j] += xgrande[r1r2[i]] * (-upbound[j - M->NumRows] / normabc);
										}

										for(j = auxx; j < aux1; j++)
										{
											wb[j] += xgrande[r1r2[i]] * (-c[j - auxx] / normabc);
										}
									}
			}

			for(i = 0; i < aux1; i++)
			{
				if(wb[i] == 0.0)
				{
					wb[i] = resif[i];	// construindo o vetor w=b-xiPi
				}
				else
				{
					wb[i] = resif[i] - wb[i];
				}
			}

			if(wb[M->NumRows + M->NumCols] == 0.0)
			{
				wb[aux1] = *resiv;	// construindo o vetor w=b-xiPi
			}
			else
			{
				wb[aux1] = *resiv - wb[aux1];
			}

			/*****************************************************************************/
			/* CONSTRUINDO A MATRIZ G                                                    */
			/*****************************************************************************/

			for(i = 0; i < pp2; i++)
			{
				valg[i] = 0;
			}

			for(i = 0; i < aux7; i++)
			{
				valg[0] += wb[i] * wb[i];	// primeiro elemento da matriz g
			}

			for(l = 0; l < p; l++)	// colocar p
			{
				for(j = 0; j < aux7; j++)
				{
					Paux1[j] = 0;  								
				}

				for(i = 0; i < l + 1; i++)
				{
					for(j = 0; j < aux7; j++)
					{
						Paux2[j] = 0;   
					}

					if(r1r2[l-i] < M->NumCols)
					{
						Paux2[aux1] = c[r1r2[l-i]] / B[r1r2[l-i]];

						for(j = Acol[r1r2[l-i]]; j < Acol[r1r2[l-i] + 1]; j++)
						{
							Paux2[M->Row[j] - 1] = val1[j];
						}

						for(k = 0; k < NumBounds; k++)
						{
							if(r1r2[l-i] == (BoundIndex[k]-1))
							{
								Paux2[k+M->NumRows] = (1 / (B[r1r2[l-i]]));
							}
						}
					}
					else
						if(r1r2[l-i] < aux)
						{
							Paux2[M->NumRows + r1r2[l-i] - M->NumCols] = 1;
						}
						else
							if(r1r2[l-i] < aux1)
							{
								Paux2[aux1] = (-b[r1r2[l-i] - aux] / B2[r1r2[l-i] - aux]);

								for(j = Acol[k = 0]; k < M->NumCols; k++)
								{
									for( ; j < Acol[k+1]; j++)
									{
										if((r1r2[l-i] - aux) == M->Row[j] - 1)
										{
											Paux2[k + auxx] = val2[j];
										}
									}
								}
							}
							else
								if(r1r2[l-i] < aux2)
								{
									Paux2[aux1] = b[r1r2[l-i] - aux1] / B2[r1r2[l-i] - aux1];

									for(j = Acol[k = 0]; k < M->NumCols; k++)
									{
										for( ; j < Acol[k+1]; j++)
										{
											if((r1r2[l-i] - aux1) == M->Row[j] - 1)
											{
												Paux2[k + auxx] = -val2[j];
											}
										}
									}
								}
								else
									if(r1r2[l-i] < aux3)
									{
										Paux2[aux1]                    = upbound[r1r2[l-i] - aux2] / normau[r1r2[l-i] - aux2];
										Paux2[auxx + r1r2[l-i] - aux2] = -1 / (normau[r1r2[l-i] - aux2]);
									}
									else
										if(r1r2[l-i] < aux4)
										{
											Paux2[auxx + r1r2[l-i] - aux3] = 1;
										}
										else
										{
											for(j = 0; j < M->NumRows; j++)
											{
												Paux2[j] = -b[j] / normabc;
											}

											for(j = M->NumRows; j < auxx; j++)
											{
												Paux2[j] = -upbound[j - M->NumRows] / normabc;
											}

											for(j = auxx; j < aux1; j++)
											{
												Paux2[j] = -c[j - auxx] / normabc;
											}
										}

					if(i == 0)
					{
						for(j = 0; j < aux7; j++)
						{
							Paux1[j] = Paux2[j];
						}
					}

					if(l == i)
					{
						valg[((l+1) * (l+2)) / 2]   = ddot(&aux7, wb, &ONE , Paux1, &ONE);
					}
				
					valg[((l+1) * (l+2)) / 2 + (l+1-i)] =  ddot(&aux7, Paux1, &ONE, Paux2, &ONE);

				}
			}

			/*****************************************************************************/
			/* FATORACAO DE CHOLESKY                                                     */
			/*****************************************************************************/

			for (i = 0; i < aux6; i++)
			{
				l1[i]   = 1;
				u[i]    = 1;
				rr1[i]  = 0;
				rr2[i]  = 0;
				rr4[i]  = 0;
				dd12[i] = 0;
				N1[i]   = 1;
			}

			v     = 1;
			N1[0] = s1;
			rr3   = 0;
			cci = 0;

			/*****************************************************************************/
			/* LAÇO WHILE RESOLVE O SUBPROBLEMA PELO METODO PONTOS INTERIORES            */
			/*****************************************************************************/

			while(zz == 1)
			{
				cci++;
				
				cholesky(aux6, pp2, p, s1, v, valgaux, valg, valg1, wwaux, l1, rr1, rr2, rr4, &rr3, u, dd12);

				/*****************************************************************************/
				/*              	RESOLUCAO DOS SISTEMAS LINEARES                      */
				/*****************************************************************************/

				resolucao_sistema(aux6, p, s1, sti, ssf, ssf2, dd12, valg1, rr4);

				c1 = (ddot(&aux6, ssf, &ONE, rr4, &ONE)-rr3) / ddot(&aux6, N1, &ONE, ssf, &ONE); //direcao de v

				mxzero  = 0;
				mxzero1 = 0;
				alpl1   = 1;
				alpu    = 1;
				alpp    = 1;

				for(i = 0; i < aux6; i++)
				{
					dl1[i]      = ssf2[i] - ssf[i] * c1;							// direcao de l1
					du[i]       = rr2[i] * (1 / l1[i]) - (u[i] / l1[i]) * dl1[i];	// direcao u
					dl1posi[i]  = i;
					dl1mzero[i] = 0;
					
					if (dl1[i] < 0.0)
					{
						mxzero++;
						dl1mzero[i] = -l1[i] / dl1[i];
					}
					
					duposi[i]  = i;
					dumzero[i] = 0;
					if(du[i] < 0.0)
					{
						mxzero1++;
						dumzero[i] = -u[i] / du[i];
					}
				}

				sortd(dl1mzero, dl1posi, aux6); 	//ordena em ordem descendente, guarda a posição ori em dl1posi				

				if(mxzero != 0)
				{
					alpl1 = dl1mzero[mxzero - 1];
				}
				
				sortd(dumzero, duposi, aux6); //ordena em ordem descendente, guarda a posição ori em duposi

				if(mxzero1 != 0)
				{
					alpu = dumzero[mxzero1 - 1];
				}

				if(alpp > 0.995 * alpl1)
				{
					alpp = 0.995 * alpl1;
				}
				
				if(alpp > 0.995 * alpu)
				{
					alpp = 0.995 * alpu;
				}

				for(i = 0; i < aux6; i++)
				{
					l1[i]  += alpp * dl1[i]; 			//ATUALZACAO DO l1
					u[i]   += alpp * du[i];		 		//atualizacao u
				}

				v += alpp * c1;    					//atualizacao
				
				wwaux[0] = 0;					
				
				i  = 0;
				p_i = 0;
				
				while(i < aux6)
				{
					wwaux[0] += valg[p_i] * l1[i]; // produto p'p*l1 a factibilidade

					i++;
					p_i += i;
				}
				
				j  = 1;
				pj = 1;
				
				while(j < aux6)
				{
					wwaux[j] = 0;
					
					pk = pj;
					
					for(k = 0; k < j; k++)
					{
						wwaux[j] += valg[pk++] * l1[k];
					}
					
					i  = j;
					p_i = i * (i + 1) / 2 + j;
					
					while(i < aux6)
					{
						wwaux[j] += valg[p_i] * l1[i];  // produto p'p*l1
						
						i++;
						p_i += i;
					}
					
					j++;
					pj += j;
				}
				

				rr1[0] = u[0] - s1 * v - wwaux[0]; // calculo para a factibilidade

				for(i = 1; i < aux6; i++)
				{
					rr1[i] = u[i] - v - wwaux[i]; // calculo para a factibilidade
				}

				v = fabs(v);
				
				if((dnrm2(&aux6, rr1, &ONE) < ep1) && (fabs(1 - ddot(&aux6, N1, &ONE, l1, &ONE)) < ep1) && ((ddot(&aux6, l1, &ONE, u, &ONE) / ( 1 + v + fabs(ddot(&aux6, l1, &ONE, wwaux, &ONE)))) < ep2))
				{
					break;
				}

				if(cci > 200)
				{
				  printf("cci =  %d \n", cci-1);
					break;
				}
			} // termina aqui o laco while

			l1[0] = 1;
			for(i = 1; i < aux6; i++)
			{
				if(u[i] > l1[i])
				{
					l1[i] = 0.0;
				}
				l1[0] -= l1[i];
			}

			l1[0] /= s1;
		
			for(i = 0; i < aux4+1; i++)
			{
				xgrande[i] *= l1[0];
			}

			for(i = 0; i < aux5; i++)
			{
				xgrande[r1r2[i]]      = l1[i + 1];
				xgrande[r1r2[i+aux5]] = l1[aux5 + i + 1];
			}

			vv1 = xgrande[aux4];

			/*****************************************************************************/
			/*           	CALCULANDO O RESIDUO EM CADA ITERAÇÃO                        */
			/*****************************************************************************/

			c11 = -vv1/normabc;

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
		}

	}

/*	for(i = 0; i < M->NumCols; i++)
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

*/	
//	printf("resif  %30.30lf \n ", resiff);
	printf("cc =  %d \n", cc-1);
	
	Free((char *)resif);
	Free((char *)Paux2);
	Free((char *)Paux1);
	Free((char *)wb);
	Free((char *)rr4);
	Free((char *)rr2);
	Free((char *)rr1);
	Free((char *)u);
	Free((char *)l1);
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
	Free((char *)wwa);
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

