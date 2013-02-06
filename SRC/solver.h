/* header file for interface to sparse symm. p.d. system solver
 *
 * PCx 1.0   08/97 
 *
 * Author: Michael Wagner
 * 
 * (C) 1997 University of Chicago. See COPYRIGHT in main directory.
 */

/**********************************************************************/
/*                                                                    */
/* For sample implementations of these specifications see the files   */
/*                                                                    */
/*        wssmp.c      and/or          NgPeyton.c                     */
/*                                                                    */
/**********************************************************************/


/**********************************************************************/
/*                                                                    */
/* Allocation and deallocation routines for FactorType                */
/*                                                                    */
/* Also allocates and deallocates storage spage behind                */ 
/* FactorType->ptr                                                    */
/*                                                                    */
/**********************************************************************/

/* these routines are solver-specific since the storage spage pointed */
/* at by FactorType->ptr might need to be allocated and initialized   */
/* (for an example of this see NgPeyton.c)                            */

FactorType *NewFactorType(/* MMTtype *A, int Ndense, int NumCols */);

void FreeFactorType(/* FactorType *Factor */);

/**********************************************************************/
/*                                                                    */
/* Ordering and Symbolic Factorization routine                        */
/*                                                                    */
/* this routine should set Factor->Perm and Factor->InvPerm           */
/*                                                                    */
/* also should set Factor->NonzerosL                                  */
/*                                                                    */
/**********************************************************************/

int Order(/* FactorType *Factor */);

/**********************************************************************/
/*                                                                    */
/* Cholesky factorization routine                                     */
/*                                                                    */
/* should perform numerical factorization (with singularity handling  */
/*                                                                    */
/**********************************************************************/

int Factorize(/* FactorType *Factor, Parameters *Inputs */);

/**********************************************************************/
/*                                                                    */
/* Sparse symmetric system solution routine                           */
/*                                                                    */
/* solves Factor->AAT * Solution = rhs                                */
/*                                                                    */
/* assumes Cholesky Factor is calculated                              */
/*                                                                    */
/**********************************************************************/

int Solve(/* FactorType *Factor, double *rhs, double *Solution */);

/**********************************************************************/ 
/*                                                                    */
/* Computes auxiliary matrices for Sherman Morrison formula to handle */
/*   dense part of system  (see manual page 9, formulas 28 and        */
/*   after (29)                                                       */
/*                                                                    */
/* The reason this is in the interface is that it depends on the fact */
/* whether routines to do a single forward- or backward subsitution   */
/* are available.                                                     */
/*                                                                    */
/* If yes, then the formula after (29) can be used and                */
/*    W = L^-1 *  P *  Aden           and                             */
/*    Lden = Cholesky factor of D^-2 + W'*W                           */
/*                                                                    */
/* (for an example of this case see NgPeyton.c)                       */
/*                                                                    */
/* If not, then (28) must be used and                                 */
/*    W = M^-1*Aden                   and                             */
/*    Lden = Cholesky factor of D^-2 + Aden'*W                        */
/*                                                                    */
/* (for an example for this case see wssmp.c)                         */
/*                                                                    */
/* W and Lden refer to Factor->W and Factor-Lden respectively         */
/**********************************************************************/

int ComputeWandLdense(/* MMTtype *Adense, FactorType *Factor, 
		      double *scale, int NumCols */);

/**********************************************************************/
/*                                                                    */
/* Solves the spd-equation (A*scale*A^t)*Solution=rhs via Sherman-    */
/* Morrison formula, using Factor->W and Factor->Lden from above.     */
/*                                                                    */
/* for examples see wssmp.c and NgPeyton.c                            */
/*                                                                    */
/**********************************************************************/

int EnhancedSolve(/* FactorType *Factor, double *rhs, double *Solution */);

/**********************************************************************/


