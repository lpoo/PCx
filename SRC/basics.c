/* basic sparse linear algebra routines 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include "memory.h"

/* FindResidualSolutionofNormalEquations.c, v1.0, 09/01/91, Sanjay Mehrotra */

/* Purpose: Finds r = ADA^T x - b.
 * 
 * KNOWN BUG: The matrix may be singular numerically but not actually.  */


/* Description of function: LDL^T=H.  It is assumed that only the lower half
 * of H (including the diagonal elements) is available and it is stored
 * columnwise.
 * 
 * The factor returned from this function is also stored columnwise.
 * 
 * LOGIC USED: The implementation uses PULL FROM BEHIND approach. It assumes
 * that the rows in the cholesky factor have been permuted in the increasing
 * order of row indices, and the structure of the transpose matrix is
 * available. */

#define ERROR(CODE)  \
  { error_code=CODE; \
  fprintf(stdout, "FindResidualSolutionofNormalEquations: error %d\n", \
        error_code); \
  Free((char *) internal_real_space_ps); \
  return error_code; }

/* extern int      CopyRealVectorxToy(); */

int             
FindResidualSolutionofNormalEquations (a_p, ira_p, ipbra_p, ipera_p, 
				       scale_p, solution_p, rhs_p,
				       residual_p, nrow_p, ncol_p)
     
     double         *a_p, *scale_p, *solution_p, *rhs_p, *residual_p;
     int            *ira_p, *ipbra_p, *ipera_p;
     int            *nrow_p, *ncol_p;
{
   double         *internal_real_space_p, *internal_real_space_ps;
   int             internal_real_space_size, error_code, i;
   double         *atsolution_p, *scaleatsolution_p;
   
   internal_real_space_size = 2 * *ncol_p;
   
   internal_real_space_p = NewDouble(internal_real_space_size,
	    "internal_real_space_p in FindResidualSolutionofNormalEquations");
   internal_real_space_ps = internal_real_space_p;
   atsolution_p = internal_real_space_p;
   scaleatsolution_p = atsolution_p + *ncol_p;
   
   if (RealSparseMatrixTransposeVectorProduct
       (a_p, ipbra_p, ipera_p, ira_p, solution_p,
	atsolution_p, nrow_p, ncol_p))
      ERROR(3);
   
   if (DiagonalMatrixTimesDiagonalMatrix(scale_p, atsolution_p, 
					 scaleatsolution_p, ncol_p))
      ERROR(4);
   
   if (CopyRealVectorxToy(rhs_p, residual_p, nrow_p))
      ERROR(5);

   if (NegateRealVector(residual_p, nrow_p))
      ERROR(6);
   
   if (RealSparseMatrixVectorProductPlusx
       (a_p, ipbra_p, ipera_p, ira_p, scaleatsolution_p,
	residual_p, nrow_p, ncol_p))
      ERROR(7);
   
   Free((char *) internal_real_space_ps);
   return 0;
}

/********************************************************************/

/* FindStructureNormalEquations.c, v1.0, 11/08/91, Sanjay Mehrotra */


/* The function in this file finds the structure of non-zeros in the
 * equations formed as H = AA^T, where A is matrix whose non-zero structure
 * is known.
 * 
 * Assumption: The matrix A is stored column-wise.
 * 
 * Input: Begining and end of each column is stored in array ipbra[] and
 * ipera[]. All the array information is passed through pointers.  These
 * pointers must point to the first element of the array.  We assume that the
 * row/column indices run from 1..m/n. As oppose to 0..m-1, which is the
 * standard C.
 * 
 * Output: The structure of equations is given provided that it does not exceed
 * a prespecified size.  If the equations exceed a pre-specified size, a flag
 * is returned and columns of equations for which the structure have been
 * found (together with appropriate pointers) are returned.  */

/* Logic...KNOWN BUG: Currently the structure of normal equations is computed
 * by using a general purpose code which gives the strucuture of H = AC,
 * where C passed as A^T.  Therefore, we do not take advantage of the
 * symmetry of AA^T, which makes the code slower, by a factor of two.  It is
 * possible to write a code that achieves this efficiency by using a logic
 * similar to the logic used in the general purpose routine, but we have not
 * done it here because the code becomes complicated when only a subset of
 * columns are needed.
 * 
 * Gaining the efficieny here would improve the overall efficiency of the LP
 * code by VERY LITTLE AMOUNT.  However, if we want to solve only one
 * least-squares problems using this package, it would be worth spending some
 * more efforts here. */
/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      TransposeStructureSparseMatrix();
extern int      FindStructureHEqualAC();

int             
FindStructureNormalEquations(ipbra_p, ipera_p, ira_p,
	                     ipbrh_p, iperh_p, irh_p, nrow_p, ncol_p, nza_p,
			     maximum_nonzeros_p, flag_maximum_nonzeros_p)
     int      *ipbra_p, *ipera_p, *ira_p, *ipbrh_p, *iperh_p, 
              *irh_p, *maximum_nonzeros_p, *flag_maximum_nonzeros_p;
     int      *ncol_p, *nrow_p, *nza_p;
{
   int            *ipbrat_p, *iperat_p, *irat_p, *ipbrat_ps, *iperat_ps,
                  *irat_ps;
   int             error_code;

   *flag_maximum_nonzeros_p = 0;
   
   /* Allocated space for arrays used internally */
   
   ipbrat_p = NewInt(*nrow_p, "ipbrat_p in FindStructureNormalEquations");
   ipbrat_ps = ipbrat_p;
   
   iperat_p = NewInt(*nrow_p, "iperat_p in FindStructureNormalEquations");
   iperat_ps = iperat_p;
   
   irat_p = NewInt(*nza_p, "irat_p in FindStructureNormalEquations");
   irat_ps = irat_p;
   
   /* Take the transpose of the matrix */
   
   
   if (TransposeStructureSparseMatrix
       (ipbra_p, ipera_p, ira_p, ipbrat_p, iperat_p, irat_p,
	nrow_p, ncol_p)) 
      {
	 error_code = 5;
	 
	 fprintf(stdout, "Error: FindStructureNormalEquations:\n");
	 fprintf(stdout, 
		 "  in TransposeStructureSparseMatrix copying vector.\n");
	 
	 Free((char *) irat_ps);
	 Free((char *) iperat_ps);
	 Free((char *) ipbrat_ps);
	 return error_code;
      }
  /* Compute the strucutre of AA^T */


   if (FindStructureHEqualAC (ipbra_p, ipera_p, ira_p, ipbrat_p, 
			      iperat_p, irat_p, ipbrh_p, iperh_p, 
			      irh_p, nrow_p, ncol_p, maximum_nonzeros_p, 
			      flag_maximum_nonzeros_p)) 
      {
	 error_code = 6;
	 
	 fprintf(stdout, "Error: FindStructureNormalEquations:\n");
	 fprintf(stdout, "  in FindStructureHEqualAC initializing array.\n");
	 
	 Free((char *) irat_ps);
	 Free((char *) iperat_ps);
	 Free((char *) ipbrat_ps);
	 return error_code;
      }
   Free((char *) irat_ps);
   Free((char *) iperat_ps);
   Free((char *) ipbrat_ps);
   return 0;
}

/********************************************************************/

/* FindNonzeroNormalEquations.c, v1.0, 11/08/91, Sanjay Mehrotra */

/* The function in this file finds the number of nonzeros in the normal
 * equations H=AA^T.  Symmetry is ignored while finding this requirement.
 * 
 * Assumption: The matrix A is stored column-wise.
 * 
 * Input: Begining and end of each column is stored in array ipbra[] and
 * ipera[]. All the array information is passed through pointers.  These
 * pointers must point to the first element of the array.  We assume that the
 * row/column indices run from 1..m/n. As oppose to 0..m-1, which is the
 * standard C.
 * 
 * Output: The structure of equations is given provided that it does not exceed
 * a prespecified size.  If the equations exceed a pre-specified size, a flag
 * is returned and columns of equations for which the structure have been
 * found (together with appropriate pointers) are returned.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      TransposeNonzeroSparseMatrix();
extern int      FindNonzeroHEqualAC();

int             
FindNonzeroNormalEquations(ipbra_p, ipera_p, ira_p, nrow_p, ncol_p, nza_p,
			   nonzero_normal_equations_p)
     int     *ipbra_p, *ipera_p, *ira_p, *nonzero_normal_equations_p;
     int     *ncol_p, *nrow_p, *nza_p;
{
   int   *ipbrat_p, *iperat_p, *irat_p, *ipbrat_ps, *iperat_ps,
         *irat_ps;
   int    error_code;
  
   /* Allocated space for arrays used internally */
   
   ipbrat_p = NewInt(*nrow_p, "ipbrat_p in FindNonzeroNormalEquations");
   ipbrat_ps = ipbrat_p;
   
   iperat_p = NewInt(*nrow_p, "iperat_p in FindNonzeroNormalEquations");
   iperat_ps = iperat_p;
   
   irat_p = NewInt(*nza_p, "irat_p in FindNonzeroNormalEquations");
   irat_ps = irat_p;
   
   /* Take the transpose of the matrix */
   
   if (TransposeStructureSparseMatrix(ipbra_p, ipera_p, ira_p,
				      ipbrat_p, iperat_p, irat_p, 
				      nrow_p, ncol_p)) 
      {
	 error_code = 5;
	 
	 fprintf(stdout, "Error: FindNonzeroNormalEquations:\n");
	 fprintf(stdout, 
		 " in TransposeNonzeroSparseMatrix copying a vector.\n");
	 
	 Free((char *) irat_ps);
	 Free((char *) iperat_ps);
	 Free((char *) ipbrat_ps);
	 return error_code;
      }
   /* Compute the strucutre of AA^T */
   
   if (FindNonzeroHEqualAC
       (ipbra_p, ipera_p, ira_p, ipbrat_p, iperat_p, irat_p,
	nrow_p, ncol_p, nonzero_normal_equations_p)) {
      error_code = 6;
      
      fprintf(stdout, "Error: FindNonzeroNormalEquations:\n");
      fprintf(stdout, "  in FindNonzeroHEqualAC initializing an array.\n");
      
      Free((char *) irat_ps);
      Free((char *) iperat_ps);
      Free((char *) ipbrat_ps);
      return error_code;
   }
   Free((char *) irat_ps);
   Free((char *) iperat_ps);
   Free((char *) ipbrat_ps);
   return 0;
}

/********************************************************************/

/* FindStructureHEqualAC.c, v1.0, 11/08/91, Sanjay Mehrotra */


/* The function in this file finds the structure of non-zeros in the
 * equations formed as H = AC, where A is matrix whose non-zero structure is
 * known.
 * 
 * Assumption: The matrix A is stored column-wise.
 * 
 * Input: Begining and end of each column is stored in array ipbra[] and
 * ipera[]. All the array information is passed through pointers.  These
 * pointers must point to the first element of the array.  We assume that the
 * row/column indices run from 1..m/n. As oppose to 0..m-1, which is the
 * standard C.
 * 
 * Output: The structure of equations is given provided that it does not exceed
 * a prespecified size.  If the equations exceed a pre-specified size, a flag
 * is returned and columns of equations for which the structure have been
 * found (together with appropriate pointers) are returned.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

/* Description of function: Permute Columns of A by (AT)T. */

extern int      ZeroIntegerSparseVector();
extern int      CopyIntegerVectorxToy();

int             
FindStructureHEqualAC(ipbra_p, ipera_p, ira_p, ipbca_p, ipeca_p, ica_p,
		      ipbrh_p, iperh_p, irh_p, nrow_p, ncol_p,
		      maximum_nonzeros_p, flag_maximum_nonzeros_p)
     int    *ipbra_p, *ipera_p, *ira_p, *ipbca_p, *ipeca_p, *ica_p, *ipbrh_p,
            *iperh_p, *irh_p, *maximum_nonzeros_p, *flag_maximum_nonzeros_p;
     int    *ncol_p, *nrow_p;
{
   int      *current_column_p, *current_column_ps, *iflag_p, *iflag_ps;
   int      *ipbra_q, *ipera_q, *ica_q, *ira_q, *iflag_q;
   int      *jbeg_p, *jend_p;
   int       i, irow, icol, kbeg, kend, nonzero_current_column, 
             nonzero_equations;
   int       error_code;

   *flag_maximum_nonzeros_p = 0;

   /* Allocated space for arrays used internally */
   
   current_column_p =
      NewInt(*nrow_p, "current_column_p in FindStructureHEqualAC");
   current_column_ps = current_column_p;
   
   iflag_p = NewInt(*nrow_p, "iflag_p in FindStructureHEqualAC");
   iflag_ps = iflag_p;
   
   /* let us now find the normal equation structure */
   
   ipbra_q = ipbra_p - 1;
   ipera_q = ipera_p - 1;
   ica_q = ica_p - 1;
   ira_q = ira_p - 1;
   *iperh_p = 0;
   iflag_q = iflag_p - 1;
   nonzero_equations = 0;
   
   if (ZeroIntegerDenseVector(iflag_p, nrow_p)) 
      {
	 error_code = 7;
	 
	 fprintf(stdout, "Error: FindStructureHEqualAC:\n");
	 fprintf(stdout, " in ZeroIntegerDenseVector \
                           initializing an array.\n");
	 
	 return error_code;
      }
   for (i = 1; i <= *nrow_p; i++) 
      {
	 /* computation for the next column in equations */
	 jbeg_p = ica_q + *ipbca_p;
	 jend_p = ica_q + *ipeca_p;
	 nonzero_current_column = 0;

	 /* added by Steve Wright 8/12/99. Force allocation of
	    space for the diagonal element of the product AC */

	 irow = i;
	 *(iflag_q + irow) = 1;
	 nonzero_current_column++;
	 *current_column_p = irow;
	 current_column_p++;
	 /*
	 printf(" adding diagonal entry for row %d \n", irow);
	 */
	 
	 /* end added segment */

	 /* start the computation for the structure of the current column */
	 for (; jbeg_p <= jend_p; jbeg_p++) 
	    {
	       icol = *jbeg_p;
	       kbeg = *(ipbra_q + icol);
	       kend = *(ipera_q + icol);

	       for (; kbeg <= kend; kbeg++) 
		  {
		     /* take the union of this column with the 
			previous columns */
		     irow = *(ira_q + kbeg);
		     if ((*(iflag_q + irow) == 0)) 
			{
			   *(iflag_q + irow) = 1;
			   nonzero_current_column++;
			   *current_column_p = irow;
			   current_column_p++;
			}
		  }
	    }
	 
	 /* reset the flag array for computing the next column */
	 
	 current_column_p -= nonzero_current_column;
	 if (ZeroIntegerSparseVector(iflag_p, nrow_p, current_column_p,
				     &nonzero_current_column)) 
	    {
	       error_code = 4;
	       
	       fprintf(stdout, "Error: FindStructureHEqualAC:\n");
	       fprintf(stdout, 
		       " in ZeroIntegerSparseVector initializing an array.\n");
	       
	       return error_code;
	    }
	 /* if there is enough space copy the column 
	    to the space for columns */
	 
	 nonzero_equations += nonzero_current_column;
	 if (nonzero_equations <= *maximum_nonzeros_p) 
	    {
	       if (CopyIntegerVectorxToy(current_column_p, irh_p,
					 &nonzero_current_column)) 
		  {
		     error_code = 5;
		     
		     fprintf(stdout, "Error: FindStructureHEqualAC:\n");
		     fprintf(stdout, 
			     " in CopyIntegerVectorxToy copying a vector.\n");
		     
		     return error_code;
		  }
	       *iperh_p = nonzero_equations;
	       iperh_p++;
	       *ipbrh_p = nonzero_equations - nonzero_current_column + 1;
	       ipbrh_p++;
	       irh_p = irh_p + nonzero_current_column;
	    } 
	 else 
	    {
	       *flag_maximum_nonzeros_p = i - 1;
	       Free((char *) iflag_ps);
	       Free((char *) current_column_ps);
	       return 0;
	    }
	 
	 /* advance pointers for the computation of next column */
	 ipbca_p++;
	 ipeca_p++;
      }
   
   ipbca_p -= *nrow_p;
   ipeca_p -= *nrow_p;
   
   Free((char *) iflag_ps);
   Free((char *) current_column_ps);
   return 0;
}

/********************************************************************/

/* FindNonzeroHEqualAC.c, v1.0, 11/08/91, Sanjay Mehrotra */

/* The function in this file finds the structure of non-zeros in the
 * equations formed as H = AC, where A is matrix whose non-zero structure is
 * known.
 * 
 * Assumption: The matrix A is stored column-wise.
 * 
 * Input: Begining and end of each column is stored in array ipbra[] and
 * ipera[]. All the array information is passed through pointers.  These
 * pointers must point to the first element of the array.  We assume that the
 * row/column indices run from 1..m/n. As oppose to 0..m-1, which is the
 * standard C.
 * 
 * Output: The structure of equations is given provided that it does not exceed
 * a prespecified size.  If the equations exceed a pre-specified size, a flag
 * is returned and columns of equations for which the structure have been
 * found (together with appropriate pointers) are returned.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      ZeroIntegerSparseVector();
extern int      CopyIntegerVectorxToy();

int             
FindNonzeroHEqualAC(ipbra_p, ipera_p, ira_p, ipbca_p, ipeca_p, ica_p, 
		    nrow_p, ncol_p, nonzero_equations_p)
     int       *ipbra_p, *ipera_p, *ira_p, *ipbca_p, *ipeca_p, *ica_p;
     int       *ncol_p, *nrow_p, *nonzero_equations_p;
{
   int    *current_column_p, *current_column_ps, *iflag_p, *iflag_ps;
   int    *ipbra_q, *ipera_q, *ica_q, *ira_q, *iflag_q;
   int    *jbeg_p, *jend_p;
   int     i, irow, icol, kbeg, kend, nonzero_current_column;
   int     error_code;

   /* Allocated space for arrays used internally */
   
   current_column_p =
      NewInt(*nrow_p, "current_column_p in FindNonzeroHEqualAC");
   current_column_ps = current_column_p;
   
   iflag_p = NewInt(*nrow_p, "iflag_p in FindNonzeroHEqualAC");
   iflag_ps = iflag_p;
   
   /* let us now find the normal equation structure */
   
   ipbra_q = ipbra_p - 1;
   ipera_q = ipera_p - 1;
   ica_q = ica_p - 1;
   ira_q = ira_p - 1;
   iflag_q = iflag_p - 1;
   *nonzero_equations_p = 0;
   
   if (ZeroIntegerDenseVector(iflag_p, nrow_p)) 
      {
	 error_code = 5;
	 fprintf(stdout, "Error: FindNonzeroHEqualAC:\n");
	 fprintf(stdout, " in ZeroIntegerDenseVector initializing array.\n");
	 return error_code;
      }
   for (i = 1; i <= *nrow_p; i++) 
      {
	 
	 /* computation for the next column in equations */
	 
	 jbeg_p = ica_q + *ipbca_p;
	 jend_p = ica_q + *ipeca_p;
	 nonzero_current_column = 0;

	 /* added by Steve Wright 8/12/99. Force allocation of
	    space for the diagonal element of the product AC */
	 
	 irow = i;
	 *(iflag_q + irow) = 1;
	 nonzero_current_column++;
	 *current_column_p = irow;
	 current_column_p++;
	 /*
	 printf(" counting diagonal entry for row %d \n", irow);
	 */
	 
	 /* end added segment */
	 
	 /* start the computation for the structure of the current column */
	 for (; jbeg_p <= jend_p; jbeg_p++) 
	    {
	       icol = *jbeg_p;
	       kbeg = *(ipbra_q + icol);
	       kend = *(ipera_q + icol);
	       for (; kbeg <= kend; kbeg++) 
		  {
		     /* take the union of this column with 
			the previous columns */
		     irow = *(ira_q + kbeg);
		     if ((*(iflag_q + irow) == 0)) 
			{
			   *(iflag_q + irow) = 1;
			   nonzero_current_column++;
			   *current_column_p = irow;
			   current_column_p++;
			}
		  }
	    }
	 
	 /* reset the flag array for computing the next column */
	 current_column_p -= nonzero_current_column;
	 if (ZeroIntegerSparseVector(iflag_p, nrow_p, current_column_p,
				     &nonzero_current_column)) 
	    {
	       error_code = 4;
	       
	       fprintf(stdout, "Error: FindNonzeroHEqualAC:\n");
	       fprintf(stdout, 
		       " in ZeroIntegerSparseVector initializing array.\n");
	       
	       return error_code;
	    }
	 /* if there is enough space copy the column 
	    to the space for columns */
	 (*nonzero_equations_p) += nonzero_current_column;
	 
	 /* advance pointers for the computation of next column */
	 ipbca_p++;
	 ipeca_p++;
      }
   
   ipbca_p -= *nrow_p;
   ipeca_p -= *nrow_p;
   
   Free((char *) iflag_ps);
   Free((char *) current_column_ps);
   return 0;
}

/* RealSparseMatrixTransposeVectorProduct.c, v1.0, 08/27/91, Sanjay Mehrotra */

/* The function in this file computes a matrix vector product.  The matrix a
 * is assumed to be stored row-wise and the given vector q is assumed to be
 * dense. The result is written in vector aq.  The begining and end of row i
 * in a is stored in ipbra[i] and ipera[i]. The column indices corresponding
 * to each row in a are stored in ira[]. All the array information is passed
 * through pointers.  The pointers point to the first element of the array.
 * The subroutine assumes that the row/column indices run from 1..m/n.  As
 * oppose to 0..m-1, which is the standard C.  All changes required are done
 * internally in a subroutine.  */

/* An alternative way of looking at this function is that it multiplies the
 * transpose of a matrix with a vector.  The matrix is assume to have been
 * stored column-wise.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

int             
RealSparseMatrixTransposeVectorProduct(a_p, ipbra_p, ipera_p, ira_p, 
				       q_p, aq_p, nrow_p, ncol_p)
     double         *a_p, *q_p, *aq_p;
     int            *ipbra_p, *ipera_p, *ira_p;
     int            *ncol_p, *nrow_p;
{
   int            *ira_ps;
   int             jbeg, jend, irow;
   double         *q_ps, *aq_ps, *a_ps;
   double          qval;
   
   /* indices start from 1 not zero so shift relevant pointers */
   a_ps = a_p - 1;
   ira_ps = ira_p - 1;
   q_ps = q_p - 1;
   
   /* now let us start the matrix vector product 
      accumulate one row at a time.  */

   aq_ps = aq_p + *ncol_p;
   for (; aq_p < aq_ps; aq_p++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 qval = 0.;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       qval += *(a_ps + jbeg) * *(q_ps + irow);
	    }
	 *aq_p = qval;
      }
   return 0;
}

/********************************************************************/

/* RealSparseMatrixTransposeVectorProductPlusx.c, v1.0, 08/27/91, Sanjay
 * Mehrotra */

/* The function in this file computes a matrix vector product.  The matrix a
 * is assumed to be stored row-wise and the given vector q is assumed to be
 * dense. The result is written in vector aq.  The begining and end of row i
 * in a is stored in ipbra[i] and ipera[i]. The column indices corresponding
 * to each row in a are stored in ira[]. All the array information is passed
 * through pointers.  The pointers point to the first element of the array.
 * The subroutine assumes that the row/column indices run from 1..m/n.  As
 * oppose to 0..m-1, which is the standard C.  All changes required are done
 * internally in a subroutine.  */

/* An alternative way of looking at this function is that it multiplies the
 * transpose of a matrix with a vector.  The matrix is assume to have been
 * stored column-wise.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

int             
RealSparseMatrixTransposeVectorProductPlusx(a_p, ipbra_p, ipera_p, ira_p, 
					    q_p, aq_p, nrow_p, ncol_p)
     double         *a_p, *q_p, *aq_p;
     int            *ipbra_p, *ipera_p, *ira_p;
     int            *ncol_p, *nrow_p;
{
   int            *ira_ps;
   int             jbeg, jend, irow;
   double         *q_ps, *aq_ps, *a_ps;
   double          qval;
   
   /* initialize *nrow_p elements of aq */
   
   /* indices start from 1 not zero so shift relevant pointers */
   a_ps = a_p - 1;
   ira_ps = ira_p - 1;
   q_ps = q_p - 1;
   
   /* now let us start the matrix vector 
      product accumulate one row at a time.  */

   aq_ps = aq_p + *ncol_p;
   for (; aq_p < aq_ps; aq_p++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 qval = 0.;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       qval += *(a_ps + jbeg) * *(q_ps + irow);
	    }
	 *aq_p += qval;
      }
   return 0;
}

/********************************************************************/


/* RealSparseMatrixVectorProduct.c, v1.0, 08/26/91, Sanjay Mehrotra */

/* The function in this file computes a matrix vector product.  The matrix a
 * is assumed to be stored column-wise and the given vector q is assumed to
 * be dense. The result is written in vector aq.  The begining and end of
 * column i in a is stored in ipbra[i] and ipera[i].  The row indices
 * corresponding to each column in a are stored in ira[]. All the array
 * information is passed through pointers.  The pointers point to the first
 * element of the array. The subroutine assumes that the row/column indices
 * run from 1..m/n. As oppose to 0..m-1, which is the standard C.  All
 * changes required are done internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      ZeroRealDenseVector();

int             
RealSparseMatrixVectorProduct(a_p, ipbra_p, ipera_p, ira_p, q_p, aq_p,
			      nrow_p, ncol_p)
     double         *a_p, *q_p, *aq_p;
     int            *ipbra_p, *ipera_p, *ira_p;
     int            *ncol_p, *nrow_p;
{
   int             error_code;
   int            *ira_ps;
   int             jbeg, jend, irow;
   double         *q_ps, *aq_ps, *a_ps;
   double          qval;
   
   /* initialize *nrow_p elements of aq */
   
   if (ZeroRealDenseVector(aq_p, nrow_p)) 
      {
	 error_code = 2;
	 fprintf(stdout, "Error: RealSparseMatrixVectorProduct: Error\n");
	 fprintf(stdout, " in ZeroRealDenseVector initializing vector.\n");
	 return error_code;
      }
   /* indices start from 1 not zero so shift relevant pointers */
   
   aq_ps = aq_p - 1;
   a_ps = a_p - 1;
   ira_ps = ira_p - 1;
   
   /* now let us start the matrix vector product accumulate one column at a
    * time.  */
   
   q_ps = q_p + *ncol_p;
   for (; q_p < q_ps; q_p++) 
      {
	 qval = *q_p;
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       *(aq_ps + irow) += *(a_ps + jbeg) * qval;
	    }
      }
   return 0;
}

/********************************************************************/

/* RealSparseMatrixVectorProductPlusx.c, v1.0, 08/26/91, Sanjay Mehrotra */

/* The function in this file computes a matrix vector product and adds it to
 * aq.  The matrix a is assumed to be stored column-wise and the given vector
 * q is assumed to be dense. The result is written in vector aq.  The
 * begining and end of column i in a is stored in ipbra[i] and ipera[i].  The
 * row indices corresponding to each column in a are stored in ira[]. All the
 * array information is passed through pointers.  The pointers point to the
 * first element of the array.  The subroutine assumes that the row/column
 * indices run from 1..m/n.  As oppose to 0..m-1, which is the standard C.
 * All changes required are done internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      ZeroRealDenseVector();

int             
RealSparseMatrixVectorProductPlusx(a_p, ipbra_p, ipera_p, ira_p, q_p, 
				   aq_p, nrow_p, ncol_p)
     double         *a_p, *q_p, *aq_p;
     int            *ipbra_p, *ipera_p, *ira_p;
     int            *ncol_p, *nrow_p;
{
   int            *ira_ps;
   int             jbeg, jend, irow;
   double         *q_ps, *aq_ps, *a_ps;
   double          qval;
   
   /* indices start from 1 not zero so shift relevant pointers */
   
   aq_ps = aq_p - 1;
   a_ps = a_p - 1;
   ira_ps = ira_p - 1;
   
   /* now let us start the matrix vector product accumulate one column at a
    * time.  */
   
   q_ps = q_p + *ncol_p;
   for (; q_p < q_ps; q_p++) 
      {
	 qval = *q_p;
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       *(aq_ps + irow) += *(a_ps + jbeg) * qval;
	    }
      }
   return 0;
}

/********************************************************************/

/* TransposeStructureSparseMatrix.c, v1.0, 08/27/91, Sanjay Mehrotra */

/* The function in this file finds the structure of the transpose of a
 * matrix. The matrix a[] is assumed to be stored column-wise. The row
 * indices corresponding to each column in a[] are stored in ira[]. Begining
 * and end of each column is stored in array ipbra[] and ipera[]. All the
 * array information is passed through pointers. These pointers must point to
 * the first element of the array.  The subroutine assumes that the
 * row/column indices run from 1..m/n.  As oppose to 0..m-1, which is the
 * standard C.  All changes required are done internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      ZeroIntegerDenseVector();

int             
TransposeStructureSparseMatrix(ipbra_p, ipera_p, ira_p, ipbrat_p,
			       iperat_p, irat_p, nrow_p, ncol_p)
     int            *ipbra_p, *ipera_p, *ira_p, *ipbrat_p, *iperat_p, *irat_p;
     int            *ncol_p, *nrow_p;
{
   int             error_code;
   int            *ipbra_ps, *ipera_ps, *ira_ps, *ipbrat_ps, *iperat_ps,
                  *irat_ps;
   int             i, irow, iloc, jbeg, jend;
   
   /* iperat is used as a temporary array in earlier part of the code */
   
   /* initialize  nrow_p elements of aq */

   if (ZeroIntegerDenseVector(iperat_p, nrow_p)) 
      {
	 error_code = 2;
	 fprintf(stdout, "Error: TransposeStructureSparseMatrix:\n");
	 fprintf(stdout, " in ZeroIntegerDenseVector initializing vector.\n");
	 
	 return error_code;
      }
   /* indices start from 1 not zero so shift relevant pointers */
   
   ira_ps = ira_p - 1;
   irat_ps = irat_p - 1;
   ipbrat_ps = ipbrat_p - 1;
   iperat_ps = iperat_p - 1;
   
   /* save certain pointers */
   
   ipbra_ps = ipbra_p;
   ipera_ps = ipera_p;
   
   /* first we find the number of nonzeros in each row of a[].  This is same
    * as the number of nonzeros in each column of at[].  */
   
   for (i = 0; i < *ncol_p; i++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       (*(iperat_ps + irow))++;
	    }
      }
   ipbra_p = ipbra_ps;
   ipera_p = ipera_ps;
   
   /* Now let us find the running sum */
   
   iperat_ps = iperat_p + *nrow_p;
   *iperat_p++;
   while (iperat_p < iperat_ps) 
      {
	 /* starting with 1 because first element remains same */
	 
	 (*iperat_p) += *(iperat_p - 1);
	 iperat_p++;
      }
   iperat_p = iperat_ps - *nrow_p;
   iperat_ps = iperat_p - 1;
   
   /* Set begining pointers for at.  These would change and would be restored
    * later.  */
   
   iperat_ps = iperat_ps + *nrow_p;
   *ipbrat_p++ = 0;
   while (iperat_p < iperat_ps) 
      {
	 /* Note: we go till  *nrow_p-1 */
	 *ipbrat_p++ = *iperat_p++;
      }
   iperat_ps = iperat_ps - *nrow_p;
   iperat_p = iperat_ps + 1;
   ipbrat_p = ipbrat_ps + 1;
   
   /* now  build the transpose structure */
   
   for (i = 1; i <= *ncol_p; i++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       iloc = *(ipbrat_ps + irow);
	       *(ipbrat_ps + irow) = ++iloc;
	       *(irat_ps + iloc) = i;
	    }
      }
   ipbra_p = ipbra_ps;
   ipera_p = ipera_ps;
   
   /* reset the pointers */
   iperat_ps = iperat_ps + *nrow_p;
   *ipbrat_p++ = 1;
   while (iperat_p < iperat_ps) 
      {
	 /* Note: we go till  *nrow_p-1 */
	 *ipbrat_p++ = *iperat_p++ + 1;
      }
   iperat_ps = iperat_ps - *nrow_p;
   iperat_p = iperat_ps + 1;
   ipbrat_p = ipbrat_ps + 1;
   
   return 0;
}

/********************************************************************/

/* TransposeSparseRealMatrix.c, v1.0, 08/27/91, Sanjay Mehrotra */

/* The function in this file finds the transpose of a matrix. The matrix a[]
 * is assumed to be stored column-wise. The row indices corresponding to each
 * column in a[] are stored in ira[]. Begining and end of each column is
 * stored in array ipbra[] and ipera[]. All the array information is passed
 * through pointers.  These pointers must point to the first element of the
 * array.  The subroutine assumes that the row/column indices run from
 * 1..m/n.  As oppose to 0..m-1, which is the standard C.  All changes
 * required are done internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

extern int      ZeroIntegerDenseVector();

int             
TransposeSparseRealMatrix(a_p, ipbra_p, ipera_p, ira_p, at_p, ipbrat_p, 
			  iperat_p, irat_p, nrow_p, ncol_p)
     double         *a_p, *at_p;
     int            *ipbra_p, *ipera_p, *ira_p, *ipbrat_p, *iperat_p, *irat_p;
     int            *ncol_p, *nrow_p;
{
   int             error_code;
   double         *a_ps, *at_ps;
   int            *ipbra_ps, *ipera_ps, *ira_ps, *ipbrat_ps, *iperat_ps,
                  *irat_ps;
   int             i, irow, iloc, jbeg, jend;
   
  /* iperat is used as a temporary array in earlier part of the code */

  /* initialize *nrow_p elements of aq */
   if (ZeroIntegerDenseVector(iperat_p, nrow_p)) 
      {
	 error_code = 2;
	 fprintf(stdout, "Error: TransposeSparseRealMatrix: in\n");
	 fprintf(stdout, " in ZeroIntegerDenseVector initializing vector.\n");
	 return error_code;
      }
   /* indices start from 1 not zero so shift relevant pointers */
   ira_ps = ira_p - 1;
   irat_ps = irat_p - 1;
   ipbrat_ps = ipbrat_p - 1;
   iperat_ps = iperat_p - 1;
   a_ps = a_p - 1;
   at_ps = at_p - 1;
   
   /* save certain pointers */
   ipbra_ps = ipbra_p;
   ipera_ps = ipera_p;
   
   /* first we find the number of nonzeros in each row of a[].  This is same
    * as the number of nonzeros in each column of at[].  */
   
   for (i = 0; i < *ncol_p; i++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       (*(iperat_ps + irow))++;
	    }
      }
   ipbra_p = ipbra_ps;
   ipera_p = ipera_ps;
   
   /* Now let us find the running sum */
   iperat_ps = iperat_p + *nrow_p;
   *iperat_p++;
   while (iperat_p < iperat_ps) 
      {
	 /* starting with 1 because first element
	  * remains same */
	 *iperat_p += *(iperat_p - 1);
	 iperat_p++;
      }
   iperat_p = iperat_ps - *nrow_p;
   iperat_ps = iperat_p - 1;
   
   /* Set begining pointers for at.  These would change and would be restored
    * later.  */
   
   iperat_ps = iperat_ps + *nrow_p;
   *ipbrat_p++ = 0;
   while (iperat_p < iperat_ps) 
      {/* Note: we go till *nrow_p-1 */
	 *ipbrat_p++ = *iperat_p++;
      }
   iperat_ps = iperat_ps - *nrow_p;
   iperat_p = iperat_ps + 1;
   ipbrat_p = ipbrat_ps + 1;
   
   /* now  build the transpose structure and matrix */
   for (i = 1; i <= *ncol_p; i++) 
      {
	 jbeg = *ipbra_p++;
	 jend = *ipera_p++;
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *(ira_ps + jbeg);
	       iloc = *(ipbrat_ps + irow);
	       *(ipbrat_ps + irow) = ++iloc;
	       *(irat_ps + iloc) = i;
	       *(at_ps + iloc) = *(a_ps + jbeg);
	    }
      }
   ipbra_p = ipbra_ps;
   ipera_p = ipera_ps;
   
   /* reset the pointers */
   iperat_ps = iperat_ps + *nrow_p;
   *ipbrat_p++ = 1;
   while (iperat_p < iperat_ps) 
      {
	 /* Note: we go till *nrow_p-1 */
	 *ipbrat_p++ = *iperat_p++ + 1;
      }
   iperat_ps = iperat_ps - *nrow_p;
   iperat_p = iperat_ps + 1;
   ipbrat_p = ipbrat_ps + 1;
   
   return 0;
}

/********************************************************************/

/* SortColumnRealSparseMatrix.c, v1.0, 11/08/91, Sanjay Mehrotra */

/* The function in this file permutes each column of a matrix a[] in
 * increasing row index order.  The matrix a[] is assumed to be stored
 * column-wise.  The row indices corresponding to each column in a[] are
 * stored in ira[]. Begining and end of each column is stored in array
 * ipbra[] and ipera[]. All the array information is passed through pointers.
 * These pointers must point to the first element of the array.  The
 * subroutine assumes that the row/column indices run from 1..m/n.  As oppose
 * to 0..m-1, which is the standard C. All changes required are done
 * internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED.  * MATRIX IS STORED IN
 * PACKED FORM WITH COLUMNS IN INCREASING ORDER AND IPERA(i)+1=IPBRA(i+1).  */

extern int      CopyIntegerVectorxToy();
extern int      TransposeStructureSparseMatrix();

int             
SortColumnRealSparseMatrix(a_p, ira_p, ipbra_p, ipera_p, nrow_p, ncol_p)
     double         *a_p;
     int            *ipbra_p, *ipera_p, *ira_p;
     int            *ncol_p, *nrow_p;
{
   
   double         *a_q, *a_qe, *q_p, *q_q, *q_ps, *a_qs;
   int            *internal_space_p, *internal_space_ps, *pira_p, *ipbrat_p,
                  *iperat_p, *irat_p;
   int            *ira_q, *pira_q, *ira_qs, *jbeg, *jend, *ipera_ps, *ira_ps,
                  *irat_ps;
   int             i, iloc, irow, nza, internal_space_size, error_code;
   
   /* Allocated space for transpose structure  */
   
   if (FindNonzeroinaMatrix(ipbra_p, ipera_p, ncol_p, &nza)) 
      {
	 /* an internal error is detected */
	 error_code = 7;
	 
	 fprintf(stdout, "Error: SortColumnRealSparseMatrix: Error in\n");
	 fprintf(stdout, "  FindNonzeroinaMatrix.\n");
	 
	 return error_code;
      }
   internal_space_size = 2 * nza + 2 * *nrow_p;
   
   internal_space_p = NewInt(internal_space_size,
			     "internal_space_p in SortColumnRealSparseMatrix");
   
   internal_space_ps = internal_space_p;
   
   ipbrat_p = internal_space_p;
   iperat_p = ipbrat_p + *nrow_p;
   irat_p = iperat_p + *nrow_p;
   pira_p = irat_p + nza;
   /* the remaining nza space is for pira_p */
   
   /* allocate memory for arrays used inside this routine */
   q_p = NewDouble(*nrow_p, "q_p in SortcolumnRealSparseMatrix");
   
   q_ps = q_p;
   
   a_q = a_p - 1;
   q_q = q_p - 1;
   ira_q = ira_p - 1;
   pira_q = pira_p - 1;
   
   /* save the original matrix */
   if (CopyIntegerVectorxToy(ira_p, pira_p, &nza)) 
      {
	 error_code = 4;
	 
	 fprintf(stdout, "Error: SortColumnRealSparseMatrix:\n");
	 fprintf(stdout, "       in CopyIntegerVectorxToy.\n");
	 
	 Free((char *) internal_space_ps);
	 Free((char *) q_ps);
	 return error_code;
      }
   /* Two calls to the transpose finding matrix would permute the matrix. */
   if (TransposeStructureSparseMatrix
       (ipbra_p, ipera_p, ira_p, ipbrat_p, iperat_p, irat_p,
	nrow_p, ncol_p)) {
      error_code = 5;
      
      fprintf(stdout, "Error: SortColumnRealSparseMatrix:\n");
      fprintf(stdout, "       in TransposeStructureSparseMatrix.\n");
      
      Free((char *) internal_space_ps);
      Free((char *) q_ps);
      return error_code;
   }
   /* now find the transpose again */
   
   /* first shift a few pointers */
   ipera_ps = ipera_p - 1;
   ira_ps = ira_p - 1;
   irat_ps = irat_p - 1;
   
   if (CopyIntegerVectorxToy(ipbra_p, ipera_p, ncol_p)) 
      {
	 error_code = 6;
	 
	 fprintf(stdout, "Error: SortColumnRealSparseMatrix:\n");
	 fprintf(stdout, "       in CopyIntegerVectorxToy.\n");
	 
	 Free((char *) internal_space_ps);
	 Free((char *) q_ps);
	 return error_code;
      }
   /* take transpose of transpose */
   for (i = 1; i <= *nrow_p; i++) 
      {
	 jbeg = irat_ps + *ipbrat_p;
	 ipbrat_p++;
	 jend = irat_ps + *iperat_p;
	 iperat_p++;
	 
	 for (; jbeg <= jend; jbeg++) 
	    {
	       irow = *jbeg;
	       iloc = *(ipera_ps + irow);
	       *(ira_ps + iloc) = i;
	       *(ipera_ps + irow) = ++iloc;
	    }
      }
   
   /* shift pointers to give them their proper value */
   ipera_ps = ipera_ps + *ncol_p;
   while (ipera_p <= ipera_ps) 
      {
	 *ipera_p = *ipera_p - 1;
	 ipera_p++;
      }
   ipera_p = ipera_ps - *ncol_p + 1;
   
   /* go through each column and permute the data also */
  for (i = 1; i <= *ncol_p; i++) 
     {
	/* first map the element to array q */
	
	ira_qs = pira_q + *ipbra_p;
	a_qs = a_q + *ipbra_p;
	a_qe = a_q + *ipera_p;
	for (; a_qs <= a_qe; a_qs++) {
      *(q_q + *ira_qs++) = *a_qs;
	}
	/* now map the elements back to the matrix */
	ira_qs = ira_q + *ipbra_p;
	a_qs = a_q + *ipbra_p;
	a_qe = a_q + *ipera_p;

    for (; a_qs <= a_qe; a_qs++) {
      *a_qs = *(q_q + *ira_qs++);
    }

    /* advance the pointers */
    ipbra_p++;
    ipera_p++;
  }

  /* free the space allocated in this routine */
  Free((char *) internal_space_ps);
  Free((char *) q_ps);

  return 0;
}

/********************************************************************/

/* FindNonzeroinaMatrix.c, v1.0, 11/08/91, Sanjay Mehrotra */

/* Given the matrix data structure, this function finds the number of
 * nonzeros in a matrix.
 * 
 * Assumption: The matrix A is stored column-wise.
 * 
 * Input: Begining and end of each column is stored in array ipbra[] and
 * ipera[]. All the array information is passed through pointers.  These
 * pointers must point to the first element of the array.  We assume that the
 * row/column indices run from 1..m/n. As oppose to 0..m-1, which is the
 * standard C.
 * 
 * Output: Number of nonzeros in the matrix.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

int             
FindNonzeroinaMatrix(ipbra_p, ipera_p, ncol_p, nza_p)
     int            *ipbra_p, *ipera_p, *ncol_p;
     int            *nza_p;
{
   int             i;
   
   *nza_p = 0;
   for (i = 0; i < *ncol_p; i++) 
      {
	 *nza_p += *ipera_p++ - *ipbra_p++ + 1;
      }
   return 0;
}

/*********************************************************************/

/* DiagonalMatrixTimesDiagonalMatrix.c, stripped, 3/13/95.  */

int             
DiagonalMatrixTimesDiagonalMatrix(x_p, y_p, z_p, nz_p)
     double         *x_p, *y_p, *z_p;
     int            *nz_p;
{
   register int    i;

   for (i = 0; i < *nz_p; i++)
      *z_p++ = *x_p++ * *y_p++;
   return 0;
}

/*********************************************************************/

/* CopyRealVectorxToy.c, stripped 3/13/96.  */

int            
CopyRealVectorxToy(x_p, y_p, nz_p)
     double         *x_p, *y_p;
     int            *nz_p;
{
   register int    i;
   
   for (i = 0; i < *nz_p; i++)
      *y_p++ = *x_p++;
   return 0;
}

/*********************************************************************/

/* NegateRealVector.c, stripped 3/16/96.  */

int             
NegateRealVector(q_p, nzq_p)
     double         *q_p;
     int            *nzq_p;
{
   double         *q_ps;
   
   q_ps = q_p + *nzq_p;
   for (; q_p < q_ps; q_p++)
      *q_p = -*q_p;
   return 0;
}

/*********************************************************************/

/* rxtx8.c, stripped,   3/16/96. */

int             
NormTwoSquareRealDenseVector(x_p, nz_p, xtx)
     double         *x_p;
     double         *xtx;
     int            *nz_p;
{
   static double  *x_ps;
   
   *xtx = 0.0;
   x_ps = x_p + *nz_p;
   for (; x_p < x_ps; x_p++)
      *xtx += *x_p * *x_p;
   return 0;
}

/*********************************************************************/

/* ZeroRealDenseVector.c, stripped,   3/16/96. */

int             
ZeroRealDenseVector(q_p, nzq_p)
     double         *q_p;
     int            *nzq_p;
{
   static double  *q_ps;
   
   q_ps = q_p + *nzq_p;
   for (; q_p < q_ps; q_p++)
      *q_p = 0.0;
   return 0;
}

/*********************************************************************/

/* ZeroIntegerDenseVector.c, stripped, 3/16/96.  */

int             
ZeroIntegerDenseVector(q_p, nzq_p)
     int            *q_p;
     int            *nzq_p;
{
   int            *q_ps;
   
   q_ps = q_p + *nzq_p;
   for (; q_p < q_ps; q_p++)
      *q_p = 0;
   return 0;
}

/*********************************************************************/

/* CopyIntegerVectorxToy.c, stripped, , 3/16/96 */

int             
CopyIntegerVectorxToy(x_p, y_p, nz_p)
     int            *x_p, *y_p;
     int            *nz_p;
{
   register int    i;
   
   for (i = 0; i < *nz_p; i++)
      *y_p++ = *x_p++;
   return 0;
}

/*********************************************************************/

/* ZeroIntegerSparseVector.c, stripped 3/16/96 */

int             
ZeroIntegerSparseVector(iq_p, nzq_p, irq_p, nzrq_p)
     int            *iq_p, *irq_p, *nzq_p, *nzrq_p;
{
   int            *irq_ps, *iq_q;
   
   iq_q = iq_p - 1;
   irq_ps = irq_p + *nzrq_p;
   for (; irq_p < irq_ps; irq_p++)
      *(iq_q + *irq_p) = 0;
   return 0;
}



