/* input and output
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include <sys/types.h>		/* for determining file size */
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "memory.h"
#include "pre.h"


/*****************************************************************/
/* contents:                                                     */
/*****************************************************************/
/*
FILE *OpenInputFile(char *infile, int *filesize, Parameters *Inputs);
 
void PrintSolution(MPStype *MPS, solution *Solution, Parameters *Inputs, 
	      char *infilename);

void ComputeAndPrintObjectives(solution *Solution, LPtype *LP, 
			       double *primal, double *dual);

void ComputeAndPrintInfeasibilities(solution *Solution, LPtype *LP);

void PrintType(int type);
*/
/*****************************************************************/

FILE *OpenInputFile(infile, filesize, Inputs)
     char           *infile;
     int            *filesize;
     Parameters     *Inputs;
     
{
   FILE           *fp;
   struct stat     buf;
   char            filename[200];
   
   fp = NULL;

   /* Always look first in the current directory */
	 
   strcpy(filename, infile);
   fp = fopen(filename, "r");
   if (fp == NULL) 
     {
       strcat(filename, ".mps");
       fp = fopen(filename, "r");
     }

   if (fp == NULL) {
     printf("File '%s' not found in current directory.\n", filename);

   /* If the specs file contains an input directory, look there */
   
     if (Inputs->InputDirectory != NULL) 
       {
         printf("Looking in the specified input directory...\n");
	 strcpy(filename, Inputs->InputDirectory);
	 strcat(filename, infile);
	 fp = fopen(filename, "r");
	 if (fp == NULL) 
	    {
	       strcat(filename, ".mps");
	       fp = fopen(filename, "r");
	    }
	 if (fp == NULL) 
	    {
	       strcpy(filename, Inputs->InputDirectory);
	       strcat(filename, "/");
	       strcat(filename, infile);
	       fp = fopen(filename, "r");
	    }
	 if (fp == NULL) 
	    {
	       strcat(filename, ".mps");
	       fp = fopen(filename, "r");
	    }
      } 
   }

   if (fp != NULL) 
      {
         /* copy actual filename to input variable */
	 strcpy(infile, filename);
         printf("Reading input from file '%s'.\n", infile);
	 stat(filename, &buf);
	 *filesize = buf.st_size;
      } 
   else 
      {
	 printf("ERROR: couldn't find file '%s'.\n", infile);
	 if (stderr != stdout) 
	    {
	       fprintf(stderr, "ERROR: couldn't find file '%s'.\n", infile);
	    }
      }
   return fp;
}

/***************************************************************************/

/* Prints the solution from "Solution" to the output file "infile".out. Also
 * writes the history file "infile".log */

void PrintSolution(MPS, Solution, Inputs, infilename)
     MPStype        *MPS;
     solution       *Solution;
     Parameters     *Inputs;
     char           *infilename;
{
   char            outfilename[200], logfilename[200], *basename,
                  *suffix, statusTxt[40], rootfilename[200];
   FILE           *outfile, *logfile;
   int             row, col, entry, i, j, len;
   int             WriteSolution, WriteHistory;
   
   WriteSolution = Inputs->WriteSolution;
   WriteHistory = Inputs->WriteHistory;
   
   if (WriteSolution || WriteHistory) 
      {
	 /* get prefix of file names by stripping off leading path and trailing
	  * .mps */
	 basename = strrchr(infilename, '/');
	 if (basename == NULL)
	    basename = infilename;
	 else
	    basename++;
	 
	 suffix = strstr(basename, ".mps");
	 if (suffix != NULL) 
	    {
	       len = strlen(basename) - strlen(suffix);
	       strncpy(rootfilename, basename, len);
	       rootfilename[len] = '\0';
	       /* strcpy(basename, outfilename); */
	    }
	 else 
	    strcpy(rootfilename, basename);
      }
   /* if requested, open the output file */
   
   if (WriteSolution) 
      {
	 strcpy(outfilename, rootfilename);
	 strcat(outfilename, ".out");
	 outfile = fopen(outfilename, "w");
	 if (outfile == NULL) 
	    {
	       printf("Unable to open outfile '%s'; solution not printed.\n",
		      outfilename);
	    }
      }
   /* if requested, open the history file */
   
   if (WriteHistory) 
      {
	 strcpy(logfilename, rootfilename);
	 strcat(logfilename, ".log");
	 logfile = fopen(logfilename, "w");
	 if (logfile == NULL) 
	    {
	       printf("Unable to open logfile '%s'; history not printed.\n",
		      logfilename);
	    }
      }
   /* determine solution status */
   
   switch (Solution->Status) 
      {
      case OPTIMAL_SOL:
	 strcpy(statusTxt, "OPTIMAL");
	 break;
      case SUBOPTIMAL_SOL:
	 strcpy(statusTxt, "SUBOPTIMAL");
	 break;
      case INFEASIBLE_SOL:
	 strcpy(statusTxt, "INFEASIBLE");
	 break;
      case UNKNOWN_SOL:
	 strcpy(statusTxt, "UNKNOWN");
	 break;
      default:
	 strcpy(statusTxt, "UNKNOWN");
      }
   
   printf("\nProblem '%s' ", MPS->ProblemName);
   /*
     printf("In original problem, have %d variables, %d constraints\n",
     MPS->NumCols, MPS->NumRows);
     printf("Iterations =  %d\n", Solution->Iterations);
     */
   printf("terminated with %s status (code %d) after %d iterations\n\n", 
	  statusTxt, Solution->Status, Solution->Iterations);
   
   printf("Primal Objective = %13.8e\n", Solution->PrimalObjective);
   printf("Dual   Objective = %13.8e\n", Solution->DualObjective);
   
   printf("\nComplementarity          = %9.2e\n", Solution->Complementarity);
   printf("Relative Complementarity = %9.2e\n",
	  Solution->RelativeComplementarity);
   
   printf("\nRelative Infeasibilities:\n");
   printf("Primal = %9.3e,    ", Solution->PrimalInfeasibility);
   printf("Dual   = %9.3e.\n", Solution->DualInfeasibility);
   /* printf("\n");    printf("Number of Orderings (for factorization) =
    * %d.\n", Solution->Factorizations); */
   printf("\nRead Time       = %.*f seconds\n", DECIMAL_WIDTH, Solution->ReadTime);
   printf("Preprocess time = %.*f seconds\n", DECIMAL_WIDTH, Solution->PreprocessTime);
   printf("Solution time   = %.*f seconds\n", DECIMAL_WIDTH, Solution->SolutionTime);
   
   if (WriteSolution && (outfile != NULL)) 
      {
	 /* If infeasible, print nothing except a note */
	 if(Solution->Status == INFEASIBLE_SOL) 
	    {
	       fprintf(outfile, "INFEASIBLE status detected by PCx()");
	    } 
	 else 
	    {
	       fprintf(outfile, "Solution for '%s'\n", MPS->ProblemName);
	       fprintf(outfile, "Variables:\n");
	       
	       fprintf(outfile, " #   Label         Value       Reduced Cost");
	       fprintf(outfile, "    Lower Bound    Upper Bound\n");
	       for (col = 0; col < MPS->NumCols; col++) 
		  {
		     fprintf(outfile, "%3d %9s  %14.7e  %14.7e  ", col,
			     MPS->ColNames[col], Solution->x[col], 
			     Solution->DualLower[col]);
		     
		     /* print lower bound */
		     switch (MPS->BoundType[col]) 
			{
			case LOWER:
			case UPPERLOWER:
			case FIX:
			   fprintf(outfile, "%14.7e ", MPS->LowBound[col]);
			   break;
			case NORMAL:
			case UPPER:
			   fprintf(outfile, "%14.7e ", 0.0);
			   break;
			case MINFTY:
			case FREE:
			   fprintf(outfile, "-Infinity      ");
			   break;
			}
		     
		     /* print upper bound */
		     switch (MPS->BoundType[col]) 
			{
			case UPPER:
			case UPPERLOWER:
			case FIX:
			   fprintf(outfile, "%14.7e ", MPS->UpBound[col]);
			   break;
			case NORMAL:
			case LOWER:
			case FREE:
			   fprintf(outfile, " Infinity");
			   break;
			case MINFTY:
			   fprintf(outfile, "%14.7e ",  MPS->UpBound[col]);
			   break;
			}
		     fprintf(outfile, "\n");
		  }				/* end col loop */
	       
	       fprintf(outfile, "\nConstraints:\n");
	       fprintf(outfile, " #   Label  Type    ");
	       fprintf(outfile, "Activity        RHS          Dual");
	       fprintf(outfile, "       Lower Bound  Upper Bound\n");
	       
	       for (row = 0; row < MPS->NumRows; row++) 
		  {
		     fprintf(outfile, "%3d %8s %c %14.7e %13.6e %13.6e",
			     row, MPS->RowNames[row], MPS->RowType[row], 
			     Solution->Activity[row], MPS->b[row], 
			     Solution->pi[row]);
		     if (toupper(MPS->RowType[row]) == 'G') 
			{
			   if (MPS->Ranges[row] != 0.0)
			      fprintf(outfile, " %13.6e %13.6e\n",
				      MPS->b[row], MPS->b[row] + 
				      fabs(MPS->Ranges[row]));
			   else
			      fprintf(outfile, " %13.6e  INFINITY\n", 
				      MPS->b[row]);
			}
		     if (toupper(MPS->RowType[row]) == 'L') 
			{
			   if (MPS->Ranges[row] != 0.0)
			      fprintf(outfile, " %13.6e %13.6e\n",
				      MPS->b[row] - fabs(MPS->Ranges[row]), 
				      MPS->b[row]);
			   else
			      fprintf(outfile, "  -INFINITY    %13.6e\n", 
				      MPS->b[row]);
			}
		     if (toupper(MPS->RowType[row]) == 'E') 
			{
			   if (MPS->Ranges[row] > 0.0)
			      fprintf(outfile, " %13.6e %13.6e\n",
				      MPS->b[row], MPS->b[row] + 
				      MPS->Ranges[row]);
			   else if (MPS->Ranges[row] < 0.0)
			      fprintf(outfile, " %13.6e %13.6e\n",
				      MPS->b[row] - MPS->Ranges[row], 
				      MPS->b[row]);
			   else			/* == 0.0 */
			      fprintf(outfile, " %13.6e %13.6e\n", 
				      MPS->b[row], MPS->b[row]);
			}
		     if (toupper(MPS->RowType[row]) == 'N')
			fprintf(outfile, "\n");
		  }
	       if (outfile != stdout)
		  fclose(outfile);
	       
	    }
      }
   /* write History file, if requested */
   
   if (WriteHistory && (logfile != NULL)) 
      {
	 fprintf(logfile, 
		 "\n******** PCx version 1.1 (Nov 1997) ************\n\n");
	 fprintf(logfile, "Problem '%s' ", MPS->ProblemName);
	 fprintf(logfile, "terminated with %s status\n", statusTxt);
	 fprintf(logfile, "Iterations=%d, Termination Code=%d\n",
		 Solution->Iterations, Solution->Status);
	 fprintf(logfile, "\nMPS formulation has %d rows, %d columns\n",
		 MPS->NumRows,  MPS->NumCols);
	 
	 fprintf(logfile, "\nPARAMETER SUMMARY\n");
	 fprintf(logfile, "=================\n\n");
	 
	 fprintf(logfile, "Maximum number of iterations: %d\n",
		 Inputs->IterationLimit);
	 fprintf(logfile, "Tolerances: Opt=%8.2e  PriFeas=%8.2e",
		 Inputs->OptTol, Inputs->PriFeasTol);
	 fprintf(logfile, " DualFeas=%8.2e\n", Inputs->DualFeasTol);

         if (Inputs->HOCorrections) {
           fprintf(logfile, "Gondzio strategy selected: ");
           fprintf(logfile, " Maximum Gondzio corrections = %d\n", Inputs->MaxCorrections);
         } else 
             fprintf(logfile, "Mehrotra predictor-corrector strategy selected\n");
	 
	 if (Inputs->Refinement)
	    {
	       fprintf(logfile, "Iterative refinement performed during");
	       fprintf(logfile, " linear system solve\n");
	    }
	 else
	    fprintf(logfile, "NO iterative refinement\n");
	 
	 if (Inputs->Preprocessing) 
	    {
	       fprintf(logfile, "Presolving was performed:\n");
	       fprintf(logfile, "   Before Presolving:  %d rows, %d columns\n",
		       Solution->PriorRows, Solution->PriorColumns);
	       fprintf(logfile, "   After  Presolving:  %d rows, %d columns",
		       Solution->ReducedRows, Solution->ReducedColumns);
	       fprintf(logfile, "  (%d %s)\n", Solution->Passes, 
		       (Solution->Passes == 1)? "pass" : "passes");
	    }
	 else
	    fprintf(logfile, "NO presolving was performed\n");
	 
	 if (Inputs->Minimize)
	    fprintf(logfile, "MINIMIZE the objective\n");
	 else
	    fprintf(logfile, "MAXIMIZE the objective\n");
	 
	 if (Inputs->WriteSolution && (outfile != NULL))
	    fprintf(logfile, "Solution written to output file %s\n",
		    outfilename);
	 else
	    fprintf(logfile, "No solution file was written\n");
	 
	 if (Inputs->ObjectiveName != NULL)
	    fprintf(logfile, "Objective Name: %s\n", Inputs->ObjectiveName);
	 
	 if (Inputs->RHSName != NULL)
	    fprintf(logfile, "RHS Name: %s\n", Inputs->RHSName);
	 
	 if (Inputs->RangeName != NULL)
	    fprintf(logfile, "Range Name: %s\n", Inputs->RangeName);
	 
	 if (Inputs->BoundName != NULL)
	    fprintf(logfile, "Bound Name: %s\n", Inputs->BoundName);
	 
	 fprintf(logfile, "\n");
	 
	 fprintf(logfile, "\nFACTORIZATION SUMMARY\n");
	 fprintf(logfile, "=====================\n\n");

	 fprintf(logfile, "code used: %s\n", Solution->FactorizationCode);

	 if(Solution->FactorizationHistory->NumDenseCols > 0)
	    fprintf(logfile, "Dense columns extracted=%d\n", 
		    Solution->FactorizationHistory->NumDenseCols);
	 fprintf(logfile, "Nonzeros in L=%d;  Density of L=%f\n",
		 Solution->FactorizationHistory->Nonzeros,
		 Solution->FactorizationHistory->Density);
		 
	 fprintf(logfile, "\nITERATION SUMMARY\n");
	 fprintf(logfile, "=================\n\n");
	 
	 
	 if(Inputs->HOCorrections && Inputs->MaxCorrections > 0) 
	    {
	       fprintf(logfile, " Iter    Primal       Dual      ");
	       fprintf(logfile, "(PriInf  DualInf)  log(mu) corr  Merit\n");
	    } 
	 else 
	    {
	       fprintf(logfile, " Iter    Primal       Dual      ");
	       fprintf(logfile, "(PriInf  DualInf)  log(mu)   Merit\n");
	    }
	 
	 
	 for (i = 0; i <= Solution->Iterations; i++) 
	    {
	       fprintf(logfile, "%3d  %11.4e  %11.4e  (%7.1e %7.1e)   %6.2f  ",
		       i, Solution->IterationHistory[i].PrimalObjective,
		       Solution->IterationHistory[i].DualObjective,
		       Solution->IterationHistory[i].PriInf,
		       Solution->IterationHistory[i].DualInf,
		       Solution->IterationHistory[i].logmu);
	       
	       if (Inputs->HOCorrections && Inputs->MaxCorrections > 0)
		  fprintf(logfile, "%2d    %7.1e\n", 
			  Solution->IterationHistory[i].NumCorrections, 
			  Solution->IterationHistory[i].phi);
	       else
		  fprintf(logfile, "  %7.1e\n", 
			  Solution->IterationHistory[i].phi);
	    }
	 
	 fprintf(logfile, "\n %d iterations\n", Solution->Iterations);
	 fprintf(logfile, "\nTerminated with status %s (code %d)\n", 
		 statusTxt, Solution->Status);
	 if (Solution->RestoredIteration != -1)
	    fprintf(logfile, "\nSolution at iteration %d:\n", 
		    Solution->RestoredIteration);
	 fprintf(logfile, "\nPrimal Objective = %13.8e\n", 
		 Solution->PrimalObjective);
	 fprintf(logfile, "Dual   Objective = %13.8e\n", 
		 Solution->DualObjective);
	 
	 fprintf(logfile, "\nComplementarity          = %9.2e\n", 
		 Solution->Complementarity);
	 fprintf(logfile, "Relative Complementarity = %9.2e\n",
		 Solution->RelativeComplementarity);
	 
	 fprintf(logfile, "\nRelative Infeasibilities:\n");
	 fprintf(logfile, "Primal = %9.3e,    ", 
		 Solution->PrimalInfeasibility);
	 fprintf(logfile, "Dual   = %9.3e.\n", Solution->DualInfeasibility);
	 


	 fprintf(logfile, "\nTIME SUMMARY\n");
	 fprintf(logfile, "============\n\n");

	 fprintf(logfile, "Time to read input file: %.*f sec\n",
		 DECIMAL_WIDTH, Solution->ReadTime);
	 if (Inputs->Preprocessing) 
	    {
	       fprintf(logfile, "Time to presolve       : %.*f sec\n\n",
		       DECIMAL_WIDTH, Solution->PreprocessTime);
	    }

#ifdef TIMING_PROFILE
	 fprintf(logfile, "InitTime               : %.*f sec ",
		 DECIMAL_WIDTH, Solution->InitTime);
	 fprintf(logfile, "(= %3.1f %% of total)\n",
		 Solution->InitTime/Solution->SolutionTime * 100);
	 
	 fprintf(logfile, "LoopTime               : %.*f sec",
		 DECIMAL_WIDTH, Solution->LoopTime);
	 fprintf(logfile, " (= %2.1f %% of total)\n",
		 Solution->LoopTime/Solution->SolutionTime * 100);
	 fprintf(logfile, "   FormADATTime               : %.*f sec",
		 DECIMAL_WIDTH, Solution->FormADATtime);
	 fprintf(logfile, " (= %4.1f %% of loop)\n",
		 Solution->FormADATtime/Solution->LoopTime * 100);
	 fprintf(logfile, "   PredictorTime              : %.*f sec",
		 DECIMAL_WIDTH, Solution->PredictorTime);
	 fprintf(logfile, " (= %4.1f %% of loop)\n",
		 Solution->PredictorTime/Solution->LoopTime * 100);
	 fprintf(logfile, "   CorrectorTime              : %.*f sec",
		 DECIMAL_WIDTH, Solution->CorrectorTime);
	 fprintf(logfile, " (= %4.1f %% of loop)\n",
		 Solution->CorrectorTime/Solution->LoopTime * 100);
	 fprintf(logfile, "   Factorization              : %.*f sec",
		 DECIMAL_WIDTH, Solution->FactorizationTime);
	 fprintf(logfile, " (= %4.1f %% of loop)\n",
		 Solution->FactorizationTime/Solution->LoopTime * 100);
#endif
	 fprintf(logfile, "Time to solve          : %.*f sec\n\n",
		 DECIMAL_WIDTH, Solution->SolutionTime);
#ifdef TIMING_PROFILE
	 fprintf(logfile, "Average time spent\n");
	 fprintf(logfile, "   for one num. factorization: %f sec\n",
		 Solution->FactorizationTime / (Solution->Iterations));
	 fprintf(logfile, "   for one SolveADAT         : %f sec\n",
		 Solution->SolveADATTime / (2*Solution->Iterations));
#endif
	 
      }
	 if ((WriteHistory) && (logfile != stdout))
	    fclose(logfile);
}

/********************************************************************/

void ComputeAndPrintObjectives(Solution, LP, primal, dual)
     solution       *Solution;
     LPtype         *LP;
     double         *primal, *dual;
{
   
   int             row, col, i;
   
   for (col = 0, *primal = 0.0; col < LP->Cols; col++)
      *primal += (LP->c[col] * Solution->x[col]);
   
   for (row = 0, *dual = 0.0; row < Solution->Rows; row++)
      *dual += (Solution->pi[row] * LP->b[row]);
   
   for (col = 0; col < LP->Cols; col++)
      if (LP->VarType[col] == UPPER)
	 *dual -= (Solution->DualUpper[col] * LP->UpBound[col]);
   
   printf(" Primal Objective = %f\n", *primal + LP->cshift);
   printf(" Dual   Objective = %f\n", *dual + LP->cshift);
}

/********************************************************************/

void ComputeAndPrintInfeasibilities(Solution, LP)
     solution       *Solution;
     LPtype         *LP;
{
   int             row, col, i, k, SparseSaxpy(), SparseSaxpyT();
   double          temp, *PrimalInf, *DualInf, complementarity, 
                   cost2norm, rhs2norm, TwoNorm2(), 
                   rel_primal, rel_dual, primal, dual;
   
   /* compute primal infeasibility */
   
   PrimalInf = NewDouble(LP->Rows, "PrimalInf in PrintInfeasibilities()");
   for (row = 0; row < LP->Rows; row++)
      PrimalInf[row] = -LP->b[row];
   SparseSaxpy(LP->A, Solution->x, PrimalInf);
   primal = TwoNorm2(PrimalInf, &(LP->Rows));
   Free((char *) PrimalInf);
   
   /* compute bound infeasibility */
   for (col = 0; col < LP->Cols; col++)
      if (LP->VarType[col] == UPPER) 
	 {
	    temp = Solution->x[col] - LP->UpBound[col];
	    if (temp > 0.0)
	       primal += temp * temp;
	 }
   /* actually, in the non-preprocessed problem there MAY BE free variables
    * --- need to test for this */
   for (col = 0; col < LP->Cols; col++)
      if (Solution->x[col] < 0.0 && LP->VarType[col] != FREE) 
	 {
	    /*
	      printf("Solution->x[%d] < 0.0: %f ", col, Solution->x[col]);
	      PrintType(LP->VarType[col]);
	      printf(" variable\n");
	      */
	    primal += Solution->x[col] * Solution->x[col];
	 }
   
   /* compute dual infeasibility */
   DualInf = NewDouble(LP->Cols, "DualInf in PrintInfeasibilities()");
   for (col = 0; col < LP->Cols; col++)
      DualInf[col] = Solution->DualLower[col] -
	 Solution->DualUpper[col] - LP->c[col];
   
   SparseSaxpyT(LP->A, Solution->pi, DualInf);
   dual = TwoNorm2(DualInf, &(LP->Cols));
   /*
     for (col = 0; col < LP->Cols; col++)
     if (fabs(DualInf[col]) > 1.0e-4) {
     if (LP->VarType[col] == NORMAL)
     printf(" Normal variable");
     else if (LP->VarType[col] == FREE)
     printf(" Free variable");
     else if (LP->VarType[col] == UPPER)
     printf(" Upper-bounded variable");
     printf(" %d dual infeasibility = %f\n", col, DualInf[col]);
     }
     */
   for (col = 0; col < LP->Cols; col++) 
      {
	 if (Solution->DualUpper[col] < 0.0) 
	    {
      /*
	printf("Solution->DualUpper[%d] < 0.0: %f\n",
	col, Solution->DualUpper[col]);
	*/
	       dual += Solution->DualUpper[col] * Solution->DualUpper[col];
	    }
	 if (Solution->DualLower[col] < 0.0) 
	    {
	       /*
		 printf("Solution->DualLower[%d] < 0.0: %f\n",
		 col, Solution->DualLower[col]);
		 */
	       dual += Solution->DualLower[col] * Solution->DualLower[col];
	    }
      }
   Free((char *) DualInf);
   
   /* compute complementarity */
   
   complementarity = 0.0;
   for (col = 0; col < LP->Cols; col++) 
      {
	 temp = Solution->x[col] * Solution->DualLower[col];
	 complementarity += temp;
	 
	 if (LP->VarType[col] == UPPER) 
	    {
	       temp = (LP->UpBound[col] - Solution->x[col]) * 
		  Solution->DualUpper[col];
	       complementarity += temp;
	    }
      }
   primal = sqrt(primal);
   dual = sqrt(dual);
   
   /* compute norm of cost vector and rhs */
   cost2norm = sqrt(TwoNorm2(LP->c, &(LP->Cols)));
   rhs2norm  = TwoNorm2(LP->b, &(LP->Rows));
   for(i=0; i<LP->NumberBounds; i++) 
      {
	 k = LP->BoundIndex[i];
	 rhs2norm += LP->UpBound[k]*LP->UpBound[k];
      }
   rhs2norm = sqrt(rhs2norm);
   
   rel_dual = dual / (1.0 + cost2norm);
   rel_primal = primal / (1.0 + rhs2norm);
   Solution->PrimalInfeasibility = rel_primal;
   Solution->DualInfeasibility = rel_dual;
   Solution->Complementarity = complementarity;
   Solution->RelativeComplementarity =
      fabs(complementarity) / (1.0 + fabs(Solution->PrimalObjective));
}

void PrintType(type)
     int             type;
{
   switch (type) 
      {
      case NORMAL:
	 printf("NORMAL");
	 break;
      case FREE:
	 printf("FREE");
	 break;
      case UPPER:
	 printf("UPPER");
	 break;
      case LOWER:
	 printf("LOWER");
	 break;
      case UPPERLOWER:
	 printf("UPPERLOWER");
	 break;
      case FIX:
	 printf("FIX");
	 break;
      case MINFTY:
	 printf("MINFTY");
	 break;
      }
}

