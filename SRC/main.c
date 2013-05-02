/* main program
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "memory.h"
#include "pre.h"

char    infile[200];
char    outfile[200];
char    sinfile[200];
char    P_name[100];

void usage(char *argv[]) {
  printf("Usage:\n");
  printf("\t%s mpsfile\n", argv[0]);
  printf("\t%s -s specification_file mpsfile\n", argv[0]);
}

main(argc, argv)
     int             argc;
     char           *argv[];
{
   
   FILE           *fp, *OpenInputFile();
   int             Preprocess(), Postprocess(), passes, PCx();
   int             CheckParameters();
   LPtype         *LP, *ReducedLP, *Convert_MPS_LP();
   MPStype        *MPS, *ReadMPS();
   Parameters     *Inputs, *NewParameters();
   ChangeStack    *Record;
   MPSchanges     *Changes;
   solution       *Solution, *NewSolution();
   void            SplitFreeVars(), UnSplitFreeVars(), PrintSolution();
   int             filesize;
   int             status;
   double          readtime = 0.0, pretime = 0.0, SysTime; 
   double          UserTime, OldSysTime, OldUserTime;
   
   extern        char            infile[200];
   extern        char            outfile[200];
   extern        char            sinfile[200];
   
  /********************************************************************
   *                                                                  *
   * Problems in the "LPtype" data structure have the following form: *
   *                                                                  *
   *    Primal:    min  c^T x                                         *
   * (pi)          s.t. A x   = b          x - primal variable        *
   * (s, r)               0 <= x <= u      "upper"  variables         *
   * (s)                  0 <= x           "normal" variables         *
   *                           x free      "free"   variables         *
   *                                                                  *
   ********************************************************************
   *                                                                  *
   *    Dual:      max  b^T pi - r^T u        pi - dual variable      *
   * ("upper"  x)  s.t. A^T pi + s - r = c     r - dual bound slack   *
   * ("normal" x)  s.t. A^T pi + s     = c     s - dual slack         *
   * ("free"   x)  s.t. A^T pi         = c                            *
   *                    pi free                                       *
   *                    (r,s) >= 0                                    *
   *                                                                  *
   ********************************************************************
   *                                                                  *
   * Before calling PCx, the free variables are split                 *
   * into positive and negative parts, to make two extra "normal"     *
   * variables.                                                       *
   *                                                                  *
   ********************************************************************/

   printf("\n******** PCx version 1.1 (Nov 1997) ************\n\n");

   /* SetFPTrap(32); */
   
   if (argc < 2) 
      {
         usage(argv);
	 exit(INVOCATION_ERROR);
      }
   /* Create the parameter data structure, insert the parameter values into
    * it, and check their validity.  */

   if (argc == 2) {
     strcpy(infile, argv[1]);
     strcpy(outfile, argv[1]);
     sinfile[0] = '\0';
   }
   else {
     if (argc == 4)
       if (strcmp("-s", argv[1]) == 0) {
         strcpy(sinfile, argv[2]);
         strcpy(infile, argv[3]);
         strcpy(outfile, argv[3]);
       }
     else {
       usage(argv);
       exit(INVOCATION_ERROR);
     }
   }
   

   /* load the default parameters */
   Inputs = NewParameters();
   
   /* read modified parameters (if any) from specs file */
   ParseSpecsFile(Inputs, infile, sinfile);
   
   /* check for errors */
   if (CheckParameters(Inputs)) 
      {
	 printf("Error return from CheckParameters\n");
	 exit(SPECS_ERROR);
      }
   /* open the input file (currently an MPS file) */
   
   fp = OpenInputFile(infile, &filesize, Inputs);
   if (fp == NULL) 
      {
	 printf("Error return from OpenInputFile\n");
	 fflush(stdout);
	 exit(INPUT_ERROR);
      }
   /* Create the MPS data structure and read the MPS input file into this
    * structure */
   
   MPS = ReadMPS(fp, filesize, Inputs, &readtime);

   if(MPS == NULL) 
      {
	 printf("\nError return from ReadMPS\n");
	 exit(INPUT_ERROR);
      }

   
   printf("\nMPS formulation: %d rows, %d columns\n", 
	  MPS->NumRows, MPS->NumCols);

   /* convert the MPStype data into LPtype data, which allows only equality
    * constraints and three kinds of x components: free, nonnegative, or
    * bounded as in 0 <= x_i <= u_i */
   
   LP = Convert_MPS_LP(MPS, &Changes);
   
   printf("LP  formulation: %d rows, %d columns\n", LP->Rows, LP->Cols);
   
   /* run the preprocessor, which converts LP to ReducedLP, and split the free
    * variables into positive and negative parts.  Double-check that there are
    * no free variables in ReducedLP.  */
   
   GetTime(&OldUserTime, &OldSysTime);
   
   if (Inputs->Preprocessing) 
      {
	 passes = Preprocess(LP, &ReducedLP, &Record, Inputs);
	 if (passes < 0) 
	    {
	       printf("Error return from Preprocessor:");
	       exit(PRESOLVE_ERROR);
	    }
      } 
   else
      ReducedLP = LP;
   
   /* PrintLP(ReducedLP); */
   
   if (Inputs->Scaling)
      ScaleLP(ReducedLP, Inputs);
   
   SplitFreeVars(ReducedLP);
   
   GetTime(&UserTime, &SysTime);
   fflush(stdout);
   pretime = UserTime - OldUserTime + SysTime - OldSysTime;
   
   /* set up a solution data structure, and put the timing information that
    * we've gathered so far into it.  */
   
   Solution = NewSolution(ReducedLP->Rows, ReducedLP->Cols,
			  Inputs->IterationLimit);
   Solution->ReadTime = readtime;
   Solution->PreprocessTime = pretime;
   
   /* solve the problem, keeping track of CPU times.  */
   
   GetTime(&OldUserTime, &OldSysTime);
   
   status = PCx(ReducedLP, Solution, Inputs);
 
   if (status != 0) 
      {
	 printf("Error in PCx(). Exiting with code %d\n", status);
	 exit(status);
      }
   GetTime(&UserTime, &SysTime);
   Solution->SolutionTime = UserTime - OldUserTime + SysTime - OldSysTime;
   
   /* recover the free variable values by combining their positive and
    * negative parts, and undo the effects of the preprocessor to recover the
    * solution to the original problem LP.  */

   UnSplitFreeVars(ReducedLP, Solution);
   
   if (Inputs->Scaling)
      UnscaleLP(ReducedLP, Solution);

   if (Inputs->Preprocessing) 
      if (Postprocess(LP, &Record, &Solution) < 0) 
	 exit(PRESOLVE_ERROR);

   /* Check the infeasibilities for the point obtained */
   ComputeAndPrintInfeasibilities(Solution, LP);

   /* Express the solution of LP in terms of the original MPS formulation */
   Solution = MPSsolution(LP, MPS, Solution, Changes, Inputs);
   DeleteChanges(Changes);
   if (ReducedLP != LP)
      DeleteLP(ReducedLP);
   DeleteLP(LP);
   /* Output the results */
   
   PrintSolution(MPS, Solution, Inputs, &outfile);

   status = Solution->Status; /* For return code. */
   FreeSolution(Solution);
   DeleteMPS(MPS);
   FreeParameters(Inputs);
   
   /* TrDump(stdout); */
   exit(status);
}

