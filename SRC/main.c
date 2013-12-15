/* main program
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 *
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */


#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "memory.h"
#include "pre.h"

#ifdef DEBUG
#include <errno.h>
#endif

char    infile[200];
char    soutfile[200];
char    sinfile[200];
char    P_name[100];
int     Eta, trocou,Etamax,aumenta,marca, param_eta;
FILE    *eigout;

void usage(char *argv[]) {
  printf("Usage:\n");
  printf("\t%s mpsfile\n", argv[0]);
  printf("\t%s -s specification_file mpsfile\n", argv[0]);
}

main(argc, argv)
     int             argc;
     char           *argv[];
{

   FILE           *fp, *outfile, *OpenInputFile();
   int             Preprocess(), Postprocess(), passes, PCx();
   int             CheckParameters();
   LPtype         *LP, *ReducedLP, *Convert_MPS_LP();
   MPStype        *MPS, *ReadMPS();
   Parameters     *Inputs, *NewParameters();
   ChangeStack    *Record;
   MPSchanges     *Changes;
   solution       *Solution, *NewSolution();
   void            SplitFreeVars(), UnSplitFreeVars(), PrintSolution();
   int             filesize, i;
   int             status;
   double          readtime = 0.0, pretime = 0.0, SysTime;
   double          UserTime, OldSysTime, OldUserTime;

   extern        char            infile[200];
   extern        char            soutfile[200];
   extern        char            sinfile[200];
   char                          eigfile[200];
   extern        char            P_name[100];
   extern        FILE            *eigout;

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

   marca=0;
   aumenta=0;

   /* Create the parameter data structure, insert the parameter values into
    * it, and check their validity.  */

   if (argc == 2) {
     strcpy(infile, argv[1]);
     strcpy(soutfile, argv[1]);
     sinfile[0] = '\0';
   }
   else {
     if (argc == 4)
       if (strcmp("-s", argv[1]) == 0) {
         strcpy(sinfile, argv[2]);
         strcpy(infile, argv[3]);
         strcpy(soutfile, argv[3]);
       }
     else {
       usage(argv);
       exit(INVOCATION_ERROR);
     }
   }

   strcpy(infile, argv[1]);
   strcpy(P_name, argv[1]);

   param_eta = 10;

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

   strcpy(eigfile, Inputs->OutputDirectory);
   // Test if output directory exist. If not, create it.
   DIR *out = NULL;
   out = opendir(eigfile);
   if (out == NULL) {
     mkdir(eigfile, S_IRWXU);  // The directory must have rwx permission
   }
   else {
     closedir(out); // Avoid corrupt the directory
   }
   // Check if need to add '/' in the output file name
   if (eigfile[strlen(eigfile)] != '/') {
     strcat(eigfile, "/");
   }
   strcat(eigfile, argv[1]);  // Need to use argv[1] because `infile` can start with `.`
   strcat(eigfile, ".out");
   eigout = fopen(eigfile,"w");
   if (eigout == NULL) {
#ifdef DEBUG
     printf("Can't open `%s` because of %d\nTo find what that means: $ grep %d /usr/include/*/errno*.h\n", eigfile, errno, errno);
#else
     printf("Can't open `%s`.\n", eigfile);
#endif
     exit(1);
   }
   fprintf(eigout, "Problem %s \n", infile);
   fflush(eigout);

  // Create a file for saving POS vetor, name= POS_file
   /*strcpy(POS_file,"./OUT/Pos_file");
   OUT_POS=fopen(POS_file,"w");
   if (OUT_POS == NULL)
     {
       printf("Nao foi possivel abrir o arquivo %s \n",POS_file);
       error(1);
     }*/

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
   //Saving b vetor
   //fprintf_vetor(ReducedLP-> b, ReducedLP->Rows, "_b");
   //Saving c vetor
   //fprintf_vetor(ReducedLP-> c, ReducedLP->Cols, "_c");

   status = PCx(ReducedLP, Solution, Inputs);

   if (status != 0)
      {
	 printf("Error in PCx(). Exiting with code %d\n", status);
	 exit(status);
      }
   GetTime(&UserTime, &SysTime);
   Solution->SolutionTime = UserTime - OldUserTime + SysTime - OldSysTime;

   fprintf(eigout, " Imprimindo o tempo de solucao \n ");
   fprintf(eigout, "\nRead Time       = %.2f seconds\n", Solution->ReadTime);
   fprintf(eigout, "Preprocess time = %.2f seconds\n", Solution->PreprocessTime);
   fprintf(eigout, "Solution time   = %.2f seconds\n", Solution->SolutionTime);

   printf( " Imprimindo o tempo de solucao \n ");
   printf("\nRead Time       = %.2f seconds\n", Solution->ReadTime);
   printf("Preprocess time = %.2f seconds\n", Solution->PreprocessTime);
   printf("Solution time   = %.2f seconds\n", Solution->SolutionTime);

   fflush(eigout);

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

   fclose(eigout);

   status = Solution->Status; /* For return code. */
   FreeSolution(Solution);
   DeleteMPS(MPS);
   FreeParameters(Inputs);

   /* TrDump(stdout); */
   exit(status);
}

