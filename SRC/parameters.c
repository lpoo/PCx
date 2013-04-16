/* allocate parameter data structure, assign defaults, read
 * specificatons file
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "memory.h"
#include "string.h"


/*=====================================================================*/
/* Define strcasecmp for Windows 32 (not supported to the best of my   */
/* knowledge (GDI 11/08/97)                                            */
/*                                                                     */
/* Added by: G. D. Irisarri                                            */
/*           Open Access Technology International                      */
/*           oatigi@winternet.com                                      */
/*                                                                     */
/*           November 1997                                             */
/*=====================================================================*/

#ifdef WIN32

static int strcasecmp(const char *, const char *);

int
strcasecmp(const char *s, const char *t)
{
	for ( ; tolower(*s) == tolower(*t); s++, t++)
		if (*s == '\0')
			return 0;
	return *s - *t;
}

#endif


/***************************************************************************/
/*   Contents:                                                             */
/***************************************************************************/
/*
Parameters  *NewParameters();

void	     FreeParameters(Parameters *Inputs);

int          CheckParameters(Parameters *ptr);

int          ParseSpecsFile(Parameters *parameters, char *infile char *sinfile);
*/
/***************************************************************************/

#define SAME 0

int             LARGEST_INT4 = 1073741824, 
                FLAG_OPERATION_COUNT = 1, 
                DEFAULT_MINIMUM_DEGREE_GRAPH_FACTOR = 3, 
                DEFAULT_MINIMUM_DEGREE_CHOLESKY_FACTOR = 5;

double          DEFAULT_CHOLESKY_PIVOT_TOLERANCE = 1.0e-14, 
                DEFAULT_1x1_PIVOT_TOLERANCE = 1.0 / (16.0 * 16.0 * 16.0 * 
						     16.0 * 16.0 * 16.0 * 
						     16.0 * 16.0 * 16.0 * 
						     16.0);

/***************************************************************************/

/* Allocate space for parameter array and set defaults */

Parameters *
NewParameters()
{
   Parameters     *ptr;
   
   ptr = (Parameters *) Malloc(sizeof(Parameters), "NewParameters");
   
   ptr->IterationLimit = 100;
   ptr->OptTol         = 1.0e-10;
   ptr->PriFeasTol     = 1.0e-8;
   ptr->DualFeasTol    = 1.0e-8;
   ptr->AlphaScale     = 0.9;
   ptr->Diagnostics    = 1;	
   ptr->ReportingLevel = 2;	
   ptr->Refinement     = NO;	/* don't do iterative refinement */
   ptr->Preprocessing  = YES;	/* do preprocessing */
   ptr->Scaling        = YES;   /* do scaling */
   ptr->HOCorrections  = YES;   /* Use Gondzio higher-order correctors */
   ptr->MaxCorrections = 0;     /* Code decides maximum number of Gondzio 
				   corrections */
   ptr->Minimize       = YES;	/* minimize, don't maximize */
   ptr->InputDirectory = NULL;	/* there's not particular place to look for
				 * the MPS input files */
   ptr->WriteSolution  = YES;	/* yes, write the solution to the file
				 * probname.out */
   ptr->WriteHistory   = YES;	/* no, don't write any history info to the
				 * file probname.log */
   ptr->ObjectiveName  = NULL;
   ptr->RHSName        = NULL;
   ptr->RangeName      = NULL;
   ptr->BoundName      = NULL;
   
   ptr->CacheSize      = 16;
   ptr->UnrollingLevel = 4;
   
   ptr->CenterExponent = 3.0;

   ptr->OrderAlg = 1;
   
   return ptr;
}


/***************************************************************************/

void		
FreeParameters(Inputs)
     Parameters		*Inputs;
{
   Free ((char*) Inputs->InputDirectory);
   Free ((char*) Inputs->RHSName);
   Free ((char*) Inputs->ObjectiveName);
   Free ((char*) Inputs->RangeName);
   Free ((char*) Inputs->BoundName);
   Free ((char*) Inputs);
}

/***************************************************************************/


/* check parameter list for outrageous values */

int             
CheckParameters(ptr)
     Parameters     *ptr;
{
   int             temp;
   
   temp = 0;
   
   if (ptr->IterationLimit <= 0 || ptr->IterationLimit > 500) 
      {
	 printf(" PARAMETER ERROR: Out-of-Range IterationLimit: %d\n",
		ptr->IterationLimit);
	 temp += 1;
      }
   if (ptr->OptTol <= 1.e-12 || ptr->OptTol > 1) 
      {
	 printf(" PARAMETER ERROR: Out-of-Range OptTol: %12.4e\n", 
		ptr->OptTol);
	 temp += 1;
      }
   if (ptr->PriFeasTol <= 1.e-12 || ptr->PriFeasTol > 1.e-1) 
      {
	 printf(" PARAMETER ERROR: Out-of-Range PriFeasTol: %12.4e\n", 
		ptr->PriFeasTol);
	 temp += 1;
      }
  if (ptr->DualFeasTol <= 1.e-12 || ptr->DualFeasTol > 1.e-1) 
     {
	printf(" PARAMETER ERROR: Out-of-Range DualFeasTol: %12.4e\n", 
	       ptr->DualFeasTol);
	temp += 1;
     }
  if (ptr->AlphaScale <= 0.0 || ptr->AlphaScale > 0.9999999) 
     {
	printf(" PARAMETER ERROR: Out-of-Range AlphaScale: %f\n",
	       ptr->AlphaScale);
	temp += 1;
     }
  if (ptr->UnrollingLevel != 1 && ptr->UnrollingLevel != 2 &&
      ptr->UnrollingLevel != 4 && ptr->UnrollingLevel != 8) 
     {
	printf(" PARAMETER ERROR: Out-of-Range UnrollingLevel: %d\n", 
	       ptr->UnrollingLevel);
	printf("   Choose 1, 2, 4, or 8\n");
	temp += 1;
     }
  if (ptr->CacheSize < 0 || ptr->CacheSize > 2048) 
     {
	printf(" PARAMETER ERROR: Out-of-Range CacheSize: %d\n", 
	       ptr->CacheSize);
	printf("   Choose CacheSize in range (0 - 2048)\n");
	temp += 1;
     }

  if (ptr->CenterExponent < 1.0 || ptr->CenterExponent > 4.0) 
     {
	printf(" PARAMETER ERROR: Out-of-Range CenterExponent: %f\n",
	       ptr->CenterExponent);
	printf("   Choose CenterExponent in range (1.0 - 4.0)\n");
	temp += 1;
     }
  
  if (ptr->MaxCorrections < 0 || ptr->MaxCorrections > 10) 
     {
	printf(" PARAMETER ERROR: Out-of-Range MaxCorrections: %d\n",
	       ptr->MaxCorrections);
	printf("   Choose MaxCorrections in range (0 - 10)\n");
	temp += 1;
     }
  
  return temp;
}

/***************************************************************************/

int 
ParseSpecsFile(parameters, infile, sinfile)
     Parameters     *parameters;
     char           *infile;
     char           *sinfile;
{
   int             match, key;
   char            line[200], rootfilename[200], filename[200];
   char           *first, *second;
   FILE           *fp;
   char           *basename, *suffix; 
   int             len;
   
   int             NumKeys = 24;
   static char    *KeyWord[] = {"max", "min", "solution", 
				"objectivename", "rhsname", "rangename",
				"boundname", "history", "presolve", 
				"preprocess", "inputdirectory", "opttol", 
				"prifeastol", "dualfeastol", "cachesize", 
				"unrollinglevel", "iterationlimit", 
				"centerexp", "refinement", "stepfactor", 
				"scaling", "hocorrections", 
				"maxcorrections", "orderalg"};
   
   /* Try finding specs file under some different names. Give priority to
    * filenames that include this particular problem name */
   
   /* get prefix of file names by stripping off leading path and trailing
    * .mps */
   
   basename = strrchr(infile, '/');
   if (basename == NULL)
      basename = infile;
   else
      basename++;
   
   suffix = strstr(basename, ".mps");
   if (suffix != NULL) 
      {
	 len = strlen(basename) - strlen(suffix);
	 strncpy(rootfilename, basename, len);
	 rootfilename[len] = '\0';
	 strcpy(basename, rootfilename);
      } 
   else 
      strcpy(rootfilename, basename);
   
   strcpy(filename, rootfilename);
   strcat(filename, ".spc");
   fp = fopen(filename, "r");
   if (fp == NULL) 
      {
	 strcpy(filename, rootfilename);
	 strcat(filename, ".specs");
	 fp = fopen(filename, "r");
      }
   if (fp == NULL && sinfile[0] != '\0')
      {
         strcpy(filename, sinfile);
	 fp = fopen(filename, "r");
      }
   if (fp == NULL) 
      {
	 strcpy(filename, "spc");
	 fp = fopen(filename, "r");
      }
   if (fp == NULL) 
      {
	 strcpy(filename, "specs");
	 fp = fopen(filename, "r");
      }
   if (fp == NULL) 
      {
	 strcpy(filename, "PCx.specs");
	 fp = fopen(filename, "r");
      }
/************** have to include this part for the sake of the DOS version ***/
   if (fp == NULL) 
      {
	 strcpy(filename, "PCx.spc");
	 fp = fopen(filename, "r");
      }
/****************************************************************************/

   /* give up at this point */
   if (fp == NULL)
    return 0;
   
   printf("Reading problem specs from file '%s'\n", filename);
   
   /* process lines in file */
   
   do 
      {
	 fgets(line, 200, fp);
      } while (line[0] == '#' && !feof(fp));
   
   while (!feof(fp)) 
      {
	 first = strtok(line, " \t\n");
	 second = strtok(NULL, " \t\n");
	 
	 if (first == NULL)
	    /* encountered blank line */
	    match = -2;
	 else {
	    for (key = 0, match = -1; key < NumKeys; key++) 
	       {
		  if (strcasecmp(KeyWord[key], first) == SAME) 
		     {
			match = key;
			break;
		     }
	       }
	 }
	 
	 switch (match) 
	    {
	    case -2:
	       break;
	    case -1:
	       printf("Unknown keyword '%s' in specs file.\n", first);
	       printf("Ignoring.\n");
	       break;
	    case 0:			/* max */
	       parameters->Minimize = NO;
	       printf("  Maximizing the objective.\n");
	       break;
	    case 1:			/* min */
	       parameters->Minimize = YES;
	       printf("  Minimizing the objective.\n");
	       break;
	    case 2:			/* write solution */
	       if (second == NULL) 
		  {
		     parameters->WriteSolution = YES;
		     printf("  Writing solution file.\n");
		  } 
	       else if ((strcasecmp("no", second) == SAME) ||
			(strcasecmp("n", second) == SAME)) 
		  {
		     parameters->WriteSolution = NO;
		     printf("  Not writing solution file.\n");
		  } 
	       else 
		  {
		     parameters->WriteSolution = YES;
		     printf("  Writing solution file.\n");
		  }
	       break;
	    case 3:			/* objectivename */
	       parameters->ObjectiveName = StrDup(second, "ObjectiveName");
	       printf("  Using objective '%s'.\n", second);
	       break;
	    case 4:			/* rhsname */
	       parameters->RHSName = StrDup(second, "RHSName");
	       printf("  Using RHS '%s'.\n", second);
	       break;
	    case 5:			/* rangename */
	       parameters->RangeName = StrDup(second, "RangeName");
	       printf("  Using Range '%s'.\n", second);
	       break;
	    case 6:			/* boundname */
	       parameters->BoundName = StrDup(second, "BoundName");
	       printf("  Using Bound '%s'.\n", second);
	       break;
	    case 7:			/* history */
	       if (second == NULL) 
		  {
		     parameters->WriteHistory = YES;
		     printf("  Writing history file.\n");
		  } 
	       else if ((strcasecmp("yes", second) == SAME) ||
			(strcasecmp("y", second) == SAME)) 
		  {
		     parameters->WriteHistory = YES;
		     printf("  Writing history file.\n");
		  }
	       else 
		  {
		     parameters->WriteHistory = NO;
		     printf("  Not writing history file.\n");
		  }
	       break;
	    case 8:
	    case 9:			/* presolve/preprocess (synonymous) */
	       if (second == NULL) 
		  {
		     parameters->Preprocessing = YES;
		     printf("  Presolving.\n");
		  } 
	       else if ((strcasecmp("no", second) == SAME) ||
			(strcasecmp("n", second) == SAME)) 
		  {
		     parameters->Preprocessing = NO;
		     printf("  NO Presolving.\n");
		  } 
	       else 
		  {
		     parameters->Preprocessing = YES;
		     printf("  Presolving.\n");
		  }
	       break;
	    case 10:			/* input directory */
	       if (second != NULL) 
		  {
		     parameters->InputDirectory = StrDup(second, 
							 "InputDirectory");
		     printf("  Look in %s if input file not found in current directory\n", second);
		  }
	       break;
	    case 11:			/* opttol */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for OptTol");
		     printf(" in parameter file.\n");
		  }
	       else 
		  {
		     parameters->OptTol = atof(second);
		     printf("  OptTol = %e\n", parameters->OptTol);
		  }
	       break;
	    case 12:			/* prifeastol */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for PriFeasTol");
		     printf(" in parameter file.\n");
		  }	       
	       else 
		  {
		     parameters->PriFeasTol = atof(second);
		     printf("  PriFeasTol = %e\n", parameters->PriFeasTol);
		  }
	       break;
	    case 13:			/* dualfeastol */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for DualFeasTol");
		     printf(" in parameter file.\n");
		  }
	       else 
		  {
		     parameters->DualFeasTol = atof(second);
		     printf("  DualFeasTol = %e\n", parameters->DualFeasTol);
		  }
	       break;
	    case 14:			/* cachesize */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for ");
		     printf("CacheSize in parameter file.\n");
		  }
	       else 
		  {
		     parameters->CacheSize = atoi(second);
		     printf("  CacheSize = %d\n", parameters->CacheSize);
		  }
	       break;
	    case 15:			/* unrollinglevel */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for ");
		     printf("UnrollingLevel in parameter file.\n");
		  }
	       else 
		  {
		     parameters->UnrollingLevel = atoi(second);
		     printf("  UnrollingLevel = %d\n", 
			    parameters->UnrollingLevel);
		  }
	       break;
	    case 16:			/* iterationlimit */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for");
		     printf(" IterationLimit in parameter file.\n");
		  }
	       else 
		  {
		     parameters->IterationLimit = atoi(second);
		     printf("  IterationLimit = %d\n", 
			    parameters->IterationLimit);
		  }
	       break;
	    case 17:                    /* centerexponent */
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for CenterExp");
		     printf(" in parameter file.\n");
		  }
	       else 
		  {
		     parameters->CenterExponent = atof(second);
		     printf("  Centering Exponent = %f\n", 
			    parameters->CenterExponent);
		  }
	       break;
	    case 18:    /* iterative refinement */
	       if (second == NULL) 
		  {
		     parameters->Refinement = YES;
		     printf("  Doing iterative refinement.\n");
		  } 
	       else if ((strcasecmp("no", second) == SAME) ||
			     (strcasecmp("n", second) == SAME)) 
		  {
		     parameters->Refinement = NO;
		     printf("  NO iterative refinement.\n");
		  } 
	       else 
		  {
		     parameters->Refinement = YES;
		     printf("  Doing iterative refinement.\n");
		  }
	       break;
	    case 19: 
	       if (second == NULL)
		  {
		     printf("  Missing numerical value for StepFactor"); 
		     printf(" in parameter file.\n");
		  }
	       else 
		  {
		     parameters->AlphaScale = atof(second);
		     printf("  StepFactor = %f\n", parameters->AlphaScale);
		  }
	       break;
	    case 20:    /* scaling */
	       if (second == NULL) 
		  {
		     parameters->Scaling = YES;
		     printf("  Doing Scaling.\n");
		  } 
	       else if ((strcasecmp("no", second) == SAME) ||
			     (strcasecmp("n", second) == SAME)) 
		  {
		     parameters->Scaling = NO;
		     printf("  NO Scaling.\n");
		  } 
	       else 
		  {
		     parameters->Scaling = YES;
		     printf("  Doing Scaling.\n");
		  }
	       break;
	    case 21:    /* Gondzio higher-order corrections */
	       if (second == NULL) 
		  {
		     parameters->HOCorrections = 1;
		     printf("  Doing Gondzio higher-order corrections.\n");
		  } 
	       else if ((strcasecmp("no", second) == SAME) ||
			(strcasecmp("n", second) == SAME)) 
		  {
		     parameters->HOCorrections = 0;
		     printf("  No Gondzio higher-order corrections.\n");
		  } 
	       else 
		  {
		     parameters->HOCorrections = atoi(second);
		     printf("  Doing Gondzio corrections.\n");
		  }
	       break;
	    case 22:    /* Gondzio higher-order corrections */
	       if (second == NULL) 
		  {
		     parameters->MaxCorrections = 0;
		     printf("  Maximum Gondzio corrections to be");
		     printf(" determined automatically by PCx.\n");
		  } 
	       else 
		  {
		     parameters->MaxCorrections = atoi(second);
		     printf("  Maximum Gondzio corrections = %d\n", 
			    parameters->MaxCorrections);
		  }
	       break;
	    case 23:    /* Order Algorithm */
	       if (second == NULL)
		  {
		     parameters->OrderAlg = 1;
		     printf("  Order Algorithm to be");
		     printf(" multiple minimum degree.\n");
		  }
	       else
		  {
                     switch (atoi(second)) {
                       case 0:
                         parameters->OrderAlg = 0;
                         printf("  Order Algorithm to be disable\n");
                         break;
                       case 2:
                         parameters->OrderAlg = 2;
                         printf("  Order Algorithm to be");
                         printf(" Reverse Cuthill-McKee.\n");
                         break;
                       case 1:
                       default:
                         parameters->OrderAlg = 1;
                         printf("  Order Algorithm to be");
                         printf(" multiple minimum degree.\n");
                         break;
                     }
		  }
	       break;
	    } /* end switch */
	 do 
	    {
	       fgets(line, 200, fp);
	    } while (line[0] == '#' && !feof(fp));
      }
   fclose(fp);
   return 0;
}

/***************************************************************************/

