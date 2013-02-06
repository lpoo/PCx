/************************************************************************

MATLAB - PCx interface

Michael Wagner
Department of Mathematics and Statistics
Old Dominion University

August 2000

 ***********************************************************************/
#include "PCx_mex.h"

/**********************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, row, col, nnz_counter, nnz, nnz_eq, nnz_ineq, n, m, m_eq, m_ineq;
  double *b_eq, *c, *b_ineq, *l, *u;
  MPStype *mps, *NewMPS();
  double *sol, *lb, *ub, *b_temp;
  double *new_sol;
  int *eq_pBeginRow, *ineq_pBeginRow, *eq_Row, *ineq_Row;
  double *eq_Value, *ineq_Value;
  solution *Solution;
  Parameters *Inputs, *NewParameters();
  
  int *matbeg, *matind, j, j1, status, dims[2];
  double *matval;
  char **field_names;

  mxArray *lambda;
  
  field_names = (char **) mxMalloc(4*sizeof(char *)); 

  for (i=0; i<4; i++)
    field_names[i] = (char *) mxMalloc(100*sizeof(char));
 
  if ((nrhs < 3) || (nrhs > 9))
    {
      char msgbuf[1024];
      sprintf(msgbuf,"Syntax: xsol = PCx_mex(c, A_ineq, b_ineq, A_eq, b_eq, lb, ub, x0, [], PCx_opts).\n\nThe first 7 arguments are mandatory. You can pass an empty structure\nif argument is not applicable.\n");
      AbortMsgTxt(msgbuf);
    }

  Inputs = NewParameters(); 
  Inputs->ReportingLevel = 0;
  Inputs->Diagnostics    = 0;

  if (nrhs > 8)
    {
      /* now read the optional options supplied by the user */
 
      Inputs = ReadOptions(prhs[8], Inputs);
      /* check for errors */
      if (CheckParameters(Inputs)) 
	{
	  char msgbuf[1024];
	  sprintf(msgbuf,"Error: PCx_opts contains nonsense.\n");
	  AbortMsgTxt(msgbuf);
	}
    }

  /* Check data type of input arguments */

  /* grab objective function coefficients */
  
  c = get_double_vec(prhs,0);
  if (c == NULL)
    {
      char msgbuf[1024];
      sprintf(msgbuf,"Error: need to specify meaningful objective function.\n");
      AbortMsgTxt(msgbuf);
    }
  n = mxGetM(prhs[0]);

  /* grab A_ineq and b_ineq */

  if (!mxIsEmpty(prhs[1]))
    {
      get_matrix(prhs, 1, &ineq_pBeginRow, &ineq_Row, &ineq_Value);
      b_ineq = get_double_vec(prhs,2);
      if (n != mxGetN(prhs[1]))
	{
	  char msgbuf[1024];
	  sprintf(msgbuf,"Error: c and A_ineq have to have the same number of columns.\n");
	  AbortMsgTxt(msgbuf);
	}
      m_ineq   = mxGetM(prhs[2]);
      if (m_ineq != mxGetM(prhs[1]))
	{
	  char msgbuf[1024];
	  printf("mxGetM(b_ineq) = %d, mxGetM(A_ineq) = %d\n", mxGetM(prhs[2]), mxGetM(prhs[1]));
	  sprintf(msgbuf,"Error: Dimensions of A_ineq and b_ineq have to be consistent.\n");
	  AbortMsgTxt(msgbuf);
	}
    }
  else
    {
      if (Inputs->Diagnostics > 0)
	printf("Matrix A_ineq is empty, ignored\n");

      ineq_pBeginRow = NULL;
      ineq_Row       = NULL;
      ineq_Value     = NULL;
      b_ineq         = NULL;
      m_ineq         = 0;
    }

  if (nrhs > 3)
    {
      /* grab A_eq and b_eq */
      
      if (!mxIsEmpty(prhs[3]))
	{
	  get_matrix(prhs, 3, &eq_pBeginRow, &eq_Row, &eq_Value);
	  if (n != mxGetN(prhs[3]))
	    {
	      char msgbuf[1024];
	      sprintf(msgbuf,"Error: Dimensions of c and A_eq have to be consistent.\n");
	      AbortMsgTxt(msgbuf);
	    }
	  
	  m_eq   = mxGetM(prhs[4]);
	  if (m_eq != mxGetM(prhs[3]))
	    {
	      char msgbuf[1024];
	      sprintf(msgbuf,"Error: Dimensions of A_eq and b_eq have to be consistent.\n");
	      AbortMsgTxt(msgbuf);
	    }

	  b_eq   = get_double_vec(prhs,4);
	}
      else
	{
	  if (Inputs->Diagnostics > 0)
	    printf("Matrix A_eq is empty, ignored.\n");
	  
	  eq_pBeginRow = NULL;
	  eq_Row = NULL;
	  eq_Value = NULL;
	  b_eq = NULL;
	  m_eq = 0;
	}
    }

  if (nrhs > 5)
    {
      /* grab bounds */
      
      lb     = get_double_vec(prhs,5);
      ub     = get_double_vec(prhs,6);
      
      if (((lb != NULL) && (n != mxGetM(prhs[5]))) || ((ub != NULL) && (n != mxGetM(prhs[6]))))
	{
	  char msgbuf[1024];
	  sprintf(msgbuf,"Error: Dimensions of bounds and c have to be consistent.\n");
	  AbortMsgTxt(msgbuf);
	}
    }

  if ((nrhs > 7) && ( !mxIsEmpty(prhs[7]) ))
    printf("Warning: starting point ignored in PCx.\n");

  
  m = m_eq + m_ineq;
  
  nnz_eq = (m_eq > 0) ? *(mxGetJc(prhs[3]) + n) : 0;
  nnz_ineq = (m_ineq > 0) ? *(mxGetJc(prhs[1]) + n) : 0;
  
  nnz = nnz_eq + nnz_ineq;
  
  /* printf("nnz = %d\n",nnz); */
  
  /* ok, start creating PCx-datastructures */
  
  mps = NewMPS(m,n,nnz);
  
  /* Now do data structure conversions*/
  mps->NumRows = m;
  mps->NumCols = n;
  
  /* easy stuff first... */
  
  /* Start setting up bounds*/
  
  for (i=0; i<m_eq; i++)
    mps->RowType[i]='E';
  
  for (i=0; i<m_ineq; i++)
    mps->RowType[m_eq+i]='L';
  
  mps->cshift = 0.0;

  for (i=0; i<n; i++)
    {
      if (lb != NULL && mxIsFinite(lb[i]))
	{
	  if (lb[i] == 0.0)
	    if (ub != NULL && mxIsFinite(ub[i]))
	      {
		mps->BoundType[i] = UPPER;
		mps->UpBound[i] = ub[i];
	      }
	    else
	      mps->BoundType[i] = NORMAL;
	  else
	    if (ub != NULL && mxIsFinite(ub[i]))
	      {
		mps->BoundType[i] = UPPERLOWER;
		mps->UpBound[i] = ub[i];
	      }
	    else
	      mps->BoundType[i] = LOWER;
	  
	  mps->LowBound[i] = lb[i];
	}
      else
	if (ub != NULL && mxIsFinite(ub[i]))
	  {
	    mps->BoundType[i] = MINFTY;
	    mps->UpBound[i] = ub[i];
	  }
	else
	  mps->BoundType[i] = FREE;
    }
  
  mps->c = c;

  for (i=0; i<m_eq; i++)
    mps->b[i] = b_eq[i];
     
  for (i=0; i<m_ineq; i++)
    mps->b[m_eq+i] = b_ineq[i];
      
  /* Now the fun begins, need to merge A_eq and A_ineq */
  
  nnz_counter = 0;
  
  if ((eq_pBeginRow != NULL) || (ineq_pBeginRow != NULL))
    for (col=0; col<n; col++)
      {
	mps->A.pBeginRow[col] = nnz_counter+1;
	if (eq_pBeginRow != NULL)
	  for (row=eq_pBeginRow[col]; row<eq_pBeginRow[col+1]; row++)
	    {
	      mps->A.Row[nnz_counter] = eq_Row[row]+1;
	      mps->A.Value[nnz_counter] = eq_Value[row];
	      nnz_counter++;
	    }
	if (ineq_pBeginRow != NULL)
	  for (row=ineq_pBeginRow[col]; row<ineq_pBeginRow[col+1]; row++)
	    {
	      mps->A.Row[nnz_counter] = ineq_Row[row] + m_eq+1;
	      mps->A.Value[nnz_counter] = ineq_Value[row];
	      nnz_counter++;
	    }
	mps->A.pEndRow[col] = nnz_counter;
      }
  
  if (nnz_counter != nnz)
    {
      char msgbuf[1024];
      sprintf(msgbuf,"Something is fishy with nnz_counter.\n");
      AbortMsgTxt(msgbuf);
    }
  
  /* one last thing before we go.... */
  
  AddMPSProblemName(mps, "from Matlab");

  /* off we go! */
  status = PCx_main(mps, Inputs, &Solution);

  return_double_vec(plhs, 0, Solution->x, Solution->Columns);

  if (nlhs > 1)
    /* return (primal) objective function value */
    plhs[1] = mxCreateDoubleScalar(Solution->PrimalObjective);
 
  if (nlhs > 2)
    /* return status indicator, should probably do more checking here.... */
    {
      switch (Solution->Status)
	{
	case UNKNOWN_SOL:
	  status = -1;
	  break;
	case OPTIMAL_SOL:
	  status = 1;
	  break;
	case INFEASIBLE_SOL:
	  status = -1;
	  break;
	case SUBOPTIMAL_SOL:
	  status = 0;
	  break;
	default:
	  status = -1;
	  break;
	}

      plhs[2] = mxCreateDoubleScalar((double) status);
    }

  if (nlhs > 3)
    {
      /* return output structure with useful (?) information */
      dims[0] = 1; dims[1] = 1;

      strcpy(field_names[0], "iterations");
      strcpy(field_names[1], "algorithm");
      strcpy(field_names[2], "cgiterations");
      plhs[3] = mxCreateStructArray(1, dims, 3, (const char**) field_names);
      
      mxSetField(plhs[3], 0, "iterations",   mxCreateDoubleScalar(Solution->Iterations));
      mxSetField(plhs[3], 0, "algorithm",    mxCreateString("PCx"));
      mxSetField(plhs[3], 0, "cgiterations", mxCreateDoubleScalar(0));
    }

  if (nlhs > 4)
    {
      /* need to return dual solution here */
      dims[0] = 1; dims[1] = 1;
   
      strcpy(field_names[0], "ineqlin");
      strcpy(field_names[1], "eqlin");
      strcpy(field_names[2], "lower");
      strcpy(field_names[3], "upper");

      plhs[4] = mxCreateStructArray(1, dims, 4, (const char **) field_names);

      lambda = mxCreateDoubleMatrix(m_ineq, 1, mxREAL);
      memcpy(mxGetPr(lambda), &(Solution->pi[m_eq]), m_ineq*sizeof(double));      
      mxSetField(plhs[4], 0, "ineqlin", lambda);

      lambda = mxCreateDoubleMatrix(m_eq, 1, mxREAL);
      memcpy(mxGetPr(lambda), Solution->pi, m_eq*sizeof(double));      
      mxSetField(plhs[4], 0, "eqlin",   lambda);

      lambda = mxCreateDoubleMatrix(n, 1, mxREAL);
      memcpy(mxGetPr(lambda), Solution->DualLower, n*sizeof(double));      
      mxSetField(plhs[4], 0, "lower",   lambda);

      lambda = mxCreateDoubleMatrix(n, 1, mxREAL);
      memcpy(mxGetPr(lambda), Solution->DualUpper, n*sizeof(double));      
      mxSetField(plhs[4], 0, "upper",   lambda);
    }

  /* Matlab claims to clean up automatically, so we don't do it */
  /* If we did the program would crash..... */
  /* Who knows why */
  
  /* I actually checked the memory usage during and after the mex-file 
     executes, and Matlab does seem to do a good job of cleaning up. 
     But I can't guarantee that nothing is leaking here.... */
  
  /*
  for (i=0; i<4; i++)
    Free((void *) field_names[i]);

  Free((void *) field_names);
  FreeSolution(Solution);
  DeleteMPS(mps);
  FreeParameters(Inputs);

  */

}

/*************************************************************/

int PCx_main(MPStype *MPS, Parameters *Inputs, solution **Solution)
{
  LPtype *LP, *ReducedLP;
  double OldUserTime, OldSysTime, pretime=0.0, UserTime, SysTime;
  double readtime = 0.0;
  int i;
  void            SplitFreeVars(), UnSplitFreeVars();
  int             PrintSolution();
  MPSchanges     *Changes;
  ChangeStack    *Record;
  solution       *NewSolution();
  int passes, status;
  
  
  /* convert the MPStype data into LPtype data, which allows only equality
   * constraints and three kinds of x components: free, nonnegative, or
   * bounded as in 0 <= x_i <= u_i */
  
  LP = Convert_MPS_LP(MPS, &Changes);
  
  if (Inputs->Diagnostics > 0)
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
  
  *Solution = NewSolution(ReducedLP->Rows, ReducedLP->Cols,
			 Inputs->IterationLimit);
  (*Solution)->ReadTime = readtime;
  (*Solution)->PreprocessTime = pretime;
  
  /* solve the problem, keeping track of CPU times.  */
  
  GetTime(&OldUserTime, &OldSysTime);
  
  status = PCx(ReducedLP, *Solution, Inputs);
  
  if (status != 0) 
    {
      char msgbuf[1024];
      sprintf(msgbuf,"Error in PCx(). Exiting with code %d\n", status);
      AbortMsgTxt(msgbuf);
    }
  
  GetTime(&UserTime, &SysTime);
  (*Solution)->SolutionTime = UserTime - OldUserTime + SysTime - OldSysTime;
  
  /* recover the free variable values by combining their positive and
   * negative parts, and undo the effects of the preprocessor to recover the
   * solution to the original problem LP.  */
  
  UnSplitFreeVars(ReducedLP, *Solution);
  
  if (Inputs->Scaling)
    UnscaleLP(ReducedLP, *Solution);
  
  if (Inputs->Preprocessing) 
    if (Postprocess(LP, &Record, Solution) < 0) 
      {
	char msgbuf[1024];
	sprintf(msgbuf,"Error in postprocessing phase of PCx(). Oops.\n");
	AbortMsgTxt(msgbuf);
      } 
  
  /* Check the infeasibilities for the point obtained */
  ComputeAndPrintInfeasibilities(*Solution, LP);
  
  /* Express the solution of LP in terms of the original MPS formulation */
  *Solution = MPSsolution(LP, MPS, *Solution, Changes, Inputs);
  
  /* The following function turns out to be a major headache, which is 
     why I completely rewrote it. Who knows why. */

   
  if (Inputs->Diagnostics > 0)
    PrintSolInfo(*Solution);
  
  /* clean up! */
  
  /*
  DeleteChangeStack(Record);
  DeleteChanges(Changes);
  if (ReducedLP != LP)
    DeleteLP(ReducedLP);
  DeleteLP(LP);
  */

  return(status);
  
  /* TrDump(stdout); */
}

/**********************************************************************/

Parameters *ReadOptions(mxArray *opts, Parameters *parameters)
{
  int n_args, i;
  int match, key;
  const char *field_name;
  char *token;
  int             NumKeys = 19;
  double field_value;
  mxArray *field;
  static char    *KeyWord[] = {"max", "presolve", "opttol", 
			       "prifeastol", "dualfeastol", "cachesize", 
			       "unrollinglevel", "iterationlimit", 
			       "centerexp", "refinement", "stepfactor", 
			       "scaling", "hocorrections", 
			       "maxcorrections", "MaxIter", "TolFun", "TolCon", "Display",
                               "Diagnostics"};
  
  parameters->ReportingLevel = 0;

  token = (char *) mxMalloc(100*sizeof(char));

  if (!mxIsStruct(opts))
    {
      printf("PCx_opts is not a structure, ignoring.....\n");
      return(parameters);
    }
  else
    {
      n_args = mxGetNumberOfFields(opts);
      
      if (n_args>0)
	for (i=0; i<n_args; i++)
	  {
	    field_name = mxGetFieldNameByNumber(opts, i);
	    /* printf("Analyzing contents of field %s\n", field_name); */
	    field     = mxGetField(opts, 0, field_name);

	    if (!mxIsEmpty(field))
	      field_value = mxGetScalar(field);
	    else
	      field_value = 0;

	    for (key = 0, match = -1; key < NumKeys; key++) 
	      {
		if (strcasecmp(KeyWord[key], field_name) == 0) 
		  {
		    match = key;
		    /* printf("Analyzing %s\n", KeyWord[key]); */
		    break;
		  }
	      }
	    
	    switch (match) 
	      {
	      case -2:
		break;
	      case -1:
		/*
		  printf("Unknown keyword '%s' in PCx_opts\n", field_name);
		  printf("Ignoring.\n");
		*/
		break;
	      case 0:			/* max */
		if (field_value > 0)
		  {
		    parameters->Minimize = NO;
		    printf("  Maximizing the objective function.\n");
		  }
		else
		  { 
		    parameters->Minimize = YES;
		    printf("  Minimizing the objective function.\n");
		  }
		break;
	      case 1:                    /* preprocessing */
		if (field_value > 0)
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
	      case 2:			/* opttol */
		parameters->OptTol = field_value;
		printf("  OptTol = %e\n", parameters->OptTol);
		break;
	      case 3:			/* prifeastol */
		parameters->PriFeasTol = field_value;
		printf("  PriFeasTol = %e\n", parameters->PriFeasTol);
		break;
	      case 4:			/* dualfeastol */
		parameters->PriFeasTol = field_value;;
		parameters->DualFeasTol = field_value;
		printf("  DualFeasTol = %e\n", parameters->DualFeasTol);
		break;
	      case 5:			/* cachesize */
		parameters->CacheSize = (int) field_value;
		printf("  CacheSize = %d\n", parameters->CacheSize);
		break;
	      case 6:			/* unrollinglevel */
		parameters->UnrollingLevel = (int) field_value;
		printf("  UnrollingLevel = %d\n", parameters->UnrollingLevel);
		break;
	      case 7:			/* iterationlimit */
		parameters->IterationLimit = (int) field_value;
		printf("  IterationLimit = %d\n", parameters->IterationLimit);
		break;
	      case 8:
		parameters->CenterExponent = field_value;
		printf("  Centering Exponent = %f\n", parameters->CenterExponent);
		break;
	      case 9:    /* iterative refinement */
		if (field_value > 0)
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
	      case 10: 
		parameters->AlphaScale = field_value;
		printf("  StepFactor = %f\n", parameters->AlphaScale);
		break;
	      case 11:    /* scaling */
		if (field_value > 0)
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
	      case 12:    /* Gondzio higher-order corrections */
		if (field_value)
		  {
		    parameters->HOCorrections = 0;
		    printf("  No Gondzio higher-order corrections.\n");
		  } 
		else 
		  {
		    parameters->HOCorrections = field_value;
		    printf("  Doing Gondzio corrections.\n");
		  }
		break;
	      case 13:    /* Gondzio higher-order corrections */
		if (field_value)
		  {
		    parameters->MaxCorrections = 0;
		    printf("  Maximum Gondzio corrections to be");
		    printf(" determined automatically by PCx.\n");
		  } 
		else 
		  {
		    parameters->MaxCorrections = (int) field_value;
		    printf("  Maximum Gondzio corrections = %d\n", 
			   parameters->MaxCorrections);
		  }
		break;
	      case 14:    /* MaxIter */	
		if (!mxIsEmpty(field))
		  {
		    parameters->IterationLimit = (int) field_value;
		    printf("  IterationLimit = %d\n", parameters->IterationLimit);
		  }
		break;
	      case 15:    /* TolFun */
		if (!mxIsEmpty(field))
		  {
		    parameters->OptTol = field_value;
		    printf("  OptTol = %e\n", parameters->OptTol);
		  }
		break;
	      case 16:    /* TolCon */
		if (!mxIsEmpty(field))
		  {
		    parameters->PriFeasTol = field_value;
		    printf("  PriFeasTol = %e\n", parameters->PriFeasTol);
		    parameters->DualFeasTol = field_value;
		    printf("  DualFeasTol = %e\n", parameters->DualFeasTol);
		  }
		break;
	      case 17:    /* Display */
		if (!mxIsEmpty(field))
		  {
		    mxGetString(field, token, 100);
		    
		    if (!strcasecmp(token, "iter"))
		      parameters->ReportingLevel = 2;
		    else if (!strcasecmp(token, "final"))
		      parameters->ReportingLevel = 1;
		    		  }
		break;
	      case 18:    /* Diagnostics */
		if (!mxIsEmpty(field))
		  {
		    mxGetString(field, token, 100);
		    
		    if (!strcasecmp(token, "on"))
		      parameters->Diagnostics = 1;
		    else if (!strcasecmp(token, "off"))
		      parameters->Diagnostics = 0;
		  }
		break;
	      } /* end switch */
	  }
    }
  /* mxFree(field_name); */
  /* printf("Exiting ReadOptions()\n"); */

  /* mxFree(token); */

  return(parameters);
}

/**********************************************************************/

/* Prints the solution from "Solution" to the output file "infile".out. Also
 * writes the history file "infile".log */
int PrintSolInfo(solution *Solution)
{
  char            statusTxt[40]; 
  int             row, col, entry, i, j, len;
  
  printf("\nProblem from Matlab terminated with ");
  
  switch (Solution->Status) 
    {
    case OPTIMAL_SOL:
      printf("OPTIMAL");
      break;
    case SUBOPTIMAL_SOL:
      printf("SUBOPTIMAL");
      break;
    case INFEASIBLE_SOL:
      printf("INFEASIBLE");
      break;
    case UNKNOWN_SOL:
      printf("UNKNOWN");
      break;
    default:
      printf("UNKNOWN");
    }
  
  printf(" status (code %d) after %d iterations\n\n", 
	 Solution->Status, Solution->Iterations);
  
  printf("Primal Objective = %13.8e\n", Solution->PrimalObjective);
  printf("Dual   Objective = %13.8e\n", Solution->DualObjective);
   
  printf("\nComplementarity          = %9.2e\n", Solution->Complementarity);
  printf("Relative Complementarity = %9.2e\n",
	 Solution->RelativeComplementarity);
  
  printf("\nRelative Infeasibilities:\n");
  printf("Primal = %9.3e,    ", Solution->PrimalInfeasibility);
  printf("Dual   = %9.3e.\n", Solution->DualInfeasibility);
  /* printf("\nRead Time       = %.2f seconds\n", Solution->ReadTime); */
  printf("\nPreprocess time = %.2f seconds\n", Solution->PreprocessTime);
  printf("Solution time   = %.2f seconds\n", Solution->SolutionTime);
  
  return(0); 
}
