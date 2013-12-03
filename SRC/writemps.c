/* writemps
 *
 * PCx  9/19/96. 
 *
 * Author: Martin Wechs
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

/*
 * utility function to dump an MPS data structure as an MPS file
 */

#include<stdio.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include"main.h" 
#include"memory.h"
#include"def.h"

void  printSTR(outfile,name,offset,length)
FILE  *outfile;
char  *name;
int   offset,length;
{
char  spaces[16];
int   i;

#if (DEBUG>3)
 printf("%d %d\n",offset,length);  
#endif

      for(i=0;i<16;i++) spaces[i]=' '; 

      spaces[offset]='\0';
      fprintf(outfile,"%s",spaces);
      spaces[offset]=' ';
      fprintf(outfile,"%s",name);
      while (*name++ !='\0') length--;
      spaces[length]='\0';
      fprintf(outfile,"%s",spaces);
}

void  printVAL(outfile,value,offset,length)
FILE  *outfile;
double value;
int    offset,length;
{
char   spaces[40];
int    i,k;
double vlog;


      for(i=0;i<16;i++) spaces[i]=' '; 

      spaces[offset]='\0';
      fprintf(outfile,"%s",spaces);
      spaces[offset]=' ';   

      vlog=0; if (value!=0) vlog=log10(fabs(value));

#if (DEBUG>3)
  printf("Value:  %4.14f %4.14f %d  %d\n",vlog,value ,offset,length); 
#endif

      if (fabs(vlog)>=length-2) {
              printf("Can not convert double to MPS format\n");
              fprintf(outfile,"Number too large");
      } else {
        i=0; if (vlog>=0) { i=fabs(vlog); i++; }

#if (DEBUG>4)
  printf("%d\n",i);
#endif
        length--;
        offset=length-i; if (value<0) offset--;
        spaces[0]='%'; 
        k=0; if (i>9) { spaces[1]='1'; i-=10; k++; } 
        spaces[k+1]=i+'0'; 
        spaces[k+2]='.'; 
        if (offset>9) { spaces[k+3]='1'; offset-=10; k++; } 
        spaces[3+k]=offset+'0';
        spaces[4+k]='f'; spaces[5+k]='\0';

#if (DEBUG>4)
  printf("%s\n",spaces);
#endif
        sprintf(&spaces[7],spaces,value);

#if (DEBUG>4)
  printf("%s\n",&spaces[7]);
#endif

        k=7; if (spaces[1]=='0') if (spaces[k++]=='-') spaces[k]='-';  
        i=k; k+=length;
	while (spaces[k]=='0') spaces[k--]=' '; 
        fprintf(outfile,"%s",&spaces[i]); 
      }
} 

void printROW(outfile,NumRows,MPSValue,RowNames,Name)
FILE     *outfile;
int      NumRows;
double   *MPSValue;
char     **RowNames;
char     *Name;
{
int      i,k=0;
double   value;

   for (i=0; i<NumRows;i++) { 
     value=MPSValue[i]; 
     if (value!=0.0) { 
       if (!(k&1)) {
         fprintf(outfile,"\n");
         printSTR(outfile,Name,4,8);
       }
       printSTR(outfile,RowNames[i],2+(k&1),8);
       printVAL(outfile,value,2,12);
       k++;
     }
   }
}

void WriteMPS(MPS,Inputs)
MPStype    *MPS;
Parameters *Inputs;
{
FILE	     *outfile; 
sparseMatrix *A;
char         c,filename[23];
char         *ProblemName,*ColName;
double       value;
int          offset1,offset2,offset3;
int          i,k;

     filename[8]='_';filename[9]='p'; filename[10]='c'; 
     filename[11]='x'; filename[12]='.'; filename[13]='m';
     filename[14]='p'; filename[15]='s'; filename[16]='\0';

     A=&(MPS->A);

#if (DEBUG>2)
  printf("Copy name\n");
#endif

     ProblemName=MPS->ProblemName;
     for (i=0;i<8;i++) filename[i]=ProblemName[i];

#if (DEBUG>2)
  printf("Test name\n"); 
#endif
   
     offset1=-1; offset2=8;
     for (i=0;i<8;i++) {
       if (filename[i]=='\0') filename[i]=' ';
       if ((filename[i]!=' ')&&(offset1==-1)) offset1=i;
       if ((filename[i]==' ')&&(offset1!=-1)) { offset2=i; i=8; }
     }     
     if (offset1<0) printf("Sorry but there is no LP-Name\n");   

#if (DEBUG>1)
  printf("Add extention %d %d\n",offset1, offset2);
#endif

     offset3=8-offset2;
     for (i=offset2-1; i>=offset1; --i) {
         c=ProblemName[i];
         filename[i+offset3]=(isupper(c) ?  tolower(c) : c);
     }

#if (DEBUG>1)
  printf("Open file %s\n",&filename[offset3]);  
#endif

     outfile=fopen(&filename[offset3],"w");

#if (DEBUG>1)
  printf("Write NAME\n");    
#endif

     fprintf(outfile,"NAME          ");
     fprintf(outfile,"%s",ProblemName);

#if (DEBUG>1)
  printf("ROWS\n");     
#endif

     filename[1]='\0';
     fprintf(outfile,"\nROWS\n");     
     for (i=0; i<MPS->NumRows;i++) { 
        filename[0]= (MPS->RowType)[i];
        fprintf(outfile," %s  ",filename);
        fprintf(outfile,"%s\n",(MPS->RowNames)[i]);   
     }
  
     fprintf(outfile," N  ");
     fprintf(outfile,"%s\n",MPS->ObjectiveName);   


#if (DEBUG>1)
  printf("COLUMNS\n");
#endif
 
     fprintf(outfile,"COLUMNS");     
     for (i=0; i<MPS->NumCols;i++) { 
       ColName=(MPS->ColNames)[i];
       for (k=(A->pBeginRow)[i]-1;k<(A->pEndRow)[i];){
           fprintf(outfile,"\n"); 
           printSTR(outfile,ColName,4,8);
           printSTR(outfile,(MPS->RowNames)[(A->Row)[k]-1],2,8);
           printVAL(outfile,(A->Value)[k++],2,12);
           if (k<A->pEndRow[i]) {
                printSTR(outfile,(MPS->RowNames)[A->Row[k]-1],3,8);
                printVAL(outfile,(A->Value)[k++],2,12);
           }  
       } 
       if (MPS->c[i]!=0) {
           k=A->pEndRow[i]-A->pBeginRow[i];
           if (k&1) {  
             fprintf(outfile,"\n"); 
             printSTR(outfile,ColName,4,8);
           }
           printSTR(outfile,MPS->ObjectiveName,3-(k&1),8);
           value=MPS->c[i];
           if (!(Inputs->Minimize)) value=-value;
           printVAL(outfile,value,2,12); 
       }
     }
  

#if (DEBUG>1)
  printf("RHS\n");
#endif

     fprintf(outfile,"\nRHS");

     printROW(outfile,MPS->NumRows,MPS->b,MPS->RowNames,MPS->RHSName);
     

#if (DEBUG>1)
  printf("SHIFT\n");
#endif

     if (MPS->cshift!=0) {
       if (k&1) {
         fprintf(outfile,"\n");
         printSTR(outfile,MPS->RHSName,4,8);
       }
       printSTR(outfile,MPS->ObjectiveName,2,8);
       printVAL(outfile,MPS->cshift,2+(k&1),12); 
     } 


#if (DEBUG>1)
  printf("RANGES\n");
#endif

     /* Here are the RANGES */
     
     if (MPS->RangeName!=NULL) {
       fprintf(outfile,"\nRANGES"); 
       printROW(outfile,MPS->NumRows,MPS->Ranges,MPS->RowNames,MPS->RangeName);
     }

#if (DEBUG>1)
  printf("BOUNDS\n");
#endif

     /* Here are the BOUNDS */

     if (MPS->BoundName!=NULL) {

       /*  filename = "PL\0 FR\0 UP\0 LO\0 FX\0 MI"; */
       filename[0] ='P'; filename[1] ='L'; filename[2] ='\0';
       filename[4] ='F'; filename[5] ='R'; filename[6] ='\0';
       filename[8] ='U'; filename[9] ='P'; filename[10]='\0';
       filename[12]='L'; filename[13]='O'; filename[14]='\0';
       filename[16]='F'; filename[17]='X'; filename[18]='\0';
       filename[20]='M'; filename[21]='I'; filename[22]='\0';
      
       fprintf(outfile,"\nBOUNDS");
       for (i=0;i<MPS->NumCols;i++) {
         k=MPS->BoundType[i];
	 if (k) {
           fprintf(outfile,"\n");
	   if (k!=UPPERLOWER) {

             offset1=4*k;

#if (DEBUG>2)
  printf("Bound: %d %d %d\n",k,offset1,i);
#endif

             printSTR(outfile,&filename[offset1],1,2);
             printSTR(outfile,MPS->BoundName,1,8);
             printSTR(outfile,MPS->ColNames[i],2,8);

#if (DEBUG>2)
  printf("Begin switch\n");
#endif

             switch (k) {
	     case LOWER:
               printVAL(outfile,MPS->LowBound[i],2,12);
               break;
             case UPPER:
               printVAL(outfile,MPS->UpBound[i],2,12);
               break;
             case FIX:
               printVAL(outfile,MPS->UpBound[i],2,12);
             }  

#if (DEBUG>2)
  printf("End switch\n");
#endif

           } else {

             printSTR(outfile,&filename[12],1,2);
             printSTR(outfile,MPS->BoundName,1,8);
             printSTR(outfile,MPS->ColNames[i],2,8); 
             printVAL(outfile,MPS->LowBound[i],2,12);

             fprintf(outfile,"\n"); 

             printSTR(outfile,&filename[8],1,2);
             printSTR(outfile,MPS->BoundName,1,8);
             printSTR(outfile,MPS->ColNames[i],2,8); 
             printVAL(outfile,MPS->UpBound[i],2,12);

           }
         }     
       } /* End FOR i=0:NumCols */
     }  /* endif               */


     fprintf(outfile,"\nENDATA\n");

#if (DEBUG>1)
  printf("Close the file\n");
#endif
    
     close(outfile);
}


MPStype *convertLP_MPS(LP,Name)
LPtype     *LP;
char       *Name;
{
MPStype    *MPS;
char       *AMPS,*ALP;
int        i,j,k;

#if (DEBUG>1)
  printf("New MPS  %s\n",Name);
#endif

  MPS = (MPStype *) Malloc(sizeof(MPStype),"MPS in convertLP_MPS");

#if (DEBUG>2)
  printf("write name\n");
#endif

  MPS->ProblemName = (char *) Malloc(9,"ProblemName in convertLP_MPS"); 
  i=0; while (i<8) { 
    (MPS->ProblemName)[i]=Name[i]; 
    if (Name[i++]=='\0') i=8;
  } 
  MPS->ProblemName[8]='\0';

#if (DEBUG>1)
  printf("Set names  %s\n",MPS->ProblemName);
#endif

  MPS->NumRows = LP->Rows; 
  MPS->NumCols = LP->Cols;
  MPS->NumEnts = LP->Ents;

  MPS->b = LP->b;
  MPS->RowType = (char *) Malloc(sizeof(char)*(LP->Rows),
                                      "RowType in convertLP_MPS");
  MPS->RowNames = (char **) Malloc(sizeof(char*)*(LP->Rows),
                                      "RowNames in convertLP_MPS");
  MPS->ColNames = (char **) Malloc(sizeof(char*)*(LP->Cols),
                                      "ColNames in convertLP_MPS");
  
#if (DEBUG>1)
  printf("RowNames\n");
#endif

  j=3; k=10;
  for (i=0;i<LP->Rows;i++) {
     (MPS->RowType)[i] = 'E';  
     if (i>=k) { j++; k*=10; }
     (MPS->RowNames)[i] = (char *) Malloc(j,"RowName in convertLP_MPS");

#if (DEBUG>3)
  printf("R%d ",i);
#endif

     sprintf(MPS->RowNames[i],"R%d",i);
  } 

#if (DEBUG>1)
  printf("ColNames\n");
#endif

  j=3; k=10;
  for (i=0;i<LP->Cols;i++) { 
     if (i>=k) { j++; k*=10; }
     (MPS->ColNames)[i] = (char *) Malloc(j,"ColName in convertLP_MPS");

#if (DEBUG>3)  
  printf("C%d ",i);
#endif

     sprintf(MPS->ColNames[i],"C%d",i);
  } 

#if (DEBUG>1)
  printf("RHS and matrix\n");
#endif

  MPS->RHSName = (char *) Malloc(4,"RHSName in convertLP_MPS");
  MPS->RHSName = "RHS";

  MPS->c = LP->c;
  MPS->ObjectiveName = (char *) Malloc(4,"ObjectiveName in convertLP_MPS");
  MPS->ObjectiveName = "Obj";

  MPS->cshift = LP->cshift;
  
  AMPS=(char *) &(MPS->A); ALP=(char *) &(LP->A);
 
  for (i=0;i<sizeof(sparseMatrix);i++) AMPS[i] = ALP[i];

#if (DEBUG>2)
  printf(" %d %d\n %d %d\n %d %d\n", (MPS->A).NumRows,(LP->A).NumRows,
                  (MPS->A).NumCols,(LP->A).NumCols,
                   (MPS->A).Nonzeros,(LP->A).Nonzeros);
#endif

  MPS->RangeName = NULL; 

#if (DEBUG>1)
  printf("Bounds\n");
#endif

  k=LP->NumberBounds+LP->NumberFree+LP->NumberSplit;

#if (DEBUG>3)
  printf("Numbers: %d %d %d %d\n",k,
      LP->NumberBounds,LP->NumberFree,LP->NumberSplit);
#endif

  if (k) {
    MPS->BoundName = (char *) Malloc(6,"BoundName in convertLP_MPS");
    MPS->BoundName = "Bound";
  
    MPS->BoundType = NewInt(LP->Cols ,"BoundType in convertLP_MPS");
    MPS->UpBound = LP->UpBound;

#if (DEBUG>2)
  printf("Set BOUND\n");
#endif


    for (i=0; i<LP->NumberBounds;i++) {

#if (DEBUG>3)
  printf("%d  %4.14f\n",i,LP->UpBound[LP->BoundIndex[i]]);
#endif

       k=LP->BoundIndex[i];
       MPS->BoundType[k] = UPPER;

    }

#if (DEBUG>1)
  printf("Set FREE\n");
#endif

    for (i=0; i<LP->NumberFree;i++) {
       k=LP->FreeIndex[i];
       MPS->BoundType[k] = FREE;
    }

#if (DEBUG>1)
  printf("Unsplit\n");
#endif

    for (i=0; i<LP->NumberSplit;i++) {
       k=LP->FreePlus[i];
       MPS->BoundType[k] = FREE;
       k=LP->FreeMinus[i];
       MPS->BoundType[k] = FIX;
       MPS->UpBound[k] = 0.0;
    }    
  } else 
       MPS->BoundName = NULL;

#if (DEBUG>1)
  printf("End convert\n");
#endif
   
   return MPS;
}


void WriteLP_MPS(LP,Name,Inputs)
LPtype      *LP;
Parameters  *Inputs;
char        *Name;
{
MPStype *MPS;

#if (DEBUG>3)
  printf("Name: %s\n",Name);
#endif

        MPS=convertLP_MPS(LP,Name);

#if (DEBUG>3)
  printf("Name: %s %d\n",MPS->ProblemName,MPS->NumCols);
#endif

        WriteMPS(MPS,Inputs);
}


