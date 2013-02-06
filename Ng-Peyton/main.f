C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  February 13, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C
C       ****************************************************************
C       ****************************************************************
C
C       MAIN.F :  DRIVER THAT SOLVES A SPARSE, SYMMETRIC, POSITIVE 
C       DEFINITE LINEAR SYSTEMS VIA SPARSE CHOLESKY FACTORIZATION.  
C       THE MOST EFFICIENT KNOWN ALGORITHMS AND ROUTINES ARE INCLUDED:  
C       THE MULTIPLE MINIMUM DEGREE ROUTINES COME FROM THE MOST RECENT 
C       RELEASE OF SPARSPAK; MOST OF THE SYMBOLIC FACTORIZATION ROUTINES 
C       AND ALL OF THE NUMERICAL FACTORIZATION ROUTINES WERE DEVELOPED 
C       AT OAK RIDGE NATIONAL LABORATORY.
C
C       ****************************************************************
C       ****************************************************************
C
C       THIS DRIVER IS NOT A PRACTICAL INSTRUMENT FOR SOLVING A USER'S 
C       SPARSE SYMMETRIC POSITIVE DEFINITE LINEAR SYSTEM.  IN 
C       PARTICULAR, THE DRIVER DOES NOT ALLOCATE AND DEALLOCATE MEMORY 
C       IN AN EFFICIENT MANNER.  WE EMPHASIZE, HOWEVER, THAT THE 
C       ROUTINES THEMSELVES ARE DESIGNED TO ENABLE EFFICIENT MEMORY 
C       MANAGEMENT, BY EITHER DYNAMIC ALLOCATION OR BY ``MANUAL' 
C       ALLOCATION (AND DE-ALLOCATION) FROM A SINGLE WORKING ARRAY.
C
C       THIS DRIVER AND ITS COMMENTS
C
C       (1) INDICATE FOR EACH ROUTINE WHICH VARIABLES AND ARRAYS
C           CONSTITUTE ITS INPUT AND WHICH CONSTITUTE ITS OUTPUT, 
C
C       (2) POINT OUT WHEN THE CONTENTS OF AN ARRAY ARE NO LONGER 
C           NEEDED, 
C
C       (3) DEMONSTRATE APPROPRIATE ERROR TRAPPING, AND 
C
C       (4) GIVE ENOUGH INFORMATION TO ENABLE EFFICIENT MEMORY 
C           MANAGEMENT:  PRIOR TO EACH SUBROUTINE CALL, THE AMOUNT OF 
C           MEMORY REQUIRED FOR EACH ARRAY ARGUMENT IS KNOWN.  
C
C       NOTE ALSO THAT THE ROUTINES USED TO CREATE THE GRAPH OF THE 
C       COEFFCIENT MATRIX, TO INPUT THE NUMERICAL VALUES, AND TO 
C       GENERATE THE RIGHT-HAND SIDE ARE NOT GENERAL-PURPOSE TOOLS FOR 
C       THESE TASKS.  THEY ARE EXTREMELY SIMPLE ROUTINES THAT THE USER 
C       WILL NEED TO REPLACE WITH ROUTINES OF HIS/HER OWN.  THE DRIVER
C       READS IN THE STRUCTURE OF THE MATRIX FROM AN INPUT FILE, WHICH 
C       CONTAINS THE ENTIRE ADJACENCY STRUCTURE STORED IN AN EXTREMELY
C       SIMPLE FORMAT.
C
C       ****************************************************************
C       ****************************************************************
C
C       THE AMOUNT OF STORAGE REQUIRED TO SOLVE A LINEAR SYSTEM
C       DEPENDS ON THE FOLLOWING PARAMETERS, EACH OF WHICH IS EITHER
C       READ IN FROM A FILE OR COMPUTED DURING THE SOLUTION PROCESS.
C
C       (1) N      -- THE NUMBER OF ROWS (COLUMNS) IN A.
C
C       (2) NSUPER -- THE NUMBER OF SUPERODES (NSUPER <= N).
C
C       (3) NNZA   -- THE NUMBER OF NONZERO ENTRIES IN A, EXCLUDING
C                     ENTRIES ON THE MAIN DIAGONAL.
C
C       (4) NSUB   -- THE NUMBER OF ROW SUBSCRIPTS NEEDED TO 
C                     REPRESENT THE ZERO-NONZERO STRUCTURE OF L. 
C
C       (5) NNZL   -- THE NUMBER OF NONZERO ENTRIES IN L, INCLUDING
C                     ENTRIES ON THE MAIN DIAGONAL.
C
C       (6) TMPSIZ -- THE SIZE OF THE FLOATING POINT WORKING STORAGE 
C                     REQUIRED BY THE FACTORIZATION ROUTINE (BLKFCT).
C
C       WITH THE EXCEPTION OF THE IWORK(*) ARRAY, THE LENGTH OF EACH 
C       ARRAY SHOULD BE N+1, NSUPER+1, OR ONE OF THE PRECEDING SIX
C       PARAMETERS. IN THIS DRIVER THE LENGTH OF EACH ARRAYS IS
C       DECLARED AT RUN TIME, WHICH IS, OF COURSE, BEFORE THE
C       APPROPRIATE LENGTH IS KNOWN.  THE FOLLOWING PARAMETERS,
C       INITIALIZED IN THE CODE BELOW, ARE UPPER BOUNDS ON THE SIX
C       VARIABLES LISTED ABOVE AND ARE USED TO DECLARE THE LENGTH OF
C       ALL ARRAYS USED BY THE PROGRAM.
C
C       (1) NMAX    -   UPPER BOUND ON N AND NSUPER.
C       (2) ANZMAX  -   UPPER BOUND ON NNZA.
C       (3) SUBMAX  -   SIMULTANEOUSLY AN UPPER BOUND ON NSUB AND NNZA.
C       (4) LNZMAX  -   UPPER BOUND ON NNZL.
C       (5) TMPMAX  -   UPPER BOUND ON TMPSIZ.
C       (6) IWMAX   -   UPPER BOUND ON IWSIZ
C
C       THE AMOUNT OF INTEGER WORK SPACE (IWORK) REQUIRED BY EACH
C       ROUTINE VARIES FROM ONE ROUTINE TO ANOTHER.  THE MAXIMUM 
C       REQUIRED IS 7*N+3.
C
C       ****************************************************************
C
C       --------------------------------------------------------------
C       PARAMETERS.
C
C       NOTE:
C           IF YOUR FORTRAN COMPILER COMPLAINS ABOUT THE DEFINITION OF
C           SUBMAX, THEN CHANGE IT TO
C                           SUBMAX      = 61648,
C       --------------------------------------------------------------
        INTEGER             IWMAX , LNZMAX, NMAX  , ANZMAX, SUBMAX, 
     &                      TMPMAX
        PARAMETER       (   NMAX        = 1824,
     &                      ANZMAX      = 61648,
     &                      SUBMAX      = MAX(ANZMAX,17876),
     &                      LNZMAX      = 112267,
     &                      TMPMAX      = 9180,
     &                      IWMAX       = 7*NMAX + 3    )
C
C       ****************************************************************
C
C       THE FOLLOWING DECLARATIONS GIVE A ROUGH OVERVIEW OF THE 
C       LOGICAL FLOW THROUGH THE DRIVER AND THE KEY OUTPUT PRODUCED AT 
C       EACH STEP.
C
C       ****************************************************************
C
C       ---------------------------------------------
C       ADJACENCY STRUCTURE READ IN FROM MATRIX FILE.
C       ---------------------------------------------
        INTEGER             N     , NNZA
        INTEGER             XADJ(NMAX+1),
     &                      ADJ(ANZMAX)
C
C       ---------------------------------------------
C       ADJACENCY STRUCTURE WITH DIAGONAL ENTRIES AND
C       NUMERICAL VALUES ADDED (CREATE).
C       ---------------------------------------------
        INTEGER             XADJD(NMAX+1),
     &                      ADJD(ANZMAX+NMAX)
        DOUBLE PRECISION    ANZ(ANZMAX+NMAX)
C
C       -------------------------------------------------
C       SOLUTION AND ASSOCIATED RIGHT-HAND SIDE (GETRHS).
C       -------------------------------------------------
        DOUBLE PRECISION    RHS(NMAX),
     &                      SOL(NMAX)
C
C       ---------------------------------------
C       PERMUTATION VECTORS (ORDNAT OR ORDMMD).
C       ---------------------------------------
        INTEGER             PERM(NMAX),
     &                      INVP(NMAX)
C
C       ----------------------------------------------
C       SUPERNODE PARTITION AND ROW AND COLUMN NONZERO 
C       COUNTS FOR L (SFINIT).
C       ----------------------------------------------
        INTEGER             NNZL  , NSUB  , NSUPER
        INTEGER             COLCNT(NMAX), 
     &                      SNODE(NMAX),
     &                      XSUPER(NMAX+1)
C
C       ---------------------------------------------
C       COMPRESSED (SUPERNODAL) REPRESENTATION OF THE 
C       ZERO-NONZERO STRUCTURE OF L (SYMFCT).
C       ---------------------------------------------
        INTEGER             XLINDX(NMAX+1),
     &                      LINDX(SUBMAX),
     &                      XLNZ(NMAX+1)
C
C       ---------------------------------------------------------
C       SIZE OF FLOATING POINT WORKING STORAGE NEEDED BY BLKFCT;
C       SPLITTING OF SUPERNODES TO BETTER EXPLOIT CACHE (BFINIT).
C       ---------------------------------------------------------
        INTEGER             TMPSIZ
        INTEGER             SPLIT(NMAX)
C
C       -----------------------------------------------------------------
C       LNZ FIRST HAS SCATTERED INTO IT THE NONZERO ENTRIES OF A (INPNV);
C       THEN IT CONTAINS THE ENTRIES OF THE CHOLESKY FACTOR (BLKFCT).
C       -----------------------------------------------------------------
        DOUBLE PRECISION    LNZ(LNZMAX)
        DOUBLE PRECISION    TMPVEC(TMPMAX)
C
C       ----------------------------------------
C       GENERAL PURPOSE INTEGER WORKING STORAGE.
C       ----------------------------------------
        INTEGER             IWORK(IWMAX)
C
C       ---------------------------------------------------
C       EXTRA WORK VECTOR FOR RIGHT-HAND SIDE AND SOLUTION.
C       ---------------------------------------------------
        DOUBLE PRECISION    NEWRHS(NMAX)
C
C       ------------------------------
C       MISCELLANEOUS LOCAL VARIABLES.
C       ------------------------------
        CHARACTER*60        STRING
        CHARACTER*60        MATFIL, OUTFIL 
        CHARACTER*40        ORDRNG(2)
        INTEGER             CACHSZ, I     , ICASE , IFLAG , IWSIZ ,
     &                      LEVEL , MATUNT, OUTUNT
        REAL                GTIMER, TIMBEG, TIMEND, TIME
        DOUBLE PRECISION    E     , ERROR
C
C       -----------------------------------
C       EXTERNAL ROUTINES PASSED TO BLKFCT.
C       -----------------------------------
        EXTERNAL            SMXPY1, SMXPY2, SMXPY4, SMXPY8
        EXTERNAL            MMPY1 , MMPY2 , MMPY4 , MMPY8
C
C       ---------------------------------------------------------
C       THE FOLLOWING STATEMENT IS NEEDED FOR THOSE MACHINES THAT
C       ALLOCATE STORAGE FROM A STACK AT RUN TIME.
C       ---------------------------------------------------------
        COMMON /BIGMEM/     ADJ   , ADJD  , COLCNT, INVP  , IWORK ,
     &                      LINDX , PERM  , SNODE , SPLIT , XADJ  ,
     &                      XADJD , XLINDX, XLNZ  , XSUPER, ANZ   ,
     &                      LNZ   , NEWRHS, RHS   , SOL   , TMPVEC
C
        DATA  ORDRNG(1)     / 'NATURAL' /,
     &        ORDRNG(2)     / 'MULTIPLE MINIMUM DEGREE' /
C
C       ****************************************************************
C       ****************************************************************
C
C       ---------------------
C       SET I/O UNIT NUMBERS.
C       ---------------------
        OUTUNT = 6
        MATUNT = 10
C
C       *****************************************
C       READ IN PROBLEM DATA FROM STANDARD INPUT.
C       *****************************************
C
C       -------------------------------------------
C       READ NAME OF OUTPUT FILE AND OPEN THE FILE.
C       -------------------------------------------
        READ *, OUTFIL
        OPEN  ( UNIT=OUTUNT, FILE=OUTFIL )
        WRITE (OUTUNT,*)  'OUTPUT IS IN ', OUTFIL
C
C       -------------------------------------------
C       READ NAME OF MATRIX FILE AND OPEN THE FILE.
C       -------------------------------------------
        READ *, MATFIL
        WRITE (OUTUNT,*)   ' '
        WRITE (OUTUNT,*)   '   MATRIX FILE: ', MATFIL
        OPEN ( UNIT=MATUNT, FILE=MATFIL )
C
C       -----------------------------------------------
C       READ CHOICE OF ORDERING AND CHECK ITS VALIDITY.
C       -----------------------------------------------
        READ *, ICASE
        IF  ( ICASE .LE. 0  .OR.  ICASE .GE. 3 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,11) ICASE
   11       FORMAT ( '*** ORDERING OPTION ICASE = ', I5,' ***' /,
     &               '*** SHOULD HAVE 1 <= ICASE <= 2 ***' )
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,*)   ' '
        WRITE (OUTUNT,22) ICASE, ORDRNG(ICASE)
   22   FORMAT ( '   ORDERING OPTION: ', I2, ' - ', A40 )
C
C       ---------------------------------------------------------
C       READ MACHINE CACHE SIZE (IN KBYTES).
C
C       THE INPUT PARAMETER CACHSZ SHOULD BE SET TO THE SIZE OF 
C       THE CACHE (IN KILOBYTES) ON THE TARGET MACHINE.  FOR MOST 
C       MACHINES (SUCH AS SUN SPARCSTATIONS), IT IS PROBABLY 32 
C       OR 64.  FOR CRAYS, CACHSZ SHOULD BE SET TO 0.
C       ---------------------------------------------------------
        READ *, CACHSZ
        IF  ( CACHSZ .LT. 0 )  THEN
            CACHSZ = 0
        ENDIF
        WRITE (OUTUNT,*)   ' '
        WRITE (OUTUNT,*)   '   CACHE SIZE (IN KBYTES): ', CACHSZ
C
C       ----------------------------------------------------
C       READ LEVEL OF LOOP UNROLLING AND CHECK ITS VALIDITY.
C       ----------------------------------------------------
        READ *, LEVEL
        IF  ( LEVEL .NE. 1  .AND.  LEVEL .NE. 2  .AND.
     &        LEVEL .NE. 4  .AND.  LEVEL .NE. 8        )  THEN
            WRITE (OUTUNT,*)   ' '
            WRITE (OUTUNT,23)  LEVEL
   23       FORMAT ('*** LOOP UNROLLING LEVEL = ', I5,' ***'/,
     &              '*** SHOULD HAVE LEVEL = 1, 2, 4, OR 8 ***' )
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,24) LEVEL
   24   FORMAT ( '   LOOP UNROLLING LEVEL: ', I3 )
C
C       *****************
C       READ MATRIX FILE.
C       *****************
C
C       -----------------------
C       GET MATRIX DESCRIPTION.
C       -----------------------
        TIMBEG = GTIMER()
        READ (MATUNT,33) STRING
   33   FORMAT ( A60 )
C
C       -------------------------------------------
C       GET MATRIX SIZE ... N AND NNZA ARE READ IN.
C       -------------------------------------------
        READ (MATUNT,44)  N, NNZA
   44   FORMAT ( I6, I8 )
C
C       --------------------------------
C       STOP IF THE MATRIX IS TOO LARGE.
C       --------------------------------
        IF  ( N .GT. NMAX  .OR.  NNZA .GT. ANZMAX )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,55) STRING
   55       FORMAT ( A60 )
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,66) 
     &          'NUMBER OF EQUATIONS                   = ', N
   66       FORMAT ( A40, I10 )
            WRITE (OUTUNT,66) 
     &          'NUMBER OF NONZEROS (INCLUDING DIAG.)  = ', N+NNZA
            WRITE (OUTUNT,66) 
     &          'NUMBER OF NONZEROS (EXCLUDING DIAG.)  = ', NNZA
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** MATRIX IS TOO LARGE ***'
            IF  ( N .GT. NMAX )  THEN
                WRITE (OUTUNT,*) ' '
                WRITE (OUTUNT,*) '*** NUMBER OF EQUATIONS = ', N
                WRITE (OUTUNT,*) '*** IS LARGER THAN NMAX = ', NMAX
            ENDIF
            IF  ( NNZA .GT. ANZMAX )  THEN
                WRITE (OUTUNT,*) ' '
                WRITE (OUTUNT,*) 
     &              '*** NUMBER OF OFF-DIAGONAL NONZEROS (IN A) = ', 
     &              NNZA
                WRITE (OUTUNT,*) '*** IS LARGER THAN ANZMAX = ', 
     &                           ANZMAX
            ENDIF
            GO TO 1700
        ENDIF
C
C       --------------------------------------------------------
C       GET MATRIX STRUCTURE ... XADJ(*) AND ADJ(*) ARE READ IN.
C       --------------------------------------------------------
        READ (MATUNT,*) (XADJ(I),I=1,N+1)
        READ (MATUNT,*) (ADJ(I),I=1,NNZA)
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,55) STRING
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,66) 
     &      'NUMBER OF EQUATIONS                   = ', N
        WRITE (OUTUNT,66) 
     &      'NUMBER OF NONZEROS (INCLUDING DIAG.)  = ', N+NNZA
        WRITE (OUTUNT,66) 
     &      'NUMBER OF NONZEROS (EXCLUDING DIAG.)  = ', NNZA
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      'TIME FOR READING THE MATRIX           = ', TIME
   77   FORMAT ( A40, F10.3 )
C
C       *****************************************************
C       COMPUTE NUMERICAL ENTRIES OF THE MATRIX;
C       COMPUTE TRUE SOLUTION AND ASSOCIATED RIGHT-HAND SIDE.
C       *****************************************************
C
C       --------------------------------------------------------
C       THE PURPOSE OF CREATE IS TO PREPARE A DATA STRUCTURE
C       FOR COMPUTING THE RIGHT-HAND SIDE FROM A GIVEN SOLUTION,
C       AND FROM WHICH IT IS VERY SIMPLE TO LOAD THE ENTRIES OF
C       A INTO L.
C
C       CREATE : ADD DIAGONAL ENTRIES TO THE STUCTURE
C                AND ADD ALL NONZERO ENTRIES, TOO.
C
C       INPUT  : N, XADJ, ADJ
C       OUTPUT : XADJD, ADJD, ANZ
C       --------------------------------------------------------
        TIMBEG = GTIMER()
        CALL  CREATE (  N     , XADJ  , ADJ   , XADJD , ADJD  ,
     &                  ANZ                                     )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      'TIME FOR CREATING FULL REPRESENTATION = ', TIME
C
C       --------------------------------------------
C       GETRHS : CONSTRUCT RIGHT HAND SIDE VECTOR
C                ASSOCIATED WITH SOLUTION IN SOL(*).
C
C       INPUT  : N, XADJD, ADJD, ANZ, SOL
C       OUTPUT : RHS
C       --------------------------------------------
        TIMBEG = GTIMER()
        DO  100  I = 1, N
            SOL(I) = I
  100   CONTINUE
        CALL  GETRHS (  N     , XADJD , ADJD  , ANZ   , SOL   ,
     &                  RHS                                     )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      'TIME FOR CONSTRUCTING RHS             = ', TIME
C
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,*) 
     &      '------------------------------------------------'
C
C       *******************
C       REORDER THE MATRIX.
C       *******************
C
        IF  ( ICASE .EQ. 1 )  THEN
C
C           --------------------------------------------------------
C           ORDNAT : NATURAL ORDERING;
C                    I.E., COMPUTE AND USE THE IDENTITY PERMUTATION.
C
C           INPUT  : N
C           OUTPUT : INVP, PERM
C           --------------------------------------------------------
            TIMBEG = GTIMER()
            CALL  ORDNAT  (  N     , PERM  , INVP   )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,77) 
     &          '   TIME FOR ORDERING                  = ', TIME
C
        ELSEIF  ( ICASE .EQ. 2 )  THEN
C
C           -------------------------------------------------------
C           COPY MATRIX STRUCTURE FROM (XADJ,ADJ) TO (XLINDX,LINDX)
C           (BECAUSE MATRIX STRUCTURE IS DESTROYED BY THE MINIMUM
C           DEGREE ORDERING ROUTINE).
C           -------------------------------------------------------
            TIMBEG = GTIMER()
            DO  200  I = 1, N+1
                XLINDX(I) = XADJ(I)
  200       CONTINUE
            DO  300  I = 1, NNZA
                LINDX(I) = ADJ(I)
  300       CONTINUE
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,77) 
     &          '   TIME FOR COPYING ADJACENCY STRUCT. = ', TIME
C
C           ------------------------------------------------------
C           ORDMMD : MULTIPLE MINIMUM DEGREE ORDERING.
C
C           INPUT  : N, XLINDX, LINDX, IWSIZ
C           OUTPUT : PERM, INVP, NSUB, IFLAG
C           WORK   : IWORK(4*N)
C           ------------------------------------------------------
            TIMBEG = GTIMER()
            IWSIZ  = 4*N
            CALL  ORDMMD (  N     , XLINDX, LINDX , INVP  , PERM  ,
     &                      IWSIZ , IWORK , NSUB  , IFLAG           )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            IF  ( IFLAG .EQ. -1 )  THEN
                WRITE (OUTUNT,*) ' '
                WRITE (OUTUNT,*) '*** SIZE OF IWORK = ',
     &                           IWSIZ
                WRITE (OUTUNT,*) '*** IS LARGER THAN IWMAX = ',
     &                           IWMAX
                GO TO 1700
            ENDIF
            WRITE (OUTUNT,77) 
     &          '   TIME FOR ORDERING                  = ', TIME
C
        ENDIF
C
C       --------------------------------------------------------------
C       SFINIT: SYMBOLIC FACTORIZATION INITIALIZATION.
C               COMPUTE SUPERNODE PARTITION AND STORAGE REQUIREMENTS 
C               FOR SYMBOLIC FACTORIZATION.  NEW ORDERING IS A 
C               POSTORDERING OF THE NODAL ELIMINATION TREE.
C
C       INPUT  : N, NNZA, XADJ, ADJ, PERM, INVP, IWSIZ
C       OUTPUT : PERM, INVP, COLCNT, NNZL, NSUB, NSUPER, SNODE, XSUPER,
C                IFLAG
C       WORK   : IWORK(7*N+3) [MAXIMUM REQUIRED BY ANY ROUTINE.]
C       --------------------------------------------------------------
        TIMBEG = GTIMER()
        IWSIZ = 7 * N + 3
        CALL  SFINIT  ( N     , NNZA  , XADJ  , ADJ   , PERM  , 
     &                  INVP  , COLCNT, NNZL  , NSUB  , NSUPER, 
     &                  SNODE , XSUPER, IWSIZ , IWORK , IFLAG   )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. -1 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** SIZE OF IWORK = ',
     &                       IWSIZ
            WRITE (OUTUNT,*) '*** IS LARGER THAN IWMAX = ',
     &                       IWMAX
            GO TO 1700
        ENDIF
        IF  ( NNZL .GT. LNZMAX )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** NUMBER OF NONZEROS IN L = ', 
     &                       NNZL
            WRITE (OUTUNT,*) '*** IS LARGER THAN LNZMAX = ', 
     &                       LNZMAX
            GO TO 1700
        ENDIF
        IF  ( NSUB .GT. SUBMAX )  THEN
            WRITE (OUTUNT,*) '*** NUMBER OF FACTOR SUBSCRIPTS = ',
     &                       NSUB
            WRITE (OUTUNT,*) '*** IS LARGER THAN SUBMAX = ',
     &                       SUBMAX
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      '   TIME FOR SYMBOLIC FACT. SETUP      = ', TIME
C
C       ---------------------------------------------------------
C       SYMFCT : PERFORM SUPERNODAL SYMBOLIC FACTORIZATION.
C
C       INPUT  : N, NNZA, XADJ, ADJ, PERM, INVP, COLCNT, NSUPER,
C                XSUPER, SNODE , NSUB, IWSIZ
C       OUTPUT : XLINDX, LINDX, XLNZ, IFLAG
C       WORK   : IWORK(NSUPER+2*N)
C       NO LONGER NEEDED: ADJ, COLCNT, XADJ
C       ---------------------------------------------------------
        TIMBEG = GTIMER()
        IWSIZ = NSUPER + 2 * N + 1
        CALL  SYMFCT ( N     , NNZA  , XADJ  , ADJ   , PERM  , 
     &                 INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                 NSUB  , XLINDX, LINDX , XLNZ  , IWSIZ , 
     &                 IWORK , IFLAG                            )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. -1 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** SIZE OF IWORK = ',
     &                       IWSIZ
            WRITE (OUTUNT,*) '*** IS LARGER THAN IWMAX = ',
     &                       IWMAX
            GO TO 1700
        ENDIF
        IF  ( IFLAG .EQ. -2 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) 
     &          '*** INCONSISTENCY IN THE INPUT TO SYMFCT ***'
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,77) 
     &      '   TIME FOR SYMBOLIC FACTORIZATION    = ', TIME
C
C       ----------------------------------------------------
C       INPNV : INPUT NUMERICAL VALUES INTO DATA STRUCTURES.
C
C       INPUT  : XADJD, ADJD, ANZ, PERM, INVP, NSUPER,
C                XSUPER, XLINDX, LINDX, XLNZ
C       OUTPUT : LNZ
C       WORK   : IWORK
C       ----------------------------------------------------
        TIMBEG = GTIMER()
        IWSIZ = N
        CALL  INPNV  (  N     , XADJD , ADJD  , ANZ   , PERM  , 
     &                  INVP  , NSUPER, XSUPER, XLINDX, LINDX , 
     &                  XLNZ  , LNZ   , IWORK                   )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      '   TIME FOR NUMERICAL INPUT           = ', TIME
C
C       --------------------------------------------------------
C       BFINIT : INITIALIZATION FOR BLOCK FACTORIZATION
C
C       INPUT  : N, NSUPER, XSUPER, SNODE, XLINDX, LINDX, CACHSZ
C       OUTPUT : TMPSIZ, SPLIT
C       --------------------------------------------------------
        TIMBEG = GTIMER()
        CALL  BFINIT (  N     , NSUPER, XSUPER, SNODE , XLINDX, 
     &                  LINDX , CACHSZ, TMPSIZ, SPLIT           )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( TMPSIZ .GT. TMPMAX )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** WORKING SPACE REQUIREMENT = ', 
     &                       TMPSIZ
            WRITE (OUTUNT,*) '*** IS LARGER THAN TMPMAX = ', 
     &                       TMPMAX
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      '   TIME FOR FACTORIZATION INIT.       = ', TIME
C
C       -----------------------------------------------------------
C       BLKFCT : NUMERICAL FACTORIZATION.
C
C       INPUT  : NSUPER, XSUPER, SNODE, SPLIT, XLINDX, LINDX, XLNZ,
C                LNZ, IWSIZ, TMPSIZ, MMPY[LVL], SMXPY[LVL]
C       OUTPUT : LNZ, IFLAG
C       WORK   : IWORK(2*N+2*NSUPER), TMPVEC(TMPSIZ)
C       -----------------------------------------------------------
        TIMBEG = GTIMER()
        IWSIZ = 2 * N + 2 * NSUPER
        IF  ( LEVEL .EQ. 1 )  THEN
            CALL  BLKFCT (  N     , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ , 
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY1 , 
     &                      SMXPY1                                  )
        ELSEIF  ( LEVEL .EQ. 2 )  THEN
            CALL  BLKFCT (  N     , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ , 
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY2 , 
     &                      SMXPY2                                  )
        ELSEIF  ( LEVEL .EQ. 4 )  THEN
            CALL  BLKFCT (  N     , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ , 
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY4 , 
     &                      SMXPY4                                  )
        ELSEIF  ( LEVEL .EQ. 8 )  THEN
            CALL  BLKFCT (  N     , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ , 
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY8 , 
     &                      SMXPY8                                  )
        ENDIF
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. -1 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) 
     &          '*** MATRIX IS NOT POSITIVE DEFINITE ***'
            GO TO 1700
        ELSEIF ( IFLAG .EQ. -2 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) 
     &          '*** INSUFFICIENT WORK STORAGE [TMPVEC(*)] ***'
            GO TO 1700
        ELSEIF ( IFLAG .EQ. -3 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** SIZE OF IWORK = ',
     &                       IWSIZ
            WRITE (OUTUNT,*) '*** IS LARGER THAN IWMAX = ',
     &                       IWMAX
            GO TO 1700
        ENDIF
        WRITE (OUTUNT,77) 
     &      '   TIME FOR NUMERICAL FACTORIZATION   = ', TIME
C
C       ------------------------------------------------------
C       BLKSLV : NUMERICAL SOLUTION.
C
C       INPUT  : NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ, RHS
C       OUTPUT : RHS
C       ------------------------------------------------------
        TIMBEG = GTIMER()
        DO  1300  I = 1, N
            NEWRHS(I) = RHS(PERM(I))
 1300   CONTINUE
        CALL  BLKSLV (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                  LNZ   , NEWRHS                          )
        DO  1400  I = 1, N
            SOL(I) = NEWRHS(INVP(I))
 1400   CONTINUE
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      '   TIME FOR TRIANGULAR SOLUTIONS      = ', TIME
C
C       ------------------------------------
C       COMPUTE AND PRINT ERROR IN SOLUTION.
C       ------------------------------------
        TIMBEG = GTIMER()
        ERROR = 0.0
        DO  1500  I = 1, N
            E = ABS(SOL(I) - I)
            IF  ( E .GT. ERROR )  ERROR = E
 1500   CONTINUE
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,77) 
     &      '   TIME FOR COMPUTING ERROR           = ', TIME
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,88) 
     &      '   MAXIMUM ABSOLUTE ERROR             = ', ERROR
        WRITE (OUTUNT,88) 
     &      '   MAXIMUM RELATIVE ERROR             = ',
     &      ERROR/N
   88   FORMAT ( A40, 1PD20.10 )
C
C       ------------------------------------------------------------
C       COMPUTE AND PRINT STATISTICS.
C
C       INPUT  : NSUPER, XSUPER, XLINDX, LINDX, XLNZ, TMPSIZ, OUTUNT
C       OUTPUT : PRINTED STATISTICS
C       ------------------------------------------------------------
        CALL  LSTATS (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                  TMPSIZ, OUTUNT                          )
C
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,*) 
     &      '------------------------------------------------'
C
C       -------------------
C       NORMAL TERMINATION.
C       -------------------
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,*) 'NORMAL TERMINATION.'
        STOP
C
C       ---------------------
C       ABNORMAL TERMINATION.
C       ---------------------
 1700   CONTINUE
        WRITE (OUTUNT,*) ' '
        WRITE (OUTUNT,*) 'ABNORMAL TERMINATION.'
        STOP
C
      END
