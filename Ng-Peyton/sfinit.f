C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  January 12, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C**************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP 
C       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION.
C
C   NOTE:
C       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E.,
C       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES).
C
C   INPUT PARAMETERS:
C       NEQNS       -   NUMBER OF EQUATIONS.
C       NNZA        -   LENGTH OF ADJACENCY STRUCTURE.
C       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
C                       TO THE ADJACENCY STRUCTURE.
C       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING
C                       THE ADJACENCY STRUCTURE.
C       PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                       POSTORDERING.
C       INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                       INVERSE OF THE POSTORDERING.
C       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE.
C
C   OUTPUT PARAMETERS:
C       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                       OF NONZEROS IN EACH COLUMN OF THE FACTOR,
C                       INCLUDING THE DIAGONAL ENTRY.
C       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING
C                       THE DIAGONAL ENTRIES.
C       NSUB        -   NUMBER OF SUBSCRIPTS.
C       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
C       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                       SUPERNODE MEMBERSHIP.
C       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE
C                       SUPERNODE PARTITIONING.
C       IFLAG(*)    -   ERROR FLAG.
C                          0: SUCCESSFUL SF INITIALIZATION.
C                         -1: INSUFFICENT WORKING STORAGE
C                             [IWORK(*)].
C
C   WORK PARAMETERS:
C       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3.
C
C   FIRST CREATED ON    NOVEMEBER 14, 1994.
C   LAST UPDATED ON     January 12, 1995.
C
C***********************************************************************
C
      SUBROUTINE  SFINIT (  NEQNS , NNZA  , XADJ  , ADJNCY, PERM  ,
     &                      INVP  , COLCNT, NNZL  , NSUB  , NSUPER,
     &                      SNODE , XSUPER, IWSIZ , IWORK , IFLAG   )
C
C       ----------- 
C       PARAMETERS.  
C       -----------
        INTEGER             IFLAG , IWSIZ , NNZA  , NEQNS , NNZL  , 
     &                      NSUB  , NSUPER
        INTEGER             ADJNCY(NNZA)    , COLCNT(NEQNS)   ,
     &                      INVP(NEQNS)     , IWORK(7*NEQNS+3),
     &                      PERM(NEQNS)     , SNODE(NEQNS)    , 
     &                      XADJ(NEQNS+1)   , XSUPER(NEQNS+1)
C
C***********************************************************************
C
C       --------------------------------------------------------
C       RETURN IF THERE IS INSUFFICIENT INTEGER WORKING STORAGE.
C       --------------------------------------------------------
        IFLAG = 0
        IF  ( IWSIZ .LT. 7*NEQNS+3 )  THEN
            IFLAG = -1
            RETURN
        ENDIF
C
C       ------------------------------------------
C       COMPUTE ELIMINATION TREE AND POSTORDERING.
C       ------------------------------------------
        CALL  ETORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  ,
     &                  IWORK(1)              ,
     &                  IWORK(NEQNS+1)        ,
     &                  IWORK(2*NEQNS+1)      ,
     &                  IWORK(3*NEQNS+1)        ) 
C
C       ---------------------------------------------
C       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS.
C       ---------------------------------------------
        CALL  FCNTHN (  NEQNS , NNZA  , XADJ  , ADJNCY, PERM  , 
     &                  INVP  , IWORK(1)      , SNODE , COLCNT, 
     &                  NNZL  ,
     &                  IWORK(NEQNS+1)        ,
     &                  IWORK(2*NEQNS+1)      ,
     &                  XSUPER                ,
     &                  IWORK(3*NEQNS+1)      ,
     &                  IWORK(4*NEQNS+2)      ,
     &                  IWORK(5*NEQNS+3)      ,
     &                  IWORK(6*NEQNS+4)        )
C
C       ---------------------------------------------------------
C       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM 
C       NUMBER OF NONZEROS IN ITS COLUMN OF L.
C       ---------------------------------------------------------
        CALL  CHORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  ,
     &                  COLCNT, 
     &                  IWORK(1)              ,
     &                  IWORK(NEQNS+1)        ,
     &                  IWORK(2*NEQNS+1)      ,
     &                  IWORK(3*NEQNS+1)        )
C
C       ----------------
C       FIND SUPERNODES.
C       ----------------
        CALL  FSUP1  (  NEQNS , IWORK(1)      , COLCNT, NSUB  , 
     &                  NSUPER, SNODE                           )
        CALL  FSUP2  (  NEQNS , NSUPER, IWORK(1)      , SNODE ,
     &                  XSUPER                                  )
C
        RETURN
      END
