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
C*************     SYMFCT ..... SYMBOLIC FACTORIZATION    **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
C       THIS ROUTINE CALLS SYMFC2 WHICH PERFORMS SUPERNODAL SYMBOLIC
C       FACTORIZATION ON A REORDERED LINEAR SYSTEM.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS
C       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST.
C       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS
C                           TO THE ADJACENCY STRUCTURE.
C       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING
C                           THE ADJACENCY STRUCTURE.
C       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           POSTORDERING.
C       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           INVERSE OF THE POSTORDERING.
C       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) NSUPER      -   NUMBER OF SUPERNODES.
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           FIRST COLUMN OF EACH SUPERNODE.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN
C                           LINDX(*).
C       (I) IWSIZ       -   SIZE OF INTEGER WORKING STORAGE.
C
C   OUTPUT PARAMETERS:
C       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS 
C                           INTO THE SUBSCRIPT VECTOR.
C       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE
C                           COMPRESSED SUBSCRIPTS.
C       (I) XLNZ        -   COLUMN POINTERS FOR L.
C       (I) FLAG        -   ERROR FLAG:
C                               0 - NO ERROR.
C                              -1 - INSUFFICIENT INTEGER WORKING SPACE.
C                              -2 - INCONSISTANCY IN THE INPUT.
C       
C   WORKING PARAMETERS:
C       (I) IWORK       -   WORKING ARRAY OF LENGTH NSUPER+2*NEQNS.
C
C***********************************************************************
C
      SUBROUTINE  SYMFCT (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  , 
     &                      INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                      NOFSUB, XLINDX, LINDX , XLNZ  , IWSIZ ,
     &                      IWORK ,
     &                      FLAG    )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             ADJLEN, FLAG  , IWSIZ , NEQNS , NOFSUB, 
     &                      NSUPER
        INTEGER             ADJNCY(ADJLEN), COLCNT(NEQNS) ,
     &                      INVP(NEQNS)   , 
     &                      IWORK(NSUPER+2*NEQNS+1),
     &                      LINDX(NOFSUB) , 
     &                      PERM(NEQNS)   , SNODE(NEQNS)  , 
     &                      XSUPER(NSUPER+1)
        INTEGER             XADJ(NEQNS+1) , XLINDX(NSUPER+1),
     &                      XLNZ(NEQNS+1)
C
C***********************************************************************
C
        FLAG = 0
        IF  ( IWSIZ .LT. NSUPER+2*NEQNS+1 )  THEN
            FLAG = -1
            RETURN
        ENDIF
        CALL  SYMFC2 (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  , 
     &                  INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                  NOFSUB, XLINDX, LINDX , XLNZ  , 
     &                  IWORK(1)              ,
     &                  IWORK(NSUPER+1)       ,
     &                  IWORK(NSUPER+NEQNS+2) ,
     &                  FLAG    )
        RETURN
      END
