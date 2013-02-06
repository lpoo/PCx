C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C****     ORDMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE CALLS LIU'S MULTIPLE MINIMUM DEGREE
C               ROUTINE.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        IWSIZ  - SIZE OF INTEGER WORKING STORAGE.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
C                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
C        IFLAG  - ERROR FLAG.
C                   0: SUCCESSFUL ORDERING
C                  -1: INSUFFICIENT WORKING STORAGE
C                      [IWORK(*)].
C
C     WORKING PARAMETERS -
C        IWORK  - INTEGER WORKSPACE OF LENGTH 4*NEQNS.
C
C***********************************************************************
C
      SUBROUTINE ORDMMD  (  NEQNS , XADJ  , ADJNCY, INVP  , PERM  ,
     1                      IWSIZ , IWORK , NOFSUB, IFLAG           )
C
C***********************************************************************
C
         INTEGER    ADJNCY(1), INVP(1)  , IWORK(1) , PERM(1)
         INTEGER    XADJ(1)
         INTEGER    DELTA , IFLAG , IWSIZ , MAXINT, NEQNS , 
     &              NOFSUB
C
C*********************************************************************
C
        IFLAG = 0
        IF  ( IWSIZ .LT. 4*NEQNS )  THEN
            IFLAG = -1
            RETURN
        ENDIF
C
C       DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C       MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                NODES.
C
        DELTA  = 0
        MAXINT = 32767
        CALL GENMMD  (  NEQNS , XADJ  , ADJNCY, INVP  , PERM  ,
     1                  DELTA , 
     1                  IWORK(1)              ,
     1                  IWORK(NEQNS+1)        ,
     1                  IWORK(2*NEQNS+1)      ,
     1                  IWORK(3*NEQNS+1)      ,
     1                  MAXINT, NOFSUB          )
         RETURN
C
      END
