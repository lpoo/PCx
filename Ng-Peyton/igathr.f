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
C******         IGATHR .... INTEGER GATHER OPERATION      **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS A STANDARD INTEGER GATHER
C               OPERATION.
C
C     INPUT PARAMETERS -
C        KLEN   - LENGTH OF THE LIST OF GLOBAL INDICES.
C        LINDX  - LIST OF GLOBAL INDICES.
C        INDMAP - INDEXED BY GLOBAL INDICES, IT CONTAINS THE
C                 REQUIRED RELATIVE INDICES.
C
C     OUTPUT PARAMETERS - 
C        RELIND - LIST RELATIVE INDICES.
C
C***********************************************************************
C
      SUBROUTINE  IGATHR ( KLEN  , LINDX, INDMAP, RELIND )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
      INTEGER             KLEN  
      INTEGER             INDMAP(*), LINDX (*), RELIND(*)
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
      INTEGER             I
C
C***********************************************************************
C
CDIR$ IVDEP
      DO  100  I = 1, KLEN  
          RELIND(I) = INDMAP(LINDX(I))
  100 CONTINUE
      RETURN
      END
