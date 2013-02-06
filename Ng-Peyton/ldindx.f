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
C******         LDINDX .... LOAD INDEX VECTOR             **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE COMPUTES THE SECOND INDEX VECTOR
C               USED TO IMPLEMENT THE DOUBLY-INDIRECT SAXPY-LIKE
C               LOOPS THAT ALLOW US TO ACCUMULATE UPDATE 
C               COLUMNS DIRECTLY INTO FACTOR STORAGE.
C
C     INPUT PARAMETERS -
C        JLEN   - LENGTH OF THE FIRST COLUMN OF THE SUPERNODE,
C                 INCLUDING THE DIAGONAL ENTRY.
C        LINDX  - THE OFF-DIAGONAL ROW INDICES OF THE SUPERNODE, 
C                 I.E., THE ROW INDICES OF THE NONZERO ENTRIES
C                 LYING BELOW THE DIAGONAL ENTRY OF THE FIRST
C                 COLUMN OF THE SUPERNODE.
C
C     OUTPUT PARAMETERS - 
C        INDMAP - THIS INDEX VECTOR MAPS EVERY GLOBAL ROW INDEX
C                 OF NONZERO ENTRIES IN THE FIRST COLUMN OF THE 
C                 SUPERNODE TO ITS POSITION IN THE INDEX LIST 
C                 RELATIVE TO THE LAST INDEX IN THE LIST.  MORE
C                 PRECISELY, IT GIVES THE DISTANCE OF EACH INDEX
C                 FROM THE LAST INDEX IN THE LIST.
C
C***********************************************************************
C
      SUBROUTINE  LDINDX ( JLEN, LINDX, INDMAP )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
      INTEGER             JLEN
      INTEGER             LINDX(*), INDMAP(*)
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
      INTEGER             CURLEN, J, JSUB
C
C***********************************************************************
C

      CURLEN = JLEN

      DO  200  J = 1, JLEN
          JSUB = LINDX(J)
          CURLEN = CURLEN - 1
          INDMAP(JSUB) = CURLEN
  200 CONTINUE
      RETURN
      END
