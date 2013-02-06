C***********************************************************************
C***********************************************************************
C
C   Written:        October 6, 1996 by SJW. Based on routine BLKSLV of
C                   Esmond G. Ng and Barry W. Peyton.
C
C   Modified:       Sept 30, 1999 to improve efficiency in the case
C                   in which the right-hand side and solution are both
C                   expected to be sparse. Happens a lot in "dense"
C                   column handling.
C
C***********************************************************************
C***********************************************************************
C*********     BLKSLF ... FORWARD TRIANGULAR SUBSTITUTION     **********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
C       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE
C       FORWARD TRIANGULAR SUBSTITUTIOn.  IT USES OUTPUT FROM BLKFCT.
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
C       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
C
C   UPDATED PARAMETERS:
C       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
C                           OUTPUT, CONTAINS THE SOLUTION.
C
C***********************************************************************
C
      SUBROUTINE  BLKSLF (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , RHS                             )
C
C***********************************************************************
C
        INTEGER             NSUPER
        INTEGER             LINDX(*)      , XSUPER(*)
        INTEGER             XLINDX(*)     , XLNZ(*)
        DOUBLE PRECISION    LNZ(*)        , RHS(*)
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        DOUBLE PRECISION    T
C
C***********************************************************************
C
        IF  ( NSUPER .LE. 0 )  RETURN
C
C       ------------------------
C       FORWARD SUBSTITUTION ...
C       ------------------------
        FJCOL = XSUPER(1)
        DO  300  JSUP = 1, NSUPER
            LJCOL  = XSUPER(JSUP+1) - 1
            IXSTRT = XLNZ(FJCOL)
            JPNT   = XLINDX(JSUP)
            DO  200  JCOL = FJCOL, LJCOL
                IXSTOP    = XLNZ(JCOL+1) - 1

                IF(RHS(JCOL) .NE. 0.D0) THEN
                   T         = RHS(JCOL)/LNZ(IXSTRT)
                   RHS(JCOL) = T
                   IPNT      = JPNT + 1
CDIR$           IVDEP
                   DO  100  IX = IXSTRT+1, IXSTOP
                      I      = LINDX(IPNT)
                      RHS(I) = RHS(I) - T*LNZ(IX)
                      IPNT   = IPNT + 1
 100               CONTINUE
                ENDIF

                IXSTRT = IXSTOP + 1
                JPNT   = JPNT + 1
  200       CONTINUE
            FJCOL = LJCOL + 1
  300   CONTINUE
C
        RETURN
      END
