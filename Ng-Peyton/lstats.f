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
C
C     -----------------------------------------
C     GATHER STATISTICS ABOUT FACTORIZATION ...
C     -----------------------------------------
C
      SUBROUTINE  LSTATS (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                      TMPSIZ, OUTUNT                          )
C
        INTEGER             NSUPER, OUTUNT, TMPSIZ
        INTEGER             XSUPER(*), XLINDX(*), LINDX(*), XLNZ(*)
C
        INTEGER             J     , JLEN  , JSIZE , JSUPER, MAXSUP,
     &                      N     , NCOLS , NOFNZ , NOFSUB, SUPSZE
        DOUBLE PRECISION    FCTOPS, SLVOPS
C
        N = XSUPER(NSUPER+1) - 1
C
        WRITE (OUTUNT,*) ' '
C       -------------------------------------------------------
C       DETERMINE THE NUMBER OF NONZEROS IN CHOLESKY FACTOR AND
C       THE NUMBER OF SUBSCRIPTS IN REPRESENTING THE SUPERNODAL
C       STRUCTURE.
C       -------------------------------------------------------
        NOFNZ = XLNZ(N+1) - 1
        NOFSUB = XLINDX(NSUPER+1) - 1
        WRITE (OUTUNT,1) 
     &      '   NUMBER OF SUPERNODES               = ', NSUPER
        WRITE (OUTUNT,1) 
     &      '   NUMBER OF NONZEROS IN L            = ', NOFNZ
        WRITE (OUTUNT,1) 
     &      '   NUMBER OF SUBSCRIPTS IN L          = ', NOFSUB
C
C       -------------------------------------------------------
C       DETERMINE THE LARGEST SUPERNODE IN THE CHOLESKY FACTOR.
C       -------------------------------------------------------
        MAXSUP = 0
        SUPSZE = 0
        DO  100  JSUPER = 1, NSUPER
C           ---------------------------------------------------
C           NCOLS IS THE NUMBER OF COLUMNS IN SUPERNODE JSUPER.
C           ---------------------------------------------------
            NCOLS = XSUPER(JSUPER+1) - XSUPER(JSUPER)
            IF  ( NCOLS .GT. MAXSUP )  MAXSUP = NCOLS
C
C           ----------------------------------------------------
C           JSIZE IS THE NUMBER OF NONZEROS IN SUPERNDOE JSUPER.
C           ----------------------------------------------------
            JLEN = XLINDX(JSUPER+1) - XLINDX(JSUPER)
            JSIZE = ((2*JLEN - NCOLS + 1)*NCOLS)/2
            IF  ( JSIZE .GT. SUPSZE )  SUPSZE = JSIZE
  100   CONTINUE
        WRITE (OUTUNT,1) 
     &      '   LARGEST SUPERNODE BY COLUMNS       = ', MAXSUP
        WRITE (OUTUNT,1) 
     &      '   LARGEST SUPERNODE BY NONZEROS      = ', SUPSZE
C
        WRITE (OUTUNT,1) 
     &      '   SIZE OF TEMPORARY WORK STORAGE     = ', TMPSIZ
C
C       ---------------------------
C       DETERMINE OPERATION COUNTS.
C       ---------------------------
        SLVOPS = 0.0
        FCTOPS = 0.0
        DO  400  J = 1, N
            JLEN = XLNZ(J+1) - XLNZ(J)
            SLVOPS = SLVOPS + 2*JLEN - 1
            FCTOPS = FCTOPS + JLEN**2 - 1
  400   CONTINUE
        SLVOPS = 2*SLVOPS
        WRITE (OUTUNT,2) 
     &      '   FACTORIZATION OPERATION COUNT      = ', FCTOPS
        WRITE (OUTUNT,2) 
     &      '   TRIANGULAR SOLN OPERATION COUNT    = ', SLVOPS
C
    1   FORMAT ( A40, I10 )
    2   FORMAT ( A40, 1PD20.10 )
        RETURN
      END
