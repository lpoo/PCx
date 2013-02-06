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
C************     ASSMB .... INDEXED ASSEMBLY OPERATION     ************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS ROUTINE PERFORMS AN INDEXED ASSEMBLY (I.E., SCATTER-ADD)
C       OPERATION, ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE
C       CHOLESKY CODES.
C
C   INPUT PARAMETERS:
C       M               -   NUMBER OF ROWS IN Y.
C       Q               -   NUMBER OF COLUMNS IN Y.
C       Y               -   BLOCK UPDATE TO BE INCORPORATED INTO FACTOR
C                           STORAGE.
C       RELIND          -   RELATIVE INDICES FOR MAPPING THE UPDATES
C                           ONTO THE TARGET COLUMNS.
C       XLNZ            -   POINTERS TO THE START OF EACH COLUMN IN THE
C                           TARGET MATRIX.
C
C   OUTPUT PARAMETERS:
C       LNZ             -   CONTAINS COLUMNS MODIFIED BY THE UPDATE
C                           MATRIX.
C
C***********************************************************************
C
      SUBROUTINE  ASSMB  (  M     , Q     , Y     , RELIND, XLNZ  ,
     &                      LNZ   , LDA                             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER             LDA   , M     , Q
        INTEGER             XLNZ(*)
        INTEGER             RELIND(*)
        DOUBLE PRECISION    LNZ(*)        , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             ICOL  , IL1   , IR    , IY1   , LBOT1 ,
     &                      YCOL  , YOFF1
C
C***********************************************************************
C
C
        YOFF1 = 0
        DO  200  ICOL = 1, Q
            YCOL = LDA - RELIND(ICOL)
            LBOT1 = XLNZ(YCOL+1) - 1
CDIR$ IVDEP
            DO  100  IR = ICOL, M
                IL1 = LBOT1 - RELIND(IR)
                IY1 = YOFF1 + IR
                LNZ(IL1) = LNZ(IL1) + Y(IY1)
                Y(IY1) = 0.0D0
  100       CONTINUE
            YOFF1 = IY1 - ICOL
  200   CONTINUE
C
      RETURN
      END
