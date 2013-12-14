c
c
c
c        Subroutine CCFeta
c
c        Compute a Controlled Cholesky Factorisation
c
      SUBROUTINE CCF( NA, NZA, A, AColPt, ARowId, NZL, L, LColPt
     t, LRowId, Eta, DimL, Eps8, LinRow, LogRow, RowCol, T, W
     u, Info )

      integer AColPt(*), ARowId(*), DimL, Eta, Info, LColPt(*)
     j, LinRow(*), LRowId(*), NA, NZA, NZL, RowCol(*)
      real*8 A(*), Eps8, L(*), T(*), W(*)
      logical*1 LogRow(*)
c
c        parameter    type      in/out  description
c
c        NA           integer   entry   order of A
c        NZA          integer   entry   number of nonzero elements of A
c        A(NZA)       real*8    entry   elements of matrix A
c        AColPt(NA+1) integer   entry   location of first entries of column of A
c        ARowId(NZA)  integer   entry   row indices of A
c        NZL          integer   exit    number of nonzero elements of L
c        L(DimL)      real*8    exit    elements of L
c        LColPt(NA+1) integer   exit    location of 1st entries of columns of L
c        LRowId(DimL) integer   exit    row indices of L
c        Eta          integer   entry   number of extra elements per column
c        DimL         integer   entry   maximum size of L
c        Eps8         real*8    entry   machine precision
c        LinRow(DimL) integer   exit    pointer for next element in the row
c        LogRow(NA)   logical*1 exit    existence of row elements
c        RowCol(3*NA) integer   exit    workspace
c        T(NA)        real*8    exit    workspace
c        W(NA)        real*8    exit    full column of L
c        Info         integer   exit    error ( <0 ) or warning ( >0 ) :
c                      0 : no error
c                     -3 : D(I,I) have to be modified more than MaxShi times
c                          where MaxShi = 15
c                      k : D(I,I) was modified k times
c
c        Local variables
c
      integer Avail, Col, FNElCo, I, J, K, LimCol, LP, LU, L1, L2, M
     j, MaxShi, NElCol, NElRow, N2
      real*8 Aux, RDiag, Sum, Shift
      logical Fail
c
c        External subprograms
c
      integer CCFind
      external CCFind, CCFsel
c
c        Intrinsic functions
c
      integer max0
      real*8 dfloat, dsqrt
c
      data MaxShi / 15 /
c



c          print*,'NA, NZA',NA,NZA
c           write(*,1),(AColPT(i),i=1,NA+1)
c           write(*,1),(ARowId(i),i=1,NZA)
c           PRINT*,'A'
c           write(*,2),(A(i),i=1,NZA)
         PRINT*,'Eta,DimL',Eta,DimL

  1        format(1x,i5)
  2        format(1pe30.16)



      fator = float(NA*Eta)/float(NZA-NA) + 1.0
      write(*,106)
      Shift = 0.0d0
      Info = 0
   10 continue
      Fail = .false.
c
c        Initialising data structures
c
      LP = 0
      LColPt(1) = 1
      N2 = 2 * NA
      do K = 1, NA
          RowCol(N2+K) = 0
          W(K) = 0.0d0
          LogRow(K) = .false.
      enddo
c
c        Performing Controlled Cholesky Factorisation
c
      do J = 1, NA
         NElRow = 0
         K = RowCol(N2+J)
         Col = J
         Sum = 0.0d0
c
c        Defining number of elements per row
c
         do while( K.gt.0 )
            Col = CCFind( K, 1, Col-1, NA, LColPt )
            NElRow = NElRow + 1
            RowCol(NA+NElRow) = Col
            T(Col) = L(K) * L(LColPt(Col))
            Sum = Sum + L(K) * T(Col)
            K = LinRow(K)
         enddo
c
c        Performing factorisation using all elements available
c


         L1 = AColPt(J)
         L2 = AColPt(J+1) - 1
         NElCol = L2 - L1
         do K = 1, NElCol
            I = ARowId(L1+K)
            RowCol(K) = I
            LogRow(I) = .true.
            W(I) = A(L1+K)
         enddo
         LU = RowCol(N2+J)
         do K = 1, NElRow
            Col = RowCol(NA+K)
            LimCol = LColPt(Col+1) - 1
            do M = LU+1, LimCol
               I = LRowId(M)
               W(I) = W(I) - L(M) * T(Col)
               if( .not.LogRow(I) ) then
                  LogRow(I) = .true.
                  NElCol = NElCol + 1
                  RowCol(NElCol) = I
               endif
            enddo
            LU = LinRow(LU)
         enddo
c
c        Selecting the largest elements in current column
c
         if( Eta.ge.0 ) then
            Avail = max0(L2-L1 + Eta,0)
         else
            Avail = max0(int(float(L2-L1)*fator),0)
         endif
         if( Avail.lt.NElCol ) then
            FNElCo = Avail
            call CCFsel( NA, FNElCo, NElCol, 2, RowCol, W )
         else
            FNElCo = NElCol
            call CCFsel( NA, FNElCo, NElCol, 1, RowCol, W )
         endif
c
c        Updating data structures
c
         LP = LP + 1
         LRowId(LP) = J
         Aux = 1.0d0 + Shift - Sum
         if( Aux.gt.Eps8 ) then
            L(LP) = Aux
         else
            Info = Info + 1
            Shift = 5.0d-4*2.0d0**dfloat(Info-1)
            Fail = .true.
            go to 20
         endif
         LinRow(LP) = 0
         RDiag = 1.0d0 / L(LP)
         do K = 1, FNElCo
            I = RowCol(K)
            LP = LP + 1
            LRowId(LP) = I
            L(LP) = W(I) * RDiag
            LinRow(LP) = RowCol(N2+I)
            RowCol(N2+I) = LP
            LogRow(I) = .false.
            W(I) = 0.0d0
         enddo
         LColPt(J+1) = LColPt(J) + FNElCo + 1
         do K = FNElCo+1, NElCol
            I = RowCol(K)
            LogRow(I) = .false.
            W(I) = 0.0d0
         enddo
      enddo
   20 continue
      if( Fail.and.Info.le.MaxShi ) go to 10
      if( Fail ) Info = -3
      NZL = LP

c
c        Substituting diagonal factor D of L.D.transpose(L) by its square root

c      do J = 1, NA
c         I = LColPt(J)
c         L(I) = dsqrt( L(I) )
c      enddo



          print*,'NZL, Info',NZL,Info
c           write(*,1),(LColpt(i),i=1,NA+1)
c           write(*,1),(LRowId(i),i=1,NZL)
c           PRINT*,'L'
c           write(*,2),(L(i),i=1,NZL)





      return
  106 format('Performing Controlled Cholesky Factorisation . . .')
      end
c
c
c
c        Function CCFind
c
c        Find an index in a vector by the binary search scheme
c
      integer function CCFind( I, LA, LB, NA, V )
c
      integer I, LA, LB, NA, V(*)
c
c        parameter    type      in/out  description
c
c        I            integer   entry   index to be found
c        LA           integer   entry   lower limit of index
c        LB           integer   entry   upper limit of index
c        NA           integer   entry   order of A or L
c        V(NA+1)      integer   entry   column of matrix L
c
      integer L, LL, UL
c
      LL = LA
      UL = LB
   10 continue
         L = ( LL + UL ) / 2
         if( V(L).lt.I ) then
            LL = L + 1
         else
            UL = L - 1
         endif
      if( V(L).gt.I.or.V(L+1).le.I ) go to 10
      CCFind = L
      return
      end
c
c
c
c        Subroutine CCFsel
c
c        Select the largest elements in the current column W of L
c        by the straight selection scheme
c
      subroutine CCFsel( NA, FNElCo, NElCol, Tp, VetRow, W )
c
      integer NA, FNElCo, NElCol, Tp, VetRow(*)
      real*8 W(*)
c
c        parameter    type      in/out  description
c
c        NA           integer   entry   order of A or L
c        FNElCo       integer   entry   final number of elements in column W
c        NElCol       integer   entry   number of elements in column W
c        Tp           integer   entry   type of selection
c        VetRow(2*NA) integer   entry   row indices
c        W(NA)        real*8    entry   current column of L
c
c        Local variables
c
      integer Iaux, Imax, Imin, J, K, M
      real*8 Aux, Max, Min
c
c        Intrinsic function
c
      real*8 dabs
c
      if( Tp.eq.2 ) then
c
c        Selecting the largest elements in W
c
         if( FNElCo.le.NElCol-FNElCo) then
            do J = 1, FNElCo
               Imax = J
               Max = dabs( W(VetRow(Imax)) )
               do K = J+1, NElCol
                  Aux = dabs( W(VetRow(K)) )
                  if( Aux.gt.Max ) then
                     Max = Aux
                     Imax = K
                  endif
               enddo
               Iaux = VetRow(J)
               VetRow(J) = VetRow(Imax)
               VetRow(Imax) = Iaux
            enddo
         else
            do J = 1, NElCol-FNElCo
               Imin = NElCol + 1 - J
               Min = dabs( W(VetRow(Imin)) )
               do K = 1, NElCol-J
                  Aux = dabs( W(VetRow(K)) )
                  if( Aux.lt.Min ) then
                     Min = Aux
                     Imin = K
                  endif
               enddo
               Iaux = VetRow(NElCol+1-J)
               VetRow(NElCol+1-J) = VetRow(Imin)
               VetRow(Imin) = Iaux
            enddo
         endif
      endif
c
c        Ordering indices of W
c
      do J = 1, FNElCo
         Imin = J
         M = VetRow(Imin)
         do K = J+1, FNElCo
            if( VetRow(K).lt.M ) then
               M = VetRow(K)
               Imin = K
            endif
         enddo
         Iaux = VetRow(J)
         VetRow(J) = VetRow(Imin)
         VetRow(Imin) = Iaux
      enddo

      return
      end
