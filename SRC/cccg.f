c
c        CCCG(eta)
c
      program CCCG_2004_0624
c
      integer DimA, DimL, DimSol, EtaMax, IteMax
      parameter ( DimA = 220000, DimSol  = 15000, EtaMax = 500
     q, DimL = DimA+EtaMax*DimSol, IteMax = 2000 )
      integer Eta, Info, IterMI(2), NA, NZA, NZL
      real TimePI(3), TimeTt
      real*8 Eps(3)
      logical CalcRv, ErroHB, PerfDS
c        arrays with dimension DimA
      integer ARowId(DimA)
      real*8 A(DimA)
c        arrays with dimension DimL
      integer LinRow(DimL), LRowId(DimL)
      real*8 L(DimL)
c        arrays with dimension DimSol
      integer AColPt(DimSol+1), LColPt(DimSol+1), RowCol(3*DimSol)
      real*8 B(DimSol), W1(DimSol), W2(DimSol), W3(DimSol), W4(DimSol)
     s, Xk(DimSol)
      logical*1 LogRow(DimSol)
c        arrays with dimension IteMax
      real*8 VAlpha(IteMax), VBeta(IteMax)
c
      integer I, J, VetEta(5)
      real SumTim
      character FilNaI*50
      data VetEta / -10, -5, 0, 10, 20/
c
c        Set linear system of equations
c



c      open(unit=5, file='arquivos_cccg', status='old' )
c   10 continue
c         read(5,105,end=20) FilNaI
c         open( unit=1, file=FilNaI, status='old' )
c         call HarBoe( NA, NZA, A, AColPt, ARowId, B, DimA, DimSol
c     d   , ErroHB, FilNaI )
         SumTim = 0.0
         PerfDS = .false.
         CalcRv = .false.
c         write(4,104)



      NA=3
      NZA=6

      A(1)=3.
      A(2)=1.
      A(3)=-1.
      A(4)=3.
      A(5)=1.
      A(6)=2.

      AColPt(1)=1
      AColPt(2)=4
      AColPt(3)=6
      AColPt(4)=7

      ARowId(1)=1
      ARowId(2)=2
      ARowId(3)=3
      ARowId(4)=2
      ARowId(5)=3
      ARowId(6)=3


      b(1)=1.
      b(2)=2.
      b(3)=3.





         do i = 1, 5
            Eta = VetEta(i)


c        call CCCG

c      ETA=10





c
            do j = 1, NA
               Xk(j) = 1.0d0
            enddo
            call CCCG( NA, NZA, A, AColPt, ARowId, B, NZL, L, LColPt
     d      , LRowId, Eta, EtaMax, DimL, Eps, IteMax, IterMI, LinRow
     e      , LogRow, RowCol, W1, W2, W3, W4, Xk, VAlpha, VBeta
     f      , TimePI, CalcRv, PerfDS, Info )
            if( Info.eq.-1.or.Info.eq.-2) then
               write(*,106) Info
               stop
            endif
            TimeTt = TimePI(1) + TimePI(2) + TimePI(3)
            SumTim = SumTim + TimeTt
            write(3,104)
            write(3,114) NZL, Eta, IterMI(2), Eps(1), Eps(2), Eps(3)
     x      , Info, TimePI(1), TimePI(2), TimePI(3), TimeTt
            write(4,114) NZL, Eta, IterMI(2), Eps(1), Eps(2), Eps(3)
     x      , Info, TimePI(1), TimePI(2), TimePI(3), TimeTt
            write(*,104)
            write(*,114) NZL, Eta, IterMI(2), Eps(1), Eps(2), Eps(3)
     x      , Info, TimePI(1), TimePI(2), TimePI(3), TimeTt
         enddo
         write(*,124) SumTim
         write(4,124) SumTim



      print*,'Xk'
      print*,(Xk(i),i=1,NA)





c      goto 10
c   20 continue
      stop
  104 format('NZero_L',2x,'Eta',1x,'Iter',3x,'Toler',3x,'CG_Res',2x
     g,'True_Res',2x,'If',2x,'DSTim',2x,'CCTim',2x,'ItTim',2x
     h,'TtTim-s')
  105 format(a50)
  106 format('* * * E R R O R : Info = ',i2)
  114 format(i7,i5,i5,1pe9.1,2(1pe9.2),i4,0pf7.2,2f7.2,f9.2)
  124 format('Total time =',f7.2,' s'/)
      end
c
c
c
c        Subroutine CCCG
c
c        Solve a system of linear equation by the Conjugate Gradient method
c        preconditioned by the Controlled Cholesky factorisation
c
      subroutine CCCG( NA, NZA, A, AColPt, ARowId, B, NZL, L, LColPt
     t, LRowId, Eta, EtaMax, DimL, Eps, IteMax, IterMI, LinRow, LogRow
     u, RowCol, W1, W2, W3, W4, Xk, VAlpha, VBeta, TimePI, CalcRv
     v, PerfDS, Info )
c
      integer AColPt(NA+1), ARowId(NZA), DimL, Eta, EtaMax, Info
     j, IterMI(2), LColPt(NA+1), LinRow(DimL), LRowId(DimL), NA, NZA
     k, NZL, RowCol(3*NA), IteMax
      real TimePI(3)
      real*8 A(NZA), B(NA), Eps(3), L(DimL), W1(NA), W2(NA)
     s, W3(NA), W4(NA), Xk(NA), VAlpha(IteMax), VBeta(IteMax)
      logical*1 LogRow(NA)
      logical CalcRv, PerfDS
c
c        parameter    type      in/out  description
c
c        NA           integer   entry   order of A
c        NZA          integer   entry   number of nonzero elements of A
c        A(NZA)       real*8    entry   elements of matrix A
c        AColPt(NA+1) integer   entry   location of first entries of column of A
c        ARowId(NZA)  integer   entry   row indices of A
c        B(NA)
c        NZL          integer   exit    number of nonzero elements of L
c        L(DimL)      real*8    exit    elements of L
c        LColPt(NA+1) integer   exit    location of 1st entries of columns of L
c        LRowId(DimL) integer   exit    row indices of L
c        Eta          integer   entry   number of extra elements per column
c        EtaMax       integer   entry   maximum value allowed for Eta
c        DimL         integer   entry   maximum size of L
c        Eps(3)       real*8            where :
c            1                  entry   tolerance for |Rk|/|Ro|
c            2                  exit    residual norm of preconditioned system
c            3                  exit    true residual norm of original system
c        IterMI(2)    integer           where :
c               1               entry   maximum number of iterations
c               2               exit    number of iterations
c        LinRow(DimL) integer   exit    pointer for next element in the row
c        LogRow(NA)   logical*1 exit    existence of row elements
c        RowCol(3*NA) integer   exit    workspace
c        W1(NA)       real*8    exit    workspace
c        W2(NA)       real*8    exit    workspace
c        W3(NA)       real*8    exit    workspace
c        W4(NA)       real*8    exit    workspace
c        Xk(NA)       real*8    exit    solution
c        TimePI(3)    real*8            where :
c               1               exit    diagonal scaling preconditioner time
c               2               exit    CCF preconditioner time
c               3               exit    CG iteration time
c        PerfDS       logical
c        CalcRv       logical   entry   calculate Ritz values
c        VAlpha(IteMax)  real*8
c        VBeta(IteMax)   real*8
c        Info         integer   exit    error ( <0 ) or warning ( >0 ) :
c                      0 : no error
c                     -1 : NA < 1 or NZA < NA or Eta > EtaMax
c                     -3 : D(I,I) have to be modified more than MaxShi times
c                          where MaxShi = 15
c                     -4 : NormD > Toler
c                     -5 : dsterf(LAPACK) failed
c                      k : D(I,I) was modified k times
c
c        Local variables
c
      real*8 Eps8
      real Ti, Tf
c
c        External subprograms
c
      real*8 d1mach
      external CCFeta, CG, diasca, d1mach, RitzVa
c
      Info = 0
      if( NA.lt.1.or.NZA.lt.NA.or.Eta.gt.EtaMax ) then
          Info = -1
          return
      endif
      if( .not.PerfDS ) then
c
c       call subroutine to perform the Diagonal Scaling
c
         call cpu_time(Ti)
         call diasca( A, AColPt, ARowId, NA, NZA, Info )
         call cpu_time(Tf)
         TimePI(1) = Tf - Ti
         PerfDS = .true.
         if( Info.eq.-2 ) return
      endif
c
c        call function to determine the machine precision
c
      Eps8 = d1mach(3)
      call cpu_time(Ti)
c
c        call subroutine to compute Controlled Cholesky Factorization
c
      call CCFeta( NA, NZA, A, AColPt, ARowId, NZL, L, LColPt, LRowId
     d, Eta, DimL, Eps8, LinRow, LogRow, RowCol, W1, W2, Info )
      call cpu_time(Tf)
      TimePI(2) = Tf - Ti
      if( Info.eq.-3 ) return
      IterMI(1) = min0(NA,IteMax)
      Eps(1) = 1.0d-10
      Eps(2) = 9.99d99
      Eps(3) = 9.99d99
      call cpu_time(Ti)
c
c        call CG solver
c
      call CG( NA, NZA, A, AColPt, ARowId, B, NZL, L, LColPt, LRowId
     d, Eps, IteMax, IterMI, W1, W2, W3, W4, Xk, VAlpha, VBeta, Info )
      call cpu_time(Tf)
      TimePI(3) = Tf - Ti
c
c        call subroutine to compute Ritz values
c
      if( CalcRv ) then
         call RitzVa( VAlpha, VBeta, IterMI(2), Info )
      endif
      return
      end
c
c
c
c        Subroutine CGpre
c                                                1/2 -1          1/2 -T
c        Perform the CCF preconditioning x = (L.D   )  . A . (L.D   )  . v
c
      subroutine CGpre( NA, NZA, A, AColPt, ARowId, NZL, L, LColPt
     t, LRowId, Work, V, X )
c
      integer AColPt(NA+1), NZL, NA, ARowId(NZA), NZA
     j, LColPt(NA+1), LRowId(NZL)
      real*8 A(NZA), L(NZL), Work(NA), V(NA)
     s, X(NA)
c
c        parameter    type      in/out  description
c
c        NA           integer   entry   order of A
c        NZA          integer   entry   number of nonzero elements of A
c        A(NZA)       real*8    entry   elements of matrix A
c        AColPt(NA+1) integer   entry   location of first entries of column of A
c        ARowId(NZA)  integer   entry   row indices of A
c        NZL          integer   entry   number of nonzero elements of L
c        L(NZL)       real*8    entry   elements of L
c        LColPt(NA+1) integer   entry   location of 1st entries of columns of L
c        LRowId(NZL)  integer   entry   row indices of L
c        Work(NA)     real*8    exit    workspace
c        V(NA)        real*8    entry   vector to be preconditioned
c        X(NA)        real*8    exit    preconditioned vector
c
c        External subprograms
c
         external dcopy, LDLtBS, LDLtFS, MaVeCo
c
c        Work <-- v
c
      call dcopy( NA, V, 1, Work, 1 )
c
c                 1/2  -T                     1/2  T
c        x = ( L.D    )  . Work   ===>   ( L.D    ) . x = Work
c
      call LDLtBS( L, NA, NZL, LColPt, LRowId, Work, X )
c
c        Work = A . X
c
      call MaVeCo( A, NA, NZA, AColPt, ARowId, X, Work )
c
c                 1/2 -1                     1/2
c        x = ( L.D   )  . Work   ===>   ( L.D   ) . y = Work
c
      call LDLtFS( L, NA, NZL, LColPt, LRowId, Work, X )
      return
      end
c
c
c
c        Subroutine CG
c
c        Solve a linear system by the Conjugate Gradient method
c        preconditioned by CCF
c
      subroutine CG( NA, NZA, A, AColPt, ARowId, B, NZL, L, LColPt
     t, LRowId, Eps, IteMax, IterMI, Rk, Wk, Pk, Tk, Xk, VAlpha
     u, VBeta, Info )
c
      integer AColPt(NA+1), ARowId(NZA), Info, IteMax, LColPt(NA+1)
     j, LRowId(NZL), NA, NZA, NZL, IterMI(2)
      real*8 A(NZA), B(NA), Eps(3), L(NZL), Pk(NA), Rk(NA), Tk(NA)
     s, Wk(NA), Xk(NA), VAlpha(IteMax), VBeta(IteMax)
c
c        parameter    type      in/out  description
c
c        NA           integer   entry   order of A
c        NZA          integer   entry   number of nonzero elements of A
c        A(NZA)       real*8    entry   elements of matrix A
c        AColPt(NA+1) integer   entry   location of first entries of column of A
c        ARowId(NZA)  integer   entry   row indices of A
c        B(NA)        real*8    entry   right hand side
c        NZL          integer   entry   number of nonzero elements of L
c        L(NZL)       real*8    entry   elements of L
c        LColPt(NA+1) integer   entry   location of 1st entries of columns of L
c        LRowId(DimL) integer   entry   row indices of L
c        Eps(3)       real*8            where :
c            1                  entry   tolerance for |Rk|/|Ro|
c            2                  exit    residual norm of preconditioned system
c            3                  exit    true residual norm of original system
c        IterMI(2)    integer           where :
c               1               entry   maximum number of iterations
c               2               exit    number of iterations
c        Pk(NA)       real*8    exit    search direction
c        Rk(NA)       real*8    exit    residual
c        Tk(NA)       real*8    exit    workspace and true residual on exit
c        Xk(NA)       real*8    exit    solution
c        Wk(NA)       real*8    exit    workspace
c        VAlpha       real*8    exit
c        VBeta        real*8    exit
c        Info         integer   exit    condition of error :
c                      0 : no error
c                     -4 : NormD > Toler
c
c        First preconditioner: diagonal scaling
c
c                  -1/2    -1/2   1/2      -1/2
c     Ax = b --> (R    .A.R    )(R   x) = R    b --> Vy = c
c
c        Second preconditioner: Controlled Cholesky Factorisation
c
c                     1/2 -1       1/2 -T      1/2 T         1/2 -1
c     Vy = c --> [(L.D   )  .V.(L.D   )  ][(L.D   ) y] = (L.D   )  c --> Hz = d
c
c
c        Local variables
c
      integer I, Iter
      real*8 Alpha, Beta, NormD, NormRk, NormR0, RDiag, Ro, Ro1, Toler
      logical Error
c
c        External subprograms
c
      real*8 ddot, dnrm2
      external CGpre, LDLtBS, LDLtFS, MaVeCo, MaVeLt
     f, daxpy, ddot, dnrm2, dscal
c
c        Intrinsic function
c
      real*8 dsqrt
c
      write(*,106)
      Toler = Eps(1)
      do I = 1, NA
         RDiag = A(AColPt(I))
c             -1/2
c        c = R    b
c
         Wk(I) = B(I) * RDiag
c             1/2
c        y = R   x
c
         Tk(I) = Xk(I) / RDiag
      enddo
c                1/2 T
c        z = (L.D   ) y, z stored in vector Xk
c
      call MaVeLt( L, NA, NZL, LColPt, LRowId, Tk, Xk )
c
c        p = Vy, Pk used as work vector
c
      call MaVeCo( A, NA, NZA, AColPt, ARowId, Tk, Pk )
c
c        p = -Vy + c
c
      call daxpy( NA, -1.0d0, Pk, 1, Wk, 1 )
c                1/2 -1
c        r = (L.D   )  (c - Vy)
c
      call LDLtFS( L, NA, NZL, LColPt, LRowId, Wk, Rk )
c
c        Pk = 0
c
      do I = 1, NA
         Pk(I) = 0.0d0
      enddo
      Ro1 = 1.0d0
      Ro = ddot( NA, Rk, 1, Rk, 1 )
      NormR0 = dsqrt( Ro )
      Iter = 0
      Error = .true.
      do while ( Error.and.Iter.lt.IterMI(1) )
         Iter = Iter + 1
         Beta = Ro / Ro1
c
c        p = beta.p + r
c
         call dscal(NA, Beta, Pk, 1)
         call daxpy( NA, 1.0d0, Rk, 1, Pk, 1 )
c                  1/2 -1       1/2 -T
c        Tk = [(L.D   )  .V.(L.D   )  ].Pk
c
         call CGpre( NA, NZA, A, AColPt, ARowId, NZL, L, LColPt, LRowId
     d   , Wk, Pk, Tk )
         Alpha = Ro / ddot( NA, Pk, 1, Tk, 1 )
c
c        x = alpha.p + x
c
         call daxpy( NA,  Alpha, Pk, 1, Xk, 1 )
c
c        r = -alpha.t + r
c
         call daxpy( NA, -Alpha, Tk, 1, Rk, 1 )
         Ro1 = Ro
         Ro = ddot( NA, Rk, 1, Rk, 1 )
         NormRk = dsqrt( Ro )
         NormD = NormRk / NormR0
         Error = NormD.gt.Toler
         VAlpha(Iter) = Alpha
         VBeta(Iter) = Beta
      enddo
      Eps(2) = NormRk
c                1/2 -T
c        y = (L.D   )  z
c
      call LDLtBS( L, NA, NZL, LColPt, LRowId, Xk, Tk )
c
c        w = Vy
c
      call MaVeCo( A, NA, NZA, AColPt, ARowId, Tk, Wk )
      do I = 1, NA
         RDiag = A(AColPt(I))
c             -1/2
c        x = R    y, solution of original system Ax = b
c
         Xk(I) = Tk(I) * RDiag
c                                        1/2
c        True Residual t = b - Ax = b - R   Vy
c
         Tk(I) = B(I) - Wk(I) / RDiag
      enddo
      Eps(3) = dnrm2( NA, Tk, 1 )
      IterMI(2) = Iter
      if( Error ) Info = -4
      return
  106 format('Solving system by CG . . .')
      end
c
c
c
c        subroutine diasca
c
c         It performs the diagonal scaling preconditioning
c
      subroutine diasca( A, AColPt, ARowId, NA, NZA, Info )
c
      integer NA, AColPt(NA+1), ARowId(NZA), Info
      real*8 A(NZA)
c
      integer I, J, J1, J2
      real*8 dsqrt, RDiag
c
c        Calculating reciprocals of the diagonal of the coefficient matrix
c
      do I = 1, NA
         if( A(AColPt(I)).gt.0.0d0 ) then
            A(AColPt(I)) = 1.0d0 / dsqrt( A(AColPt(I)) )
         else
            Info = -2
            return
         endif
      enddo
c
c        Performing diagonal scaling
c
      do I = 1, NA
         RDiag = A(AColPt(I))
         J1 = AColPt(I) + 1
         J2 = AColPt(I+1) - 1
         do J = J1, J2
            A(J) = A(J) * RDiag * A(AColPt(ARowId(J)))
         enddo
      enddo
      return
      end
c
c
c
c        Subroutine RitzVa
c
      subroutine RitzVa( VAlpha, VBeta, Iter, Info)
c
      integer Iter, Info
      real*8 VAlpha(Iter), VBeta(Iter)
c
      integer I, InfoRv
      real*8 Alpha, Beta, dsqrt, R
c
      Alpha = VAlpha(1)
      R = 1.0d0 / Alpha
      VAlpha(1) = R
      do I = 2, Iter
         Alpha = VAlpha(I)
         Beta = VBeta(I)
         VAlpha(I) = 1.0d0 / Alpha + Beta * R
         VBeta(I-1) =  -dsqrt( Beta ) * R
         R = 1.0d0 / Alpha
      enddo
      call dsterf( Iter, VAlpha, VBeta, InfoRv )
      if( InfoRv.eq.0 ) then
         write(3,103) VAlpha(Iter), VAlpha(1), VAlpha(Iter)/VAlpha(1)
         write(*,103) VAlpha(Iter), VAlpha(1), VAlpha(Iter)/VAlpha(1)
      else
         do I = 1, Iter
            VAlpha(I) = 0.0d0
         enddo
         Info = -5
      endif
  103 format('Ritzvalues: max =',1pe12.5,5x,'min =',1pe12.5,5x
     g,'ratio = ',1pe12.5)
      return
      end
c
c
c
c        subroutine LDLtBS
c                                                    1/2  T
c         It performs the backward substitution ( L.D    ) x = v
c
      subroutine LDLtBS( L, NA, NZA, ColPtr, RowInd, V, X )
c
      integer NA, NZA, ColPtr(NA+1), RowInd(NZA)
      real*8 L(NZA), V(NA), X(NA)
c
      integer I, J, LL, LU
      real*8 Sum
c
      do J = NA, 1, -1
         LL = ColPtr(J)   + 1
         LU = ColPtr(J+1) - 1
         Sum = 0.0d0
         do I = LL, LU
            Sum = Sum + L(I) * X(RowInd(I))
         enddo
         X(J) = V(J) / L(ColPtr(J)) - Sum
      enddo
      return
      end
c
c
c
c        subroutine LDLtFS
c                                                 1/2
c        It performs the forward substituion ( L.D   )x = v
c
      subroutine LDLtFS( L, NA, NZA, ColPtr, RowInd, V, X )
c
      integer NA, NZA, ColPtr(NA+1), RowInd(NZA)
      real*8 L(NZA), V(NA), X(NA)
c
      integer I, J, K, LL, LU
      external dcopy
c
c        X <-- V
c
      call dcopy( NA, V, 1, X, 1 )
      do J = 1, NA
         LL = ColPtr(J)   + 1
         LU = ColPtr(J+1) - 1
         do K = LL, LU
            I = RowInd(K)
            X(I) = X(I) - L(K) * X(J)
         enddo
         X(J) = X(J) / L(ColPtr(J))
      enddo
      return
      end
c
c
c
c        subroutine MaVeCo
c
c        It performs the product of a scaled symmetric matrix by a vector Mv = x
c
      subroutine MaVeCo( A, NA, NZA, ColPtr, RowInd, V, X )
c
      integer NA, NZA, ColPtr(NA+1), RowInd(NZA)
      real*8 A(NZA), V(NA), X(NA)
c
      integer I, J, L, LL, LU
c
      do J = NA, 1, -1
         X(J) = V(J)
         LL = ColPtr(J) + 1
         LU = ColPtr(J+1) - 1
         do L = LL, LU
            I = RowInd(L)
            X(I) = X(I) + A(L) * V(J)
            X(J) = X(J) + A(L) * V(I)
         enddo
      enddo
      return
      end
c
c
c
c        subroutine MaVeLt
c                                                1/2 T
c        It performs the multiplication x = ( L.D   ) v
c
      subroutine MaVeLt( L, NA, NZL, ColPtr, RowInd, V, X )
c
      integer NA, NZL, ColPtr(NA+1), RowInd(NZL)
      real*8 L(NZL), V(NA), X(NA)
c
      integer I, J, LL, LU
c
      do I = 1, NA
         LL = ColPtr(I)   + 1
         LU = ColPtr(I+1) - 1
         X(I) = V(I)
         do J = LL, LU
            X(I) = X(I) + L(J) * V(RowInd(J))
         enddo
         X(I) = X(I) * L(LL-1)
      enddo
      return
      end
c
      include 'harboe.f'
      include 'ccf.f'
      include 'blas.f'
      include 'eigenv.f'
