c    subroutine diasca
c
c         It performs the diagonal scaling preconditioning
c
      subroutine dscale(A, AColPt, ARowId, NA, NZA,  Info )
c
      integer NA, AColPt(*), ARowId(*), Info
      real*8 A(*)
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
