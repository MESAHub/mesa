module test_newuoa

   use num_def
   use num_lib
   use const_def, only: dp

   implicit none

   integer :: nfcn

contains

   subroutine do_test_newuoa

!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
!     with NPT = 2N+1.
!
      real(dp), dimension(10) :: X
      real(dp), dimension(10000) :: W
      real(dp), parameter :: max_valid_value = 1d99
      real(dp) :: f, RHOBEG, RHOend
      integer :: IPRINT, I, N, NPT, MAXFUN
      include 'formats'
      IPRINT = 0
      MAXFUN = 5000
      RHOend = 1.0D-6
      do N = 2, 6, 2
         nfcn = 0
         NPT = 2*N + 1
         do I = 1, N
            X(I) = DBLE(I)/DBLE(N + 1)
         end do
         RHOBEG = 0.2D0*X(1)
         write (*, '(4X,A,I2,A,I3)') 'test NEWUOA with N =', N, ' and NPT =', NPT
         call NEWUOA(N, NPT, X, RHOBEG, RHOend, IPRINT, MAXFUN, W, CALFUN, max_valid_value)
         call calfun(n, x, f)
         !write(*,2) 'f', nfcn, f
         if (abs(f) > 1d-10) write (*, *) 'failed in test of newuoa: min f', f
      end do
   end subroutine do_test_newuoa

   subroutine calfun(n, x, f)
      use const_def, only: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f

      integer :: I, J, IW, NP
      real(dp) :: Y(10, 10), sum
      nfcn = nfcn + 1
      do J = 1, N
         Y(1, J) = 1.0D0
         Y(2, J) = 2.0D0*X(J) - 1.0D0
      end do
      do I = 2, N
      do J = 1, N
         Y(I + 1, J) = 2.0D0*Y(2, J)*Y(I, J) - Y(I - 1, J)
      end do
      end do
      F = 0.0D0
      NP = N + 1
      IW = 1
      do I = 1, NP
         SUM = 0.0D0
         do J = 1, N
            SUM = SUM + Y(I, J)
         end do
         SUM = SUM/DBLE(N)
         IF (IW > 0) SUM = SUM + 1.0D0/DBLE(I*I - 2*I)
         IW = -IW
         F = F + SUM*SUM
      end do
      return
   end subroutine CALFUN

end module test_newuoa
