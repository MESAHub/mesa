module test_bobyqa

   use num_def
   use num_lib

   implicit none

   integer :: nfcn

contains

   subroutine do_test_bobyqa
      real(dp), dimension(100) :: X, XL, XU
      real(dp), dimension(10000) :: W
      real(dp), parameter :: max_valid_value = 1d99
      real(dp), parameter :: BDL = -1.0d0
      real(dp), parameter :: BDU = 1.0d0
      real(dp) :: f, RHOBEG, RHOend
      integer :: I, IPRINT, N, MAXFUN, NPT
      include 'formats'
      IPRINT = 0
      MAXFUN = 5000
      RHOend = 1.0D-6
      do N = 2, 6, 2
         nfcn = 0
         NPT = 2*N + 1
         do I = 1, N
            XL(I) = BDL
            XU(I) = BDU
            X(I) = DBLE(I)/DBLE(N + 1)
         end do
         RHOBEG = 0.2D0*X(1)
         write (*, '(4X,A,I2,A,I3)') 'test BOBYQA with N =', N, ' and NPT =', NPT
         call BOBYQA(N, NPT, X, XL, XU, RHOBEG, RHOend, IPRINT, MAXFUN, W, CALFUN, max_valid_value)
         call calfun(n, x, f)
         !write(*,2) 'f', nfcn, f
         if (abs(f) > 1d-10) write (*, *) 'failed in test of BOBYQA: min f', f
      end do
   end subroutine do_test_bobyqa

   subroutine calfun(n, x, f)
      use const_def, only: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f
      integer :: i, j, iw, np
      real(dp) :: sum

      real(dp) :: Y(10, 10)
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
      RETURN
   end subroutine CALFUN

   subroutine xCALFUN(N, X, F)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f
      integer :: i, j
      real(dp) :: temp
      F = 0.0D0
      do I = 4, N, 2
      do J = 2, I - 2, 2
         TEMP = (X(I - 1) - X(J - 1))**2 + (X(I) - X(J))**2
         TEMP = DMAX1(TEMP, 1.0D-6)
         F = F + 1.0D0/DSQRT(TEMP)
      end do
      end do
      RETURN
   end subroutine xCALFUN

end module test_bobyqa
