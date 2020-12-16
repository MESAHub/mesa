      module test_newuoa
      
      use num_def
      use num_lib
      
      integer :: nfcn

      contains
      
      subroutine do_test_newuoa      

!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
!     with NPT = 2N+1.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(10),W(10000)
      real*8, parameter :: max_valid_value = 1d99
      include 'formats'
      IPRINT=0
      MAXFUN=5000
      RHOEND=1.0D-6
      DO 30 N=2,6,2
      nfcn = 0
      NPT=2*N+1
      DO 10 I=1,N
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)
      PRINT 20, N,NPT
   20 FORMAT (4X,'test NEWUOA with N =',I2,' and NPT =',I3)
      CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN,max_valid_value)
      call calfun(n,x,f)
      !write(*,2) 'f', nfcn, f
      if (abs(f) > 1d-10) write(*,*) 'failed in test of newuoa: min f', f
   30 CONTINUE
      END subroutine do_test_newuoa



      subroutine calfun(n,x,f)
      use const_def, only: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f

      real(dp) :: Y(10,10)
      nfcn = nfcn + 1
      DO 10 J=1,N
      Y(1,J)=1.0D0
   10 Y(2,J)=2.0D0*X(J)-1.0D0
      DO 20 I=2,N
      DO 20 J=1,N
   20 Y(I+1,J)=2.0D0*Y(2,J)*Y(I,J)-Y(I-1,J)
      F=0.0D0
      NP=N+1
      IW=1
      DO 40 I=1,NP
      SUM=0.0D0
      DO 30 J=1,N
   30 SUM=SUM+Y(I,J)
      SUM=SUM/DFLOAT(N)
      IF (IW .GT. 0) SUM=SUM+1.0D0/DFLOAT(I*I-2*I)
      IW=-IW
   40 F=F+SUM*SUM
      RETURN
      END SUBROUTINE CALFUN      
      
      end module test_newuoa
