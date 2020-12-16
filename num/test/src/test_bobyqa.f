      module test_bobyqa
      
      use num_def
      use num_lib
      
      integer :: nfcn

      contains
      
      subroutine do_test_bobyqa      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(100),XL(100),XU(100),W(10000)
      real*8, parameter :: max_valid_value = 1d99
      include 'formats'
      BDL=-1.0D0
      BDU=1.0D0
      IPRINT=0
      MAXFUN=5000
      RHOEND=1.0D-6
      DO 30 N=2,6,2
      nfcn = 0
      NPT=2*N+1
      DO 10 I=1,N
      XL(I)=BDL
      XU(I)=BDU
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)
      PRINT 20, N,NPT
   20 FORMAT (4X,'test BOBYQA with N =',I2,' and NPT =',I3)
      CALL BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN,max_valid_value)
      call calfun(n,x,f)
      !write(*,2) 'f', nfcn, f
      if (abs(f) > 1d-10) write(*,*) 'failed in test of BOBYQA: min f', f
   30 CONTINUE
      END subroutine do_test_bobyqa



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


      SUBROUTINE xCALFUN (N,X,F)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f
      integer :: i, j
      real(dp) :: temp
      F=0.0D0
      DO 10 I=4,N,2
      DO 10 J=2,I-2,2
      TEMP=(X(I-1)-X(J-1))**2+(X(I)-X(J))**2
      TEMP=DMAX1(TEMP,1.0D-6)
   10 F=F+1.0D0/DSQRT(TEMP)
      RETURN
      END SUBROUTINE xCALFUN
      
      
      end module test_bobyqa
