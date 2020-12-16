      module mod_newuoa
      use const_def, only: dp
      use math_lib

      contains
      
      
      

      SUBROUTINE do_newuoa (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,
     >         CALFUN,max_valid_value)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
      interface
#include "num_newuoa_proc.dek"
      end interface
      real*8, intent(in) :: max_valid_value
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     There can be some freedom in the interpolation conditions, which is
!     taken up by minimizing the Frobenius norm of the change to the second
!     derivative of the quadratic model, beginning with a zero matrix. The
!     arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in the
!       interval [N+2,(N+1)(N+2)/2].
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
!
!     The above settings provide a partition of W for subroutine NEWUOB.
!     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
!     W plus the space that is needed by the last array of NEWUOB.
!
      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW),CALFUN,max_valid_value)
   20 RETURN
      END subroutine do_newuoa


      SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
     1  XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,
     1  CALFUN,max_valid_value)
      IMPLICIT REAL*8 (A-H,O-Z)
      interface
#include "num_newuoa_proc.dek"
      end interface
      real*8, intent(in) :: max_valid_value
      logical :: do_replace
      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),FVAL(*),
     1  GQ(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),VLAG(*),W(*)
!
!     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
!       to the corresponding arguments in SUBROUTINE NEWUOA.
!     XBASE will hold a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     FVAL will hold the values of F at the interpolation points.
!     GQ will hold the gradient of the quadratic model at XBASE.
!     HQ will hold the explicit second derivatives of the quadratic model.
!     PQ will contain the parameters of the implicit second derivatives of
!       the quadratic model.
!     BMAT will hold the last N columns of H.
!     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
!       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
!       the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     D is reserved for trial steps from XOPT.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     The array W will be used for working space. Its length must be at least
!       10*NDIM = 10*(NPT+N).
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
!
!     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,NH
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPTM
   40 ZMAT(K,J)=ZERO
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF,.).
!
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              XPT(NF,NFM)=RHOBEG
          ELSE IF (NFM .GT. N) THEN
              XPT(NF,NFMM)=-RHOBEG
          END IF
      ELSE
          ITEMP=(NFMM-1)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XIPT=RHOBEG
          IF (FVAL(IPT+NP) .LT. FVAL(IPT+1)) XIPT=-XIPT
          XJPT=RHOBEG
          IF (FVAL(JPT+NP) .LT. FVAL(JPT+1)) XJPT=-XJPT
          XPT(NF,IPT)=XIPT
          XPT(NF,JPT)=XJPT
      END IF
!
!     Calculate the next value of F, label 70 being reached immediately
!     after this calculation. The least function value so far and its index
!     are required.
!
      DO 60 J=1,N
   60 X(J)=XPT(NF,J)+XBASE(J)
      GOTO 310
   70 FVAL(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          FOPT=F
          KOPT=1
      ELSE IF (F .LT. FOPT) THEN
          FOPT=F
          KOPT=NF
      END IF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in
!     the cases when NF is at most 2*N+1.
!
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              GQ(NFM)=(F-FBEG)/RHOBEG
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHOBEG
                  BMAT(NF,NFM)=ONE/RHOBEG
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM .GT. N) THEN
              BMAT(NF-N,NFMM)=HALF/RHOBEG
              BMAT(NF,NFMM)=-HALF/RHOBEG
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              TEMP=(FBEG-F)/RHOBEG
              HQ(IH)=(GQ(NFMM)-TEMP)/RHOBEG
              GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
          END IF
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          IF (XIPT .LT. ZERO) IPT=IPT+N
          IF (XJPT .LT. ZERO) JPT=JPT+N
          ZMAT(1,NFMM)=RECIP
          ZMAT(NF,NFMM)=RECIP
          ZMAT(IPT+1,NFMM)=-RECIP
          ZMAT(JPT+1,NFMM)=-RECIP
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/(XIPT*XJPT)
      END IF
      IF (NF .LT. NPT) GOTO 50
!
!     Begin the iterative procedure, because the initial model is complete.
!
      RHO=RHOBEG
      DELTA=RHO
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      XOPTSQ=ZERO
      DO 80 I=1,N
      XOPT(I)=XPT(KOPT,I)
   80 XOPTSQ=XOPTSQ+XOPT(I)**2
   90 NFSAV=NF
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve the model.
!
  100 KNEW=0
      CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
      DSQ=ZERO
      DO 110 I=1,N
  110 DSQ=DSQ+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          IF (NF .LE. NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!     to BMAT that do not depend on ZMAT.
!
  120 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          TEMPQ=0.25D0*XOPTSQ
          DO 140 K=1,NPT
          SUM=ZERO
          DO 130 I=1,N
  130     SUM=SUM+XPT(K,I)*XOPT(I)
          TEMP=PQ(K)*SUM
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          DO 140 I=1,N
          GQ(I)=GQ(I)+TEMP*XPT(K,I)
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          VLAG(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
          IP=NPT+I
          DO 140 J=1,I
  140     BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
          DO 180 K=1,NPTM
          SUMZ=ZERO
          DO 150 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
  150     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 170 J=1,N
          SUM=TEMPQ*SUMZ*XOPT(J)
          DO 160 I=1,NPT
  160     SUM=SUM+W(I)*XPT(I,J)
          VLAG(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 170 I=1,NPT
  170     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 180 I=1,N
          IP=I+NPT
          TEMP=VLAG(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 180 J=1,I
  180     BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
!
!     The following instructions complete the shift of XBASE, including
!     the changes to the parameters of the quadratic model.
!
          IH=0
          DO 200 J=1,N
          W(J)=ZERO
          DO 190 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  190     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 200 I=1,J
          IH=IH+1
          IF (I .LT. J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
          GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  200     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 210 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
  210     XOPT(J)=ZERO
          XOPTSQ=ZERO
      END IF
!
!     Pick the model step if KNEW is positive. A different choice of D
!     may be made later, if the choice of D by BIGLAG causes substantial
!     cancellation in DENOM.
!
      IF (KNEW .GT. 0) THEN
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP,
     1      D,ALPHA,VLAG,VLAG(NPT+1),W,W(NP),W(NP+N))
      END IF
!
!     Calculate VLAG and BETA for the current choice of D. The first NPT
!     components of W_check will be held in W.
!
      DO 230 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 220 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  220 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
  230 VLAG(K)=SUM
      BETA=ZERO
      DO 250 K=1,NPTM
      SUM=ZERO
      DO 240 I=1,NPT
  240 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 250 I=1,NPT
  250 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      DO 280 J=1,N
      SUM=ZERO
      DO 260 I=1,NPT
  260 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 270 K=1,N
  270 SUM=SUM+BMAT(JP,K)*D(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  280 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     If KNEW is positive and if the cancellation in DENOM is unacceptable,
!     then BIGDEN calculates an alternative model step, XNEW being used for
!     working space.
!
      IF (KNEW .GT. 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) .LE. 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,D,W,VLAG,BETA,XNEW,W(NDIM+1),W(6*NDIM+1))
          END IF
      END IF
!
!     Calculate the next value of the objective function.
!
  290 DO 300 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  300 X(I)=XBASE(I)+XNEW(I)
      NF=NF+1
  310 IF (NF .GT. NFTEST) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) PRINT 320
  320     FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 530
      END IF
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 330, NF,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .LE. NPT) GOTO 70
      IF (KNEW .EQ. -1) GOTO 530
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and set DIFF to the error of this prediction.
!
      VQUAD=ZERO
      IH=0
      DO 340 J=1,N
      VQUAD=VQUAD+D(J)*GQ(J)
      DO 340 I=1,J
      IH=IH+1
      TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) TEMP=HALF*TEMP
  340 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 350 K=1,NPT
  350 VQUAD=VQUAD+PQ(K)*W(K)
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. The branch when KNEW is positive occurs if D is not
!     a trust region step.
!
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          FOPT=F
          XOPTSQ=ZERO
          DO 360 I=1,N
          XOPT(I)=XNEW(I)
  360     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 410
!
!     Pick the next value of DELTA after a trust region step.
!
      IF (VQUAD .GE. ZERO) THEN
          IF (IPRINT .GT. 0) PRINT 370
  370     FORMAT (/4X,'Return from NEWUOA because a trust',
     1      ' region step has failed to reduce Q.')
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO .LE. TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 400 K=1,NPT
      HDIAG=ZERO
      DO 380 J=1,NPTM
      TEMP=ONE
      IF (J .LT. IDZ) TEMP=-ONE
  380 HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
      TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
      DISTSQ=ZERO
      DO 390 J=1,N
  390 DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
      IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)*(DISTSQ/RHOSQ)*(DISTSQ/RHOSQ)
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          KNEW=K
      END IF
  400 CONTINUE
      IF (KNEW .EQ. 0) GOTO 460
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!     can be moved. Begin the updating of the quadratic model, starting
!     with the explicit second derivative term.
!
  410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      FVAL(KNEW)=F
      IH=0
      DO 420 I=1,N
      TEMP=PQ(KNEW)*XPT(KNEW,I)
      DO 420 J=1,I
      IH=IH+1
  420 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      PQ(KNEW)=ZERO
!
!     Update the other second derivative parameters, and then the gradient
!     vector of the model. Also include the new interpolation point.
!
      DO 440 J=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 440 K=1,NPT
  440 PQ(K)=PQ(K)+TEMP*ZMAT(K,J)
      GQSQ=ZERO
      DO 450 I=1,N
      GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
      GQSQ=GQSQ+GQ(I)**2
  450 XPT(KNEW,I)=XNEW(I)
!
!     If a trust region step makes a small change to the objective function,
!     then calculate the gradient of the least Frobenius norm interpolant at
!     XBASE, and store it in W, using VLAG for a vector of right hand sides.
!
      IF (KSAVE .EQ. 0 .AND. DELTA .EQ. RHO) THEN
          IF (DABS(RATIO) .GT. 1.0D-2) THEN
              ITEST=0
          ELSE
              DO 700 K=1,NPT
  700         VLAG(K)=FVAL(K)-FVAL(KOPT)
              GISQ=ZERO
              DO 720 I=1,N
              SUM=ZERO
              DO 710 K=1,NPT
  710         SUM=SUM+BMAT(K,I)*VLAG(K)
              GISQ=GISQ+SUM*SUM
  720         W(I)=SUM
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
              ITEST=ITEST+1
              IF (GQSQ .LT. 1.0D2*GISQ) ITEST=0
              do_replace = (ITEST .GE. 3)
              if (.not. do_replace) then ! check for "invalid" value
                do k=1,npt
                  if (fval(k) > max_valid_value) then
                     do_replace = .true.
                     !write(*,*) 'newuoa value > max valid', fval(k), max_valid_value
                     exit
                  end if
                end do
              end if 
              IF (do_replace) THEN    
                  DO 730 I=1,N
  730             GQ(I)=W(I)
                  DO 740 IH=1,NH
  740             HQ(IH)=ZERO
                  DO 760 J=1,NPTM
                  W(J)=ZERO
                  DO 750 K=1,NPT
  750             W(J)=W(J)+VLAG(K)*ZMAT(K,J)
  760             IF (J .LT. IDZ) W(J)=-W(J)
                  DO 770 K=1,NPT
                  PQ(K)=ZERO
                  DO 770 J=1,NPTM
  770             PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
                  ITEST=0
              END IF
          END IF
      END IF
      IF (F .LT. FSAVE) KOPT=KNEW
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case KSAVE>0 occurs
!     when the new function value was calculated by a model step.
!
      IF (F .LE. FSAVE+TENTH*VQUAD) GOTO 100
      IF (KSAVE .GT. 0) GOTO 100
!
!     Alternatively, find out if the interpolation points are close enough
!     to the best point so far.
!
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO 480 K=1,NPT
      SUM=ZERO
      DO 470 J=1,N
  470 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  480 CONTINUE
!
!     If KNEW is positive, then set DSTEP, and branch back for the next
!     iteration, which will generate a "model step".
!
      IF (KNEW .GT. 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO .GT. ZERO) GOTO 100
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 100
!
!     The calculations with the current value of RHO are complete. Pick the
!     next values of RHO and DELTA.
!
  490 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 500
  500         FORMAT (5X)
              PRINT 510, RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 90
      END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
      IF (KNEW .EQ. -1) GOTO 290
  530 IF (FOPT .LE. F) THEN
          DO 540 I=1,N
  540     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 550, NF
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
          PRINT 520, F,(X(I),I=1,N)
      END IF
      RETURN
      END SUBROUTINE NEWUOB

      SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,D,W,VLAG,BETA,S,WVEC,PROD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  W(*),VLAG(*),S(*),WVEC(NDIM,*),PROD(NDIM,*)
      DIMENSION DEN(9),DENEX(9),PAR(9)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     D will be set to the step from XOPT to the new point, and on entry it
!       should be the D that was calculated by the last call of BIGLAG. The
!       length of the initial D provides a trust region bound on the final D.
!     W will be set to Wcheck for the final choice of D.
!     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
!     BETA will be set to the value that will occur in the updating formula
!       when the KNEW-th interpolation point is moved to its new position.
!     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
!       for working space.
!
!     D is calculated in a way that should provide a denominator with a large
!     modulus in the updating formula when the KNEW-th interpolation point is
!     shifted to the new position XOPT+D.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      QUART=0.25D0
      TWO=2.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*atan(ONE)
      NPTM=NPT-N-1
!
!     Store the first NPT elements of the KNEW-th column of H in W(N+1)
!     to W(N+NPT).
!
      DO 10 K=1,NPT
   10 W(N+K)=ZERO
      DO 20 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 20 K=1,NPT
   20 W(N+K)=W(N+K)+TEMP*ZMAT(K,J)
      ALPHA=W(N+KNEW)
!
!     The initial search direction D is taken from the last call of BIGLAG,
!     and the initial S is set below, usually to the direction from X_OPT
!     to X_KNEW, but a different direction to an interpolation point may
!     be chosen, in order to prevent S from being nearly parallel to D.
!
      DD=ZERO
      DS=ZERO
      SS=ZERO
      XOPTSQ=ZERO
      DO 30 I=1,N
      DD=DD+D(I)**2
      S(I)=XPT(KNEW,I)-XOPT(I)
      DS=DS+D(I)*S(I)
      SS=SS+S(I)**2
   30 XOPTSQ=XOPTSQ+XOPT(I)**2
      IF (DS*DS .GT. 0.99D0*DD*SS) THEN
          KSAV=KNEW
          DTEST=DS*DS/SS
          DO 50 K=1,NPT
          IF (K .NE. KOPT) THEN
              DSTEMP=ZERO
              SSTEMP=ZERO
              DO 40 I=1,N
              DIFF=XPT(K,I)-XOPT(I)
              DSTEMP=DSTEMP+D(I)*DIFF
   40         SSTEMP=SSTEMP+DIFF*DIFF
              IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
                  KSAV=K
                  DTEST=DSTEMP*DSTEMP/SSTEMP
                  DS=DSTEMP
                  SS=SSTEMP
              END IF
          END IF
   50     CONTINUE
          DO 60 I=1,N
   60     S(I)=XPT(KSAV,I)-XOPT(I)
      END IF
      SSDEN=DD*SS-DS*DS
      ITERC=0
      DENSAV=ZERO
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction.
!
   70 ITERC=ITERC+1
      TEMP=ONE/DSQRT(SSDEN)
      XOPTD=ZERO
      XOPTS=ZERO
      DO 80 I=1,N
      S(I)=TEMP*(DD*S(I)-DS*D(I))
      XOPTD=XOPTD+XOPT(I)*D(I)
   80 XOPTS=XOPTS+XOPT(I)*S(I)
!
!     Set the coefficients of the first two terms of BETA.
!
      TEMPA=HALF*XOPTD*XOPTD
      TEMPB=HALF*XOPTS*XOPTS
      DEN(1)=DD*(XOPTSQ+HALF*DD)+TEMPA+TEMPB
      DEN(2)=TWO*XOPTD*DD
      DEN(3)=TWO*XOPTS*DD
      DEN(4)=TEMPA-TEMPB
      DEN(5)=XOPTD*XOPTS
      DO 90 I=6,9
   90 DEN(I)=ZERO
!
!     Put the coefficients of Wcheck in WVEC.
!
      DO 110 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      TEMPC=ZERO
      DO 100 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
      TEMPB=TEMPB+XPT(K,I)*S(I)
  100 TEMPC=TEMPC+XPT(K,I)*XOPT(I)
      WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
      WVEC(K,2)=TEMPA*TEMPC
      WVEC(K,3)=TEMPB*TEMPC
      WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
  110 WVEC(K,5)=HALF*TEMPA*TEMPB
      DO 120 I=1,N
      IP=I+NPT
      WVEC(IP,1)=ZERO
      WVEC(IP,2)=D(I)
      WVEC(IP,3)=S(I)
      WVEC(IP,4)=ZERO
  120 WVEC(IP,5)=ZERO
!
!     Put the coefficents of THETA*Wcheck in PROD.
!
      DO 190 JC=1,5
      NW=NPT
      IF (JC .EQ. 2 .OR. JC .EQ. 3) NW=NDIM
      DO 130 K=1,NPT
  130 PROD(K,JC)=ZERO
      DO 150 J=1,NPTM
      SUM=ZERO
      DO 140 K=1,NPT
  140 SUM=SUM+ZMAT(K,J)*WVEC(K,JC)
      IF (J .LT. IDZ) SUM=-SUM
      DO 150 K=1,NPT
  150 PROD(K,JC)=PROD(K,JC)+SUM*ZMAT(K,J)
      IF (NW .EQ. NDIM) THEN
          DO 170 K=1,NPT
          SUM=ZERO
          DO 160 J=1,N
  160     SUM=SUM+BMAT(K,J)*WVEC(NPT+J,JC)
  170     PROD(K,JC)=PROD(K,JC)+SUM
      END IF
      DO 190 J=1,N
      SUM=ZERO
      DO 180 I=1,NW
  180 SUM=SUM+BMAT(I,J)*WVEC(I,JC)
  190 PROD(NPT+J,JC)=SUM
!
!     Include in DEN the part of BETA that depends on THETA.
!
      DO 210 K=1,NDIM
      SUM=ZERO
      DO 200 I=1,5
      PAR(I)=HALF*PROD(K,I)*WVEC(K,I)
  200 SUM=SUM+PAR(I)
      DEN(1)=DEN(1)-PAR(1)-SUM
      TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
      DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
      DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
      DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
      DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
      DEN(4)=DEN(4)-TEMPA-PAR(2)+PAR(3)
      TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
      DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
      DEN(8)=DEN(8)-PAR(4)+PAR(5)
      TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
  210 DEN(9)=DEN(9)-HALF*TEMPA
!
!     Extend DEN so that it holds all the coefficients of DENOM.
!
      SUM=ZERO
      DO 220 I=1,5
      PAR(I)=HALF*PROD(KNEW,I)**2
  220 SUM=SUM+PAR(I)
      DENEX(1)=ALPHA*DEN(1)+PAR(1)+SUM
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
      DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
      DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
      DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
      DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
      DENEX(4)=ALPHA*DEN(4)+TEMPA+PAR(2)-PAR(3)
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
      DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
      DENEX(8)=ALPHA*DEN(8)+PAR(4)-PAR(5)
      DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
!
!     Seek the value of the angle that maximizes the modulus of DENOM.
!
      SUM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
      DENOLD=SUM
      DENMAX=SUM
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      PAR(1)=ONE
      DO 250 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      PAR(2)=cos(ANGLE)
      PAR(3)=sin(ANGLE)
      DO 230 J=4,8,2
      PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
  230 PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
      SUMOLD=SUM
      SUM=ZERO
      DO 240 J=1,9
  240 SUM=SUM+DENEX(J)*PAR(J)
      IF (DABS(SUM) .GT. DABS(DENMAX)) THEN
          DENMAX=SUM
          ISAVE=I
          TEMPA=SUMOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=SUM
      END IF
  250 CONTINUE
      IF (ISAVE .EQ. 0) TEMPA=SUM
      IF (ISAVE .EQ. IU) TEMPB=DENOLD
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-DENMAX
          TEMPB=TEMPB-DENMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
!
!     Calculate the new parameters of the denominator, the new VLAG vector
!     and the new D. Then test for convergence.
!
      PAR(2)=cos(ANGLE)
      PAR(3)=sin(ANGLE)
      DO 260 J=4,8,2
      PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
  260 PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
      BETA=ZERO
      DENMAX=ZERO
      DO 270 J=1,9
      BETA=BETA+DEN(J)*PAR(J)
  270 DENMAX=DENMAX+DENEX(J)*PAR(J)
      DO 280 K=1,NDIM
      VLAG(K)=ZERO
      DO 280 J=1,5
  280 VLAG(K)=VLAG(K)+PROD(K,J)*PAR(J)
      TAU=VLAG(KNEW)
      DD=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      DO 290 I=1,N
      D(I)=PAR(2)*D(I)+PAR(3)*S(I)
      W(I)=XOPT(I)+D(I)
      DD=DD+D(I)**2
      TEMPA=TEMPA+D(I)*W(I)
  290 TEMPB=TEMPB+W(I)*W(I)
      IF (ITERC .GE. N) GOTO 340
      IF (ITERC .GT. 1) DENSAV=DMAX1(DENSAV,DENOLD)
      IF (DABS(DENMAX) .LE. 1.1D0*DABS(DENSAV)) GOTO 340
      DENSAV=DENMAX
!
!     Set S to half the gradient of the denominator with respect to D.
!     Then branch for the next iteration.
!
      DO 300 I=1,N
      TEMP=TEMPA*XOPT(I)+TEMPB*D(I)-VLAG(NPT+I)
  300 S(I)=TAU*BMAT(KNEW,I)+ALPHA*TEMP
      DO 320 K=1,NPT
      SUM=ZERO
      DO 310 J=1,N
  310 SUM=SUM+XPT(K,J)*W(J)
      TEMP=(TAU*W(N+K)-ALPHA*VLAG(K))*SUM
      DO 320 I=1,N
  320 S(I)=S(I)+TEMP*XPT(K,I)
      SS=ZERO
      DS=ZERO
      DO 330 I=1,N
      SS=SS+S(I)**2
  330 DS=DS+D(I)*S(I)
      SSDEN=DD*SS-DS*DS
      IF (SSDEN .GE. 1.0D-8*DD*SS) GOTO 70
!
!     Set the vector W before the RETURN from the subroutine.
!
  340 DO 350 K=1,NDIM
      W(K)=ZERO
      DO 350 J=1,5
  350 W(K)=W(K)+WVEC(K,J)*PAR(J)
      VLAG(KOPT)=VLAG(KOPT)+ONE
      RETURN
      END SUBROUTINE BIGDEN
      
      
      SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,
     1  DELTA,D,ALPHA,HCOL,GC,GD,S,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  HCOL(*),GC(*),GD(*),S(*),W(*)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DELTA is the current trust region bound.
!     D will be set to the step from XOPT to the new point.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     HCOL, GC, GD, S and W will be used for working space.
!
!     The step D is calculated in a way that attempts to maximize the modulus
!     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
!     the KNEW-th Lagrange function.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*atan(ONE)
      DELSQ=DELTA*DELTA
      NPTM=NPT-N-1
!
!     Set the first NPT components of HCOL to the leading elements of the
!     KNEW-th column of H.
!
      ITERC=0
      DO 10 K=1,NPT
   10 HCOL(K)=ZERO
      DO 20 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 20 K=1,NPT
   20 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
!
!     Set the unscaled initial direction D. Form the gradient of LFUNC at
!     XOPT, and multiply D by the second derivative matrix of LFUNC.
!
      DD=ZERO
      DO 30 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
      GC(I)=BMAT(KNEW,I)
      GD(I)=ZERO
   30 DD=DD+D(I)**2
      DO 50 K=1,NPT
      TEMP=ZERO
      SUM=ZERO
      DO 40 J=1,N
      TEMP=TEMP+XPT(K,J)*XOPT(J)
   40 SUM=SUM+XPT(K,J)*D(J)
      TEMP=HCOL(K)*TEMP
      SUM=HCOL(K)*SUM
      DO 50 I=1,N
      GC(I)=GC(I)+TEMP*XPT(K,I)
   50 GD(I)=GD(I)+SUM*XPT(K,I)
!
!     Scale D and GD, with a sign change if required. Set S to another
!     vector in the initial two dimensional subspace.
!
      GG=ZERO
      SP=ZERO
      DHD=ZERO
      DO 60 I=1,N
      GG=GG+GC(I)**2
      SP=SP+D(I)*GC(I)
   60 DHD=DHD+D(I)*GD(I)
      SCALE=DELTA/DSQRT(DD)
      IF (SP*DHD .LT. ZERO) SCALE=-SCALE
      TEMP=ZERO
      IF (SP*SP .GT. 0.99D0*DD*GG) TEMP=ONE
      TAU=SCALE*(DABS(SP)+HALF*SCALE*DABS(DHD))
      IF (GG*DELSQ .LT. 0.01D0*TAU*TAU) TEMP=ONE
      DO 70 I=1,N
      D(I)=SCALE*D(I)
      GD(I)=SCALE*GD(I)
   70 S(I)=GC(I)+TEMP*GD(I)
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction, except that termination occurs if
!     the given D and S are nearly parallel.
!
   80 ITERC=ITERC+1
      DD=ZERO
      SP=ZERO
      SS=ZERO
      DO 90 I=1,N
      DD=DD+D(I)**2
      SP=SP+D(I)*S(I)
   90 SS=SS+S(I)**2
      TEMP=DD*SS-SP*SP
      IF (TEMP .LE. 1.0D-8*DD*SS) GOTO 160
      DENOM=DSQRT(TEMP)
      DO 100 I=1,N
      S(I)=(DD*S(I)-SP*D(I))/DENOM
  100 W(I)=ZERO
!
!     Calculate the coefficients of the objective function on the circle,
!     beginning with the multiplication of S by the second derivative matrix.
!
      DO 120 K=1,NPT
      SUM=ZERO
      DO 110 J=1,N
  110 SUM=SUM+XPT(K,J)*S(J)
      SUM=HCOL(K)*SUM
      DO 120 I=1,N
  120 W(I)=W(I)+SUM*XPT(K,I)
      CF1=ZERO
      CF2=ZERO
      CF3=ZERO
      CF4=ZERO
      CF5=ZERO
      DO 130 I=1,N
      CF1=CF1+S(I)*W(I)
      CF2=CF2+D(I)*GC(I)
      CF3=CF3+S(I)*GC(I)
      CF4=CF4+D(I)*GD(I)
  130 CF5=CF5+S(I)*GD(I)
      CF1=HALF*CF1
      CF4=HALF*CF4-CF1
!
!     Seek the value of the angle that maximizes the modulus of TAU.
!
      TAUBEG=CF1+CF2+CF4
      TAUMAX=TAUBEG
      TAUOLD=TAUBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=cos(ANGLE)
      STH=sin(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      IF (DABS(TAU) .GT. DABS(TAUMAX)) THEN
          TAUMAX=TAU
          ISAVE=I
          TEMPA=TAUOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=TAU
      END IF
  140 TAUOLD=TAU
      IF (ISAVE .EQ. 0) TEMPA=TAU
      IF (ISAVE .EQ. IU) TEMPB=TAUBEG
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-TAUMAX
          TEMPB=TEMPB-TAUMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
!
!     Calculate the new D and GD. Then test for convergence.
!
      CTH=cos(ANGLE)
      STH=sin(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      DO 150 I=1,N
      D(I)=CTH*D(I)+STH*S(I)
      GD(I)=CTH*GD(I)+STH*W(I)
  150 S(I)=GC(I)+GD(I)
      IF (DABS(TAU) .LE. 1.1D0*DABS(TAUBEG)) GOTO 160
      IF (ITERC .LT. N) GOTO 80
  160 RETURN
      END SUBROUTINE BIGLAG
      
      SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  D,G,HD,HS,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*),
     1  D(*),G(*),HD(*),HS(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
!       in order to define the current quadratic model Q.
!     DELTA is the trust region radius, and has to be positive.
!     STEP will be set to the calculated trial step.
!     The arrays D, G, HD and HS will be used for working space.
!     CRVMIN will be set to the least curvature of H along the conjugate
!       directions that occur, except that it is set to zero if STEP goes
!       all the way to the trust region boundary.
!
!     The calculation of STEP begins with the truncated conjugate gradient
!     method. If the boundary of the trust region is reached, then further
!     changes to STEP may be made, each one being in the 2D space spanned
!     by the current STEP and the corresponding gradient of Q. Thus STEP
!     should provide a substantial reduction to Q within the trust region.
!
!     Initialization, which includes setting HD to H times XOPT.
!
      HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*atan(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      DO 10 I=1,N
   10 D(I)=XOPT(I)
      GOTO 170
!
!     Prepare for the first line search.
!
   20 QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I)=GQ(I)+HD(I)
      D(I)=-G(I)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
!
!     Calculate the step to the trust region boundary and the product HD.
!
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
!
!     Update CRVMIN and set the step-length ALPHA.
!
      ALPHA=BSTEP
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
!
!     Update STEP and HS.
!
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
!
!     Begin another conjugate direction iteration if required.
!
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 0.01D0*QRED) GOTO 160
          IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
          IF (ITERC .EQ. ITERMAX) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (DS .LE. ZERO) GOTO 160
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
!
!     Test whether an alternative iteration is required.
!
   90 IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99D0) GOTO 160
!
!     Begin the alternative iteration by calculating D and HD and some
!     scalar products.
!
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
!
!     Seek the value of the angle that minimizes Q.
!
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=cos(ANGLE)
      STH=sin(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
!
!     Calculate the new STEP and HS. Then test for convergence.
!
      CTH=cos(ANGLE)
      STH=sin(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
  160 RETURN
!
!     The following instructions act as a subroutine for setting the vector
!     HD to the vector D multiplied by the second derivative matrix of Q.
!     They are called from three different places, which are distinguished
!     by the value of ITERC.
!
  170 DO 180 I=1,N
  180 HD(I)=ZERO
      DO 200 K=1,NPT
      TEMP=ZERO
      DO 190 J=1,N
  190 TEMP=TEMP+XPT(K,J)*D(J)
      TEMP=TEMP*PQ(K)
      DO 200 I=1,N
  200 HD(I)=HD(I)+TEMP*XPT(K,I)
      IH=0
      DO 210 J=1,N
      DO 210 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
  210 HD(I)=HD(I)+HQ(IH)*D(J)
      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
      END SUBROUTINE TRSAPP

      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
!
!     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
!     interpolation point that has index KNEW. On entry, VLAG contains the
!     components of the vector Theta*Wcheck+e_b of the updating formula
!     (6.11), and BETA holds the value of the parameter that has this name.
!     The vector W is used for working space.
!
!     Set some constants.
!
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
      JL=1
      DO 20 J=2,NPTM
      IF (J .EQ. IDZ) THEN
          JL=IDZ
      ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN
          TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,JL)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 10 I=1,NPT
          TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
   10     ZMAT(I,JL)=TEMP
          ZMAT(KNEW,J)=ZERO
      END IF
   20 CONTINUE
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
      TEMPA=ZMAT(KNEW,1)
      IF (IDZ .GE. 2) TEMPA=-TEMPA
      IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL)
      DO 30 I=1,NPT
      W(I)=TEMPA*ZMAT(I,1)
      IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
   30 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
!
!     Complete the updating of ZMAT when there is only one nonzero element
!     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
!     then the first column of ZMAT will be exchanged with another one later.
!
      IFLAG=0
      IF (JL .EQ. 1) THEN
          TEMP=DSQRT(DABS(DENOM))
          TEMPB=TEMPA/TEMP
          TEMPA=TAU/TEMP
          DO 40 I=1,NPT
   40     ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
          IF (IDZ .EQ. 1 .AND. TEMP .LT. ZERO) IDZ=2
          IF (IDZ .GE. 2 .AND. TEMP .GE. ZERO) IFLAG=1
      ELSE
!
!     Complete the updating of ZMAT in the alternative case.
!
          JA=1
          IF (BETA .GE. ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*DSQRT(DABS(DENOM))
          DO 50 I=1,NPT
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
   50     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          IF (DENOM .LE. ZERO) THEN
              IF (BETA .LT. ZERO) IDZ=IDZ+1
              IF (BETA .GE. ZERO) IFLAG=1
          END IF
      END IF
!
!     IDZ is reduced in the following case, and usually the first column
!     of ZMAT is exchanged with a later one.
!
      IF (IFLAG .EQ. 1) THEN
          IDZ=IDZ-1
          DO 60 I=1,NPT
          TEMP=ZMAT(I,1)
          ZMAT(I,1)=ZMAT(I,IDZ)
   60     ZMAT(I,IDZ)=TEMP
      END IF
!
!     Finally, update the matrix BMAT.
!
      DO 70 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 70 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   70 CONTINUE
      RETURN
      END SUBROUTINE UPDATE
      
      

      end module mod_newuoa

