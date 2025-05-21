      module mod_bobyqa
      use const_def, only: dp

      contains

      subroutine do_BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN,max_valid_value)
      implicit REAL(dp) (A-H,O-Z)
      integer :: N, NPT, IPRINT, MAXFUN
      dimension X(:),XL(:),XU(:),W(*)
      interface
#include "num_bobyqa_proc.dek"
      end interface
      real(dp), intent(in) :: max_valid_value

!     This subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. There is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the Frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. The values of the variables are constrained by upper and
!     lower bounds. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!
!     subroutine CALFUN (N,X,F) has to be provided by the user. It must set
!     F to the value of the objective function for the current values of the
!     variables X(1),X(2),...,X(N), which are generated automatically in a
!     way that satisfies the bounds given in XL and XU.
!
!     Return if the value of NPT is unacceptable.
!
      NP=N+1
      if (NPT < N+2 .or. NPT > ((N+2)*NP)/2) then
          PRINT 10
   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in the required interval')
          GOTO 40
      end if
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
      ZERO=0.0D0
      do J=1,N
        TEMP=XU(J)-XL(J)
        if (TEMP < RHOBEG+RHOBEG) then
           PRINT 20
   20      FORMAT (/4X,'Return from BOBYQA because one of the differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
           GOTO 40
        end if
        JSL=ISL+J-1
        JSU=JSL+N
        W(JSL)=XL(J)-X(J)
        W(JSU)=XU(J)-X(J)
        if (W(JSL) >= -RHOBEG) then
            if (W(JSL) >= ZERO) then
                X(J)=XL(J)
                W(JSL)=ZERO
                W(JSU)=TEMP
            else
                X(J)=XL(J)+RHOBEG
                W(JSL)=-RHOBEG
                W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
            end if
        else if (W(JSU) <= RHOBEG) then
            if (W(JSU) <= ZERO) then
                X(J)=XU(J)
                W(JSL)=-TEMP
                W(JSU)=ZERO
            else
                X(J)=XU(J)-RHOBEG
                W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
                W(JSU)=RHOBEG
            end if
        end if
      end do
!
!     Make the call of BOBYQB.
!
      CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     &  W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     &  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW),CALFUN,
     &  max_valid_value)
   40 RETURN
      end subroutine do_BOBYQA


      subroutine BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     &  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     &  SL,SU,XNEW,XALT,D,VLAG,W,CALFUN,max_valid_value)
      implicit real(dp) (A-H,O-Z)
      integer :: N, NPT, IPRINT, MAXFUN, NDIM
      dimension X(:),XL(:),XU(:),XBASE(*),XPT(NPT,*),FVAL(*),
     &  XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     &  SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
      interface
#include "num_bobyqa_proc.dek"
      end interface
      real(dp), intent(in) :: max_valid_value

      logical :: do_replace
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in subroutine BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiteness.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the components of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by subroutine TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentioned.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a one-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      TEN=10.0D0
      TENTH=0.1D0
      TWO=2.0D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
      CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,
     &  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,CALFUN)
      XOPTSQ=ZERO
      do I=1,N
         XOPT(I)=XPT(KOPT,I)
         XOPTSQ=XOPTSQ+XOPT(I)**2
      end do
      FSAVE=FVAL(1)
      if (NF < NPT) then
          if (IPRINT > 0) PRINT 390
          GOTO 720
      end if
      KBASE=1
!
!     Complete the settings that are required for the iterative procedure.
!
      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      NFSAV=NF
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
   20 if (KOPT /= KBASE) then
          IH=0
          do J=1,N
            do I=1,J
                IH=IH+1
                if (I < J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
                GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
            end do
          end do
          if (NF > NPT) then
              do K=1,NPT
                TEMP=ZERO
                do J=1,N
                    TEMP=TEMP+XPT(K,J)*XOPT(J)
                end do
                TEMP=PQ(K)*TEMP
                do I=1,N
                    GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
                end do
              end do
          end if
      end if
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
     &  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      if (DNORM < HALF*RHO) then
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          if (NF <= NFSAV+2) GOTO 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          if (CRVMIN > ZERO .and. ERRBIG > FRHOSQ*CRVMIN) GOTO 650
          BDTOL=ERRBIG/RHO
          do J=1,N
            BDTEST=BDTOL
            if (XNEW(J) == SL(J)) BDTEST=W(J)
            if (XNEW(J) == SU(J)) BDTEST=-W(J)
            if (BDTEST < BDTOL) then
                CURV=HQ((J+J*J)/2)
                do K=1,NPT
                    CURV=CURV+PQ(K)*XPT(K,J)**2
                end do
                BDTEST=BDTEST+HALF*CURV*RHO
                if (BDTEST < BDTOL) GOTO 650
            end if
          end do
          GOTO 680
      end if
      NTRITS=NTRITS+1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
   90 if (DSQ <= 1.0D-3*XOPTSQ) then
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          do K=1,NPT
            SUMPQ=SUMPQ+PQ(K)
            SUM=-HALF*XOPTSQ
            do I=1,N
                SUM=SUM+XPT(K,I)*XOPT(I)
            end do
            W(NPT+K)=SUM
            TEMP=FRACSQ-HALF*SUM
              do I=1,N
                 W(I)=BMAT(K,I)
                 VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
                 IP=NPT+I
              do J=1,I
                 BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
              end do
            end do
          end do
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
          do JJ=1,NPTM
            SUMZ=ZERO
            SUMW=ZERO
            do K=1,NPT
                SUMZ=SUMZ+ZMAT(K,JJ)
                VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
                SUMW=SUMW+VLAG(K)
            end do
            do J=1,N
                SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
                do K=1,NPT
                    SUM=SUM+VLAG(K)*XPT(K,J)
                end do
                W(J)=SUM
                do K=1,NPT
                    BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
                end do
            end do
            do I=1,N
                IP=I+NPT
                TEMP=W(I)
                do J=1,I
                    BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
                end do
            end do
        end do
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
          IH=0
          do J=1,N
          W(J)=-HALF*SUMPQ*XOPT(J)
          do K=1,NPT
             W(J)=W(J)+PQ(K)*XPT(K,J)
             XPT(K,J)=XPT(K,J)-XOPT(J)
          end do
            do I=1,J
            IH=IH+1
            HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
            BMAT(NPT+I,J)=BMAT(NPT+J,I)
            end do
          end do
          do I=1,N
            XBASE(I)=XBASE(I)+XOPT(I)
            XNEW(I)=XNEW(I)-XOPT(I)
            SL(I)=SL(I)-XOPT(I)
            SU(I)=SU(I)-XOPT(I)
            XOPT(I)=ZERO
          end do
          XOPTSQ=ZERO
      end if
      if (NTRITS == 0) GOTO 210
      GOTO 230
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
  190 NFSAV=NF
      KBASE=KOPT
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,
     &  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,
     &  VLAG,W,W(N+NP),W(NDIM+NP),CALFUN)
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
      XOPTSQ=ZERO
      if (KOPT /= KBASE) then
          do I=1,N
             XOPT(I)=XPT(KOPT,I)
             XOPTSQ=XOPTSQ+XOPT(I)**2
          end do
      end if
      if (NF < 0) then
          NF=MAXFUN
          if (IPRINT > 0) PRINT 390
          GOTO 720
      end if
      NRESC=NF
      if (NFSAV < NF) then
          NFSAV=NF
          GOTO 20
      end if
      if (NTRITS > 0) GOTO 60
!
!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     &  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      do I=1,N
         D(I)=XNEW(I)-XOPT(I)
      end do
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
  230 do K=1,NPT
        SUMA=ZERO
        SUMB=ZERO
        SUM=ZERO
        do J=1,N
            SUMA=SUMA+XPT(K,J)*D(J)
            SUMB=SUMB+XPT(K,J)*XOPT(J)
            SUM=SUM+BMAT(K,J)*D(J)
        end do
        W(K)=SUMA*(HALF*SUMA+SUMB)
        VLAG(K)=SUM
        W(NPT+K)=SUMA
      end do
      BETA=ZERO
      do JJ=1,NPTM
        SUM=ZERO
        do K=1,NPT
            SUM=SUM+ZMAT(K,JJ)*W(K)
        end do
        BETA=BETA-SUM*SUM
        do K=1,NPT
            VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
        end do
      end do
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      do J=1,N
      DSQ=DSQ+D(J)**2
      SUM=ZERO
      do K=1,NPT
         SUM=SUM+W(K)*BMAT(K,J)
      end do
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      do I=1,N
         SUM=SUM+BMAT(JP,I)*D(I)
      end do
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
      DX=DX+D(J)*XOPT(J)
      end do
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
      if (NTRITS == 0) then
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          if (DENOM < CAUCHY .and. CAUCHY > ZERO) then
              do I=1,N
                XNEW(I)=XALT(I)
                D(I)=XNEW(I)-XOPT(I)
              end do
              CAUCHY=ZERO
              GOTO 230
          end if
          if (DENOM <= HALF*VLAG(KNEW)**2) then
              if (NF > NRESC) GOTO 190
              if (IPRINT > 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA because of much cancellation in a denominator.')
              GOTO 720
          end if
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
      else
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          do K=1,NPT
            if (K == KOPT) CYCLE
            HDIAG=ZERO
            do JJ=1,NPTM
               HDIAG=HDIAG+ZMAT(K,JJ)**2
            end do
            DEN=BETA*HDIAG+VLAG(K)**2
            DISTSQ=ZERO
            do J=1,N
                DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
            end do
            TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
            if (TEMP*DEN > SCADEN) then
                SCADEN=TEMP*DEN
                KNEW=K
                DENOM=DEN
            end if
            BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
          end do
          if (SCADEN <= HALF*BIGLSQ) then
              if (NF > NRESC) GOTO 190
              if (IPRINT > 0) PRINT 320
              GOTO 720
          end if
      end if
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
  360 do I=1,N
        X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
        if (XNEW(I) == SL(I)) X(I)=XL(I)
        if (XNEW(I) == SU(I)) X(I)=XU(I)
      end do
      if (NF >= MAXFUN) then
          if (IPRINT > 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been called MAXFUN times.')
          GOTO 720
      end if
      NF=NF+1
      CALL CALFUN (N,X,F)
      if (IPRINT == 3) then
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     &       '    The corresponding X is:'/(2X,5D15.6))
      end if
      if (NTRITS == -1) then
          FSAVE=F
          GOTO 720
      end if
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      do J=1,N
        VQUAD=VQUAD+D(J)*GOPT(J)
            do I=1,J
            IH=IH+1
            TEMP=D(I)*D(J)
            if (I == J) TEMP=HALF*TEMP
            VQUAD=VQUAD+HQ(IH)*TEMP
            end do
      end do
      do K=1,NPT
         VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
      end do
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      if (DNORM > RHO) NFSAV=NF
!
!     Pick the next value of DELTA after a trust region step.
!
      if (NTRITS > 0) then
          if (VQUAD >= ZERO) then
              if (IPRINT > 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA because a trust',
     &          ' region step has failed to reduce Q.')
              GOTO 720
          end if
          RATIO=(F-FOPT)/VQUAD
          if (RATIO <= TENTH) then
              DELTA=DMIN1(HALF*DELTA,DNORM)
          else if (RATIO. LE. 0.7D0) then
              DELTA=DMAX1(HALF*DELTA,DNORM)
          else
              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
          end if
          if (DELTA <= 1.5D0*RHO) DELTA=RHO
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
          if (F < FOPT) then
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              do K=1,NPT
                HDIAG=ZERO
                do JJ=1,NPTM
                    HDIAG=HDIAG+ZMAT(K,JJ)**2
                end do
                DEN=BETA*HDIAG+VLAG(K)**2
                DISTSQ=ZERO
                do J=1,N
                    DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
                end do
                TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
                if (TEMP*DEN > SCADEN) then
                    SCADEN=TEMP*DEN
                    KNEW=K
                    DENOM=DEN
                end if
                BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              end do
              if (SCADEN <= HALF*BIGLSQ) then
                  KNEW=KSAV
                  DENOM=DENSAV
              end if
          end if
      end if
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      IH=0
      PQOLD=PQ(KNEW)
      PQ(KNEW)=ZERO
      do I=1,N
         TEMP=PQOLD*XPT(KNEW,I)
         do J=1,I
         IH=IH+1
           HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
        end do
      end do
      do JJ=1,NPTM
        TEMP=DIFF*ZMAT(KNEW,JJ)
        do K=1,NPT
           PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
        end do
      end do
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
      FVAL(KNEW)=F
      do I=1,N
        XPT(KNEW,I)=XNEW(I)
        W(I)=BMAT(KNEW,I)
      end do
      do K=1,NPT
        SUMA=ZERO
        do JJ=1,NPTM
            SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
        end do
        SUMB=ZERO
        do J=1,N
           SUMB=SUMB+XPT(K,J)*XOPT(J)
        end do
        TEMP=SUMA*SUMB
        do I=1,N
            W(I)=W(I)+TEMP*XPT(K,I)
        end do
      end do
      do I=1,N
         GOPT(I)=GOPT(I)+DIFF*W(I)
      end do
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
      if (F < FOPT) then
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          do J=1,N
            XOPT(J)=XNEW(J)
            XOPTSQ=XOPTSQ+XOPT(J)**2
            do I=1,J
                IH=IH+1
                if (I < J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
                GOPT(I)=GOPT(I)+HQ(IH)*D(J)
            end do
          end do
          do K=1,NPT
            TEMP=ZERO
            do J=1,N
                TEMP=TEMP+XPT(K,J)*D(J)
            end do
            TEMP=PQ(K)*TEMP
            do I=1,N
                GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
            end do
          end do
      end if
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
      if (NTRITS > 0) then
          do K=1,NPT
            VLAG(K)=FVAL(K)-FVAL(KOPT)
            W(K)=ZERO
          end do
          do J=1,NPTM
            SUM=ZERO
            do K=1,NPT
                SUM=SUM+ZMAT(K,J)*VLAG(K)
            end do
          end do
          do K=1,NPT
             W(K)=W(K)+SUM*ZMAT(K,J)
          end do
          do K=1,NPT
            SUM=ZERO
            do J=1,N
                SUM=SUM+XPT(K,J)*XOPT(J)
            end do
            W(K+NPT)=W(K)
            W(K)=SUM*W(K)
          end do
          GQSQ=ZERO
          GISQ=ZERO
          do I=1,N
            SUM=ZERO
            do K=1,NPT
                SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
            end do
            if (XOPT(I) == SL(I)) then
                GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
                GISQ=GISQ+DMIN1(ZERO,SUM)**2
            else if (XOPT(I) == SU(I)) then
                GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
                GISQ=GISQ+DMAX1(ZERO,SUM)**2
            else
                GQSQ=GQSQ+GOPT(I)**2
                GISQ=GISQ+SUM*SUM
            end if
            VLAG(NPT+I)=SUM
          end do
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
          ITEST=ITEST+1
          if (GQSQ < TEN*GISQ) ITEST=0
          do_replace = (ITEST >= 3)
          if (.not. do_replace) then ! check for "invalid" value
            do k=1,npt
               if (fval(k) > max_valid_value) then
                  do_replace = .true.
                  !write(*,*) 'bobyqa value > max valid', fval(k), max_valid_value
                  exit
               end if
            end do
          end if
          if (do_replace) then
              !stop 'bobyqa: do_replace'
              do I=1,MAX0(NPT,NH)
                if (I <= N) GOPT(I)=VLAG(NPT+I)
                if (I <= NPT) PQ(I)=W(NPT+I)
                if (I <= NH) HQ(I)=ZERO
                ITEST=0
              end do
          end if
      end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
      if (NTRITS == 0) GOTO 60
      if (F <= FOPT+TENTH*VQUAD) GOTO 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      do K=1,NPT
        SUM=ZERO
        do J=1,N
            SUM=SUM+(XPT(K,J)-XOPT(J))**2
        end do
        if (SUM > DISTSQ) then
            KNEW=K
            DISTSQ=SUM
        end if
      end do
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
      if (KNEW > 0) then
          DIST=DSQRT(DISTSQ)
          if (NTRITS == -1) then
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              if (DELTA <= 1.5D0*RHO) DELTA=RHO
          end if
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      end if
      if (NTRITS == -1) GOTO 680
      if (RATIO > ZERO) GOTO 60
      if (DMAX1(DELTA,DNORM) > RHO) GOTO 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
  680 if (RHO > RHOEND) then
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          if (RATIO <= 16.0D0) then
              RHO=RHOEND
          else if (RATIO <= 250.0D0) then
              RHO=DSQRT(RATIO)*RHOEND
          else
              RHO=TENTH*RHO
          end if
          DELTA=DMAX1(DELTA,RHO)
          if (IPRINT >= 2) then
              if (IPRINT >= 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     &          ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     &          'The corresponding X is:'/(2X,5D15.6))
          end if
          NTRITS=0
          NFSAV=NF
          GOTO 60
      end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
      if (NTRITS == -1) GOTO 360
  720 if (FVAL(KOPT) <= FSAVE) then
          do I=1,N
            X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
            if (XOPT(I) == SL(I)) X(I)=XL(I)
            if (XOPT(I) == SU(I)) X(I)=XU(I)
          end do
          F=FVAL(KOPT)
      end if
      if (IPRINT >= 1) then
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
      end if
      RETURN
      end subroutine BOBYQB


      subroutine ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     &  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      implicit real(dp) (A-H,O-Z)
      integer :: N, NPT, NDIM, KOPT, KNEW
      dimension XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
     &  SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      CONST=ONE+DSQRT(2.0D0)
      do K=1,NPT
         HCOL(K)=ZERO
      end do
      do J=1,NPT-N-1
        TEMP=ZMAT(KNEW,J)
        do K=1,NPT
            HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
        end do
      end do
      ALPHA=HCOL(KNEW)
      HA=HALF*ALPHA
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
      do I=1,N
         GLAG(I)=BMAT(KNEW,I)
      end do
      do K=1,NPT
        TEMP=ZERO
        do J=1,N
            TEMP=TEMP+XPT(K,J)*XOPT(J)
        end do
        TEMP=HCOL(K)*TEMP
        do I=1,N
           GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
        end do
      end do
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
      PRESAV=ZERO
      do K=1,NPT
        if (K == KOPT) CYCLE
        DDERIV=ZERO
        DISTSQ=ZERO
        do I=1,N
            TEMP=XPT(K,I)-XOPT(I)
            DDERIV=DDERIV+GLAG(I)*TEMP
            DISTSQ=DISTSQ+TEMP*TEMP
        end do
        SUBD=ADELT/DSQRT(DISTSQ)
        SLBD=-SUBD
        ILBD=0
        IUBD=0
        SUMIN=DMIN1(ONE,SUBD)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
        do I=1,N
            TEMP=XPT(K,I)-XOPT(I)
            if (TEMP > ZERO) then
                if (SLBD*TEMP < SL(I)-XOPT(I)) then
                    SLBD=(SL(I)-XOPT(I))/TEMP
                    ILBD=-I
                end if
                if (SUBD*TEMP > SU(I)-XOPT(I)) then
                    SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
                    IUBD=I
                end if
            else if (TEMP < ZERO) then
                if (SLBD*TEMP > SU(I)-XOPT(I)) then
                    SLBD=(SU(I)-XOPT(I))/TEMP
                    ILBD=I
                end if
                if (SUBD*TEMP < SL(I)-XOPT(I)) then
                    SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
                    IUBD=-I
                end if
            end if
        end do
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
        if (K == KNEW) then
            DIFF=DDERIV-ONE
            STEP=SLBD
            VLAG=SLBD*(DDERIV-SLBD*DIFF)
            ISBD=ILBD
            TEMP=SUBD*(DDERIV-SUBD*DIFF)
            if (DABS(TEMP) > DABS(VLAG)) then
                STEP=SUBD
                VLAG=TEMP
                ISBD=IUBD
            end if
            TEMPD=HALF*DDERIV
            TEMPA=TEMPD-DIFF*SLBD
            TEMPB=TEMPD-DIFF*SUBD
            if (TEMPA*TEMPB < ZERO) then
                TEMP=TEMPD*TEMPD/DIFF
                if (DABS(TEMP) > DABS(VLAG)) then
                    STEP=TEMPD/DIFF
                    VLAG=TEMP
                    ISBD=0
                end if
            end if
!
!     Search along each of the other lines through XOPT and another point.
!
        else
            STEP=SLBD
            VLAG=SLBD*(ONE-SLBD)
            ISBD=ILBD
            TEMP=SUBD*(ONE-SUBD)
            if (DABS(TEMP) > DABS(VLAG)) then
                STEP=SUBD
                VLAG=TEMP
                ISBD=IUBD
            end if
            if (SUBD > HALF) then
                if (DABS(VLAG) < 0.25D0) then
                    STEP=HALF
                    VLAG=0.25D0
                    ISBD=0
                end if
            end if
            VLAG=VLAG*DDERIV
        end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
        TEMP=STEP*(ONE-STEP)*DISTSQ
        PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
        if (PREDSQ > PRESAV) then
            PRESAV=PREDSQ
            KSAV=K
            STPSAV=STEP
            IBDSAV=ISBD
        end if
      end do
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
      do I=1,N
         TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
         XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
      end do
      if (IBDSAV < 0) XNEW(-IBDSAV)=SL(-IBDSAV)
      if (IBDSAV > 0) XNEW(IBDSAV)=SU(IBDSAV)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
      BIGSTP=ADELT+ADELT
      IFLAG=0
  100 WFIXSQ=ZERO
      GGFREE=ZERO
      do I=1,N
        W(I)=ZERO
        TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
        TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
        if (TEMPA > ZERO .or. TEMPB < ZERO) then
            W(I)=BIGSTP
            GGFREE=GGFREE+GLAG(I)**2
        end if
      end do
      if (GGFREE == ZERO) then
          CAUCHY=ZERO
          GOTO 200
      end if
!
!     Investigate whether more components of W can be fixed.
!
  120 TEMP=ADELT*ADELT-WFIXSQ
      if (TEMP > ZERO) then
          WSQSAV=WFIXSQ
          STEP=DSQRT(TEMP/GGFREE)
          GGFREE=ZERO
          do I=1,N
            if (W(I) == BIGSTP) then
                TEMP=XOPT(I)-STEP*GLAG(I)
                if (TEMP <= SL(I)) then
                    W(I)=SL(I)-XOPT(I)
                    WFIXSQ=WFIXSQ+W(I)**2
                else if (TEMP >= SU(I)) then
                    W(I)=SU(I)-XOPT(I)
                    WFIXSQ=WFIXSQ+W(I)**2
                else
                    GGFREE=GGFREE+GLAG(I)**2
                end if
            end if
          end do
          if (WFIXSQ > WSQSAV .and. GGFREE > ZERO) GOTO 120
      end if
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
      GW=ZERO
      do I=1,N
        if (W(I) == BIGSTP) then
            W(I)=-STEP*GLAG(I)
            XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
        else if (W(I) == ZERO) then
            XALT(I)=XOPT(I)
        else if (GLAG(I) > ZERO) then
            XALT(I)=SL(I)
        else
            XALT(I)=SU(I)
        end if
        GW=GW+GLAG(I)*W(I)
      end do
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
      CURV=ZERO
      do K=1,NPT
        TEMP=ZERO
        do J=1,N
            TEMP=TEMP+XPT(K,J)*W(J)
        end do
        CURV=CURV+HCOL(K)*TEMP*TEMP
      end do
      if (IFLAG == 1) CURV=-CURV
      if (CURV > -GW .and. CURV < -CONST*GW) then
          SCALE=-GW/CURV
          do I=1,N
             TEMP=XOPT(I)+SCALE*W(I)
             XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
          end do
          CAUCHY=(HALF*GW*SCALE)**2
      else
          CAUCHY=(GW+HALF*CURV)**2
      end if
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
      if (IFLAG == 0) then
          do I=1,N
            GLAG(I)=-GLAG(I)
            W(N+I)=XALT(I)
          end do
          CSAVE=CAUCHY
          IFLAG=1
          GOTO 100
      end if
      if (CSAVE > CAUCHY) then
          do I=1,N
             XALT(I)=W(N+I)
          end do
          CAUCHY=CSAVE
      end if
  200 RETURN
      end subroutine ALTMOV


      subroutine PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,
     &  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,CALFUN)
      implicit real(dp) (A-H,O-Z)
      dimension X(:),XL(:),XU(:),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
     &  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
      interface
#include "num_bobyqa_proc.dek"
      end interface
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in subroutine BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintained as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     subroutine PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
      do J=1,N
        XBASE(J)=X(J)
        do K=1,NPT
            XPT(K,J)=ZERO
        end do
        do I=1,NDIM
            BMAT(I,J)=ZERO
        end do
      end do
      do IH=1,(N*NP)/2
         HQ(IH)=ZERO
      end do
      do K=1,NPT
        PQ(K)=ZERO
        do J=1,NPT-NP
            ZMAT(K,J)=ZERO
        end do
      end do
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      if (NFM <= 2*N) then
          if (NFM >= 1 .and. NFM <= N) then
              STEPA=RHOBEG
              if (SU(NFM) == ZERO) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          else if (NFM > N) then
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              if (SL(NFX) == ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
              if (SU(NFX) == ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          end if
      else
          ITEMP=(NFM-NP)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          if (IPT > N) then
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          end if
          XPT(NF,IPT)=XPT(IPT+1,IPT)
          XPT(NF,JPT)=XPT(JPT+1,JPT)
      end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
      do J=1,N
        X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
        if (XPT(NF,J) == SL(J)) X(J)=XL(J)
        if (XPT(NF,J) == SU(J)) X(J)=XU(J)
      end do
      CALL CALFUN (N,X,F)
      if (IPRINT == 3) then
         PRINT 70, NF,F,(X(I),I=1,N)
   70    FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))
      end if
      FVAL(NF)=F
      if (NF == 1) then
         FBEG=F
         KOPT=1
      else if (F < FVAL(KOPT)) then
         KOPT=NF
      end if
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
      if (NF <= 2*N+1) then
          if (NF >= 2 .and. NF <= N+1) then
              GOPT(NFM)=(F-FBEG)/STEPA
              if (NPT < NF+N) then
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              end if
          else if (NF >= N+2) then
              IH=(NFX*(NFX+1))/2
              TEMP=(F-FBEG)/STEPB
              DIFF=STEPB-STEPA
              HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
              GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
              if (STEPA*STEPB < ZERO) then
                  if (F < FVAL(NF-N)) then
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      if (KOPT == NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  end if
              end if
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          end if
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
      else
          IH=(IPT*(IPT-1))/2+JPT
          ZMAT(1,NFX)=RECIP
          ZMAT(NF,NFX)=RECIP
          ZMAT(IPT+1,NFX)=-RECIP
          ZMAT(JPT+1,NFX)=-RECIP
          TEMP=XPT(NF,IPT)*XPT(NF,JPT)
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
      end if
      if (NF < NPT .and. NF < MAXFUN) GOTO 50
      RETURN
      end subroutine PRELIM


      subroutine RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,
     &  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,
     &  KOPT,VLAG,PTSAUX,PTSID,W,CALFUN)
      implicit real(dp) (A-H,O-Z)
      dimension XL(:),XU(:),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
     &  GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
     &  VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
      interface
#include "num_bobyqa_proc.dek"
      end interface
!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NP=N+1
      SFRAC=HALF/DBLE(NP)
      NPTM=NPT-NP
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
      SUMPQ=ZERO
      WINC=ZERO
      do K=1,NPT
        DISTSQ=ZERO
        do J=1,N
            XPT(K,J)=XPT(K,J)-XOPT(J)
            DISTSQ=DISTSQ+XPT(K,J)**2
        end do
        SUMPQ=SUMPQ+PQ(K)
        W(NDIM+K)=DISTSQ
        WINC=DMAX1(WINC,DISTSQ)
         do J=1,NPTM
            ZMAT(K,J)=ZERO
         end do
      end do
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
      IH=0
      do J=1,N
        W(J)=HALF*SUMPQ*XOPT(J)
        do K=1,NPT
            W(J)=W(J)+PQ(K)*XPT(K,J)
        end do
        do I=1,J
            IH=IH+1
            HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
        end do
      end do
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
      do J=1,N
        XBASE(J)=XBASE(J)+XOPT(J)
        SL(J)=SL(J)-XOPT(J)
        SU(J)=SU(J)-XOPT(J)
        XOPT(J)=ZERO
        PTSAUX(1,J)=DMIN1(DELTA,SU(J))
        PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
        if (PTSAUX(1,J)+PTSAUX(2,J) < ZERO) then
            TEMP=PTSAUX(1,J)
            PTSAUX(1,J)=PTSAUX(2,J)
            PTSAUX(2,J)=TEMP
        end if
        if (DABS(PTSAUX(2,J)) < HALF*DABS(PTSAUX(1,J))) then
            PTSAUX(2,J)=HALF*PTSAUX(1,J)
        end if
        do I=1,NDIM
            BMAT(I,J)=ZERO
        end do
      end do
      FBASE=FVAL(KOPT)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
      PTSID(1)=SFRAC
      do J=1,N
        JP=J+1
        JPN=JP+N
        PTSID(JP)=DBLE(J)+SFRAC
        if (JPN <= NPT) then
            PTSID(JPN)=DBLE(J)/DBLE(NP)+SFRAC
            TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
            BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
            BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
            BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
            ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
            ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
            ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
        else
            BMAT(1,J)=-ONE/PTSAUX(1,J)
            BMAT(JP,J)=ONE/PTSAUX(1,J)
            BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
        end if
      end do
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
      if (NPT >= N+NP) then
          do K=2*NP,NPT
            IW=(DBLE(K-NP)-HALF)/DBLE(N)
            IP=K-NP-IW*N
            IQ=IP+IW
            if (IQ > N) IQ=IQ-N
            PTSID(K)=DBLE(IP)+DBLE(IQ)/DBLE(NP)+SFRAC
            TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
            ZMAT(1,K-NP)=TEMP
            ZMAT(IP+1,K-NP)=-TEMP
            ZMAT(IQ+1,K-NP)=-TEMP
          ZMAT(K,K-NP)=TEMP
          end do
      end if
      NREM=NPT
      KOLD=1
      KNEW=KOPT
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
   80 do J=1,N
        TEMP=BMAT(KOLD,J)
        BMAT(KOLD,J)=BMAT(KNEW,J)
        BMAT(KNEW,J)=TEMP
      end do
      do J=1,NPTM
        TEMP=ZMAT(KOLD,J)
        ZMAT(KOLD,J)=ZMAT(KNEW,J)
        ZMAT(KNEW,J)=TEMP
      end do
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      if (KNEW /= KOPT) then
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          if (NREM == 0) GOTO 350
          do K=1,NPT
             W(NDIM+K)=DABS(W(NDIM+K))
          end do
      end if
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
  120 DSQMIN=ZERO
      do K=1,NPT
        if (W(NDIM+K) > ZERO) then
            if (DSQMIN == ZERO .or. W(NDIM+K) < DSQMIN) then
                KNEW=K
                DSQMIN=W(NDIM+K)
            end if
        end if
      end do
      if (DSQMIN == ZERO) GOTO 260
!
!     Form the W-vector of the chosen original interpolation point.
!
      do J=1,N
         W(NPT+J)=XPT(KNEW,J)
      end do
      do K=1,NPT
      SUM=ZERO
      if (K == KOPT) then
        ! do nothing
      else if (PTSID(K) == ZERO) then
          do J=1,N
             SUM=SUM+W(NPT+J)*XPT(K,J)
          end do
      else
          IP=PTSID(K)
          if (IP > 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DBLE(NP)*PTSID(K)-DBLE(IP*NP)
          if (IQ > 0) then
              IW=1
              if (IP == 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          end if
      end if
         W(K)=HALF*SUM*SUM
      end do
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
      do K=1,NPT
        SUM=ZERO
        do J=1,N
            SUM=SUM+BMAT(K,J)*W(NPT+J)
        end do
        VLAG(K)=SUM
      end do
      BETA=ZERO
      do J=1,NPTM
        SUM=ZERO
        do K=1,NPT
            SUM=SUM+ZMAT(K,J)*W(K)
        end do
        BETA=BETA-SUM*SUM
        do K=1,NPT
           VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
        end do
      end do
      BSUM=ZERO
      DISTSQ=ZERO
      do J=1,N
        SUM=ZERO
        do K=1,NPT
            SUM=SUM+BMAT(K,J)*W(K)
        end do
        JP=J+NPT
        BSUM=BSUM+SUM*W(JP)
        do IP=NPT+1,NDIM
            SUM=SUM+BMAT(IP,J)*W(IP)
        end do
        BSUM=BSUM+SUM*W(JP)
        VLAG(JP)=SUM
        DISTSQ=DISTSQ+XPT(KNEW,J)**2
      end do
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
      DENOM=ZERO
      VLMXSQ=ZERO
      do K=1,NPT
        if (PTSID(K) /= ZERO) then
            HDIAG=ZERO
            do J=1,NPTM
                HDIAG=HDIAG+ZMAT(K,J)**2
            end do
            DEN=BETA*HDIAG+VLAG(K)**2
            if (DEN > DENOM) then
                KOLD=K
                DENOM=DEN
            end if
        end if
        VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      end do
      if (DENOM <= 1.0D-2*VLMXSQ) then
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          GOTO 120
      end if
      GOTO 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
  260 do KPT=1,NPT
        if (PTSID(KPT) == ZERO) CYCLE
        if (NF >= MAXFUN) then
            NF=-1
            GOTO 350
        end if
        IH=0
        do J=1,N
            W(J)=XPT(KPT,J)
            XPT(KPT,J)=ZERO
            TEMP=PQ(KPT)*W(J)
            do I=1,J
                IH=IH+1
                HQ(IH)=HQ(IH)+TEMP*W(I)
            end do
        end do
        PQ(KPT)=ZERO
        IP=PTSID(KPT)
        IQ=DBLE(NP)*PTSID(KPT)-DBLE(IP*NP)
        if (IP > 0) then
            XP=PTSAUX(1,IP)
            XPT(KPT,IP)=XP
        end if
        if (IQ > 0) then
            XQ=PTSAUX(1,IQ)
            if (IP == 0) XQ=PTSAUX(2,IQ)
            XPT(KPT,IQ)=XQ
        end if
!
!     Set VQUAD to the value of the current model at the new point.
!
        VQUAD=FBASE
        if (IP > 0) then
            IHP=(IP+IP*IP)/2
            VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
        end if
        if (IQ > 0) then
            IHQ=(IQ+IQ*IQ)/2
            VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
            if (IP > 0) then
                IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
                VQUAD=VQUAD+XP*XQ*HQ(IW)
            end if
        end if
        do K=1,NPT
            TEMP=ZERO
            if (IP > 0) TEMP=TEMP+XP*XPT(K,IP)
            if (IQ > 0) TEMP=TEMP+XQ*XPT(K,IQ)
            VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
        end do
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
        do I=1,N
            W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
            if (XPT(KPT,I) == SL(I)) W(I)=XL(I)
            if (XPT(KPT,I) == SU(I)) W(I)=XU(I)
        end do
        NF=NF+1
        CALL CALFUN (N,W,F)
        if (IPRINT == 3) then
            PRINT 300, NF,F,(W(I),I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))
        end if
        FVAL(KPT)=F
        if (F < FVAL(KOPT)) KOPT=KPT
        DIFF=F-VQUAD
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
        do I=1,N
            GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
        end do
        do K=1,NPT
            SUM=ZERO
            do J=1,NPTM
            SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
            end do
            TEMP=DIFF*SUM
            if (PTSID(K) == ZERO) then
                PQ(K)=PQ(K)+TEMP
            else
                IP=PTSID(K)
                IQ=DBLE(NP)*PTSID(K)-DBLE(IP*NP)
                IHQ=(IQ*IQ+IQ)/2
                if (IP == 0) then
                    HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
                else
                    IHP=(IP*IP+IP)/2
                    HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
                    if (IQ > 0) then
                        HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                        IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                        HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
                    end if
                end if
            end if
        end do
        PTSID(KPT)=ZERO
      end do
  350 RETURN
      end subroutine RESCUE


      subroutine TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,
     &  XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
      implicit real(dp) (A-H,O-Z)
      dimension XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
     &  XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are zero.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0D0 is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      ONEMIN=-1.0D0
      ZERO=0.0D0
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
      ITERC=0
      NACT=0
      SQSTP=ZERO
      do I=1,N
        XBDI(I)=ZERO
        if (XOPT(I) <= SL(I)) then
            if (GOPT(I) >= ZERO) XBDI(I)=ONEMIN
        else if (XOPT(I) >= SU(I)) then
            if (GOPT(I) <= ZERO) XBDI(I)=ONE
        end if
        if (XBDI(I) /= ZERO) NACT=NACT+1
        D(I)=ZERO
        GNEW(I)=GOPT(I)
      end do
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
   20 BETA=ZERO
   30 STEPSQ=ZERO
      do I=1,N
        if (XBDI(I) /= ZERO) then
            S(I)=ZERO
        else if (BETA == ZERO) then
            S(I)=-GNEW(I)
        else
            S(I)=BETA*S(I)-GNEW(I)
        end if
        STEPSQ=STEPSQ+S(I)**2
      end do
      if (STEPSQ == ZERO) GOTO 190
      if (BETA == ZERO) then
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      end if
      if (GREDSQ*DELSQ <= 1.0D-4*QRED*QRED) GOTO 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
      GOTO 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      do I=1,N
        if (XBDI(I) == ZERO) then
            RESID=RESID-D(I)**2
            DS=DS+S(I)*D(I)
            SHS=SHS+S(I)*HS(I)
        end if
      end do
      if (RESID <= ZERO) GOTO 90
      TEMP=DSQRT(STEPSQ*RESID+DS*DS)
      if (DS < ZERO) then
          BLEN=(TEMP-DS)/STEPSQ
      else
          BLEN=RESID/(TEMP+DS)
      end if
      STPLEN=BLEN
      if (SHS > ZERO) then
          STPLEN=DMIN1(BLEN,GREDSQ/SHS)
      end if

!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
      IACT=0
      do I=1,N
        if (S(I) /= ZERO) then
            XSUM=XOPT(I)+D(I)
            if (S(I) > ZERO) then
                TEMP=(SU(I)-XSUM)/S(I)
            else
                TEMP=(SL(I)-XSUM)/S(I)
            end if
            if (TEMP < STPLEN) then
                STPLEN=TEMP
                IACT=I
            end if
        end if
      end do
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
      SDEC=ZERO
      if (STPLEN > ZERO) then
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          if (IACT == 0 .and. TEMP > ZERO) then
              CRVMIN=DMIN1(CRVMIN,TEMP)
              if (CRVMIN == ONEMIN) CRVMIN=TEMP
          end if
          GGSAV=GREDSQ
          GREDSQ=ZERO
          do I=1,N
            GNEW(I)=GNEW(I)+STPLEN*HS(I)
            if (XBDI(I) == ZERO) GREDSQ=GREDSQ+GNEW(I)**2
            D(I)=D(I)+STPLEN*S(I)
          end do
          SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      end if
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
      if (IACT > 0) then
          NACT=NACT+1
          XBDI(IACT)=ONE
          if (S(IACT) < ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          if (DELSQ <= ZERO) GOTO 90
          GOTO 20
      end if
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
      if (STPLEN < BLEN) then
          if (ITERC == ITERMAX) GOTO 190
          if (SDEC <= 0.01D0*QRED) GOTO 190
          BETA=GREDSQ/GGSAV
          GOTO 30
      end if
   90 CRVMIN=ZERO
!
!     Prepare for the alternative iteration by calculating some scalars
!     and by multiplying the reduced D by the second derivative matrix of
!     Q, where S holds the reduced D in the call of GGMULT.
!
  100 if (NACT >= N-1) GOTO 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      do I=1,N
        if (XBDI(I) == ZERO) then
            DREDSQ=DREDSQ+D(I)**2
            DREDG=DREDG+D(I)*GNEW(I)
            GREDSQ=GREDSQ+GNEW(I)**2
            S(I)=D(I)
        else
            S(I)=ZERO
        end if
      end do
      ITCSAV=ITERC
      GOTO 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      if (TEMP <= 1.0D-4*QRED*QRED) GOTO 190
      TEMP=DSQRT(TEMP)
      do I=1,N
        if (XBDI(I) == ZERO) then
            S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
        else
            S(I)=ZERO
        end if
      end do
      SREDG=-TEMP
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
      ANGBD=ONE
      IACT=0
      do I=1,N
        if (XBDI(I) == ZERO) then
            TEMPA=XOPT(I)+D(I)-SL(I)
            TEMPB=SU(I)-XOPT(I)-D(I)
            if (TEMPA <= ZERO) then
                NACT=NACT+1
                XBDI(I)=ONEMIN
                GOTO 100
            else if (TEMPB <= ZERO) then
                NACT=NACT+1
                XBDI(I)=ONE
                GOTO 100
            end if
            RATIO=ONE
            SSQ=D(I)**2+S(I)**2
            TEMP=SSQ-(XOPT(I)-SL(I))**2
            if (TEMP > ZERO) then
                TEMP=DSQRT(TEMP)-S(I)
                if (ANGBD*TEMP > TEMPA) then
                    ANGBD=TEMPA/TEMP
                    IACT=I
                    XSAV=ONEMIN
                end if
            end if
            TEMP=SSQ-(SU(I)-XOPT(I))**2
            if (TEMP > ZERO) then
                TEMP=DSQRT(TEMP)+S(I)
                if (ANGBD*TEMP > TEMPB) then
                    ANGBD=TEMPB/TEMP
                    IACT=I
                    XSAV=ONE
                end if
            end if
        end if
      end do
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
      GOTO 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      do I=1,N
        if (XBDI(I) == ZERO) then
            SHS=SHS+S(I)*HS(I)
            DHS=DHS+D(I)*HS(I)
            DHD=DHD+D(I)*HRED(I)
        end if
      end do
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0D0*ANGBD+3.1D0
      do I=1,IU
        ANGT=ANGBD*DBLE(I)/DBLE(IU)
        STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
        TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
        REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
        if (REDNEW > REDMAX) then
            REDMAX=REDNEW
            ISAV=I
            RDPREV=REDSAV
        else if (I == ISAV+1) then
            RDNEXT=REDNEW
        end if
        REDSAV=REDNEW
      end do
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
      if (ISAV == 0) GOTO 190
      if (ISAV < IU) then
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DBLE(ISAV)+HALF*TEMP)/DBLE(IU)
      end if
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      if (SDEC <= ZERO) GOTO 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
      DREDG=ZERO
      GREDSQ=ZERO
      do I=1,N
        GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
        if (XBDI(I) == ZERO) then
            D(I)=CTH*D(I)+STH*S(I)
            DREDG=DREDG+D(I)*GNEW(I)
            GREDSQ=GREDSQ+GNEW(I)**2
        end if
        HRED(I)=CTH*HRED(I)+STH*HS(I)
      end do
      QRED=QRED+SDEC
      if (IACT > 0 .and. ISAV == IU) then
          NACT=NACT+1
          XBDI(IACT)=XSAV
          GOTO 100
      end if
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
      if (SDEC > 0.01D0*QRED) GOTO 120
  190 DSQ=ZERO
      do I=1,N
        XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
        if (XBDI(I) == ONEMIN) XNEW(I)=SL(I)
        if (XBDI(I) == ONE) XNEW(I)=SU(I)
        D(I)=XNEW(I)-XOPT(I)
        DSQ=DSQ+D(I)**2
      end do
      RETURN

!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
  210 IH=0
      do J=1,N
        HS(J)=ZERO
        do I=1,J
            IH=IH+1
            if (I < J) HS(J)=HS(J)+HQ(IH)*S(I)
            HS(I)=HS(I)+HQ(IH)*S(J)
        end do
      end do
      do K=1,NPT
        if (PQ(K) /= ZERO) then
            TEMP=ZERO
            do J=1,N
               TEMP=TEMP+XPT(K,J)*S(J)
            end do
            TEMP=TEMP*PQ(K)
            do I=1,N
               HS(I)=HS(I)+TEMP*XPT(K,I)
            end do
        end if
      end do
      if (CRVMIN /= ZERO) GOTO 50
      if (ITERC > ITCSAV) GOTO 150
      do I=1,N
         HRED(I)=HS(I)
      end do
      GOTO 120
      end subroutine TRSBOX


      subroutine UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      implicit real(dp) (A-H,O-Z)
      dimension BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
!
!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      ZTEST=ZERO
      do K=1,NPT
        do J=1,NPTM
            ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
        end do
      end do
      ZTEST=1.0D-20*ZTEST
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
      JL=1
      do J=2,NPTM
        if (DABS(ZMAT(KNEW,J)) > ZTEST) then
            TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
            TEMPA=ZMAT(KNEW,1)/TEMP
            TEMPB=ZMAT(KNEW,J)/TEMP
            do I=1,NPT
                TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
                ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
                ZMAT(I,1)=TEMP
            end do
        end if
        ZMAT(KNEW,J)=ZERO
      end do

!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.

      do I=1,NPT
         W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
      end do
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE

!     Complete the updating of ZMAT.

      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      do I=1,NPT
         ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
      end do

!     Finally, update the matrix BMAT.

      do J=1,N
        JP=NPT+J
        W(JP)=BMAT(KNEW,J)
        TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
        TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
        do I=1,JP
            BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
            if (I > NPT) BMAT(JP,I-NPT)=BMAT(I,J)
        end do
      end do
      RETURN
      end subroutine UPDATE

      end module mod_bobyqa
