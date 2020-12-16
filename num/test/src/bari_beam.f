c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Beam (ODE case)
c        ODE of dimension 80
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/beam.f
c
c     This is revision
c     $Id: beam.F,v 1.2 2006/10/02 10:29:13 testset Exp $
c
c-----------------------------------------------------------------------
      subroutine beam_init(neqn,y,yprime,consis)
      integer neqn
      double precision y(neqn),yprime(neqn)
      logical consis

      integer i

      do 10 i=1,neqn
         y(i) = 0d0
   10 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine beam_feval(nvar,t,th,df,ierr,rpar,ipar)
      use math_lib
      IMPLICIT REAL*8 (A-H,O-Z)
      integer ierr,ipar(*)
      integer, parameter :: N=40, NN=2*N, NCOM=N, NSQ=N*N, NQUATR=NSQ*NSQ, NNCOM=NN
      double precision rpar(*), an, deltas
        DIMENSION DF(NN),TH(150),U(150),V(150),W(150)
        DIMENSION ALPHA(150),BETA(150),STH(150),CTH(150)
C --- SET DEFAULT VALUES
      if (nvar /= nn) stop 'bad nvar for beam_feval'
      AN=N
      DELTAS=1.0D+0/AN
C ----- CALCUL DES TH(I) ET DES SIN ET COS -------------
        DO 22 I=2,N
        THDIFF=TH(I)-TH(I-1)
        STH(I)=sin(THDIFF)
   22   CTH(I)=cos(THDIFF)
C -------- CALCUL DU COTE DROIT DU SYSTEME LINEAIRE -----
        IF(T.GT.3.14159265358979324D0)THEN
C --------- CASE T GREATER PI ------------
C ---------- I=1 ------------
            TERM1=(-3.D0*TH(1)+TH(2))*NQUATR
            V(1)=TERM1
C -------- I=2,..,N-1 -----------
            DO 32 I=2,N-1
            TERM1=(TH(I-1)-2.D0*TH(I)+TH(I+1))*NQUATR
  32        V(I)=TERM1
C ----------- I=N -------------
            TERM1=(TH(N-1)-TH(N))*NQUATR
            V(N)=TERM1
        ELSE
C --------- CASE T LESS EQUAL PI ------------
            FABS=1.5D0*sin(T)*sin(T)
            FX=-FABS
            FY= FABS
C ---------- I=1 ------------
            TERM1=(-3.D0*TH(1)+TH(2))*NQUATR
            TERM2=NSQ*(FY*cos(TH(1))-FX*sin(TH(1)))
            V(1)=TERM1+TERM2
C -------- I=2,..,N-1 -----------
            DO 34 I=2,N-1
            TERM1=(TH(I-1)-2.D0*TH(I)+TH(I+1))*NQUATR
            TERM2=NSQ*(FY*cos(TH(I))-FX*sin(TH(I)))
  34        V(I)=TERM1+TERM2
C ----------- I=N -------------
            TERM1=(TH(N-1)-TH(N))*NQUATR
            TERM2=NSQ*(FY*cos(TH(N))-FX*sin(TH(N)))
            V(N)=TERM1+TERM2
        END IF
C -------- COMPUTE PRODUCT DV=W ------------
        W(1)=STH(2)*V(2)
        DO 43 I=2,N-1
   43   W(I)=-STH(I)*V(I-1)+STH(I+1)*V(I+1)
        W(N)=-STH(N)*V(N-1)
C -------- TERME 3 -----------------
        DO 435 I=1,N
  435   W(I)=W(I)+TH(N+I)*TH(N+I)
C ------- SOLVE SYSTEM CW=W ---------
        ALPHA(1)=1.D0
        DO 44 I=2,N
        ALPHA(I)=2.D0
   44   BETA(I-1)=-CTH(I)
        ALPHA(N)=3.D0
        DO 45 I=N-1,1,-1
        Q=BETA(I)/ALPHA(I+1)
        W(I)=W(I)-W(I+1)*Q
   45   ALPHA(I)=ALPHA(I)-BETA(I)*Q
        W(1)=W(1)/ALPHA(1)
        DO 46 I=2,N
   46   W(I)=(W(I)-BETA(I-1)*W(I-1))/ALPHA(I)
C -------- COMPUTE U=CV+DW ---------
        U(1)=V(1)-CTH(2)*V(2)+STH(2)*W(2)
        DO 47 I=2,N-1
   47   U(I)=2.D0*V(I)-CTH(I)*V(I-1)-CTH(I+1)*V(I+1)
     &                -STH(I)*W(I-1)+STH(I+1)*W(I+1)
        U(N)=3.D0*V(N)-CTH(N)*V(N-1)-STH(N)*W(N-1)
C -------- PUT  DERIVATIVES IN RIGHT PLACE -------------
        DO 54 I=1,N
        DF(I)=TH(N+I)
  54    DF(N+I)=U(I)
        RETURN
        END
c-----------------------------------------------------------------------
      subroutine beam_jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)
c
c     dummy subroutine
c
      return
      end
c-----------------------------------------------------------------------
      subroutine beam_solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c
c computed using double precision RADAU on an 
c     Alphaserver DS20E, with a 667 MHz EV67 processor.
c          
c          uround = 1.01d-19
c          rtol = atol = h0 = 1.1d-18
c
c
      y(  1) = -0.5792366591285007D-002             
      y(  2) = -0.1695298550721735D-001              
      y(  3) = -0.2769103312973094D-001             
      y(  4) = -0.3800815655879158D-001             
      y(  5) = -0.4790616859743763D-001            
      y(  6) = -0.5738710435274594D-001            
      y(  7) = -0.6645327313454617D-001             
      y(  8) = -0.7510730581979037D-001             
      y(  9) = -0.8335219765414992D-001      
      y( 10) = -0.9119134654647947D-001          
      y( 11) = -0.9862858700132091D-001            
      y( 12) = -0.1056682200378002D+000            
      y( 13) = -0.1123150395409595D+000           
      y( 14) = -0.1185743552727245D+000             
      y( 15) = -0.1244520128755561D+000           
      y( 16) = -0.1299544113264161D+000           
      y( 17) = -0.1350885180610398D+000             
      y( 18) = -0.1398618819194422D+000             
      y( 19) = -0.1442826441015242D+000             
      y( 20) = -0.1483595472463012D+000           
      y( 21) = -0.1521019429001447D+000            
      y( 22) = -0.1555197978061129D+000             
      y( 23) = -0.1586236993420229D+000             
      y( 24) = -0.1614248603702127D+000           
      y( 25) = -0.1639351238193223D+000             
      y( 26) = -0.1661669673440852D+000           
      y( 27) = -0.1681335081778558D+000              
      y( 28) = -0.1698485080602243D+000              
      y( 29) = -0.1713263782440888D+000              
      y( 30) = -0.1725821847462537D+000            
      y( 31) = -0.1736316537975654D+000              
      y( 32) = -0.1744911773840049D+000             
      y( 33) = -0.1751778187863392D+000             
      y( 34) = -0.1757093178712902D+000             
      y( 35) = -0.1761040960228576D+000           
      y( 36) = -0.1763812607175549D+000              
      y( 37) = -0.1765606097564671D+000            
      y( 38) = -0.1766626352260565D+000              
      y( 39) = -0.1767085270807460D+000            
      y( 40) = -0.1767201761075488D+000              
      y( 41) =  0.3747362681329794D-001              
      y( 42) =  0.1099117880217593D+000            
      y( 43) =  0.1798360474312799D+000            
      y( 44) =  0.2472427305442391D+000             
      y( 45) =  0.3121293820596567D+000            
      y( 46) =  0.3744947377019500D+000           
      y( 47) =  0.4343386073492798D+000           
      y( 48) =  0.4916620354601748D+000           
      y( 49) =  0.5464677854586807D+000            
      y( 50) =  0.5987609702624270D+000            
      y( 51) =  0.6485493611110740D+000           
      y( 52) =  0.6958435169132503D+000            
      y( 53) =  0.7406572668520808D+000            
      y( 54) =  0.7830081747813241D+000            
      y( 55) =  0.8229176659201515D+000            
      y( 56) =  0.8604110305190560D+000             
      y( 57) =  0.8955175502377805D+000            
      y( 58) =  0.9282708263127953D+000           
      y( 59) =  0.9587089334522034D+000       
      y( 60) =  0.9868747821728363D+000             
      y( 61) =  0.1012816579961883D+001              
      y( 62) =  0.1036587736679858D+001            
      y( 63) =  0.1058246826481355D+001            
      y( 64) =  0.1077857811433353D+001            
      y( 65) =  0.1095490222005369D+001            
      y( 66) =  0.1111219164294898D+001             
      y( 67) =  0.1125125269286501D+001              
      y( 68) =  0.1137294526609229D+001            
      y( 69) =  0.1147818025153607D+001             
      y( 70) =  0.1156792132004482D+001           
      y( 71) =  0.1164318845130183D+001            
      y( 72) =  0.1170505992596124D+001             
      y( 73) =  0.1175467424299550D+001            
      y( 74) =  0.1179323003228859D+001             
      y( 75) =  0.1182198586299667D+001           
      y( 76) =  0.1184226111223146D+001             
      y( 77) =  0.1185543909805575D+001             
      y( 78) =  0.1186297084203716D+001            
      y( 79) =  0.1186637618908127D+001            
      y( 80) =  0.1186724615113034D+001              
 


      return
      end
