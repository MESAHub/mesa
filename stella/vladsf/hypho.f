      SUBROUTINE HYPHO(ZED,N,LMIN,LMAX,NFREQ,FREQ,ANL)
C
C          GENERATE HYDROGENIC PHOTOIONIZATION CROSS SECTIONS ANL
C          INITIALIZE BY SETTING ZED.LE.0.0D0
C
C          INPUT:
C
C               ZED = NUCLEAR CHARGE
C                 N = PRINCIPAL QUANTUM NUMBER OF LEVEL
C         LMIN,LMAX = MINIMUM AND MAXIMUM OF RANGE OF L TO BE COMPUTED
C                     WHERE L IS ANGULAR MOMENTUM QUANTUM NUMBER OF LEVEL.
C                     VALUES OF L OUTSIDE THE NATURAL RANGE 0 TO N-1 WILL
C                     NOT BE COMPUTED.
C             NFREQ = NUMBER OF FREQUENCIES IN ARRAY FREQ
C          FREQ(IF) = ARRAY OF FREQUENCIES IN UNITS OF RYD*ZED**2 AT WHICH 
C                     ANL IS TO BE COMPUTED, IF=1,NFREQ.  THESE FREQUENCIES 
C                     SHOULD BE IN ASCENDING ORDER, STARTING WITH THE 
C                     THRESHOLD FREQUENCY.
C
C          OUTPUT:
C
C         ANL(IL,IF) = ARRAY OF CROSS SECTIONS, IL=LMIN+1,LMAX+1 AND IF=1,NFREQ
C   
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      DIMENSION G2(100,2),G3(100,2),MM(100),AB(100),ALO(1500),FAL(1500),
     ~
     .FREQ(100),ANL(50,100)


       
C
C          INITIALIZATION (ZED.LE.0.0D0)
C
      IF(ZED.GT.0.0D0) GOTO 10
C
C         DOUBLE PRECISION CONSTANTS
C
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      FOUR=4.0D0
      TEN=1.0D1
      HALF=ONE/TWO
      CON1=8.5594D-19
C
C          EVALUATE LOG FACTORIALS
C
      FAL(1) = ZERO
      DO  I = 2,1500
           AI = dble(I - 1)
           ALO(I-1)=DLOG10(AI)
           FAL(I) = ALO(I-1) + FAL(I-1)
      ENDDO
      INIT=999
      RETURN
C
C          CALCULATION OF HYDROGENIC PHOTOIONIZATION CROSS SECTION
C
10    IF(INIT.NE.999) THEN
           WRITE(6,15)
15         FORMAT(16H NOT INITIALIZED//)
           RETURN
      ENDIF
C
C         INITIALIZATION FOR N
C
      FN =dble(N)
      SN = FN*FN
      SN4 = FOUR*SN
      CON2=CON1*(FN/ZED)**2
      FTH=ONE/SN
      GN0=2.3052328943D0-2.302585093D0*FAL(N+N)-FN*0.61370563888D0+
     .ALO(N)*(FN+ONE)*2.30258093D0
      LMAX1=MIN(LMAX,N-1)
      ILMAX=N-LMIN
C
C         INITIALIZE G'S
C
      DO  I = 1,NFREQ
           MM(I) = 0
           DO  J = 1,2
                G2(I,J) = ZERO
                G3(I,J) = ZERO
           ENDDO
      ENDDO
C
C          L LOOP
C
      DO IL=1,ILMAX
           L=N-IL
           LIN=L-LMIN+1
           M=0
           AL = dble(L)
           K = N - L - 2
           CON3=CON2/(TWO*AL+ONE)
C
C          FREQUENCY LOOP (FREQ UNITS ARE RYD*ZED**2)
C
           DO 70 IF=1,NFREQ
                IF(FREQ(IF).LT.FTH) THEN
                     IF(L.LE.LMAX1) ANL(LIN,IF)=ZERO
                     GOTO 70
                ENDIF
                M=M+1
                SE=FREQ(IF)-FTH
                E=DSQRT(SE)
                X=ONE+SN*SE
                IF(K.LT.0) THEN
                     IF (E.GE.0.314D0) THEN
                          EE=6.2831853D0/E
                          P=0.5D0*DLOG10(EXP1(-EE))
                     ELSE
                          P = ZERO
                     ENDIF
                     IF(E.GT.1.0D-6) THEN
                          A=TWO*(FN-DATAN(FN*E)/E)
                     ELSE
                          A=ZERO
                     ENDIF
                     AB(M)=(GN0+A)/2.302585D0-P-(FN+TWO)*DLOG10(X)
                     GNE = 0.1D0
                     GN1E = X*GNE/(FN + FN)
                     G3(M,2) = GNE
                     G3(M,1) = GN1E
                     G2(M,2) = GNE*FN*X*(FN+FN-ONE)
                     G2(M,1) = GN1E*(FN+FN-ONE)*(FOUR+(FN-ONE)*X)
                ENDIF
                G22 = G2(M,2)
                G32 = G3(M,2)
                G21 = G2(M,1)
                G31 = G3(M,1)
                MULS = MM(M)
                IF(K)20,30,25
C
C                     L.LT.N-2
C
25              IF(K.GT.1) GOTO 28
                LL = N - 1
                LM = N - 2
28              SL = dble(LL*LL)
                SL4 = FOUR*SL
                FLL = dble(LL)
                G12 = (SN4 - SL4 + (TWO*SL - FLL)*X)*G22 -
     .          SN4*(SN - SL)*(ONE + (FLL + ONE)**2*SE)*G32
                IF (L.EQ.0) GO TO 12
                SM = dble(LM*LM)
                SM4 = FOUR*SM
                FLM = dble(LM)
                G11 = (SN4 - SM4 + (TWO*SM + FLM)*X)*G21 -
     .          SN4*(SN - (FLM + ONE)**2)*(ONE + SM*SE)*G31
                G31 = G21
                G21 = G11
12              G32 = G22
                G22 = G12
                IF (M.NE.NFREQ) GO TO 31
                LL = LL - 1
                IF (L.EQ.0) GO TO 31
                LM = LM - 1
31              CONTINUE
                IF (G12.LT.1.D20) GO TO 40
                MULS = MULS + 35
                G22 = G22*1.D-35
                G32 = G32*1.D-35
                G12 = G12*1.D-35
                IF (L.EQ.0) GO TO 40
                G11 = G11*1.D-35
                G21 = G21*1.D-35
                G31 = G31*1.D-35
                GO TO 40
C
C                    L.EQ.N-1
C
20              G11 = G31
                IF (L.EQ.0) G11 = ZERO
                G12 = G32
                GO TO 40
C
C                    L.EQ.N-2
C
30              G11 = G21
                IF (L.EQ.0) G11 = ZERO
                G12 = G22
40              CONTINUE
                MM(M) = MULS
                G2(M,2) = G22
                G3(M,2) = G32
                G2(M,1) = G21
                G3(M,1) = G31
                ALFAC = FAL(N+L+1) - FAL(N-L) + TWO*(AL-FN)*ALO(2*N)
                P1 = ONE
                LLL = L + 1
                LLM = L - 1
                MULR = 0
                IF (LLM.LT.1) GO TO 46
                DO  I = 1,LLM
                     AI = dble(I)
                     P1 = P1*(ONE + AI*AI*SE)
                     IF (P1.GE.1.D20) THEN
                          P1 = P1*1.D-10
                          MULR = MULR + 10
                     ENDIF
                ENDDO
46              P2 = P1
                LLK = LLM + 1
                IF (LLK.LT.1) LLK = 1
                DO  I = LLK,LLL
                     AI = dble(I)
                     P2 = P2*(ONE + AI*AI*SE)
                ENDDO
                MULP = 0
17             IF (G12.LT.ONE) GO TO 117
               MULP = MULP - 10
               G12 = G12*1.D-10
               IF (L.EQ.0) GO TO 17
               G11 = G11*1.D-10
               GO TO 17
117            SUM=ALFAC+dble(MULR)+TWO*(AB(M)+dble(MULS-MULP+1))
               FAC = ZERO
               IF(DABS(SUM).LT.50.D0) FAC=TEN**SUM
71             IF (L.EQ.0) GO TO 18
               G11 = G11*P1*G11
18             G12 = G12*P2*G12
               IF(L.LE.LMAX1) ANL(LIN,IF)=CON3*(G11 *AL + G12 *(AL+ONE))
     .         *X*FAC
70         ENDDO
      ENDDO
      RETURN
      END
      DOUBLE PRECISION FUNCTION EXP1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DX=DABS(X)
      IF(DX.LT.1.0D-9) GOTO1
      IF(DX.LT.1.0D-5) GOTO2
      IF(DX.LT.1.0D-3) GOTO3
      EXP1=1.0D0-DEXP(X)
      RETURN
1     EXP1=-X
      RETURN
2     EXP1=((-X*0.5D0)-1.0D0)*X
      RETURN
3     EXP1=(((-X*0.1666666667D0)-0.5D0)*X-1.0D0)*X
      RETURN
      END
