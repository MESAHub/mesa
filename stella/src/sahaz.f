      SUBROUTINEURSOSZ
      IMPlICITREAL*8(A-H,O-Z)
      PARAMETER(NVARS=3)
      include '../obj/nfreq_and_mzone.inc'
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
      Parameter(Is=5)
      PARAMETER(NZ=3000000)
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
      LogicalLSYSTEM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+2
     *3,AMbrun=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-10,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551
     *000000D+10,ULGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP
     *=1.0000000000D+05,URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-05,ERGEV=1.6021764864D-12,GRADeV=1.16045052
     *85D+04,RADC=7.5657680191D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
      LogicalRadP
      COMMON/RadP/RadP
      COMMON/ARG/TP,PL,CHEM,LST,KENTR,JURS
      COMMON/RESULT/P,Egas,Sgas,ENG,CAPPA,PT,ET,ST,ENGT,CAPT,NZR
      COMMON/ABUND/XYZA,Yat(Natom)
      COMMON/AZ/AS,ZS,SCN
      COMMON/STR/PPL,EPL,SPL,ENGPL,CAPPL,CP,GAM,DA,DPE,DSE,DSP,BETgas
      COMMON/XELECT/XE,XET,XEPL,PE,Ycomp
      COMMON/URScap/Tpsqrt,Psicap,Scap,ScapT,ScapPl,ZMean,YZMean,ZMT,ZMPl,YZMT,YZMPl
      Parameter(tret=1.d0/3.d0)
      Parameter(tret2=2.d0/3.d0)
      Parameter(epsurs=1.d-8,Xeflo=1.d-20)
      Parameter(epslow=1.d-20)
      DimensionYal(Natom),Yaldz(Natom),al(Natom),aldT(Natom),sigYal(Natom),AD1(Natom),AD2(Natom),AR1(Natom),AR2(Natom)
      Parameter(Nstagez=2)
      DimensionXion(Nstagez,Natom)
      LogicalLOW,HIGH
      COMMON/NIT/Xelow,Xion,Nit,LOW,HIGH
      Parameter(CPSI=3.0176945495D-04,CPRESS=1.0035847556D-04,CRGAS=8.3141678120D-01,CAvD=6.0220450000D+17,CSaha=1.5270993016D+23,DL
     *OW=1.0000000000D+05,DUP=1.0000000000D+08,CWW=1.4476482730D-01,CLOW=1.6666666667D+00,PIONUP=6.3151191931D+00)
      Parameter(Zf=1.2000000000D+01,Tpzb=3.8868074375D+00,Tpzf=1.4821399533D+02)
      DIMENSIONPION(Nstagez,Natom)
      DIMENSIONPsum(Nstagez,Natom)
      DIMENSIONGWL(Nstagez,Natom)
      DIMENSIONRGWR(Nstagez,Natom)
      Character*2Chames(Natom)
      DATAPION/1.578025D+00,0.000000D+00,2.853397D+00,6.315119D+00,1.306704D+00,2.829607D+00,1.686647D+00,3.435147D+00,1.580346D+00,
     *4.075154D+00,2.502467D+00,4.753573D+00,5.963725D-01,5.487462D+00,8.873057D-01,1.744787D+00,6.946654D-01,2.184958D+00,9.460262D
     *-01,1.896926D+00,1.828806D+00,3.206300D+00,5.037659D-01,3.670029D+00,7.094036D-01,1.377610D+00,9.167560D-01,1.878537D+00,8.862
     *941D-01,2.108423D+00/
      DATAPsum/1.578025D+00,1.578025D+00,2.853397D+00,9.168516D+00,1.306704D+00,4.136311D+00,1.686647D+00,5.121794D+00,1.580346D+00,
     *5.655500D+00,2.502467D+00,7.256039D+00,5.963725D-01,6.083835D+00,8.873057D-01,2.632093D+00,6.946654D-01,2.879624D+00,9.460262D
     *-01,2.842953D+00,1.828806D+00,5.035106D+00,5.037659D-01,4.173795D+00,7.094036D-01,2.087014D+00,9.167560D-01,2.795293D+00,8.862
     *941D-01,2.994717D+00/
      DATAGWL/-6.931472D-01,0.000000D+00,6.931472D-01,-6.931472D-01,-4.054651D-01,-1.791759D+00,8.109302D-01,-4.054651D-01,-8.109302
     *D-01,8.109302D-01,1.791759D+00,4.054651D-01,-6.931472D-01,1.791759D+00,6.931472D-01,-6.931472D-01,-1.791759D+00,6.931472D-01,-
     *4.054651D-01,-1.791759D+00,1.791759D+00,4.054651D-01,-6.931472D-01,1.791759D+00,6.931472D-01,-6.931472D-01,1.823216D-01,-1.823
     *216D-01,-7.419373e-01,7.419373D-01/
      DATARGWR/2.000000D+00,1.000000D+00,5.000000D-01,2.000000D+00,1.500000D+00,6.000000D+00,4.444444D-01,1.500000D+00,2.250000D+00,
     *4.444444D-01,1.666667D-01,6.666667D-01,2.000000D+00,1.666667D-01,5.000000D-01,2.000000D+00,6.000000D+00,5.000000D-01,1.500000D
     *+00,6.000000D+00,1.666667D-01,6.666667D-01,2.000000D+00,1.666667D-01,5.000000D-01,2.000000D+00,8.333333D-01,1.200000D+00,2.100
     *000d+00,4.761905d-01/
      DataChames/'H ','He','C ','N ','O ','Ne','Na','Mg','Al','Si','Ar','K ','Ca','Fe','Ni'/
      DataAD2(1)/0.d0/
      GOTO09993
09999 CONTINUE
      Dsaha=log(Csaha*Tpsqrt*Rnbar)
      DLPOT=DSaha+1.d0
      If((Tp*DLPOT).GT.Tpzb)then
      If((Tp*DLPOT)*(1.D0+Epsurs).LT.Tpzf)then
      ZMean=Zf*((Tp*DLPOT)-Tpzb)/(Tpzf-Tpzb)
      PotZM=(Tpzf-Tpzb)*ZMean/Zf+Tpzb
      PotInt=(Tpzf-Tpzb)*ZMean**2/(2.d0*Zf)+Tpzb*ZMean
      If(LST.NE.0)Then
      ZMT=Zf*(DLPOT+1.5D0)/(Tpzf-Tpzb)
      ZMPl=-Zf*Tp/(Tpzf-Tpzb)
      WZMT=ZMT/(Zf-ZMean)
      WZMPl=ZMPl/(Zf-ZMean)
      XeZMT=YZ*ZMT
      XeZMPl=ZMPl*YZ
      endif
      else
      ZMean=Zf
      PotInt=(0.5d0*Zf)*(Tpzf+Tpzb)
      ZMT=0.D0
      WZMT=0.D0
      ZMPl=0.D0
      WZMPl=0.D0
      XeZMT=0.D0
      XeZMPl=0.D0
      endif
      XeZM=ZMean*YZ
      Egas=CRgas*YZ*PotInt
      If(LST.NE.0)Then
      PT=0.D0
      ET=CRgas*XeZMT*PotZM
      PPl=0.D0
      EPl=CRgas*XeZMPl*PotZM
      endif
      YZMean=XeZM/Zf
      else
      ZMean=0.D0
      YZMean=0.D0
      XeZM=0.D0
      Egas=0.D0
      If(LST.NE.0)Then
      XeZMT=0.D0
      ZMT=0.D0
      WZMT=0.D0
      XeZMPl=0.D0
      ZMPl=0.D0
      WZMPl=0.D0
      PT=0.D0
      ET=0.D0
      PPl=0.D0
      EPl=0.D0
      endif
      endif
      If(Zf-ZMean.GT.Epsurs*Zf)then
      If(exp(Dsaha-Pionup*RTp)*epsurs.GT.1.d0)then
      LOW=.FALSE.
      HIGH=.TRUE.
      Xe=XeZM+Yat(1)
      DO09990i=2,Natom
      Xe=Xe+dble(Nstagez)*Yat(i)
09990 CONTINUE
      else
      Asum=0.d0
      DO09987i=1,Natom
      AD1(i)=exp(Dsaha-Pion(1,i)*RTp+GWl(1,i))
      Asum=Asum+AD1(i)*Yat(i)
09987 CONTINUE
      If(Asum.LT.EpsLOW)Then
      LOW=.TRUE.
      HIGH=.FALSE.
      Xe=Sqrt(Asum)
      XeLOW=Sqrt(Asum)
      else
      Xe=CHEM
      LOW=.FALSE.
      HIGH=.FALSE.
      DO09984i=2,Natom
      AD2(i)=AD1(i)*exp(Dsaha-Pion(2,i)*RTp+GWl(2,i))
09984 CONTINUE
      endif
      RXe=1.d0/(max(Xe,Xeflo))
      endif
      else
      LOW=.false.
      HIGH=.true.
      Xe=Yat(1)+XeZM
      DO09981i=2,Natom
      Xe=Xe+2.D0*Yat(i)
09981 CONTINUE
      Scap=Yat(1)*Rgwr(1,1)+Yat(2)*Rgwr(2,2)*4.d0
      ScapT=0.D0
      ScapPl=0.D0
      endif
      NIT=0
      If(.NOT.LOW.AND..NOT.HIGH)Then
      RXeold=0.d0
09978 IF(.NOT.(ABS(Rxeold-Rxe).GT.epsurs*Rxe.AND.NIT.LT.64))GOTO09977
      RXeold=RXe
      AR1(1)=AD1(1)*Rxe
      alpha=1.d0/(1.d0+AR1(1))
      Yal(1)=Yat(1)*alpha
      al(1)=alpha
      eta=AR1(1)*alpha
      U=Yat(1)*eta
      Yaldz(1)=U
      V=U-U*eta
      DO09975i=2,Natom
      AR1(i)=AD1(i)*RXe
      AR2(i)=(AD2(i)*RXe)*RXe
      alpha=1.d0/(1.d0+AR2(i)+AR1(i))
      al(i)=alpha
      Yal(i)=Yat(i)*alpha
      dzeta=2.D0*AR2(i)+AR1(i)
      eta=2.D0*AR2(i)+dzeta
      Yaldz(i)=Yal(i)*dzeta
      U=U+Yaldz(i)
      V=V+Yal(i)*(eta-dzeta*(dzeta*alpha))
09975 CONTINUE
      RXe=abs((Rxeold*V+1.D0)/(V+U+XeZM))
      NIT=NIT+1
      GOTO09978
09977 CONTINUE
      If(NIT.GE.64)then
      Write(6,*)' Iterations in SAHAZ fail to converge'
      Write(6,*)' Tp=',Tp,'  Pl=',Pl,'   Rxeold=',Rxeold,'   Rxe=',Rxe,'  epsurs=',epsurs
      STOP
      endif
      Xe=1.d0/Rxe
      Endif
      Xe=max(Xe,Xeflo)
      If(HIGH)Then
      Xion(1,1)=1.D0
      Xion(1,2)=0.D0
      Xion(2,2)=1.D0
      Scap=Yat(1)*Rgwr(1,1)+Yat(2)*Rgwr(2,2)*4.d0
      ScapZZ=0.D0
      DO09972i=3,Natom
      ScapZZ=ScapZZ+Yat(i)*(1.d0-ZMean/Zf)*Rgwr(2,i)*4.d0
      Xion(1,i)=0.D0
      Xion(2,i)=1.D0
09972 CONTINUE
      Scap=Scap+ScapZZ
      elseif(LOW)Then
      Xion(1,1)=AD1(1)*RXe
      Scap=0.D0
      DO09969i=1,Natom
      Scap=Scap+Yat(i)*AD1(i)*Rgwr(1,i)
      Xion(1,i)=AD1(i)*RXe
      Xion(2,i)=0.D0
09969 CONTINUE
      Scap=Scap*RXe
      else
      Xion(1,1)=Ar1(1)*al(1)
      Xion(1,2)=Ar1(2)*al(2)
      Xion(2,2)=Ar2(2)*al(2)
      Scap=Yal(1)*Ar1(1)*Rgwr(1,1)+Yal(2)*(Ar1(2)*Rgwr(1,2)+Ar2(2)*Rgwr(2,2)*4.d0)
      ScapZZ=0.D0
      DO09966i=3,Natom
      ScapZZ=ScapZZ+Yal(i)*(1.d0-ZMean/Zf)*(Ar1(i)*Rgwr(1,i)+Ar2(i)*Rgwr(2,i)*4.d0)
      Xion(1,i)=Ar1(i)*al(i)
      Xion(2,i)=Ar2(i)*al(i)
09966 CONTINUE
      Scap=Scap+ScapZZ
      endif
      P=CRgas*(XYZA+Xe)*Pl*Tp
      If(LOW)Then
      WEgas=0.d0
      DO09963i=1,Natom
      WEgas=WEgas+AD1(i)*Psum(1,i)*Yat(i)*RXe
09963 CONTINUE
      Elseif(HIGH)Then
      WEgas=Psum(1,1)*Yat(1)+Psum(Nstagez,2)*Yat(2)
      EZZ=0.D0
      DO09960i=3,Natom
      EZZ=EZZ+Psum(Nstagez,i)*Yat(i)
09960 CONTINUE
      WEgas=WEgas+EZZ
      Else
      WEgas=Psum(1,1)*AR1(1)*Yal(1)
      sigYal(1)=WEgas
      DO09957i=2,Natom
      W=(Psum(2,i)*AR2(i)+Psum(1,i)*AR1(I))*Yal(i)
      WEgas=WEgas+W
      sigYal(i)=W
09957 CONTINUE
      Endif
      Egas=Egas+CRgas*(1.5d0*(XYZA+Xe)*Tp+WEgas)
      IF(LST.NE.0)THEN
      If(HIGH)Then
      XeT=XeZMT
      XePl=XeZMPl
      ET=ET+CRgas*(1.5D0*(XYZA+Xe+Tp*XeT))
      EPl=EPl+CRgas*(1.5D0*Tp*XePl)
      ScapT=-ScapZZ*WZMT
      ScapPl=-ScapZZ*WZMPl
      elseif(LOW)Then
      XeT=0.D0
      DO09954i=1,Natom
      XeT=XeT+Yat(i)*(1.5D0+Psum(1,i)*RTp)*AD1(i)
09954 CONTINUE
      XeT=0.5d0*XeT*(RXe*RTp)
      XePl=0.5d0*Asum*RXe-Xe
      XeTl=XeT*RXe
      DET=0.D0
      DePl=0.d0
      ScapT=0.D0
      ScapPl=0.D0
      DO09951i=1,Natom
      DET=DET+Yat(i)*(AD1(i)*RXe)*Psum(1,i)*((1.5D0+Psum(1,i)*RTp)*RTp-XeTl)
      DEPl=DEPl+Yat(i)*(AD1(i)*RXe)*Psum(1,i)
      ScapT=ScapT+Yat(i)*(AD1(i)*RXe)*RGWR(1,i)*((1.5D0+Psum(1,i)*RTp)*RTp-XeTl)
      ScapPl=ScapPl-Yat(i)*(AD1(i)*RXe)*RGWR(1,i)
09951 CONTINUE
      ET=ET+CRgas*(1.5D0*(XYZA+Xe+Tp*XeT)+DET)
      EPl=EPl+CRgas*(1.5D0*Tp*XePl-DePl*(XePl*RXe+1.d0))
      ScapPl=ScapPl*(XePl*RXe+1.D0)
      else
      alpT=(1.5D0+Psum(1,1)*RTp)*AR1(1)*RTp
      DGT=(Yat(1)-Yaldz(1))*alpT*al(1)+XeZMT
      aldT(1)=alpT
      DO09948i=2,Natom
      alpT=(((3.D0+Psum(2,i)*RTp)*AR2(i))+((1.5D0+Psum(1,i)*RTp)*AR1(i)))*RTp
      dzT=(2.D0*((3.D0+Psum(2,i)*RTp)*AR2(i))+((1.5D0+Psum(1,i)*RTp)*AR1(i)))*RTp
      DGT=DGT+Yal(i)*dzT-Yaldz(i)*alpT*al(i)
      aldT(i)=alpT
09948 CONTINUE
      XeT=DGT/(V*RXe+1.D0)
      XePl=(U+XeZM+XeZMPl)/(V*RXe+1.D0)-Xe
      XeTl=XeT*RXe
      sigmT=(aldT(1)-XeTl*AR1(1))*Psum(1,1)
      DET=Yal(1)*sigmT-sigYal(1)*(aldT(1)-XeTl*AR1(1))*al(1)
      signe=AR1(1)*Psum(1,1)
      DEPl=Yal(1)*signe-sigYal(1)*AR1(1)*al(1)
      ScapT=Yal(1)*RGWR(1,1)*(aldT(1)-AR1(1)*((aldT(1)-XeTl*AR1(1))*al(1)+XeTl))
      ScapPl=Yal(1)*RGWR(1,1)*AR1(1)*(AR1(1)*al(1)-1.D0)
      DO09945i=2,Natom
      sigmT=(Psum(2,i)*AR2(i))*((3.D0+Psum(2,i)*RTp)*RTp-2.d0*XeTl)+(Psum(1,i)*AR1(i))*((1.5D0+Psum(1,i)*RTp)*RTp-XeTl)
      ScapT=ScapT+Yal(i)*((((3.D0+Psum(2,i)*RTp)*AR2(i))*RGWR(2,i)*4.D0+((1.5D0+Psum(1,i)*RTp)*AR1(i))*RGWR(1,i))*RTp-(Ar1(i)*Rgwr(1
     *,i)+Ar2(i)*Rgwr(2,i)*4.d0)*(aldT(i)-XeTl*(2.D0*AR2(i)+AR1(i)))*al(i)-XeTl*(Ar1(i)*Rgwr(1,i)+Ar2(i)*Rgwr(2,i)*8.d0))
      DET=DET+Yal(i)*sigmT-sigYal(i)*(aldT(i)-XeTl*(2.D0*AR2(i)+AR1(i)))*al(i)
      signe=2.d0*AR2(i)*Psum(2,i)+AR1(i)*Psum(1,i)
      DEPl=DEPl+Yal(i)*signe-sigYal(i)*(2.D0*AR2(i)+AR1(i))*al(i)
      ScapPl=ScapPl+Yal(i)*((2.D0*AR2(i)+AR1(i))*al(i)*(Ar1(i)*Rgwr(1,i)+Ar2(i)*Rgwr(2,i)*4.d0)-(Ar1(i)*Rgwr(1,i)+Ar2(i)*Rgwr(2,i)*8
     *.d0))
09945 CONTINUE
      ET=ET+CRgas*(1.5D0*(XYZA+Xe+Tp*XeT)+DET)
      EPl=EPl+CRgas*(1.5D0*Tp*XePl-DePl*(XePl*RXe+1.d0))
      ScapT=ScapT-ScapZZ*WZMT
      ScapPl=ScapPl*(XePl*RXe+1.D0)-ScapZZ*WZMPl
      endif
      PPl=PPl+CRgas*(XYZA+Xe+XePl)*Tp
      PT=PT+CRgas*(XYZA+Xe+Tp*XeT)*Pl
      YZMT=XeZMT/Zf
      YZMPl=XeZMPl/Zf
      ENDIF
      CHEM=Xe
      GOTO(09997,09994),I6R001
09993 CONTINUE
      GOTO09942
09998 CONTINUE
      ZMean=Zf
      XeZM=ZMean*YZ
      YZMean=YZ
      Scap=Yat(1)*Rgwr(1,1)+Yat(2)*Rgwr(2,2)*4.d0
      Ye=Yat(1)+2.D0*Yat(2)+XeZM
      dY=(Pl*Ye)**tret2
      PSI=Cpsi*Dy*RTp
      P=CPRESS*(Pl*Ye)*dY+CRgas*(XYZA+Ye/(1.D0+0.4D0*PSI))*Pl*Tp
      PotInt=(0.5d0*Zf)*(Tpzf+Tpzb)
      Egas=1.5D0*P*RPl
      If(LST.NE.0)Then
      PPl=(5.D0*tret)*CPRESS*Ye*dY+CRgas*(XYZA+Ye/(1.D0+0.4D0*PSI))*Tp-CRgas*Ye/(1.D0+0.4D0*PSI)**2*(0.4D0*tret2)*PSI*Tp
      PT=CRgas*(XYZA+Ye/(1.D0+0.4D0*PSI))*Pl+CRgas*Ye/(1.D0+0.4D0*PSI)**2*0.4D0*PSI*Pl
      EPl=1.5D0*PPl-Egas
      ET=1.5D0*PT*RPl
      ZMT=0.D0
      ZMPl=0.D0
      YZMT=0.D0
      YZMPl=0.D0
      ScapT=0.D0
      ScapPl=0.D0
      Endif
      Egas=Egas+CRgas*(YZ*PotInt+Psum(1,1)*Yat(1)+Psum(Nstagez,2)*Yat(2))
      Xion(1,1)=1.D0
      DO09939i=2,Natom
      Xion(1,i)=0.D0
      Xion(2,i)=1.D0
09939 CONTINUE
      GOTO(09996,09995),I6R002
09942 CONTINUE
      RPl=1.d0/Pl
      RTp=1.d0/Tp
      Tpsqrt=Tp*sqrt(Tp)
      Barnum=Pl*CAvD
      Rnbar=RPl*(1.d0/CAvD)
      YZ=0.D0
      DO09936I=3,Natom
      YZ=YZ+Yat(i)
09936 CONTINUE
      XYZA=YZ+Yat(1)+Yat(2)
      IF(.NOT.(Pl.LE.DLOW))GOTO09933
      I6R001=1
      GOTO09999
09997 CONTINUE
      Psicap=0.d0
      GOTO09931
09933 CONTINUE
      IF(.NOT.(Pl.GE.DUP))GOTO09932
      I6R002=1
      GOTO09998
09996 CONTINUE
      Xe=Ye
      Psicap=Psi
      GOTO09931
09932 CONTINUE
      WW=CWW*LOG(Pl)-CLOW
      WD=CWW*RPl
      WP=WW**3*(4.D0-3.D0*WW)
      I6R002=2
      GOTO09998
09995 CONTINUE
      Pdeg=P
      Edeg=Egas
      Scapdg=Scap
      IF(LST.NE.0)THEN
      PTdeg=PT
      ETdeg=ET
      PPldeg=PPl
      EPldeg=EPl
      ENDIF
      I6R001=2
      GOTO09999
09994 CONTINUE
      IF(LST.NE.0)THEN
      PT=PTdeg*WP+PT*(1.D0-WP)+(ETdeg-ET)*Pl**2*(12.D0*WW**2*(1.D0-WW)*WD)
      ET=ETdeg*WP+ET*(1.D0-WP)
      PPl=PPLdeg*WP+PPl*(1.D0-WP)+(Pdeg-P)*(12.D0*WW**2*(1.D0-WW)*WD)+(EPLdeg-EPl)*Pl*(12.D0*WW**2*(1.D0-WW)*WD)+(Edeg-Egas)*(Pl**2*
     *12.D0*WW*(2.D0-3.D0*WW)*WD**2-12.D0*CWW*WW**2*(1.D0-WW)+24.D0*Pl*WW**2*(1.D0-WW)*WD)
      EPl=EPLdeg*WP+EPl*(1.D0-WP)+(Edeg-Egas)*(Pl*12.D0*WW**2*(1.D0-WW)*WD)
      ScapPl=ScapPl*(1.D0-WP)+(Scapdg-Scap)*(Pl*12.D0*WW**2*(1.D0-WW)*WD)
      ScapT=ScapT*(1.D0-WP)
      ZMPl=ZMPl*(1.D0-WP)+(zF-ZMEaN)*(Pl*12.D0*WW**2*(1.D0-WW)*WD)
      ZMT=ZMT*(1.D0-WP)
      YZMPl=YZMPl*(1.D0-WP)+(YZ-YZMEaN)*(Pl*12.D0*WW**2*(1.D0-WW)*WD)
      YZMT=YZMT*(1.D0-WP)
      ENDIF
      P=Pdeg*WP+P*(1.D0-WP)+(Edeg-Egas)*Pl**2*(12.D0*WW**2*(1.D0-WW)*WD)
      Egas=Edeg*WP+Egas*(1.D0-WP)
      Scap=Scapdg*WP+Scap*(1.D0-WP)
      Zmean=Zf*WP+Zmean*(1.D0-WP)
      YZmean=YZ*WP+YZmean*(1.D0-WP)
      Xe=Ye*WP+Xe*(1.D0-WP)
      Psicap=Psi*WP
09931 CONTINUE
      If(RadP)Then
      P=P+(RadC*Tp**4)*tret
      Egas=Egas+(RadC*Tp**4)*RPl
      If(LST.NE.0)Then
      PT=PT+4.D0*(RadC*Tp**3)*tret
      ET=ET+4.D0*(RadC*Tp**3)*RPl
      EPl=EPl-(RadC*Tp**4)*RPl
      Endif
      Endif
      If(LST.NE.0)Then
      DPE=(EPl+Tp*PT*RPl)*(Pl/P)-1.D0
      W=(PT/ET)*RPl
      GAM=(Pl*PPl+W*(Tp*PT))/P
      DA=W/GAM
      Endif
      RETURN
      END
