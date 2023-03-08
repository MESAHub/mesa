! ***********************************************************************
!
!   Copyright (C) 2019  Radek Smolec
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module rsp_lina
      use star_def, only: star_info
      use utils_lib, only: is_bad, mesa_error
      use const_def, only: dp, crad
      use rsp_def
      
      implicit none
      
      private
      public :: mesa_eos_kap, SORT, do_LINA

      
      contains
      

      subroutine do_LINA(s, L0, NZN, NMODES, VEL, PERS, ETO, &
         M, DM, DM_BAR, R, Vol, T, Et, Lr, ierr)
!     LINear Analysis (LINA) USING LAPACK/EISPACK
!     SETS THE MATRIX FOR THE LINEAR EIGENVALUE PROBLEM

      type (star_info), pointer :: s
      real(dp), intent(in) :: L0
      integer, intent(in) :: NMODES, NZN
      real(dp), intent(out) :: VEL(:,:) ! (NZN+1,15)
      real(dp), intent(out), dimension(15) :: PERS, ETO
      real(dp), intent(inout), dimension(:) :: &
         M, DM, DM_BAR, R, Vol, T, Et, Lr
      integer, intent(out) :: ierr
      
      real(dp), dimension(15) :: OMEG, EK
      real(dp), allocatable, dimension(:) :: &
         E, P, Lc, Hp_face, Y_face, K, CPS, QQS, &
         MX10,MX00,MX01,MY00,MY01, &
         MU10,MU00,MU01,MZ00,MZ01, &
         EVUU0,EVUUM, &
         EX20,EX10,EX00,EX01, &
         EY10,EY00,EY01, &
         EU10,EU00,EZ10,EZ00,EZ01, &
         CX10,CX00,CX01,CY00,CY01, &
         CZ00,CX20,CY10,CU00,CU10, &
         CZ01,CZ10, &
         dC_dr_in, dC_dr_00, dC_dr_out, dC_dT_out, &
         dC_dT_00, dC_dT_in,dC_dw_00,dC_dr_in2, &
         DPV,dP_dT_00,DEV,dE_dT_00,dK_dV_00, &
         dK_dT_00,DVR,DVRM, &
         ELRP,ELRM,ELR,ELTP,ELT, &
         ELZ0,ELZP, &
         CPV,dCp_dT_00,QQV,dQQ_dT_00 , &
         dP_dr_00, dCp_dr_00, dQQ_dr_00, &
         dP_dr_in,dCp_dr_in,dQQ_dr_in, &
         dHp_dr_in,dHp_dr_out,dHp_dr_00,dHp_dT_00,dHp_dT_out, &
         dY_dr_in,dY_dr_00,dY_dr_out,dY_dT_00, &
         dY_dT_out, &
         PII,dPII_dr_in,dPII_dr_00,dPII_dr_out, &
         dPII_dT_00,dPII_dT_out,DPIIZ0, &
         dsrc_dr_in,dsrc_dr_00,dsrc_dr_out,dsrc_dT_00, &
         dsrc_dT_out,dsrc_dw_00,source, &
         d_damp_dr_in,d_damp_dr_00,d_damp_dr_out,d_damp_dT_00, &
         d_damp_dT_out,d_damp_dw_00,DAMP, &
         dsrc_dr_in2,dsrc_dT_in,d_damp_dr_in2,d_damp_dT_in, &
         d_dampR_dr_in,d_dampR_dr_00,d_dampR_dr_out,d_dampR_dT_00, &
         d_dampR_dT_out,d_dampR_dw_00,DAMPR, &
         d_dampR_dr_in2,d_dampR_dT_in, &
         DLCXM,DLCX0,DLCXP,DLCY0, &
         DLCYP,DLCZ0,DLCZP, &
         DLTXM,DLTX0,DLTXP,DLTY0, &
         DLTYP,DLTZ0,DLTZP,Lt, &
         PTURB,dPtrb_dw_00,dPtrb_dr_00,dPtrb_dr_in,XXL
      real(dp), allocatable, dimension(:,:) :: &
         QWK,QWKEV,PHR,PHT,PHL,PHU,PHC,ZZZ
      complex(8), allocatable, dimension(:,:) :: &
         VRR,VRT,VRL,VRC,DVRR,DVRT,DVRL,DVRC,VRU
      
      real(dp) :: FFXM,FFX0,FFXP,FFY0,FFYP,FF,POM1
      real(dp) :: IGR1,IGR1XM,IGR1X0,IGR1XP,IGR1Y0,IGR1YP
      real(dp) :: IGR2,IGR2XM,IGR2X0,IGR2XP,IGR2Y0,IGR2YP
      real(dp) :: GGXM,GGX0,GGXP,GGY0,GGYP,GG,GPF
      real(dp) :: DFCX0,DFCXM,DFCXP,DFCY0,DFCYP,DFCZ0,DFCZP
      real(dp) :: DFCMX0,DFCMXM,DFCMXP,DFCMY0,DFCMYP,DFCMZ0,DFCMZP
      real(dp) :: FLIM,FLD
      real(dp) :: T1,DLR,DLRP,DLRM,DLT,DLTP,DLMR,DLMRP,DLMRM,DLMT,DLMTP
      real(dp) :: W_00,W_out,BW,BK,T2,T3,DLK,DLKP,MINI
      real(dp) :: T4,T5,P2,P3,P5,P6,P7,P8,P9,P10,POM3,POM2,POM,POM4
      integer :: I,J,NZN3,IG,IE,IR,IC,INFO,IMI,LD_VL,LD_VR,n,op_err
      real(dp) :: RES(12),VRRS(15),Q(15)
      character (len=250) FILENAME
      character (len=1) NUMER1
      character (len=2) NUMER2
      complex(8):: DP_0,DV_0,VTTS(15),SCALE(15),CPOM,DPEV,dP_dT_00URB
      real(dp) :: SGRP,SGRM
      real(dp) :: PSIG,TEMI,TEMM,TEM1
      real(dp) :: NORMC
      real(dp) :: QCHECK(15),ETOIEV(15),QWKPT(1000,15),ETOIPT(15),SGR(15)
      real(dp) :: TT4,TT4P,EFL02   
      real(dp) :: ETOI(15)
      
      !write(*,'(a55,i12,99(1pd26.16))') 'start LINA s% w(2)**2', &
      !   2, s% w(2)**2, Et(NZN-1)
      
      if (.false.) then
      write(*,*) 'L0', L0
      do I=1,NZN
         write(*,*) 'M', I, M(I)
         write(*,*) 'DM', I, DM(I)
         write(*,*) 'DM_BAR', I, DM_BAR(I)
         write(*,*) 'R', I, R(I)
         write(*,*) 'R**3', I, R(I)**3
         write(*,*) 'Vol', I, Vol(I)
         write(*,*) 'T', I, T(I)
         write(*,*) 'Lr', I, Lr(I)
         write(*,*) 'Et', I, Et(I)
      end do
      call mesa_error(__FILE__,__LINE__,'do_LINA')
      end if
      
      ierr = 0
      EFL02 = EFL0*EFL0
      n = NZN+1
      allocate( &
         E(n), P(n), Lc(n), Hp_face(n), Y_face(n), K(n), CPS(n), QQS(n), &
         MX10(n),MX00(n),MX01(n),MY00(n),MY01(n), &
         MU10(n),MU00(n),MU01(n),MZ00(n),MZ01(n), &
         EVUU0(n),EVUUM(n), &
         EX20(n),EX10(n),EX00(n),EX01(n), &
         EY10(n),EY00(n),EY01(n), &
         EU10(n),EU00(n),EZ10(n),EZ00(n),EZ01(n), &
         CX10(n),CX00(n),CX01(n),CY00(n),CY01(n), &
         CZ00(n),CX20(n),CY10(n),CU00(n),CU10(n), &
         CZ01(n),CZ10(n), &
         dC_dr_in(n), dC_dr_00(n), dC_dr_out(n), dC_dT_out(n), &
         dC_dT_00(n), dC_dT_in(n),dC_dw_00(n),dC_dr_in2(n), &
         DPV(n),dP_dT_00(n),DEV(n),dE_dT_00(n),dK_dV_00(n), &
         dK_dT_00(n),DVR(n),DVRM(n), &
         ELRP(n),ELRM(n),ELR(n),ELTP(n),ELT(n), &
         ELZ0(n),ELZP(n), &
         CPV(n),dCp_dT_00(n),QQV(n),dQQ_dT_00(n) , &
         dP_dr_00(n), dCp_dr_00(n), dQQ_dr_00(n), &
         dP_dr_in(n),dCp_dr_in(n),dQQ_dr_in(n), &
         dHp_dr_in(n),dHp_dr_out(n),dHp_dr_00(n),dHp_dT_00(n),dHp_dT_out(n), &
         dY_dr_in(n),dY_dr_00(n),dY_dr_out(n),dY_dT_00(n), &
         dY_dT_out(n), &
         PII(n),dPII_dr_in(n),dPII_dr_00(n),dPII_dr_out(n), &
         dPII_dT_00(n),dPII_dT_out(n),DPIIZ0(n), &
         dsrc_dr_in(n),dsrc_dr_00(n),dsrc_dr_out(n),dsrc_dT_00(n), &
         dsrc_dT_out(n),dsrc_dw_00(n),source(n), &
         d_damp_dr_in(n),d_damp_dr_00(n),d_damp_dr_out(n),d_damp_dT_00(n), &
         d_damp_dT_out(n),d_damp_dw_00(n),DAMP(n), &
         dsrc_dr_in2(n),dsrc_dT_in(n),d_damp_dr_in2(n),d_damp_dT_in(n), &
         d_dampR_dr_in(n),d_dampR_dr_00(n),d_dampR_dr_out(n),d_dampR_dT_00(n), &
         d_dampR_dT_out(n),d_dampR_dw_00(n),DAMPR(n), &
         d_dampR_dr_in2(n),d_dampR_dT_in(n), &
         DLCXM(n),DLCX0(n),DLCXP(n),DLCY0(n), &
         DLCYP(n),DLCZ0(n),DLCZP(n), &
         DLTXM(n),DLTX0(n),DLTXP(n),DLTY0(n), &
         DLTYP(n),DLTZ0(n),DLTZP(n),Lt(n), &
         PTURB(n),dPtrb_dw_00(n),dPtrb_dr_00(n),dPtrb_dr_in(n), &
         QWK(n,15),QWKEV(n,15),VRU(n,15), &
         VRR(n,15), VRT(n,15), VRL(n,15), VRC(n,15), &
         DVRR(n,15),DVRT(n,15),DVRL(n,15),DVRC(n,15), &
         PHR(n,15), PHT(n,15), PHL(n,15), PHU(n,15), PHC(n,15), &
         ZZZ(4*n,4*n),XXL(4*n))

      NZN3=4*NZN
      LD_VL = LD_LLL
      LD_VR = LD_LLL

!     SET ALL NECESSARY DERIVATIVES AND VALUES EQUAL 0
!     FOR FROZEN-IN APPROXIMATION, OR IN CASE ALFA=0 (RADIATIVE)
      if(ALFA.eq.0.d0)then
         do I=1,NZN
            Et(I)=0.d0
            EFL02=0.d0
            DLCX0(I)=0.d0
            DLCXM(I)=0.d0
            DLCXP(I)=0.d0
            DLCY0(I)=0.d0
            DLCYP(I)=0.d0
            DLCZ0(I)=0.d0
            DLCZP(I)=0.d0
            DLTX0(I)=0.d0
            DLTXM(I)=0.d0
            DLTXP(I)=0.d0
            DLTY0(I)=0.d0
            DLTYP(I)=0.d0
            DLTZ0(I)=0.d0
            DLTZP(I)=0.d0
            EVUU0(I)= 0.d0
            EVUUM(I)= 0.d0
            PTURB(I)= 0.d0
            dPtrb_dr_00(I)= 0.d0
            dPtrb_dr_in(I)= 0.d0
            dPtrb_dw_00(I)= 0.d0
            dC_dr_00(I) = 0.d0
            dC_dr_out(I) = 0.d0
            dC_dr_in(I) = 0.d0
            dC_dr_in2(I)= 0.d0
            dC_dT_in(I) = 0.d0
            dC_dT_00(I) = 0.d0
            dC_dT_out(I) = 0.d0
            dC_dw_00(I) = 0.d0
         enddo
      endif

!     LOOP 1 .. EOS
 
      !$OMP PARALLEL DO PRIVATE(I,T1,op_err) SCHEDULE(dynamic,2)
      do 1 I=1,NZN

         if(Et(I).le.EFL02) Et(I)=EFL02

         call mesa_eos_kap(s,0, &
              T(I),Vol(I),P(I),DPV(I),dP_dT_00(I), &
              E(I),DEV(I),dE_dT_00(I),CPS(I),CPV(I),dCp_dT_00(I), &
              QQS(I),QQV(I),dQQ_dT_00(I),K(I),dK_dV_00(I),dK_dT_00(I),op_err)
         if (op_err /= 0) ierr = op_err
         if (ierr /= 0) cycle

         T1=P43/dm(I)
         DVR(I)=3.d0*T1*R(I)**2    
         if(I.eq.1) goto 2        
         DVRM(I)=-3.d0*T1*R(max(1,I-1))**2 
             ! bp: max(1,i-1) to prevent bogus warning from gfortran
         goto 3
 2       DVRM(I)=-3.d0*T1*s% R_center**2
 3       continue
         dP_dr_00(I) =DPV(I)*DVR(I) 
         dP_dr_in(I)=DPV(I)*DVRM(I)
         dCp_dr_00(I) =CPV(I)*DVR(I)
         dCp_dr_in(I)=CPV(I)*DVRM(I)
         dQQ_dr_00(I) =QQV(I)*DVR(I)
         dQQ_dr_in(I)=QQV(I)*DVRM(I)
         
 1    continue
      !$OMP END PARALLEL DO
      if (ierr /= 0) return

!     SKIP ALL DERIVATIVE CALCULATIONS IN CASE OF FROZEN-IN
!     APPROXIMATION OR ALFA=0 (RADIATIVE CASE)
      if(ALFA.eq.0.d0) goto 999

!     SET E_T (NOW =w) BELOW AND ABOVE BOUNDARIES
      Et(NZN) = 0.d0
      do I=1,IBOTOM
         Et(I) = 0.d0
      enddo

      do I=1,NZN-1
         POM=  (R(I)**2)/(2.d0*G*M(I))
         POM2=POM*(P(I)*Vol(I)+P(I+1)*Vol(I+1))
         Hp_face(I)=POM2
         dHp_dr_in(I)=POM*(P(I)*DVRM(I)+Vol(I)*dP_dr_in(I))
         dHp_dr_00(I)=2.d0*Hp_face(I)/R(I)+POM*(P(I)*DVR(I)+Vol(I)*dP_dr_00(I) &
                   +P(I+1)*DVRM(I+1)+Vol(I+1)*dP_dr_in(I+1)) 
         dHp_dr_out(I)=POM*(P(I+1)*DVR(I+1)+Vol(I+1)*dP_dr_00(I+1))
         dHp_dT_00(I)=POM*Vol(I)*dP_dT_00(I)
         dHp_dT_out(I)=POM*Vol(I+1)*dP_dT_00(I+1)
      enddo
      POM=(R(NZN)**2)/(2.d0*G*M(NZN))
      Hp_face(NZN)=POM*P(NZN)*Vol(NZN)
      dHp_dr_in(NZN)=POM*(P(NZN)*DVRM(NZN)+Vol(NZN)*dP_dr_in(NZN))     
      dHp_dr_00(NZN)=2.d0*Hp_face(NZN)/R(NZN)+POM* &
      (P(NZN)*DVR(NZN)+Vol(NZN)*dP_dr_00(NZN))
      dHp_dr_out(NZN)=0.d0
      dHp_dT_00(NZN)=POM*Vol(NZN)*dP_dT_00(NZN)
      dHp_dT_out(NZN)=0.d0
      
      !call mesa_error(__FILE__,__LINE__,'Hp_face')

      do I=1,NZN-1
         POM=0.5d0*(QQS(I)/CPS(I)+QQS(I+1)/CPS(I+1))
         POM2=0.5d0*(P(I+1)-P(I))
         IGR1=0.5d0*(QQS(I)/CPS(I)+QQS(I+1)/CPS(I+1))*(P(I+1)-P(I)) &
              -(dlog(T(I+1))-dlog(T(I)))
         IGR1X0=POM2*((dQQ_dr_00(I)-QQS(I)/CPS(I)*dCp_dr_00(I))/CPS(I)+ &
                 (dQQ_dr_in(I+1)-QQS(I+1)/CPS(I+1)*dCp_dr_in(I+1))/CPS(I+1)) &
               +POM*(dP_dr_in(I+1)-dP_dr_00(I))
         IGR1XM=POM2*(dQQ_dr_in(I)-QQS(I)/CPS(I)*dCp_dr_in(I))/CPS(I) &
                +POM*(-dP_dr_in(I))
         IGR1XP=POM2*(dQQ_dr_00(I+1)-QQS(I+1)/CPS(I+1)*dCp_dr_00(I+1))/CPS(I+1) &
                +POM*(dP_dr_00(I+1))
         IGR1Y0=POM2*(dQQ_dT_00(I)-QQS(I)/CPS(I)*dCp_dT_00(I))/CPS(I) &
                +POM*(-dP_dT_00(I))+1.d0/T(I)
         IGR1YP=POM2*(dQQ_dT_00(I+1)-QQS(I+1)/CPS(I+1)*dCp_dT_00(I+1))/CPS(I+1) &
                +POM*(dP_dT_00(I+1))-1.d0/T(I+1)

         POM=2.d0/(Vol(I)+Vol(I+1))
         POM2=8.d0*PI*(R(I)**2)/dm_bar(I)*Hp_face(I)
         IGR2=4.d0*PI*(R(I)**2)*Hp_face(I)*POM/dm_bar(I)
         IGR2X0=2.d0*IGR2/R(I)+IGR2/Hp_face(I)*dHp_dr_00(I) &
               -POM2/(Vol(I)+Vol(I+1))**2*(DVR(I)+DVRM(I+1))
         IGR2XM=-POM2/(Vol(I)+Vol(I+1))**2*DVRM(I)    &
               +IGR2/Hp_face(I)*dHp_dr_in(I)
         IGR2XP=-POM2/(Vol(I)+Vol(I+1))**2*DVR(I+1) &
               +IGR2/Hp_face(I)*dHp_dr_out(I)
         IGR2Y0=IGR2/Hp_face(I)*dHp_dT_00(I)
         IGR2YP=IGR2/Hp_face(I)*dHp_dT_out(I)

         Y_face(I)=IGR1*IGR2
         dY_dr_00(I)=IGR1*IGR2X0+IGR2*IGR1X0
         dY_dr_in(I)=IGR1*IGR2XM+IGR2*IGR1XM
         dY_dr_out(I)=IGR1*IGR2XP+IGR2*IGR1XP
         dY_dT_00(I)=IGR1*IGR2Y0+IGR2*IGR1Y0
         dY_dT_out(I)=IGR1*IGR2YP+IGR2*IGR1YP
      enddo
      
      !call mesa_error(__FILE__,__LINE__,'IGRS')

      do I=1,NZN-1
         POM=sqrt(2.d0/3.d0)*0.5d0
         FF=POM*((E(I)  +P(I)  *Vol(I)  )/T(I) &
              +(E(I+1)+P(I+1)*Vol(I+1))/T(I+1))
         FFY0=POM*(dE_dT_00(I)+dP_dT_00(I)*Vol(I)-(E(I)+P(I)*Vol(I))/T(I))/T(I)
         FFYP=POM*(dE_dT_00(I+1)+dP_dT_00(I+1)*Vol(I+1)-(E(I+1)+P(I+1)*Vol(I+1)) &
            /T(I+1))/T(I+1)
         FFXP=POM* ((DEV(I+1)+P(I+1))*DVR(I+1)+Vol(I+1)*dP_dr_00(I+1))/T(I+1)
         FFX0=POM*(((DEV(I)  +P(I)  )*DVR(I)  +Vol(I)  *dP_dr_00(I)  )/T(I)+ &
              ((DEV(I+1)+P(I+1))*DVRM(I+1)+Vol(I+1)*dP_dr_in(I+1)) &
              /T(I+1))
         FFXM=POM* ((DEV(I)  +P(I)  )*DVRM(I) +Vol(I)*dP_dr_in(I))/T(I)
         
         POM=ALFAS*ALFA
         POM2=0.5d0*(CPS(I)+CPS(I+1))
         GG=POM*POM2*Y_face(I)
         GGXM=POM*(POM2*dY_dr_in(I)+Y_face(I)*0.5d0*dCp_dr_in(I))
         GGX0=POM*(POM2*dY_dr_00(I)+Y_face(I)*0.5d0*(dCp_dr_00(I)+dCp_dr_in(I+1)))
         GGXP=POM*(POM2*dY_dr_out(I)+Y_face(I)*0.5d0*dCp_dr_00(I+1))
         GGY0=POM*(POM2*dY_dT_00(I)+Y_face(I)*0.5d0*dCp_dT_00(I))
         GGYP=POM*(POM2*dY_dT_out(I)+Y_face(I)*0.5d0*dCp_dT_00(I+1))
         GPF=GG/FF

!        corelation PI defined without e_t
         POM=1.d0

         PII(I)=POM*GG
         DPIIZ0(I)=0.d0
         dPII_dr_00(I)=POM*GGX0
         dPII_dr_in(I)=POM*GGXM 
         dPII_dr_out(I)=POM*GGXP
         dPII_dT_00(I)=POM*GGY0
         dPII_dT_out(I)=POM*GGYP
      enddo
 
      !call mesa_error(__FILE__,__LINE__,'LINA')

      do I=IBOTOM+1,NZN-1
!        SOURCE TERM
         POM=0.5d0*(PII(I)/Hp_face(I)+PII(I-1)/Hp_face(I-1))
         POM2=T(I)*P(I)*QQS(I)/CPS(I)
         POM3=sqrt(Et(I))
         SOURCE(I)=POM*POM2*POM3
         TEM1=POM2*POM3*0.5d0
         TEMI=-PII(I)/Hp_face(I)**2
         TEMM=-PII(I-1)/Hp_face(I-1)**2

         dsrc_dr_out(I)= TEM1*(dPII_dr_out(I)/Hp_face(I) &
                         +TEMI*dHp_dr_out(I))
         dsrc_dr_00(I)= TEM1*(dPII_dr_00(I)/Hp_face(I)+dPII_dr_out(I-1)/Hp_face(I-1) &
                         +TEMI*dHp_dr_00(I)+TEMM*dHp_dr_out(I-1))
         dsrc_dr_in(I)= TEM1*(dPII_dr_in(I)/Hp_face(I)+dPII_dr_00(I-1)/Hp_face(I-1) &
                         +TEMI*dHp_dr_in(I)+TEMM*dHp_dr_00(I-1))
         dsrc_dr_in2(I)=TEM1*(dPII_dr_in(I-1)/Hp_face(I-1) &
                         +TEMM*dHp_dr_in(I-1))

         dsrc_dT_out(I)= TEM1*(dPII_dT_out(I)/Hp_face(I) &
                         +TEMI*dHp_dT_out(I))
         dsrc_dT_00(I)= TEM1*(dPII_dT_00(I)/Hp_face(I)+dPII_dT_out(I-1)/Hp_face(I-1) &
                         +TEMI*dHp_dT_00(I)+TEMM*dHp_dT_out(I-1))
         dsrc_dT_in(I)= TEM1*(dPII_dT_00(I-1)/Hp_face(I-1) &
                         +TEMM*dHp_dT_00(I-1))

         dsrc_dw_00(I)= POM*POM2*0.5d0/sqrt(Et(I))

         POM=POM*POM3

         dsrc_dT_00(I)=dsrc_dT_00(I)+ &
            POM/CPS(I)*(P(I)*QQS(I)+T(I)*QQS(I)*dP_dT_00(I)+T(I)*P(I)*dQQ_dT_00(I) &
                            -T(I)*P(I)*QQS(I)/CPS(I)*dCp_dT_00(I))
         dsrc_dr_00(I)=dsrc_dr_00(I)+ &
               POM*T(I)/CPS(I)*(QQS(I)*dP_dr_00(I)+P(I)*dQQ_dr_00(I) &
                           -P(I)*QQS(I)/CPS(I)*dCp_dr_00(I))
         dsrc_dr_in(I)=dsrc_dr_in(I)+ &
                POM*T(I)/CPS(I)*(QQS(I)*dP_dr_in(I)+P(I)*dQQ_dr_in(I) &
                            -P(I)*QQS(I)/CPS(I)*dCp_dr_in(I))

!       SOURCE TERM ALWAYS POSITIVE
!       THIS IS FORMULATION USED IN BUDAPEST-FLORIDA CODE
!       IT REQUIRES E_T AS MAIN VARIABLE - OTHERWISE CONVERGENCE
!       IS EXTREMELY SLOW, AND POSSIBLE ONLY IF LOW ACCURACY
!       FOR CONVERGENCE CONDITION, DXXC.LT. 1d-2 - 1.d-4, IS SET
!        if(SOURCE(I).le.0.d0)then
!           SOURCE(I) = 0.d0
!           dsrc_dr_out(I) = 0.d0
!           dsrc_dr_00(I) = 0.d0
!           dsrc_dr_in(I) = 0.d0
!           dsrc_dr_in2(I)= 0.d0
!           dsrc_dT_out(I) = 0.d0
!           dsrc_dT_in(I) = 0.d0
!           dsrc_dT_00(I) = 0.d0
!           dsrc_dw_00(I) = 0.d0
!        endif

!        DAMP TERM
         POM=(CEDE/ALFA)*(Et(I)**1.5d0-EFL02**1.5d0)
         POM2=0.5d0*(Hp_face(I)+Hp_face(I-1))
         DAMP(I)=POM/POM2
         TEM1=-0.5d0*POM/POM2**2
         d_damp_dr_out(I) =TEM1*dHp_dr_out(I)
         d_damp_dr_00(I) =TEM1*(dHp_dr_00(I)+dHp_dr_out(I-1))
         d_damp_dr_in(I) =TEM1*(dHp_dr_in(I)+dHp_dr_00(I-1))
         d_damp_dr_in2(I)=TEM1*(dHp_dr_in(I-1))
         d_damp_dT_00(I) =TEM1*(dHp_dT_00(I)+dHp_dT_out(I-1))
         d_damp_dT_out(I) =TEM1*dHp_dT_out(I)
         d_damp_dT_in(I) =TEM1*dHp_dT_00(I-1)
         d_damp_dw_00(I) =1.5d0*(CEDE/ALFA)/POM2*sqrt(Et(I))

!        RADIATIVE DAMP TERM
         if(GAMMAR.eq.0.d0)then
            DAMPR(I)   = 0.d0
            d_dampR_dr_out(I) = 0.d0
            d_dampR_dr_00(I) = 0.d0
            d_dampR_dr_in(I) = 0.d0
            d_dampR_dr_in2(I)= 0.d0
            d_dampR_dT_out(I) = 0.d0
            d_dampR_dT_in(I) = 0.d0
            d_dampR_dT_00(I) = 0.d0
            d_dampR_dw_00(I) = 0.d0
         else
            POM=(GAMMAR**2)/(ALFA**2)*4.d0*SIG
            POM2=T(I)**3*Vol(I)**2/(CPS(I)*K(I))
            POM3=Et(I)
            POM4=0.5d0*(Hp_face(I)**2+Hp_face(I-1)**2)
            DAMPR(I)=POM*POM2*POM3/POM4
            TEM1=-DAMPR(I)/POM4
            d_dampR_dr_out(I) =TEM1*(Hp_face(I)*dHp_dr_out(I))
            d_dampR_dr_00(I) =TEM1*(Hp_face(I)*dHp_dr_00(I) &
                             +Hp_face(I-1)*dHp_dr_out(I-1))
            d_dampR_dr_in(I) =TEM1*(Hp_face(I)*dHp_dr_in(I) &
                             +Hp_face(I-1)*dHp_dr_00(I-1))
            d_dampR_dr_in2(I)=TEM1*(Hp_face(I-1)*dHp_dr_in(I-1))
            d_dampR_dT_out(I) =TEM1*(Hp_face(I)*dHp_dT_out(I))
            d_dampR_dT_00(I) =TEM1*(Hp_face(I)*dHp_dT_00(I) &
                             +Hp_face(I-1)*dHp_dT_out(I-1))
            d_dampR_dT_in(I) =TEM1*(Hp_face(I-1)*dHp_dT_00(I-1))

            d_dampR_dw_00(I)=POM*POM2/POM4

            TEM1=POM*POM3/POM4
            d_dampR_dr_00(I)=d_dampR_dr_00(I) &
                     +TEM1*T(I)**3*(2.d0*Vol(I)*DVR(I) &
                    -Vol(I)**2*( 1.d0/CPS(I)*dCp_dr_00(I) &
                                 +1.d0/K(I)*dK_dV_00(I)*DVR(I))) &
                      /(CPS(I)*K(I))              

            d_dampR_dr_in(I)=d_dampR_dr_in(I) &
                    +TEM1*T(I)**3*(2.d0*Vol(I)*DVRM(I) &
                    -Vol(I)**2*( 1.d0/CPS(I)*dCp_dr_in(I) &
                                 +1.d0/K(I)*dK_dV_00(I)*DVRM(I))) &
                      /(CPS(I)*K(I))

            d_dampR_dT_00(I)=d_dampR_dT_00(I) &
                          +TEM1*Vol(I)**2*(3.d0*T(I)**2 &
                          -T(I)**3*( 1.d0/CPS(I)*dCp_dT_00(I) &
                                       +1.d0/K(I)*dK_dT_00(I))) &
                         /(CPS(I)*K(I))

         endif
  
         dC_dr_00(I) =dsrc_dr_00(I) -d_damp_dr_00(I) -d_dampR_dr_00(I)
         dC_dr_out(I) =dsrc_dr_out(I) -d_damp_dr_out(I) -d_dampR_dr_out(I)
         dC_dr_in(I) =dsrc_dr_in(I) -d_damp_dr_in(I) -d_dampR_dr_in(I)
         dC_dr_in2(I)=dsrc_dr_in2(I)-d_damp_dr_in2(I)-d_dampR_dr_in2(I)
         dC_dT_in(I) =dsrc_dT_in(I) -d_damp_dT_in(I) -d_dampR_dT_in(I)
         dC_dT_00(I) =dsrc_dT_00(I) -d_damp_dT_00(I) -d_dampR_dT_00(I)
         dC_dT_out(I) =dsrc_dT_out(I) -d_damp_dT_out(I) -d_dampR_dT_out(I)
         dC_dw_00(I) =dsrc_dw_00(I) -d_damp_dw_00(I) -d_dampR_dw_00(I)
         
      enddo
 
       !call mesa_error(__FILE__,__LINE__,'LINA')

      do I=IBOTOM,NZN-1
!        CONVECTIVE LUMINOSITY
         POM=4.d0*PI*(R(I)**2)*(T(I)/Vol(I)+T(I+1)/Vol(I+1))*0.5d0* &
             (ALFAC/ALFAS)
         POM3=(sqrt(Et(I))+sqrt(Et(I+1)))*0.5d0
         Lc(I)=POM*PII(I)*POM3
         DLCZ0(I)=0d0
         DLCZP(I)=0d0
         if (Et(I) /= 0d0) DLCZ0(I)=POM*PII(I)*0.25d0/sqrt(Et(I))
         if (Et(I+1) /= 0d0) DLCZP(I)=POM*PII(I)*0.25d0/sqrt(Et(I+1))
         POM2=4.d0*PI*(R(I)**2)*PII(I)*0.5d0*(ALFAC/ALFAS)*POM3
         POM=POM*POM3
         DLCX0(I)=POM*dPII_dr_00(I)+2.d0*Lc(I)/R(I) &
                            -POM2*(T(I)  /(Vol(I)**2  )*DVR(I)+ &
                                   T(I+1)/(Vol(I+1)**2)*DVRM(I+1))
         DLCXM(I)=POM*dPII_dr_in(I)-POM2*T(I)  /(Vol(I)**2  )*DVRM(I)
         DLCXP(I)=POM*dPII_dr_out(I)-POM2*T(I+1)/(Vol(I+1)**2)*DVR(I+1)
         DLCY0(I)=POM*dPII_dT_00(I)+POM2/Vol(I)
         DLCYP(I)=POM*dPII_dT_out(I)+POM2/Vol(I+1)

         if(I.eq.IBOTOM) DLCZ0(I)=0.d0
         if(I.eq.NZN-1)  DLCZP(I)=0.d0
         if(Et(I).lt.EFL02*1d-10)then
            Lc(I)=0.d0
            DLCX0(I)=0.d0
            DLCXM(I)=0.d0
            DLCXP(I)=0.d0
            DLCY0(I)=0.d0
            DLCYP(I)=0.d0
            DLCZ0(I)=0.d0
            DLCZP(I)=0.d0
         endif

!         if(PII(I).lt.0.d0.or.ALFA.eq.0.d0)then
!            Lc(I)=0.d0
!            DLCX0(I)=0.d0
!            DLCXM(I)=0.d0
!            DLCXP(I)=0.d0
!            DLCY0(I)=0.d0
!            DLCYP(I)=0.d0
!            DLCZ0(I)=0.d0
!            DLCZP(I)=0.d0
!         endif

!        TURBULENT LUMINOSITY
         if(ALFAT.eq.0.d0.or.ALFA.eq.0.d0)then
            Lt(I)=0.d0
            DLTX0(I)=0.d0
            DLTXM(I)=0.d0
            DLTXP(I)=0.d0
            DLTY0(I)=0.d0
            DLTYP(I)=0.d0
            DLTZ0(I)=0.d0
            DLTZP(I)=0.d0
         else
            POM=-2.d0/3.d0*ALFA*ALFAT*(4.d0*PI*(R(I)**2))**2
            POM2=Hp_face(I)*(1.d0/Vol(I)**2+1.d0/Vol(I+1)**2)*0.5d0
            POM3=(Et(I+1)**1.5d0-Et(I)**1.5d0)/dm_bar(I)
            Lt(I)=POM*POM2*POM3
            DLTX0(I)=4.d0*Lt(I)/R(I) &
                 +Lt(I)/POM2*Hp_face(I)*(-1.d0/Vol(I)**3*DVR(I) &
                                        -1.d0/Vol(I+1)**3*DVRM(I+1))  &
                 +Lt(I)/Hp_face(I)*dHp_dr_00(I)

            DLTXM(I)=Lt(I)/POM2*Hp_face(I)*(-1.d0/Vol(I)**3*DVRM(I)) &
                    +Lt(I)/Hp_face(I)*dHp_dr_in(I)
            DLTXP(I)=Lt(I)/POM2*Hp_face(I)*(-1.d0/Vol(I+1)**3*DVR(I+1)) &
                    +Lt(I)/Hp_face(I)*dHp_dr_out(I)
            DLTY0(I)=Lt(I)/Hp_face(I)*dHp_dT_00(I)
            DLTYP(I)=Lt(I)/Hp_face(I)*dHp_dT_out(I)
            DLTZ0(I)=-POM*POM2*1.5d0*sqrt(Et(I)  )/dm_bar(I) 
            DLTZP(I)= POM*POM2*1.5d0*sqrt(Et(I+1))/dm_bar(I) 
         endif
      enddo

!     TURBULENT PRESSURE (ZONE)
      do I=IBOTOM+1,NZN-1
         if(ALFAP.eq.0.d0)then
            PTURB(I) = 0.d0
            dPtrb_dw_00(I) = 0.d0
            dPtrb_dr_00(I) = 0.d0
            dPtrb_dr_in(I) = 0.d0
         else
            PTURB(I) =  ALFAP*Et(I)/Vol(I)
            dPtrb_dw_00(I) =  ALFAP/Vol(I)
            TEM1=-ALFAP*Et(I)/Vol(I)**2
            dPtrb_dr_00(I) = TEM1*DVR(I)
            dPtrb_dr_in(I) = TEM1*DVRM(I)
         endif
      enddo

!     EDDY VISCOSITY
      do I=IBOTOM+1,NZN-1
         if(ALFAM.ge.0d0) then
!           Kuhfuss (1986) tensor EDDY VISCOSITY
            POM=16.d0/3.d0*PI*ALFA*ALFAM*sqrt(Et(I))
            POM1=1.d0/Vol(I)**2/dm(I)
            POM2=(R(I)**6+R(I-1)**6)*(Hp_face(I)+Hp_face(I-1))*0.25d0
            EVUU0(I)= POM*POM1*POM2/R(I)
            EVUUM(I)=-POM*POM1*POM2/R(I-1)
         else
!           Kollath et al. 2002  EDDY VISCOSITY pressure 
            POM=-(16.d0/3.d0)*PI*ALFA*abs(ALFAM)*sqrt(Et(I))
            POM1=1.d0/Vol(I)**2/dm(I)
            POM2=(R(I)**3+R(I-1)**3)*(Hp_face(I)+Hp_face(I-1))*0.25d0
            EVUU0(I)= POM*POM1*POM2/R(I)
            EVUUM(I)=-POM*POM1*POM2/R(I-1) 
          endif
      enddo

      do I=1,IBOTOM
         Lc(I)= 0.d0
         Lt(I)= 0.d0
         PTURB(I)= 0.d0
         dPtrb_dw_00(I)= 0.d0
         dPtrb_dr_00(I)= 0.d0
         dPtrb_dr_in(I)= 0.d0
         dC_dr_00(I) = 0.d0
         dC_dr_out(I) = 0.d0
         dC_dr_in(I) = 0.d0
         dC_dr_in2(I)= 0.d0
         dC_dT_in(I) = 0.d0
         dC_dT_00(I) = 0.d0
         dC_dT_out(I) = 0.d0
         dC_dw_00(I) = 0.d0
         EVUU0(I) = 0.d0
         EVUUM(I) = 0.d0
      enddo
      do I=1,IBOTOM-1
         DLCX0(I)=0.d0
         DLCXM(I)=0.d0
         DLCXP(I)=0.d0
         DLCY0(I)=0.d0
         DLCYP(I)=0.d0
         DLCZ0(I)=0.d0
         DLCZP(I)=0.d0
         DLTX0(I)=0.d0
         DLTXM(I)=0.d0
         DLTXP(I)=0.d0
         DLTY0(I)=0.d0
         DLTYP(I)=0.d0
         DLTZ0(I)=0.d0
         DLTZP(I)=0.d0
      enddo
      Lc(NZN)= 0.d0
      Lt(NZN)= 0.d0
      PTURB(NZN)= 0.d0
      dPtrb_dw_00(NZN)= 0.d0
      dPtrb_dr_00(NZN)= 0.d0
      dPtrb_dr_in(NZN)= 0.d0
      DLCX0(NZN)= 0.d0
      DLCXM(NZN)= 0.d0
      DLCXP(NZN)= 0.d0
      DLCY0(NZN)= 0.d0
      DLCYP(NZN)= 0.d0
      DLCZ0(NZN)= 0.d0
      DLCZP(NZN)= 0.d0
      DLTX0(NZN)=0.d0
      DLTXM(NZN)=0.d0
      DLTXP(NZN)=0.d0
      DLTY0(NZN)=0.d0
      DLTYP(NZN)=0.d0
      DLTZ0(NZN)=0.d0
      DLTZP(NZN)=0.d0
      dC_dr_00(NZN) = 0.d0
      dC_dr_out(NZN) = 0.d0
      dC_dr_in(NZN) = 0.d0
      dC_dr_in2(NZN)= 0.d0
      dC_dT_in(NZN) = 0.d0
      dC_dT_00(NZN) = 0.d0
      dC_dT_out(NZN) = 0.d0
      dC_dw_00(NZN) = 0.d0 
      EVUU0(NZN)= 0.d0
      EVUUM(NZN)= 0.d0

      DLCZP(NZN-1)=0.d0
      DLTZP(NZN-1)=0.d0

 999  continue

!     LOOP 2 .. LUM PLUGS
      DLR  =  0.d0
      DLRP =  0.d0
      DLRM =  0.d0
      DLT  =  0.d0
      DLTP =  0.d0
      DLR  =  0.d0!-1.d0

      DFCX0 = 0.d0
      DFCXM = 0.d0
      DFCXP = 0.d0
      DFCY0 = 0.d0
      DFCYP = 0.d0
      DFCZ0 = 0.d0
      DFCZP = 0.d0
      do 5 I=1,NZN
!        SET LUM(I-1)
         DLMR  = DLR
         DLMRP = DLRP
         DLMRM = DLRM
         DLMT  = DLT
         DLMTP = DLTP

         DFCMX0 = DFCX0
         DFCMXM = DFCXM
         DFCMXP = DFCXP
         DFCMY0 = DFCY0
         DFCMYP = DFCYP
         DFCMZ0 = DFCZ0
         DFCMZP = DFCZP
         if(I.eq.NZN) goto 6
!        Lr(I)=Eq. A.4, Stellingwerf 1975, Appendix A
!        CALC LUM(I)
         W_00=T(I)**4
         W_out=T(I+1)**4
         BW=dlog(W_out/W_00)
         BK=dlog(K(I+1)/K(I))
         T1=-CL*R(I)**4/dm_bar(I)
         T2=(W_out/K(I+1)-W_00/K(I))/(1.d0-BK/BW)
         T3=T1/(BW-BK)
         DLK=  (T3/K(I))  *(W_00*BW/K(I)  -T2) !dL(i)/dK(i)
         DLKP=-(T3/K(I+1))*(W_out*BW/K(I+1)-T2) !dL(i)/dK(i+1)
         DLRP= DLKP*dK_dV_00(I+1)*DVR(I+1)
         DLRM= DLK *dK_dV_00(I)  *DVRM(I)         
         DLR= 4.d0*T1*T2/R(I)+DLK*dK_dV_00(I)*DVR(I)+DLKP*dK_dV_00(I+1)*DVRM(I+1)
         DLTP=4.d0*(T3/T(I+1))*(W_out*BW/K(I+1)-T2*BK/BW)+DLKP*dK_dT_00(I+1)
         DLT=-4.d0*(T3/T(I))*(W_00*BW/K(I)-T2*BK/BW)+DLK*dK_dT_00(I)
!-----------------------------
         DFCX0=DLCX0(I)
         DFCXM=DLCXM(I)
         DFCXP=DLCXP(I)
         DFCY0=DLCY0(I)
         DFCYP=DLCYP(I)
         DFCZ0=DLCZ0(I)
         DFCZP=DLCZP(I)
         goto 7
!        OUTER LUM BOUNDARY CONDITION
 6       continue
         DLT = 4.d0*L0/T(I)  !L=4piR^2sigT^4
         DLR = 2.d0*L0/R(I)  !L=4piR^2sigT^4
         DLRM = 0.d0
         DLRP = 0.d0
         DLTP = 0.d0
         DFCX0 = 0.d0
         DFCXM = 0.d0
         DFCXP = 0.d0
         DFCY0 = 0.d0
         DFCYP = 0.d0
         DFCZ0 = 0.d0
         DFCZP = 0.d0
 7       continue
         ELRP(I) = DLRP  +DFCXP +DLTXP(I)
         ELRM(I) = DLRM  +DFCXM +DLTXM(I)
         ELR(I)  = DLR   +DFCX0 +DLTX0(I)
         ELTP(I) = DLTP  +DFCYP +DLTYP(I)
         ELT(I)  = DLT   +DFCY0 +DLTY0(I)
         ELZ0(I) =        DFCZ0 +DLTZ0(I)
         ELZP(I) =        DFCZP +DLTZP(I) 

!        CALC ENERGY EQUATION(I)
         T2=P4/dm(I)
         T3=T2*(P(I)+DEV(I))
         T4=1.d0/dm(I)
         EX20(I) = -T4*(    -DLMRM) -dC_dr_in2(I) -T4*(     -DFCMXM)
         EX10(I) = -T4*(DLRM-DLMR ) -dC_dr_in(I)  -T4*(DFCXM-DFCMX0)
         EX00(I) = -T4*(DLR -DLMRP) -dC_dr_00(I)  -T4*(DFCX0-DFCMXP)
         EX01(I) = -T4*(DLRP      ) -dC_dr_out(I)  -T4*(DFCXP       )
         EY10(I) = -T4*(    -DLMT ) -dC_dT_in(I)  -T4*(     -DFCMY0)
         EY00(I) = -T4*(DLT -DLMTP) -dC_dT_00(I)  -T4*(DFCY0-DFCMYP)
         EY01(I) = -T4*(DLTP      ) -dC_dT_out(I)  -T4*(DFCYP       )
         EZ10(I) =                             -T4*(     -DFCMZ0)
         EZ00(I) =                  -dC_dw_00(I)  -T4*(DFCZ0-DFCMZP)
         EZ01(I) =                             -T4*(DFCZP       )
         if(I.eq.1) then
            EU10(I) =  T3*s% R_center**2
         else
            EU10(I) =  T3*R(max(1,I-1))**2
             ! bp: max(1,I-1) to prevent bogus warning from gfortran
         endif
         EU00(I) = -T3*R(I)**2

!        CALC CONVECTIVE ENERGY EQUATION(I)
         T3=P4/dm(I)*PTURB(I)
         T4=1.d0/dm(I)
         if (I == 1) then
            CX20(I) = dC_dr_in2(I) -T4*(        -0)
            CX10(I) = dC_dr_in(I)  -T4*(DLTXM(I)-0)
            CX00(I) = dC_dr_00(I)  -T4*(DLTX0(I)-0)
            CY10(I) = dC_dT_in(I)  -T4*(        -0)
            CY00(I) = dC_dT_00(I)  -T4*(DLTY0(I)-0)
            CZ10(I) =           -T4*(        -0)
            CZ00(I) = dC_dw_00(I)  -T4*(DLTZ0(I)-0)
         else
            CX20(I) = dC_dr_in2(I) -T4*(        -DLTXM(I-1))
            CX10(I) = dC_dr_in(I)  -T4*(DLTXM(I)-DLTX0(I-1))
            CX00(I) = dC_dr_00(I)  -T4*(DLTX0(I)-DLTXP(I-1))
            CY10(I) = dC_dT_in(I)  -T4*(        -DLTY0(I-1))
            CY00(I) = dC_dT_00(I)  -T4*(DLTY0(I)-DLTYP(I-1))
            CZ10(I) =           -T4*(        -DLTZ0(I-1))
            CZ00(I) = dC_dw_00(I)  -T4*(DLTZ0(I)-DLTZP(I-1))
         end if
         CY01(I) = dC_dT_out(I)  -T4*(DLTYP(I)           )
         CX01(I) = dC_dr_out(I)  -T4*(DLTXP(I)           )
         CZ01(I) =           -T4*(DLTZP(I)           )
         if(I.eq.1) then
            CU10(I) =  T3*s% R_center**2
         else
            CU10(I) =  T3*R(max(1,I-1))**2
             ! bp: max(1,I-1) to prevent bogus warning from gfortran
         endif
         CU00(I) = -T3*R(I)**2

!        CALC MOMENTUM EQUATION(I)
         T1=P4*R(I)**2/dm_bar(I)
         if(ALFAM.ge.0.d0) T4=P4/(dm_bar(I)*R(I))
         if(ALFAM.lt.0.d0) T4=-T1
         MU10(I) =  T4*(-EVUUM(I))
         MZ00(I) = -T1*(-dPtrb_dw_00(I))
         MX10(I) = -T1*(-dP_dr_in(I)-dPtrb_dr_in(I))
         MY00(I) = -T1*(-dP_dT_00(I))
         if(I.ne.NZN)then
            MX00(I) =  4.d0*G*M(I)/R(I)**3 &
                      -T1*(dP_dr_in(I+1)-dP_dr_00(I)+dPtrb_dr_in(I+1)-dPtrb_dr_00(I))
            MX01(I) = -T1*(dP_dr_00(I+1)        +dPtrb_dr_00(I+1))
            MY01(I) = -T1*(dP_dT_00(I+1))
            MU00(I) =  T4*(EVUUM(I+1)-EVUU0(I))
            MU01(I) =  T4*(EVUU0(I+1))
            MZ01(I) = -T1*(dPtrb_dw_00(I+1))
         else
            MX00(I) = 4.d0*G*M(I)/R(I)**3 &
                     -T1*(-dP_dr_00(I))
            MX01(I) = 0.d0
            MY01(I) = 0.d0
            MU00(I) = T4*(-EVUU0(I))
            MU01(I) = 0.d0
            MZ01(I) = 0.d0
         endif
      
  5    continue

      do I=1,NZN3
         do J=1,NZN3
            LLL(I,J)=0.d0
         enddo
      enddo

!
! VELOCITY DEFINITION
! (dR/dT) = u
! MOMENTUM EQUATION
! (dU/dt) = - (4 PI R^2)(dp/dm) - (GM)/R^2
! ENERGY EQUATION
! (c_v)(dT/dt) = - (p+(de/dV)_T)(dV/dR)U - (dLr/dm)
! TURBULENT ENERGY EQUATION
! (de_t/dt) = 
! SEE CODE DOCUMENTATION FOR LINEARIZATION
!
! EIGENVALUE PROBLEM: LLL*X=SIGMA*X
! X={dR1,dU1,dT1,dw1,...,dRN,dUN,dTN,dwN}
! VELOCITY DEF. EQ. IN ROWS IG
! MOMENTUM EQ. IN ROWS IR
! ENERGY EQ. IN ROWS IE
! TURBULENT ENERGY EQ. IN ROWS IC
! NAMING SCHEME FOR ELEMENTS OF LLL:
! EX00(I) - d(energy.eq.)/dR_i
! EX01(I) - d(energy.eq.)/dR_i+1
! EX10(I) - d(energy.eq.)/dR_i-1
! EY00(I) - d(energy.eq.)/dT_i
! EU00(I) - d(energy.eq.)/dU_i
! MY00(I) - d(moment.eq.)/dT_i
! e.t.c.
! VELOCITY DEF. AND MOMENTUM EQ. AT THE INTERFACE 1....NZN
! ENERGY EQUATIONS IN THE ZONE 1....NZN
!

      do I=1,NZN
         IG=4*I-3
         IR=4*I-2
         IE=4*I-1
         IC=4*I

         if(IG+1.le.NZN3) LLL(IG,IG+1) = 1.d0
!---------------------------------------
         if(IR-4.ge.1)    LLL(IR,IR-4)  = MU10(I)
                          LLL(IR,IR)    = MU00(I)
         if(IR+4.le.NZN3) LLL(IR,IR+4)  = MU01(I)

         if(IR-5.ge.1)    LLL(IR,IR-5)  = MX10(I)
         if(IR-1.ge.1)    LLL(IR,IR-1)  = MX00(I)
         if(IR+3.le.NZN3) LLL(IR,IR+3)  = MX01(I)

         if(IR+1.le.NZN3) LLL(IR,IR+1)  = MY00(I)
         if(IR+5.le.NZN3) LLL(IR,IR+5)  = MY01(I)

         if(IR+2.le.NZN3) LLL(IR,IR+2)  = MZ00(I)
         if(IR+6.le.NZN3) LLL(IR,IR+6)  = MZ01(I)
!---------------------------------------
         if(IE-5.ge.1)    LLL(IE,IE-5)  = EU10(I)/dE_dT_00(I)
         if(IE-1.ge.1)    LLL(IE,IE-1)  = EU00(I)/dE_dT_00(I)

         if(IE-10.ge.1)   LLL(IE,IE-10) = EX20(I)/dE_dT_00(I)
         if(IE-6.ge.1)    LLL(IE,IE-6)  = EX10(I)/dE_dT_00(I)
         if(IE-2.ge.1)    LLL(IE,IE-2)  = EX00(I)/dE_dT_00(I)
         if(IE+2.le.NZN3) LLL(IE,IE+2)  = EX01(I)/dE_dT_00(I)

         if(IE-4.ge.1)    LLL(IE,IE-4)  = EY10(I)/dE_dT_00(I)
                          LLL(IE,IE)    = EY00(I)/dE_dT_00(I)
         if(IE+4.le.NZN3) LLL(IE,IE+4)  = EY01(I)/dE_dT_00(I)

         if(IE-3.ge.1)    LLL(IE,IE-3)  = EZ10(I)/dE_dT_00(I)
         if(IE+1.le.NZN3) LLL(IE,IE+1)  = EZ00(I)/dE_dT_00(I)
         if(IE+5.le.NZN3) LLL(IE,IE+5)  = EZ01(I)/dE_dT_00(I)
!---------------------------------------
         if(IC-11.ge.1)   LLL(IC,IC-11) = CX20(I)
         if(IC-7.ge.1)    LLL(IC,IC-7)  = CX10(I)
         if(IC-3.ge.1)    LLL(IC,IC-3)  = CX00(I)
         if(IC+1.le.NZN3) LLL(IC,IC+1)  = CX01(I)
  
         if(IC-6.ge.1)    LLL(IC,IC-6)  = CU10(I) 
         if(IC-2.ge.1)    LLL(IC,IC-2)  = CU00(I) 

         if(IC-5.ge.1)    LLL(IC,IC-5)  = CY10(I) 
         if(IC-1.ge.1)    LLL(IC,IC-1)  = CY00(I) 
         if(IC+3.le.NZN3) LLL(IC,IC+3)  = CY01(I)

         if(IC-4.ge.1)    LLL(IC,IC-4)  = CZ10(I) 
                          LLL(IC,IC)    = CZ00(I)
         if(IC+4.le.NZN3) LLL(IC,IC+4)  = CZ01(I)
         
         IF (IE+4 <= NZN3) then
            if(LLL(IE,IE+4).lt.0.d0)then
               !write(*,*) 'rerrrrrrrrrrrrrrrrrrrrrrrr',i
            endif
         endif
      enddo
      
      if (s% RSP_trace_RSP_build_model) &
         write(*,*) 'waiting for DGEEV to solve eigenvalue problem....'
      call DGEEV('n','v',NZN3,LLL,LD_LLL,WRx,WIx,VLx,LD_VL,VRx,LD_VR, &
                 WORKx,4*NZN3,INFO)         
      if(INFO.ne.0)then
         write(*,*) 'FAILED!'
         write(*,*) 'LAPACK/DGEEV error, ier= ',INFO
         ierr = -1
         return
         stop
      endif

! THIS IS A MODIFIED "NUMERICAL RECIPIES" SORTING SUBROUTINE
! IT SORTS THE VECTOR WI (DIMENSION NZN3) IN ASCENDING ORDER
! AND IN THE SAME WAY REARRANGES THE ELEMENTS OF WR (WR IS NOT SORTED)
! THE INTEGER VECTOR ISORT CONTAINS INFORMATION ABOUT CHANGES DONE
! NAMELY ELEMENT I OF SORTED WI (AND REARRANGED WR) WAS ISORT(I) ELEMENT
! OF UNSORTED WI (AND NOT REARANGED WR)
      call SORT(NZN3,WIx,WRx,ISORTx)

!     THIS LOOP FINDS THE SMALLEST BUT POSITIVE VALUE OF VECTOR WI (MINI)
!     AND RETURNS THE CORRESPONDING INDEX (IMI)
      FILENAME=trim(s% log_directory) // '/' // 'LINA_period_growth_info.data'
      !open(15,file=trim(FILENAME),status='unknown')
      MINI=5.d0
      IMI=0
      !write(15,'(a)') '     I                R              PERIOD'
      do I=1,NZN3
!         if(WIx(I).lt.MINI.and.WIx(I).gt.1d-9)then
         if((WIx(I).lt.MINI) .and.(WIx(I).gt.1.d-9))then
            if(P4*WRx(I)/WIx(I).gt.-.3d+1)then
               MINI=WIx(I)
               IMI=I
            end if
         endif
         if (abs(WIx(I)) > 1d-50) then
            !write(15,'(2(d14.8,tr2),f11.6)') &
            !   WIx(I),WRx(I),2.d0*PI/WIx(I)/86400d0
         else
            !write(15,'(2(d14.8,tr2),f11.6)') &
            !   0d0,WRx(I),0d0
         end if
      enddo
      !close(15)

      do J=1,NMODES
 444     continue
         if(P4*WRx(IMI+J-1)/WIx(IMI+J-1).lt.-.5d+1)then
            IMI=IMI+1
            goto 444
         endif
!        VRRS(J) IS THE MODULI OF THE SUTFACE R-EIGENVECTOR OF THE MODE J
         VRRS(J)=sqrt(VRx(4*NZN-3,ISORTx(IMI+J-1))**2+ &
                       VRx(4*NZN-3,ISORTx(IMI+J-1)+1)**2) !surface value
         VTTS(J)=VRx(4*NZN-3,ISORTx(IMI+J-1))+(0.d0,1.d0)* &
                       VRx(4*NZN-3,ISORTx(IMI+J-1)+1)
 
         SCALE(J)=R(NZN)/VTTS(J)
!        PERS(J) IS THE PERIOD OF THE MODE J (IN SECONDS)
         OMEG(J)=WIx(IMI+J-1)
         PERS(J)=2.d0*PI/OMEG(J)
         Q(J)=PERS(J)*SQRT((M(NZN)/SUNM)*(SUNR/R(NZN))**3)/86400.d0
!        ETO(J) IS THE GROWTH RATE OF MODE J
         ETO(J)= P4*WRx(IMI+J-1)/OMEG(J)
         EK(J)=0.d0
         SGRP=0.d0
         SGRM=0.d0
         NORMC=-1.d50
         do I=1,NZN
!           SEE LAPACK USERS GUIDE FOR CONSTRUCTION OF EIGENVECTORS

!           VRR(I,J) IS THE dR EIGENVECTOR OF THE MODE J
!           dR_i 
            VRR(I,J)=VRx(4*I-3,ISORTx(IMI+J-1))+(0.d0,1.d0)* &
                     VRx(4*I-3,ISORTx(IMI+J-1)+1)
!           DVRR(I,J) IS SCALED dR/R EIGENVECTOR OF THE MODE (J)
!           (dR_i/R_i)/(dR_NZN/R_NZN)
            DVRR(I,J)=VRR(I,J)/R(I)*SCALE(J)
            PHR(I,J)=datan2(aimag(DVRR(I,J)),dble(DVRR(I,J)))

!           VRT(I,J) IS THE dT EIGENVECTOR OF THE MODE J 
!           dT_i
            VRT(I,J)=VRx(4*I-1,ISORTx(IMI+J-1))+(0.d0,1.d0)* &
                     VRx(4*I-1,ISORTx(IMI+J-1)+1)
!           DVRT(I,J) IS SCALED dT/T EIGENVECTOR OF THE MODE (J)
!           (dT_i/T_i)/(dR_NZN/R_NZN)
            DVRT(I,J)=VRT(I,J)/T(I)*SCALE(J)
            PHT(I,J)=datan2(aimag(DVRT(I,J)),dble(DVRT(I,J)))

!           VRC(I,J) IS THE dT EIGENVECTOR OF THE MODE J 
!           dOMEGA_i
            VRC(I,J)=VRx(4*I,ISORTx(IMI+J-1))+(0.d0,1.d0)* &
                     VRx(4*I,ISORTx(IMI+J-1)+1)
!           DVRC(I,J) IS SCALED dOMEGA/OMEGA EIGENVECTOR OF THE MODE (J)
!           (dOMEGA_i/OMEGA_i)/(dR_NZN/R_NZN)
            DVRC(I,J)=VRC(I,J)
            if(NORMC.lt.abs(DVRC(I,J))) NORMC=abs(DVRC(I,J))
            PHC(I,J)=datan2(aimag(DVRC(I,J)),dble(DVRC(I,J)))

!           VRU(I,J) IS THE dU EIGENVECTOR OF THE MODE J 
!           dU_i
            VRU(I,J)=VRx(4*I-2,ISORTx(IMI+J-1))+(0.d0,1.d0)* &
                     VRx(4*I-2,ISORTx(IMI+J-1)+1)
            PHU(I,J)=datan2(aimag(VRU(I,J)),dble(VRU(I,J)))

!           KINETIC ENERGY OF THE MODE
            EK(J)=EK(J)+0.5d0*dm_bar(I)*(OMEG(J)*abs(VRR(I,J)*SCALE(J)))**2

!           WORK DONE IN ZONE I (FOR THE MODE J)
            if(I.eq.1) then
               DV_0=DVR(I)*VRR(I,J)*SCALE(J)
            else
               DV_0=(DVR(I)*VRR(I,J)+DVRM(I)*VRR(I-1,J))*SCALE(J)
            end if
            DP_0=(DPV(I)*DV_0+dP_dT_00(I)*VRT(I,J)*SCALE(J))

            DPEV=EVUU0(I)*VRU(I,J)
            if(I.ge.2)DPEV=DPEV+EVUUM(I)*VRU(I-1,J)
            DPEV=DPEV*SCALE(J)

            dP_dT_00URB=dPtrb_dr_00(I)*VRR(I,J)
            if(I.ge.2) dP_dT_00URB=dP_dT_00URB+dPtrb_dr_in(I)*VRR(I-1,J)
            dP_dT_00URB=dP_dT_00URB+dPtrb_dw_00(I)*VRC(I,J)
            dP_dT_00URB=dP_dT_00URB*SCALE(J)

            QWK(I,J)=-PI*dm(I)*aimag(conjg(DP_0)*DV_0) 
            if(ALFAM.lt.0.d0)QWKEV(I,J)=-PI*dm(I)*aimag(conjg(DPEV)*DV_0) 
            QWKPT(I,J)=-PI*dm(I)*aimag(conjg(dP_dT_00URB)*DV_0) 
            if(ALFAM.gt.0.d0)then
               QWKEV(I,J)=PI*dm(I)*aimag(conjg(DPEV)*(DV_0/R(I)**3- &
                        3.d0*Vol(I)/R(I)**4*VRR(I,J)*SCALE(J)))
            endif

            if(QWK(I,J)+QWKEV(I,J)+QWKPT(I,J).ge.0.d0) &
              SGRP=SGRP+QWK(I,J)+QWKEV(I,J)+QWKPT(I,J)
            if(QWK(I,J)+QWKEV(I,J)+QWKPT(I,J).lt.0.d0) &
              SGRM=SGRM+abs(QWK(I,J)+QWKEV(I,J)+QWKPT(I,J))
        
            VEL(I,J)=abs(VRR(I,J))/VRRS(J)
            if(abs(PHR(I,J)).gt.1.57d0)VEL(I,J)=-VEL(I,J)
         enddo

!        WRITE WORK-INTEGRALS INTO FILE
         if(J.le.9)then
            write(NUMER1,'(I1)') J
            FILENAME=trim(s% log_directory) // '/' // 'LINA_work'//NUMER1//'.data'
         endif
         if(J.ge.10)then
            write(NUMER2,'(I2)') J
            FILENAME=trim(s% log_directory) // '/' // 'LINA_work'//NUMER2//'.data'
         endif

         open(57,file=trim(FILENAME),status='unknown')
         write(57,'(a)') &
            '#ZONE  log(T)         X              WORK(P)' // &
            '         WORK(P_NU)      WORK(P_T)        CWORK(P)' // &
            '        CWORK(P_NU)        CWORK(P_T)'
         ETOI(J)=0.d0
         ETOIEV(J)=0.d0
         ETOIPT(J)=0.d0
         do I=1,NZN
            ETOI(J)  =ETOI(J)  +QWK(I,J)
            ETOIEV(J)=ETOIEV(J)+QWKEV(I,J)
            ETOIPT(J)=ETOIPT(J)+QWKPT(I,J)
            write(57,'(I3,2(tr1,e15.8),6(tr1,e15.8))') &
              NZN+1-I,dlog10(T(I)),R(I)/R(NZN), &
                QWK(I,J)/EK(J), QWKEV(I,J)/EK(J), QWKPT(I,J)/EK(J), &
                ETOI(J)/EK(J),ETOIEV(J)/EK(J),ETOIPT(J)/EK(J)
         enddo
         write(57,'("#KINETIC ENERGY:",e15.8)') EK(J)
         close(57)

!        CALCULATE LUMINOSITY EIGENFUNCTIONS
         do I=1,NZN
!           VRL(I,J) IS THE dL EIGENVECTOR OF THE MODE J 
!           dL_i
            VRL(I,J)= ELT(I)*VRT(I,J) &
                     +ELR(I)*VRR(I,J) &
                     +ELZ0(I)*VRC(I,J)
            if(I.ne.1)   VRL(I,J)=VRL(I,J)+ELRM(I)*VRR(I-1,J)
            if(I.ne.NZN) VRL(I,J)=VRL(I,J)+ELTP(I)*VRT(I+1,J) &
                                          +ELRP(I)*VRR(I+1,J) &
                                          +ELZP(I)*VRC(I+1,J)
!           DVRL(I,J) IS SCALED dL/L EIGENVECTOR OF THE MODE (J)
!           (dL_i/L0)/(dR_NZN/R_NZN)
            DVRL(I,J)=VRL(I,J)/L0*SCALE(J)
            PHL(I,J)=datan2(aimag(DVRL(I,J)),dble(DVRL(I,J)))
         enddo

!        WRITE EIGENVECTORS TO FILE
         if(J.le.9)then
            write(NUMER1,'(I1)') J
            FILENAME=trim(s% log_directory) // '/' // 'LINA_eigen'//NUMER1//'.data'
         endif
         if(J.ge.10)then
            write(NUMER2,'(I2)') J
            FILENAME=trim(s% log_directory) // '/' // 'LINA_eigen'//NUMER2//'.data'
         endif
         open(56,file=trim(FILENAME),status='unknown')
         write(56,'(a)') &
            '#ZONE  TEMP.           FRAC. RADIUS  ABS(dR/R)       ' // &
            'PH(dR/R)        ABS(dT/T)       PH(dT/T)        ' // &
            'ABS(dL/L)       PH(dL/L)        ABS(dE_T)       PH(dE_T)'
         do I=1,NZN
            write(56,'(I3,tr1,e15.6,tr1,e15.6,8(tr1,e15.8))') &
                NZN+1-i,T(I),R(I)/R(NZN),abs(DVRR(I,J)),PHR(I,J), &
                abs(DVRT(I,J)),PHT(I,J), &
                abs(DVRL(I,J)),PHL(I,J), &
                abs(DVRC(I,J))/NORMC,PHC(I,J)
         end do
         close(56)

         SGR(J)=(SGRP-SGRM)/(SGRP+SGRM)

         PSIG=M(NZN)/(P43*R(NZN)**3)
         PSIG=sqrt(P4*G*PSIG)

         QCHECK(J)=100.d0*(ETO(J)-(ETOI(J)+ETOIEV(J)+ETOIPT(J))/EK(J)) &
                          /ETO(J)
 19      format('nonadiabatic solution, MODE: ',I2)
 20      format('eigenvalue: (',1P,D15.8,' , ',D15.8,')')
 21      format(14X,'freq.[Hz]: ',D11.5,' , sigma: ',D12.6)
 22      format(14X,'period[d]:',F11.6,1X,' ,     Q:',F9.6)
 23      format(9X,'KE growte rate:',F13.9,2X,',    KE: ',D12.6)
 24      format('control: KE growte rate:',F13.9)
 25      format('         f: (',1P,D15.8,' , ',D15.8,'), abs(f):',D15.8)
 26      format('control: f: (',1P,D15.8,' , ',D15.8,'), abs(f):',D15.8)

      enddo
      
      !write(*,'(a55,i12,99(1pd26.16))') 'end LINA s% w(2)**2 Et', &
      !   2, s% w(2)**2, Et(NZN-1)

      return
      end subroutine do_LINA
     
     
      SUBROUTINE mesa_eos_kap (s,k,G,H, &   !input: temp,volume  &
               P,PV,PT,E,EV,ET,CP,CPV,dCp_dT_00, &
               Q,QV,QT,OP,OPV,OPT,ierr) 
      use rsp_eval_eos_and_kap, only : eval_mesa_eos_and_kap
      implicit none
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      integer :: k, j
      real(8) :: G,H,P,PV,PT, &
         E,EV,ET,CP,CPV,dCp_dT_00, &
         Q,QV,QT,OP,OPV,OPT,cs, &
         Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
         egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT
      include 'formats'
      if (k <= 0 .or. k > NZN) then
         j = 0
      else
         j = NZN+1-k
      end if
      if (is_bad(G+H)) then
         write(*,2) 'LINA mesa_eos_kap G H', k, G, H
         call mesa_error(__FILE__,__LINE__,'mesa_eos_kap')
      end if
      call eval_mesa_eos_and_kap(s,j,G,H, &
               Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT,&
               egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
               cs,CP,CPV,dCp_dT_00, &
               Q,QV,QT,OP,OPV,OPT,ierr) 
      if (ierr /= 0) return
      E = egas + erad
      EV = d_egas_dV + d_erad_dV
      ET = d_egas_dT + d_erad_dT
      P = Pgas + Prad
      PV = d_Pg_dV
      PT = d_Pg_dT + d_Pr_dT
      end SUBROUTINE mesa_eos_kap


! THIS IS A MODIFIED "NUMERICAL RECIPIES" SORTING SUBROUTINE
! IT SORTS THE VECTOR RA (DIMENSION N) IN ASCENDING ORDER
! AND IN THE SAME WAY REARRANGES THE ELEMENTS OF RB (RB IS NOT SORTED)
! THE INTEGER VECTOR ISORT CONTAINS INFORMATION ABOUT CHANGES DONE
! NAMELY ELEMENT I OF SORTED RA (AND REARRANGED RB) WAS ISORT(I) ELEMENT
! OF UNSORTED RA (AND NOT REARANGED RB)
      subroutine SORT(N,RA,RB,ISORT)
      implicit none
      integer :: N,ISORT(N)
      real(dp) :: RA(N),RB(N)
      integer :: L,IR,I,J,RRI
      real(dp) ::  RRA,RRB
   
      do I=1,N
         ISORT(I)=I
      enddo

      L=N/2+1
      IR=N
10    continue
        if(L.gt.1)then
           L=L-1
           RRA=RA(L)
           RRB=RB(L)
           RRI=ISORT(L)
        else
           RRA=RA(IR)
           RRB=RB(IR)
           RRI=ISORT(IR)
           RA(IR)=RA(1)
           RB(IR)=RB(1)
           ISORT(IR)=ISORT(1)
           IR=IR-1
           if(IR.eq.1)then
              RA(1)=RRA
              RB(1)=RRB
              ISORT(1)=RRI
              return
           endif
        endif
        I=L
        J=L+L
 20     if(J.le.IR)then
           if(J.lt.IR)then
              if(RA(J).lt.RA(J+1))J=J+1
           endif
           if(RRA.lt.RA(J))then
              RA(I)=RA(J)
              RB(I)=RB(J)
              ISORT(I)=ISORT(J)
              I=J
              J=J+J
           else
              J=IR+1
           endif
           goto 20
        endif
        RA(I)=RRA
        RB(I)=RRB
        ISORT(I)=RRI
      goto 10
      end subroutine SORT

            
      end module rsp_lina
