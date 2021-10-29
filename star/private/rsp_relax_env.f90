! ***********************************************************************
!
!   Copyright (C) 2018-2019  Radek Smolec & The MESA Team
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

      module rsp_relax_env
      use star_def, only: star_info
      use utils_lib, only: is_bad, mesa_error
      use const_def, only: dp, crad
      use eos_lib, only: Radiation_Pressure
      use rsp_def
      use rsp_eval_eos_and_kap, only: X, Y, Z
      use rsp_lina, only: mesa_eos_kap
      
      implicit none
      
      private
      public :: EOP, RELAX_ENV
      
      contains


      subroutine EOP(s,k,T,P,V,E,CP,QQ,SVEL,OP,ierr)
      use rsp_eval_eos_and_kap, only: &
         eval1_mesa_Rho_given_PT, eval1_gamma_PT_getRho
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      real(dp) :: T,V,P1,DPV,E,OP,F2,P,U,SVEL,CP,QQ
      real(dp) :: DU1,DU2,DU3,DU4,DU5,DU6,DU7,DU8,DU9,DU10,DU11,DU12
      integer :: k,I
      real(dp) :: RES(12)

      real(dp) :: R,A,PRECEQ,rho_guess,rho,Prad,P_test
      data R,A/8.317d7,7.5641d-15/
      data PRECEQ/5.d-13/

      if(P.le.0.d0) goto 100 !negative pressure, stop
      
      ierr = 0
      Prad = Radiation_Pressure(T)
      if (P <= Prad) then
         write(*,*) 'P <= Prad', k, P, Prad, T
         ierr = -1
         return
         call mesa_error(__FILE__,__LINE__,'EOP')
      end if
      rho_guess = eval1_gamma_PT_getRho(s, 0, P, T, ierr)
      if (ierr /= 0) then
         write(*,*) 'eval1_gamma_PT_getRho failed', k
         write(*,*) 'P', P
         write(*,*) 'T', T
         ierr = -1
         return
         call mesa_error(__FILE__,__LINE__,'EOP')
      end if
      call eval1_mesa_Rho_given_PT(s, 0, P, T, rho_guess, rho, ierr)
      if (ierr /= 0) then
         write(*,*) 'eval1_mesa_Rho_given_PT failed', k
         write(*,*) 'P', P
         write(*,*) 'T', T
         write(*,*) 'rho_guess', rho_guess
         ierr = -1
         return
         call mesa_error(__FILE__,__LINE__,'EOP')
      end if
      V = 1d0/rho ! initial guess to be improved below
      
         !write(*,*) 'eval1_mesa_Rho_given_PT', k
         !write(*,*) 'P', P
         !write(*,*) 'T', T
         !write(*,*) 'rho_guess', rho_guess
         !write(*,*) 'V', 0, V
      
      I=0
      !     NEWTON-RAPHSON ITERATION (TO MAKE P1->P)
 1    I=I+1
      call mesa_eos_kap(s,0, &
        T,V,P1,DPV,DU1,E,DU2,DU3,CP,DU4,DU5,QQ,DU6,DU7,OP &
           ,DU8,DU9,ierr)
      if (ierr /= 0) return
      if(I.gt.25) goto 2 !no convergence, stop
      F2=P1/P-1.d0
      V=V-F2*P/DPV
      if (is_bad(V)) then
         write(*,*) 'V', V
         write(*,*) 'F2', F2
         write(*,*) 'DPV', DPV
         write(*,*) 'P', P
         write(*,*) 'EOP I', I
         ierr = -1
         return
         stop
      end if
      
         !write(*,*) 'V', I, V
         
      if(abs(F2).lt.PRECEQ) goto 3
      goto 1
 3    continue
       if (abs(P - P_test) > 1d-10*P) then
       end if
      return

 2    write(*,*) 'NO CONVERGENCE IN EOP; T,P= ',T,P
         ierr = -1
         return
      stop
 100  write(*,*) 'NEGATIVE OR ZERO PRESSURE IN EOP'
         ierr = -1
         return
      stop
      end subroutine EOP


      subroutine RELAX_ENV(s, L0, TH0, TE, NZT, NZN, &         
         M, DM, DM_BAR, R, Vol, T, Et, ierr)
      use star_utils, only: rand
      type (star_info), pointer :: s
      integer, intent(in) :: NZN, NZT
      real(dp), intent(in) :: L0, TH0, TE
      real(dp), intent(inout), dimension(:) :: &
         M, DM, DM_BAR, R, Vol, T, Et
      integer, intent(out) :: ierr
      
      logical, parameter :: RSP_eddi = .true. ! use Eddington approx at surface
      
      real(dp), allocatable, dimension(:) :: &
         E, P, Lr, Lc, Hp_face, Y_face, K, CPS, QQS, &
         dC_dr_in,dC_dr_00,dC_dr_out,dC_dT_00,dC_dT_out, &
         dC_dw_00,dC_dr_in2,dC_dT_in, &
         DPV,dP_dT_00,DEV,dE_dT_00,dK_dV_00, &
         dK_dT_00,DVR,DVRM, &
         CPV, dCp_dT_00, QQV,dQQ_dT_00, & 
         dP_dr_00, dCp_dr_00, dQQ_dr_00, &
         dP_dr_in,dCp_dr_in,dQQ_dr_in, &
         dHp_dr_in,dHp_dr_out,dHp_dr_00,dHp_dT_00,dHp_dT_out, &
         dY_dr_in,dY_dr_00,dY_dr_out,dY_dT_00, &
         dY_dT_out, &
         PII,dPII_dr_in,dPII_dr_00,dPII_dr_out, &
         dPII_dT_00,dPII_dT_out,DPIIZ0, &
         dsrc_dr_in,dsrc_dr_00,dsrc_dr_out,dsrc_dT_00, &
         dsrc_dT_out,dsrc_dw_00,SOURCE, &
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
         PTURB,dPtrb_dw_00,dPtrb_dr_00,dPtrb_dr_in
      
      real(dp) :: T1,DLR,DLRP,DLRM,DLT,DLTP,DLMR,DLMRP,DLMRM,DLMT,DLMTP
      real(dp) :: W_00,W_out,BW,BK,T2,T3,DLK,DLKP,SVEL,T_0
      real(dp) :: DU1,DU2,POM3,POM2,POM,MAXFLUXD
      integer :: I,J,IW,IR,IC, IMAXR, IMAXW, IMAXC
      real(dp) :: MAXR,MAXW,MAXC,PB,PPB,ELB,ELMB,RES(12)
      
      integer :: INFO,II,ITROUBC,ITROUBT,IZIP
      real(dp) :: XXR,XXC,XXT,EZH,DXH,DXXC,DXXT,DXKC,DXKT, &
         IGR1,IGR1XM,IGR1X0,IGR1XP,IGR1Y0,IGR1YP, &
         IGR2,IGR2XM,IGR2X0,IGR2XP,IGR2Y0,IGR2YP
      
      real(dp) :: FFXM,FFX0,FFXP,FFY0,FFYP,FF
      real(dp) :: GGXM,GGX0,GGXP,GGY0,GGYP,GG,GPF
      real(dp) :: FLIM,FLD,EFL02

      real(dp) :: POM4,DXXR,TEM1,TEMI,TEMM

      real(dp) :: TT,dmN,dmNL,TNL,DDT,MSTAR,PRECR
      real(dp) :: RINNER,ROUTER,TINNER,MENVEL
      integer :: IOP,NEGFLU
      logical ending

      real(dp) :: SUMM,AONE,HAHA
      real(dp) :: AAA(100),AAP(100),AAT(100)
      real(dp) :: AALFA, AALFAP, AALFAT
      integer :: ICAA, ICAP, ICAT
      integer :: NDIVAA, NDIVAP, NDIVAT, dmN_cnt, max_dmN_cnt
      logical :: ok_to_adjust_mass, ok_to_adjust_Tsurf
      real(dp) :: EDFAC, Psurf, CFIDDLE, FSUB, EMR, ELR
      real(dp) :: PREC1, SOL, DAMPS,DAMPRS,DELTA,SOURS

      real(dp) :: XX_max, XX_max_val, XX_max_dx
      integer :: i_XX_max, var_XX_max, n, op_err

      !write(*,*) 'RELAX_ENV'
      EFL02 = EFL0*EFL0
      FSUB = s% RSP_dq_1_factor         
      EMR = s% RSP_mass
      ELR = s% RSP_L

      n = NZN+1
      allocate( &
         E(n), P(n), Lr(n), Lc(n), Hp_face(n), Y_face(n), K(n), CPS(n), QQS(n), &
         dC_dr_in(n),dC_dr_00(n),dC_dr_out(n),dC_dT_00(n),dC_dT_out(n), &
         dC_dw_00(n),dC_dr_in2(n),dC_dT_in(n), &
         DPV(n),dP_dT_00(n),DEV(n),dE_dT_00(n),dK_dV_00(n), &
         dK_dT_00(n),DVR(n),DVRM(n), &
         CPV(n), dCp_dT_00(n), QQV(n),dQQ_dT_00(n), &
         dP_dr_00(n), dCp_dr_00(n), dQQ_dr_00(n), &
         dP_dr_in(n),dCp_dr_in(n),dQQ_dr_in(n), &
         dHp_dr_in(n),dHp_dr_out(n),dHp_dr_00(n),dHp_dT_00(n),dHp_dT_out(n), &
         dY_dr_in(n),dY_dr_00(n),dY_dr_out(n),dY_dT_00(n), &
         dY_dT_out(n), &
         PII(n),dPII_dr_in(n),dPII_dr_00(n),dPII_dr_out(n), &
         dPII_dT_00(n),dPII_dT_out(n),DPIIZ0(n), &
         dsrc_dr_in(n),dsrc_dr_00(n),dsrc_dr_out(n),dsrc_dT_00(n), &
         dsrc_dT_out(n),dsrc_dw_00(n),SOURCE(n), &
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
         PTURB(n),dPtrb_dw_00(n),dPtrb_dr_00(n),dPtrb_dr_in(n))

!     STORE SOME VALUES FOR FURTHER COMPARISON
      MSTAR  = M(NZN)
      TINNER = T(1)
      RINNER = s% R_center ! R0
      ROUTER = R(NZN)
      MENVEL = MSTAR-M(1)+dm(1)

      TNL = 0d0 

      if (s% RSP_use_Prad_for_Psurf) then
         if(.not.RSP_eddi) then !     EXACT GREY RELATION
            T_0= pow(sqrt(3.d0)/4.d0,0.25d0)*TE !0.811194802d0*TE
         else !     EDDINGTON APPROXIMATION
            T_0= pow(0.5d0, 0.25d0)*TE ! T_0= pow(1.0d0/2.d0,0.25d0)*TE 
         endif      
         Psurf = crad*T_0*T_0*T_0*T_0/3d0
      else
         Psurf = 0d0
      end if

      ending=.false.
      if (s% RSP_trace_RSP_build_model) &
         write(*,*) '*** relax envelope ***'  
      
      NEGFLU=0
      NDIVAA = 20 ! s% RSP_NDIVAA
      NDIVAP = 20 ! s% RSP_NDIVAP
      NDIVAT = 20 ! s% RSP_NDIVAT
      ok_to_adjust_mass = .true. ! s% RSP_ok_to_adjust_mass
      ok_to_adjust_Tsurf = .true. ! s% RSP_ok_to_adjust_Tsurf
      dmN_cnt = 0
      max_dmN_cnt = s% RSP_relax_max_tries
      PREC1 = s% RSP_relax_dm_tolerance

!     SUMM IS ZONE OF THE ENVELOPE UP TO ANCHOR 
!     (NOT CHANGED IN THE ITERATIONS)
      SUMM=0.d0
      do I=1,NZN-NZT+1
         SUMM=SUMM+dm(I)
      enddo

!     IF BOTH ALFAP NAD ALFAT .ne. 0, IT IS NECESSARY TO ITERATE WITHOUT
!     TURBULENT FLUX, WITH TURBULENT PRESSURE ONLY, AND AFTER
!     CONVERGENCE TURBULENT FLUX IS ITERATED
      AALFA = -1 ! turn off ALFA relax
      if (AALFA <= 0d0) AALFA = ALFA
      AALFAT = ALFAT
      AALFAP = ALFAP
      ALFAT  = 0.d0
      ALFAP  = 0.d0

!     SET ALFA TO ITERATE
      do I=1,NDIVAA
         AAA(I)=ALFA+(AALFA-ALFA)*I/dble(NDIVAA)
      enddo
      ICAA=1
      
!     SET ITERATIONS FOR ALFAP
      if(AALFAP.ne.0.d0)then
         do I=1,NDIVAP
            AAP(I)=AALFAP*I/dble(NDIVAP)
         enddo
         ICAP=1
      endif

!     SET ITERATIONS FOR ALFAT
      if(AALFAT.ne.0.d0)then
         do I=1,NDIVAT
            AAT(I)=AALFAT*I/dble(NDIVAT)
         enddo
         ICAT=1
      endif
!-
      PRECR = 1.d-10  !PRECISION FOR NEWTON-RHAPSON ITERAIONS
      DXH = 0.01d0   !UNDERCORRECTION FOR NEWTON-RHAPSON CORRECTIONS
      DDT = -dm(NZN)/1000.d0 !INITIAL CHANGE IN OUTER ZONE MASS

      IOP = 0           

 999  continue

      IZIP = 0
      if(IOP.eq.0) goto 100

!------ MASS TRICKS ---
!     dmN = dm(NZN) ! THIS IS NOT WORKING IF FSUB.ne.1
      dmN = dm(NZN-1) !THIS IS WORKING FOR ALL FSUB VALUES
      
      TT  = T(NZN-NZT+1)-TH0

      if(IOP.eq.1) goto 24

      DDT = TT*(dmN-dmNL)/(TT-TNL)
      !if (is_bad(DDT)) then
      !   write(*,*) 'DDT', DDT
      !   write(*,*) 'TT-TNL', TT-TNL
      !   write(*,*) 'dmN-dmNL', dmN-dmNL
      !   write(*,*) 'TT', TT
      !   write(*,*) 'TNL', TNL
      !   write(*,*) 'dmNL', dmNL
      !   write(*,*) 'dmN', dmN
      !   call mesa_error(__FILE__,__LINE__,'failed in RELAX_ENV')
      !end if
      CFIDDLE=0.1d0
      if(abs(DDT/dmN).gt.CFIDDLE) DDT=CFIDDLE*dmN*(DDT/abs(DDT))

!     CHECK IF ALFA ITERATION IS FINISHED 
      if(abs(DDT/dmN).lt.PREC1.and.ALFA.ne.AALFA) then
         ALFA=AAA(ICAA)
         ICAA=ICAA+1
!        APPLY ARTIFICIAL ZONE MASS CHANGE TO ALLOW ITERATIONS
         DDT = -dm(NZN)/100000.d0
         write(*,*) 'ALFA =',ALFA
         goto 24
      endif
      if(abs(DDT/dmN).lt.PREC1.and.ICAA.eq.NDIVAA+1)then
         write(*,*) '***** ALFA ITERATION FINISHED *****'
         write(*,*) 'final ALFA =', ALFA, s% RSP_alfa
         dmN_cnt = 0
         NDIVAA=-99
      endif
      
      if (s% RSP_relax_alfap_before_alfat) then
      
         ! CHECK IF ALFAP ITERATION IS FINISHED 
         if(abs(DDT/dmN).lt.PREC1.and.ALFAP.ne.AALFAP) then
            ALFAP=AAP(ICAP)
            ICAP=ICAP+1
            ! APPLY ARTIFICIAL ZONE MASS CHANGE TO ALLOW ITERATIONS
            DDT = -dm(NZN)/100000.d0
            write(*,*) 'ALFAP=',ALFAP
            goto 24
         endif
         if(abs(DDT/dmN).lt.PREC1.and.ICAP.eq.NDIVAP+1)then
            write(*,*) '***** ALFAP ITERATION FINISHED *****'
            dmN_cnt = 0
            NDIVAP=-99
         endif

         ! CHECK IF ALFAT ITERATION IS FINISHED 
         if(abs(DDT/dmN).lt.PREC1.and.ALFAT.ne.AALFAT) then
            ALFAT=AAT(ICAT)
            ICAT=ICAT+1
            ! APPLY ARTIFICIAL ZONE MASS CHANGE TO ALLOW ITERATIONS
            DDT = -dm(NZN)/100000.d0
            write(*,*) 'ALFAT =',ALFAT
            goto 24
         endif
         if(abs(DDT/dmN).lt.PREC1.and.ICAT.eq.NDIVAT+1)then
            write(*,*) '***** ALFAT ITERATION FINISHED *****'
            dmN_cnt = 0
            NDIVAT=-99
         endif
      
      else

         ! CHECK IF ALFAT ITERATION IS FINISHED 
         if(abs(DDT/dmN).lt.PREC1.and.ALFAT.ne.AALFAT) then
            ALFAT=AAT(ICAT)
            ICAT=ICAT+1
            ! APPLY ARTIFICIAL ZONE MASS CHANGE TO ALLOW ITERATIONS
            DDT = -dm(NZN)/100000.d0
            write(*,*) 'ALFAT =',ALFAT
            goto 24
         endif
         if(abs(DDT/dmN).lt.PREC1.and.ICAT.eq.NDIVAT+1)then
            write(*,*) '***** ALFAT ITERATION FINISHED *****'
            dmN_cnt = 0
            NDIVAT=-99
         endif
      
         ! CHECK IF ALFAP ITERATION IS FINISHED 
         if(abs(DDT/dmN).lt.PREC1.and.ALFAP.ne.AALFAP) then
            ALFAP=AAP(ICAP)
            ICAP=ICAP+1
            ! APPLY ARTIFICIAL ZONE MASS CHANGE TO ALLOW ITERATIONS
            DDT = -dm(NZN)/100000.d0
            write(*,*) 'ALFAP=',ALFAP
            goto 24
         endif
         if(abs(DDT/dmN).lt.PREC1.and.ICAP.eq.NDIVAP+1)then
            write(*,*) '***** ALFAP ITERATION FINISHED *****'
            dmN_cnt = 0
            NDIVAP=-99
         endif
      
      end if

      if(abs(DDT/dmN).lt.PREC1) then
         if(NEGFLU.ge.1)then
            if (NEGFLU.eq.2)then
               !write(*,*) '*** ENVELOPE IS RELAXED ***'
               ending=.true.
               goto 100
            endif
            !write(*,*) 'NEGATIVE FLUX IS ON'
            DDT = -dm(NZN)/100000.d0
            NEGFLU=2
            goto 24
         endif
         !write(*,*) 'NEGATIVE FLUX WITHOUT Z-DERIV. IS ON'
         DDT = -dm(NZN)/100000.d0
         NEGFLU=1
         goto 24
      endif
      
      !write(*,*) 'abs(DDT/dmN), PREC1', DDT, dmN, abs(DDT/dmN), PREC1

      if(abs(DDT/dmN).lt.PREC1) then
         !write(*,*) '*** ENVELOPE IS RELAXED ***'
         ending=.true.
         goto 100
      endif
      
      if (dmN_cnt >= max_dmN_cnt) then
         write(*,*) 'RELAX_ENV has reached max num allowed tries for outer dm', max_dmN_cnt
         ierr = -1
         return
         call mesa_error(__FILE__,__LINE__)
      end if

 24   TNL  = TT
      dmNL = dmN
      dmN_cnt = dmN_cnt + 1
      
      if (s% RSP_relax_adjust_inner_mass_distribution) then
      
         dmN  = dmN-DDT

         if(IOP.ge.1) then
   !        CALCULATE HAHA - ZONE MASS RATIO BELOW ANCHOR
   !        TOTAL MASS BELOW ANCHOR: SUMM, FIRST ZONE MASS: AONE
            AONE=dmN
            HAHA=1.01d0
            do I=1,100
               POM=1.d0/(NZN-NZT+1)*dlog10(1.d0-SUMM/AONE*(1.d0-HAHA))
               if(dabs(HAHA-10.d0**POM).lt.1d-10) goto 22
               HAHA=10.d0**POM         
            enddo
            write(*,*) 'NO CONVERGENCE IN RELAX_ENV ITERATION FOR H'
            stop
 22      continue
            HAHA=10.d0**POM  
   !        SET NEW MASS DISTRIBUTION
            do I=NZN,1,-1
               if(I.eq.NZN) then
   !              SPECIAL DEFINITIONS FOR THE OUTER ZONE
                  M(I)=MSTAR
                  dm(I)=dmN*FSUB
                  dm_bar(I)=dm(I)*0.5d0
               else
                  if(I.ge.NZN-NZT+1) dm(I)=dmN !DEFINITION DOWN TO ANCHOR
                  if(I.lt.NZN-NZT+1.and.ok_to_adjust_mass)  &
                                     dm(I)=AONE*pow(HAHA,(NZN-NZT+1)-I)
                  dm_bar(I)=(dm(I)+dm(I+1))*0.5d0
                  M(I)=M(I+1)-dm(I+1)             
               endif
            enddo
         endif
         
      end if

 100  continue

      do 102 II=1,8000

!     LOOP 1 .. EOS
      !$OMP PARALLEL DO PRIVATE(I,T1,POM,EDFAC,op_err) SCHEDULE(dynamic,2)
      do 1 I=2,NZN
         T1=P43/dm(I)
         Vol(I)=T1*(R(I)**3-R(I-1)**3)
         DVRM(I) = -3.d0*T1*R(I-1)**2
         DVR(I)  =  3.d0*T1*R(I)**2

         if(I.eq.NZN.and.ok_to_adjust_Tsurf)then
            if(RSP_eddi)then
               EDFAC=1.d0/2.d0
            else
               EDFAC=sqrt(3.d0)/4.d0
            end if
            POM=(EDFAC*L0/(4.d0*PI*SIG*R(NZN)**2))**0.25d0
            T(NZN)=POM
         endif

         call mesa_eos_kap(s,0, &
           T(I),Vol(I),P(I),DPV(I),dP_dT_00(I),E(I), &
           DEV(I),dE_dT_00(I),CPS(I),CPV(I),dCp_dT_00(I), &
           QQS(I),QQV(I),dQQ_dT_00(I),K(I),dK_dV_00(I),dK_dT_00(I),op_err)
         if (op_err /= 0) ierr = op_err
         if (ierr /= 0) cycle
         if (is_bad(dQQ_dT_00(I))) then
            write(*,*) 'dQQ_dT_00(I)', i, dQQ_dT_00(I)
            write(*,*) 'QQS(I)', i, QQS(I)
            write(*,*) 'T(I)', i, T(I)
            write(*,*) 'Vol(I)', i, Vol(I)
            stop
         end if
         
         dP_dr_00(I)  = DPV(I)*DVR(I)
         dP_dr_in(I) = DPV(I)*DVRM(I)
         dCp_dr_00(I)  = CPV(I)*DVR(I)
         dCp_dr_in(I) = CPV(I)*DVRM(I)
         dQQ_dr_00(I)  = QQV(I)*DVR(I)
         dQQ_dr_in(I) = QQV(I)*DVRM(I) 
 1    continue
      !$OMP END PARALLEL DO
      if (ierr /= 0) return

!     FIRST ZONE, CALCULATION OF R0 AND EOS
      P(1)=P(2)+G*M(1)*dm_bar(1)/(P4*R(1)**4)
!     WITH GIVEN T AND P ITERATE FOR V AND CALCULATE ITS DERIVATIVES
      call EOP(s,0,T(1),P(1),Vol(1),E(1),CPS(1),QQS(1),SVEL,K(1),ierr)
      if (ierr /= 0) return
      s% R_center=pow(R(1)**3-3.d0*Vol(1)*dm(1)/P4,1.d0/3.d0)
      T1=P43/dm(1)
      DVRM(1) = -3.d0*T1*s% R_center**2
      DVR(1)  =  3.d0*T1*R(1)**2
!     CALCULATE EOS FOR THE FIRST ZONE TO OBTAIN DERIVATIVES
      call mesa_eos_kap(s,0, &
           T(1),Vol(1),P(1),DPV(1),dP_dT_00(1),E(1), &
           DEV(1),dE_dT_00(1),CPS(1),CPV(1),dCp_dT_00(1), &
           QQS(1),QQV(1),dQQ_dT_00(1),K(1),dK_dV_00(1),dK_dT_00(1),ierr)
      if (ierr /= 0) return
       !write(*,*) '1st zone T V P',T(1),Vol(1),P(1)

      dP_dr_00(1)  = DPV(1)*DVR(1)
      dP_dr_in(1) = DPV(1)*DVRM(1)
      dCp_dr_00(1)  = CPV(1)*DVR(1)
      dCp_dr_in(1) = CPV(1)*DVRM(1)
      dQQ_dr_00(1)  = QQV(1)*DVR(1)
      dQQ_dr_in(1) = QQV(1)*DVRM(1)

!     SET E_T (NOW =w) BELOW AND ABOVE BOUNDARIES
!     w = EFL02 BELOW IBOTOM, INSREAD OF =0 BETTER, BECAUSE OF
!     TURBULENT FLUX, NEAR THE IBOTOM
      Et(NZN) = 0.d0
      do I=1,IBOTOM
         Et(I) = 0.d0
      enddo

      do I=1,NZN-1
         POM=  (R(I)**2)/(2.d0*G*M(I))
         Hp_face(I)=POM*(P(I)*Vol(I)+P(I+1)*Vol(I+1))
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

!     PROTECT FROM 'ALWAYS ZERO' SOLUTION WHEN (Y>0)
!     (SEEMS UNNECESSARY)
      do I=IBOTOM+1,NZN-1
!         if((Y_face(I)+Y_face(I-1)).gt.0.d0.and.Et(I).eq.0.d0) 
!     x                    Et(I)=1.d+6!1.d-6
         if(ALFA.eq.0.d0) Et(I)=0.d0
      enddo

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

      do I=IBOTOM+1,NZN-1
!        SOURCE TERM
         POM=0.5d0*(PII(I)/Hp_face(I)+PII(I-1)/Hp_face(I-1))
         POM2=T(I)*P(I)*QQS(I)/CPS(I)
         POM3=sqrt(Et(I))
         SOURCE(I)=POM*POM2*POM3
         TEM1=POM2*POM3*0.5d0
         TEMI=-PII(I)/Hp_face(I)**2
         TEMM=-PII(I-1)/Hp_face(I-1)**2
         
         ! X -> w
         ! Y -> T
         ! Z -> R

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

      do I=1,NZN
!        CONVECTIVE LUMINOSITY
         if(I.lt.IBOTOM.or.I.ge.NZN.or.alfa.eq.0d0) then
            Lc(I)=0.d0
            DLCX0(I)=0.d0
            DLCXM(I)=0.d0
            DLCXP(I)=0.d0
            DLCY0(I)=0.d0
            DLCYP(I)=0.d0
            DLCZ0(I)=0.d0
            DLCZP(I)=0.d0
         else
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

            if(NEGFLU.eq.0.and.PII(I).lt.0.d0)then
               Lc(I)=0.d0
               DLCX0(I)=0.d0
               DLCXM(I)=0.d0
               DLCXP(I)=0.d0
               DLCY0(I)=0.d0
               DLCYP(I)=0.d0
               DLCZ0(I)=0.d0
               DLCZP(I)=0.d0
            endif
            if(NEGFLU.eq.1.and.PII(I).lt.0.d0)then
               if(II.gt.300)then
                  DLCZ0(I)=0.d0
                  DLCZP(I)=0.d0
               endif
            endif
         end if
 
!        TURBULENT LUMINOSITY
         if(I.lt.IBOTOM.or.I.ge.NZN.or. &
            ALFAT.eq.0.d0.or.ALFA.eq.0.d0)then
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
            if (is_bad(Lt(I))) then
               write(*,*) 'Lt(I)', I, Lt(I)
               write(*,*) 'POM', I, POM
               write(*,*) 'POM2', I, POM2
               write(*,*) 'POM3', I, POM3
               stop
            end if
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

      do I=1,IBOTOM
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
      enddo
      do I=1,IBOTOM-1
         DLCX0(I)= 0.d0
         DLCXM(I)= 0.d0
         DLCXP(I)= 0.d0
         DLCY0(I)= 0.d0
         DLCYP(I)= 0.d0
         DLCZ0(I)= 0.d0
         DLCZP(I)= 0.d0
         DLTX0(I)= 0.d0
         DLTXM(I)= 0.d0
         DLTXP(I)= 0.d0
         DLTY0(I)= 0.d0
         DLTYP(I)= 0.d0
         DLTZ0(I)= 0.d0
         DLTZP(I)= 0.d0
      enddo
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
      DLTX0(NZN)= 0.d0
      DLTXM(NZN)= 0.d0
      DLTXP(NZN)= 0.d0
      DLTY0(NZN)= 0.d0
      DLTYP(NZN)= 0.d0
      DLTZ0(NZN)= 0.d0
      DLTZP(NZN)= 0.d0
      dC_dr_00(NZN) = 0.d0
      dC_dr_out(NZN) = 0.d0
      dC_dr_in(NZN) = 0.d0
      dC_dr_in2(NZN)= 0.d0
      dC_dT_in(NZN) = 0.d0
      dC_dT_00(NZN) = 0.d0
      dC_dT_out(NZN) = 0.d0
      dC_dw_00(NZN) = 0.d0
      
      if(ALFA.eq.0.d0) then     !RADIATIVE CASE
         SOURCE(I) = 0.d0
         DAMP(I)   = 0.d0
         DAMPR(I)  = 0.d0
         dC_dr_00(I)   = 0.d0
         dC_dr_out(I)   = 0.d0
         dC_dr_in(I)   = 0.d0
         dC_dr_in2(I)  = 0.d0
         dC_dT_in(I)   = 0.d0
         dC_dT_00(I)   = 0.d0
         dC_dT_out(I)   = 0.d0
         dC_dw_00(I)   = 0.d0
      endif


!     INITIALIZE HD(11,3*NZN)
      do I=1,3*NZN       
         do J=1,11
            HD(J,I)=0.d0
         enddo
      enddo

!     LOOP 2 .. LUM PLUGS
      DLR  =  0.d0
      DLRP =  0.d0
      DLRM =  0.d0
      DLT  =  0.d0
      DLTP =  0.d0
      DLR  =  0.d0 !!-1.d0
      do 5 I=1,NZN
         IR = 3+3*(I-1)
         IC = 2+3*(I-1)
         IW = 1+3*(I-1)
!        SET LUM(I-1)
         DLMR  = DLR
         DLMRP = DLRP
         DLMRM = DLRM
         DLMT  = DLT
         DLMTP = DLTP
         if(I.eq.NZN) goto 6
!        Lr(I)=Eq. A.4, Stellingwerf 1975, Appendix A
!        CALC LUM(I)
         W_00=T(I)**4
         W_out=T(I+1)**4
         BW=dlog(W_out/W_00)
         BK=dlog(K(I+1)/K(I))
         T1=-CL*R(I)**4/dm_bar(I)
         T2=(W_out/K(I+1)-W_00/K(I))/(1.d0-BK/BW)
         Lr(I)=T1*T2
         T3=T1/(BW-BK)
         DLK=  (T3/K(I))  *(W_00*BW/K(I)  -T2) !dL(i)/dK(i)
         DLKP=-(T3/K(I+1))*(W_out*BW/K(I+1)-T2) !dL(i)/dK(i+1)
         DLRP= DLKP*dK_dV_00(I+1)*DVR(I+1)
         DLRM= DLK *dK_dV_00(I)  *DVRM(I)         
         DLR= 4.d0*T1*T2/R(I)+DLK*dK_dV_00(I)*DVR(I)+DLKP*dK_dV_00(I+1)*DVRM(I+1)
         DLTP=4.d0*(T3/T(I+1))*(W_out*BW/K(I+1)-T2*BK/BW)+DLKP*dK_dT_00(I+1)
         DLT=-4.d0*(T3/T(I))*(W_00*BW/K(I)-T2*BK/BW)+DLK*dK_dT_00(I)
         goto 7
!        OUTER LUM BOUNDARY CONDITION
 6       continue
         Lr(I)=L0
         DLT  = 4.d0*L0/T(I)  !L=4piR^2sigT^4
         DLR  = 2.d0*L0/R(I)  !L=4piR^2sigT^4
         DLRM = 0.d0
         DLRP = 0.d0
         DLTP = 0.d0
 7       continue
!        CALC ENERGY EQUATION(I)
         HR(IW)    = -(Lr(I)/L0+Lc(I)/L0+Lt(I)/L0-1.d0)
         !write(*,*) 'energy', I, HR(IW)
         if (is_bad(HR(IW))) then
            write(*,*) 'L0', I, L0
            write(*,*) 'Lt(I)', I, Lt(I)
            write(*,*) 'Lc(I)', I, Lc(I)
            write(*,*) 'Lr(I)', I, Lr(I)
            stop
         end if
         HD(6,IW)  = (DLT +DLCY0(I)+DLTY0(I))/L0     !Y(i)
         HD(9,IW)  = (DLTP+DLCYP(I)+DLTYP(I))/L0     !Y(i+1)
         HD(5,IW)  = (DLRM+DLCXM(I)+DLTXM(I))/L0     !X(i-1)
         HD(8,IW)  = (DLR +DLCX0(I)+DLTX0(I))/L0     !X(i)
         HD(11,IW) = (DLRP+DLCXP(I)+DLTXP(I))/L0     !X(i+1)
         HD(7,IW)  = (     DLCZ0(I)+DLTZ0(I))/L0     !Z(i)
         HD(10,IW) = (     DLCZP(I)+DLTZP(I))/L0     !Z(i+1) 
!        CALC MOMEMTUM EQUATION(I)
         T1=-P4*(R(I)**2)/dm_bar(I)
         if(I.eq.NZN) then
            HR(IR) = -T1*(Psurf-P(I)-PTURB(I))+G*M(I)/R(I)**2
         else
            HR(IR)=-T1*(P(I+1)-P(I)+PTURB(I+1)-PTURB(I))+G*M(I)/R(I)**2
         end if
         !write(*,*) 'momentum', I, HR(IR)
         HD(3,IR) = +T1*(-dP_dr_in(I)-dPtrb_dr_in(I))  !X(i-1)
         HD(4,IR) = +T1*(-dP_dT_00(I))            !Y(i)

         HD(2,IR) = 0.d0                     !Z(i-1)
         HD(5,IR) = +T1*(-dPtrb_dw_00(I))          !Z(i)
         HD(8,IR) = +T1*(dPtrb_dw_00(I+1))         !Z(i+1)

         if(I.eq.NZN) goto 111
         HD(6,IR) = +T1*(dP_dr_in(I+1)-dP_dr_00(I)+dPtrb_dr_in(I+1)-dPtrb_dr_00(I))+ &
                    4.d0*G*M(I)/(R(I)**3)    !X(i)
         HD(7,IR) = +T1*(dP_dT_00(I+1))           !Y(i+1)
         HD(9,IR) = +T1*(dP_dr_00(I+1)+dPtrb_dr_00(I+1))!X(i+1)
         goto 112
 111     HD(7,IR)=  0.d0
         HD(9,IR)=  0.d0
         HD(6,IR)= +T1*(-dP_dr_00(I))+4.d0*G*M(I)/(R(I)**3)
 112     continue
 
!        CALC CONVECTIVE ENERGY EQUATION
         if(I.le.IBOTOM.or.I.eq.NZN.or.ALFA.eq.0.d0)then
            do J=1,11
               HD(J,IC)=0.d0
            enddo
            HD(6,IC)=1.d0
            HR(IC)=0.d0
         else
            T1=-1.d0/dm(I)
            HR(IC)= -(T1*(Lt(I)-Lt(I-1))+SOURCE(I)-DAMP(I)-DAMPR(I))
            !write(*,*) 'conv energy', I, HR(IC)
            if (is_bad(HR(IC))) then
               write(*,*) 'Lt(I', I, Lt(I)
               write(*,*) 'Lt(I-1)', I-1, Lt(I-1)
               write(*,*) 'SOURCE(I)', I, SOURCE(I)
               write(*,*) 'DAMP(I)', I, DAMP(I)
               write(*,*) 'DAMPR(I)', I, DAMPR(I)
               stop
            end if
            HD(3,IC)  =          +T1*(        -DLTZ0(I-1))  !Z(i-1)
            HD(6,IC)  = dC_dw_00(I) +T1*(DLTZ0(I)-DLTZP(I-1))  !Z(i)
            HD(9,IC)  =          +T1*(DLTZP(I))             !Z(i+1)
            HD(1,IC)  = dC_dr_in2(I)+T1*(        -DLTXM(I-1))  !X(i-2)
            HD(4,IC)  = dC_dr_in(I) +T1*(DLTXM(I)-DLTX0(I-1))  !X(i-1)
            HD(7,IC)  = dC_dr_00(I) +T1*(DLTX0(I)-DLTXP(I-1))  !X(i)
            HD(10,IC) = dC_dr_out(I) +T1*(DLTXP(I))             !X(i+1)
            HD(2,IC)  = dC_dT_in(I) +T1*(        -DLTY0(I-1))  !Y(i-1)
            HD(5,IC)  = dC_dT_00(I) +T1*(DLTY0(I)-DLTYP(I-1))  !Y(i)
            HD(8,IC)  = dC_dT_out(I) +T1*(DLTYP(I))             !Y(i+1)
         end if

!        SET ANCHOR T(NZN)=CONST
         if(I.eq.NZN) then
            do j=1,11
               HD(J,IW) = 0.d0 
            enddo
!           SET dT(NZN)=0
            HD(6,IW) = 1.d0 
            HR(IW)   = 0.d0 
!           SET DERIVATIVES d(*)/dT(NZN)=0 IN ZONE/ITERFACE NZN
            HD(4,IR) = 0.d0 
            HD(5,IC) = 0.d0 
         endif   
         if(I.eq.NZN-1)then
!           SET DERIVATIVES d(*)/dT(NZN)=0 IN ZONE/ITERFACE NZN-1
            HD(9,IW)  = 0.d0
            HD(7,IR)  = 0.d0
            HD(8,IC)  = 0.d0
         endif

  5   continue

      if(ending) goto 889

      do J=1,5     !translate hd into band storage of LAPACK/LINPACK
         do I=1,3*NZN+1
            ABB(J,I)=0.0d0
         enddo
      enddo      
      do J=1,5      
         do I=1,3*NZN+1-J
            ABB(11-J,I+J)=HD(6+J,I) !upper diagonals
            ABB(11+J,I)=HD(6-J,I+J) !lower diagonals
         enddo  
      enddo
      do I=1,3*NZN+1
         ABB(11,I)=HD(6,I)
      enddo  

!-    LAPACK
      call DGBTRF(3*NZN,3*NZN,5,5,ABB,LD_ABB,IPVT,INFO)
      if(INFO.ne.0) then
         write(*,*) 'hyd: LAPACK/dgbtrf problem no., iter.',INFO,II
         open(61,file='x.dat',status='unknown')
         write(61,'(F6.3,tr2,f8.2,tr2,f7.2,tr2,d9.3)')EMR,ELR,TE-0.5d0, &
              gam
         close(61) 
         open(61,file='logic.dat',status='unknown')
         write(61,'(I1)') 3
         close(61)
         stop
      endif
      call DGBTRS('n',3*NZN,5,5,1,ABB,LD_ABB,IPVT,HR,3*NZN,INFO)
      if(INFO.ne.0) then
         write(*,*) 'hyd: LAPACK/dgbtrs problem no., iter.',INFO,II
         stop
      endif            
!-
      do I=1,3*NZN,1
         DX(I) = HR(I)
      enddo

      EZH = 1.d0
      XXR = 0.d0
      XXT = 0.d0
      XXC = 0.d0
         XX_max = 0
         XX_max_val = 0
         XX_max_dx = 0
         i_XX_max = 0
         var_XX_max = 0
      do 14 I=2,NZN
         IR = 3+3*(I-1)
         IC = IR-1
         IW = IR-2
         if(Et(I).gt.(1.d-20)*EFL02) XXC=dabs(DX(IC)/Et(I))
         XXR=((DX(IR-3)-DX(IR))/(R(I)-R(I-1)))/DXH
         XXT=dabs(DX(IW)/T(I))/DXH
         
         if (XXC > XX_max) then
            XX_max = XXC
            XX_max_val = Et(I)
            XX_max_dx = DX(IC)
            i_XX_max = i
            var_XX_max = 1
         end if
         if (XXR > XX_max) then
            XX_max = XXR
            XX_max_val = R(I)
            XX_max_dx = DX(IR)
            i_XX_max = i
            var_XX_max = 2
         end if
         if (XXT > XX_max) then
            XX_max = XXT
            XX_max_val = T(I)
            XX_max_dx = DX(IW)
            i_XX_max = i
            var_XX_max = 3
         end if
            
         EZH=1.d0/dmax1(1.d0/EZH,XXR,XXT,XXC)
   14 continue

      DXXR = 0.d0
      DXXT = 0.d0
      DXXC = 0.d0
      ITROUBT = 0
      ITROUBC = 0
!     It seems that for BL Her models fixed undercorrection factor works best
      do 103 I=1,NZN
         IW=1+3*(I-1)
         IR=IW+2
         IC=IW+1  
         T(I)      = T(I)     +EZH*DX(IW)/2d0
         R(I)      = R(I)     +EZH*DX(IR)/2d0
         Et(I) = Et(I)+EZH*DX(IC)/2d0
         if(Et(I).le.0.d0) then
            Et(I)=EFL02*1d-4*rand(s)
         end if
         DXKT=DXXT
         DXKC=DXXC
         DXXR=dmax1(DXXR,dabs(DX(IR)/R(I)))
         DXXT=dmax1(DXXT,dabs(DX(IW)/T(I)))
         if(Et(I).gt.(1.d-20)*EFL02)  &
            DXXC=dmax1(DXXC,dabs(DX(IC)/Et(I)))
         if(DXXC.gt.DXKC) ITROUBC=I
         if(DXXT.gt.DXKT) ITROUBT=I
 103  continue
      !write(*,*) 'apply corr', II,dabs(DXXC),dabs(DXXR),dabs(DXXT),ITROUBC, &
      !      Et(ITROUBC),EZH
      !call mesa_error(__FILE__,__LINE__,'debug')
      !write(*,*) 'DXXC/PRECR', ITROUBC, DXXC/PRECR, DXXC
      if(dabs(DXXC).lt.PRECR.and.dabs(DXXR).lt.PRECR.and. &
         dabs(DXXT).lt.PRECR) then
         IZIP=1
         IOP=IOP+1
         goto 999
      endif

 102  continue
 
      write(*,*)
      write(*,*) ' NO CONVERGENCE IN RELAX_ENV, ITERATION: ',II
      write(*,*) ' try increasing RSP_relax_dm_tolerance', s% RSP_relax_dm_tolerance
      write(*,*)
      ierr = -1
      return
      stop

 889  continue


      !MAXFLUXD=-2d10
      !do I=1,NZN
      !   if(abs(HR(1+3*(I-1))).gt.MAXFLUXD) &
      !      MAXFLUXD=abs(HR(1+3*(I-1)))
      !enddo
      !write(*,*) 'MAXFLUXD= ',MAXFLUXD

      POM=MSTAR-M(1)+dm(1)
      if (s% RSP_trace_RSP_build_model) then
         write(*,*) ' inner temper. rel. diff',(T(1)-TINNER)/TINNER
         write(*,*) '  inner radius rel. diff',(R(1)-RINNER)/RINNER
         write(*,*) '  outer radius rel. diff',(R(NZN)-ROUTER)/ROUTER
         write(*,*) ' envelope mass rel. diff',(POM-MENVEL)/MENVEL
      end if
 222  format(A,d14.8)
      if (.false.) then
         MAXR=-1d0
         MAXW=-1d0
         MAXC=-1d0
         write(*,*) 'static structure check'
         !open(15,file='ss_lin.dat',status='unknown')
         do J=1,NZN
            IR=4+3*(J-1)
            IW=2+3*(J-1)            
            PB=P(J)+PTURB(J)
            if (J == NZN) then
               PPB = 0
            else
               PPB=P(J+1)+PTURB(J+1)
            end if
            T1=P4*(R(J)**2)
            T2=(PPB-PB)/dm_bar(J)
            T3=G*M(J)/(R(J)**2)
            if (dabs((T3+T1*T2)/T3).gt.MAXR)then
               MAXR=dabs((T3+T1*T2)/T3)
               IMAXR=J
            endif
            ELB=Lc(J)+Lr(J)+Lt(J)
            if (J == 1) then
               ELMB=0
            else
               ELMB=Lc(J-1)+Lr(J-1)+Lt(J-1)
            end if
            if (dabs((ELB-L0)/L0).gt.MAXW)then
               MAXW=dabs((ELB-L0)/L0)
               IMAXW=J
            endif
            if (J > IBOTOM .and. J < NZN) then
               ELB=Lt(J)
               ELMB=Lt(J-1) 
               T1=(ELB-ELMB)/dm(J)            
               T2=SOURCE(J)-DAMP(J)-DAMPR(J)
               if ((T2.ne.0.d0).and.(T1.ne.0.d0))then
                  if (dabs(T1-T2)/max(T1,T2).gt.MAXC)then
                     MAXC=dabs(T1-T2)/max(T1,T2)
                     IMAXC=J
                  endif
               endif
            end if
         enddo
         if (MAXR/=-1d0)write(*,*) 'MAX DIFFERENCE R: ',MAXR,' ZONE: ',IMAXR
         if (MAXW/=-1d0)write(*,*) 'MAX DIFFERENCE W: ',MAXW,' ZONE: ',IMAXW
         if (MAXC/=-1d0)write(*,*) 'MAX DIFFERENCE C: ',MAXC,' ZONE: ',IMAXC
      end if
      return 
      end subroutine RELAX_ENV

            
      end module rsp_relax_env
