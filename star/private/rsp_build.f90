! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton, Radek Smolec & The MESA Team
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

      module rsp_build
      use star_def, only: star_info
      use utils_lib, only: is_bad
      use const_def, only: dp, crad
      use eos_lib, only: Radiation_Pressure
      use rsp_def
      use rsp_eval_eos_and_kap, only: X, Y, Z
      use rsp_lina, only: mesa_eos_kap, do_LINA
      use rsp_relax_env, only: EOP, RELAX_ENV
      
      implicit none
      
      private
      public :: do_rsp_build

      real(dp) :: PREC,FSUB,TIN,CFIDDLE,ALF, &
         HHFAC,DdmFAC,SVEL,EFL02,EMR,ELR, &
         E_0,E_1,T_0,T_1,V_0,V_1,P_0,P_1,QQ_0,QQ_1, &
         CP_0,CP_1,OP_0,OP_1,R_1,M_0,dm_bar_0
      real(dp), dimension(15) :: PERS,ETO
      real(dp), pointer, dimension(:) :: &
         M, DM, DM_BAR, R, Vol, T, w, Et, E, P, Lr, Lc, Hp_face, Y_face, K, CPS, QQS
      logical, parameter :: RSP_eddi = .true. ! use Eddington approx at surface
      
      contains
      
      
      subroutine do_rsp_build(s,ierr)
      !use rsp_create_env, only: do_rsp_create_env
      type (star_info), pointer :: s    
      integer, intent(out) :: ierr
      integer :: NZT    ! warst. z joniz./il. warstw 
      real(dp) :: TH0                ! temp. warstwy jonizacji              
      real(dp) :: Mass,L
      real(dp) :: H,dmN
      integer :: NMODES            ! ilosc modow rozwazanych (N=NMODES)
      integer :: NDIM1,NDIM2       ! maks. ilosc modow/warstw 
      real(dp) :: VEL0(15)
      character (len=10) PARA
      character (len=30) HEAD
      character (len=72) HEAD2      
      integer :: I,J,kk,NSEQ,NIT
      real(dp) ::  GEFF,MBOL,DU,TET
      integer :: IARG
      character (len=8) II1,II2,II3,II4,II5
      character (len=15) PROGNAME
      character (len=250) FILENAME
      integer :: IFROZEN,IRELAX,ICASTOR
      
      integer :: IO,II,IX,iter
      real(dp) :: SS,AA,BB,XX
      real(dp), allocatable :: TA(:), VEL(:,:), TEMP(:)
      real(dp) :: TAUTEFF,TAUATTEFF
      logical RELAX
      complex(8) xC
      real(dp) :: X1,X2,Y1,Y2,amix1,amix2
      integer :: ISTAT
      type (star_info), target :: copy_info
      type (star_info), pointer :: c, prv
      
      ierr = 0
   
      if (.not. associated(M)) then
         allocate( &
            M(NZN+1), DM(NZN+1), DM_BAR(NZN+1), T(NZN+1), E(NZN+1), w(NZN+1), &
            Et(NZN+1), P(NZN+1), R(NZN+1), Vol(NZN+1), &
            Lr(NZN+1), Lc(NZN+1), Hp_face(NZN+1), Y_face(NZN+1), K(NZN+1), &
            CPS(NZN+1), QQS(NZN+1))
      end if
         
      allocate(TA(NZN+1),VEL(NZN+1,15),TEMP(NZN))

      TH0 = s% RSP_T_anchor
      TIN = s% RSP_T_inner
      NZT = s% RSP_nz_outer
      FSUB = s% RSP_dq_1_factor         
      ddmfac = 1d0 ! s% RSP_ddmfac
      hhfac = 1.02d0 ! s% RSP_hhfac
      nmodes = s% RSP_nmodes
      RELAX = s% RSP_relax_initial_model
      EFL02 = EFL0*EFL0

!     START MAIN CYCLE
   
!     INITIALIZE ARRAYS
      NDIM1=15 !max. number of modes
      NDIM2=NZN+1
      do J=1,NDIM1
         do I=1,NDIM2
            VEL(I,J)=0.d0
         enddo
      enddo
      do I=1,NDIM1
         PERS(I)= 0.d0
         VEL0(I)= 0.d0
         ETO(I) = 0.d0
         !SGR(I) = 0.d0
      enddo
      do I=1,NDIM2
         R(I) = 0.d0
         P(I) = 0.d0
         Vol(I) = 0.d0
         E(I) = 0.d0
         T(I) = 0.d0
         K(I)= 0.d0
         M(I)= 0.d0
         dm(I)= 0.d0
         dm_bar(I)= 0.d0
         w(I)=0.d0
         Hp_face(I)    = 0.d0
         Y_face(I)   = 0.d0
         CPS(I)    = 0.d0
         QQS(I)    = 0.d0
         Lc(I) = 0.d0
         Lr(I) = 0.d0
         w(I) = 0.d0
      enddo

!     SET STANDARD PARAMETERS
!     NSEQ = SEQUENCE NUMBER (DUMMY)
      NSEQ    = 1

      CFIDDLE = 0.02d0 !0.02
      ALF     = 1.0d-6
!     PRECISIONS
      PREC    = 1d-10
      
      EMR = s% RSP_mass
      ELR = s% RSP_L
      TE = s% RSP_Teff
      Mass=EMR*SUNM
      L=ELR*SUNL       
      
      if (s% RSP_trace_RSP_build_model) then
         write(*,*) '*** build initial model ***'
         write(*,'(a9,f15.5)') 'M/Msun', EMR
         write(*,'(a9,f15.5)') 'L/Lsun', ELR
         write(*,'(a9,f15.5)') '  Teff', TE
      end if
            
      call STAH(s,Mass,L,TE,H,dmN,TH0,NZT,NZN,ierr)
      if (ierr /= 0) return

!     INITIAL GUESS FOR w=E_T (NOW DEFINED AT THE ZONE)
      TEMP(1)=0.d0
      do I=2,NZN
         !TEMP(I)=0.5d0*(w(I)**2+w(I-1)**2)
         !TEMP(I)=sqrt(w(I)**2*w(I-1)**2)
         TEMP(I)=(w(I)**2*dm(I)+w(I-1)**2*dm(I-1))/ &
                 (dm(I)+dm(I-1))
      enddo
      do I=1,NZN
         w(I)=TEMP(I)
!        LINE BELOW MUST BE PRESENT IF MAIN VARIABLE IS E_T
!        (DERIVATIVES ~1/sqrt(E_T) AND WITH w=0, LAPACK
!         PROBLEM ARISES; NO SUCH PROBLEM WHEN E_T^2 IS USED
!         AND LINE BELOW IS NOT NECESSARY THEN)
         if(w(I).le.EFL02) w(I)=EFL02
         !write(*,*) NZN-I+1, sqrt(w(i))
      enddo
      
      ! NOTE: w(I) now holds Et = w**2
      ! watch out

      if(RELAX) then
         call RELAX_ENV(s, L, TH0, TE, NZT, NZN, &
            M, DM, DM_BAR, R, Vol, T, w, ierr)
         if (ierr /= 0) return
      end if
      
      ! optical depth.  just for output.
      TA(NZN) = s% tau_factor*s% tau_base
      SS=TA(NZN)
      do I=NZN,2,-1
         SS=SS+K(I)*(R(I)-R(I-1))/Vol(I)
         TA(I-1)=SS
      enddo
      II=0; IO=0
      do I=NZN,1,-1
         if(TA(I).gt.2.d0/3.d0)then
            II=I
            IO=I+1
            goto 77
         endif
      enddo
 77   continue
      if (IO == 0) then
         write(*,*) 'failed to find photosphere'
         ierr = -1
         return
      end if
      AA=(T(IO)-T(II))/(TA(IO)-TA(II))
      BB=T(IO)-AA*TA(IO)
      TAUTEFF=AA*2.d0/3.d0+BB

      do I=NZN,1,-1
         if(T(I).gt.TE)then
            II=I
            IO=I+1
            goto 78
         endif
      enddo
 78   continue
      AA=(T(IO)-T(II))/(TA(IO)-TA(II))
      BB=T(IO)-AA*TA(IO)
      TAUATTEFF=(TE-BB)/AA
      
      call cleanup_for_LINA(s, M, DM, DM_BAR, R, Vol, T, w, P, ierr)

 1    format(1X,1P,5D26.16)

      GEFF=G*Mass/R(NZN)**2
      MBOL=-2.5d0*dlog10(ELR)+4.79d0    
      
      if(NMODES.eq.0) goto 11 ! jesli masz liczyc tylko static envelope
      
      if (.not. (s% use_RSP_new_start_scheme .or. s% use_other_RSP_linear_analysis)) then          
         if (s% RSP_trace_RSP_build_model) write(*,*) '*** linear analysis ***'
         do I=1,NZN ! LINA changes Et, so make a work copy for it
            Et(I) = w(I)
         end do
         call do_LINA(s, L, NZN, NMODES, VEL, PERS, ETO, &
            M, DM, DM_BAR, R, Vol, T, Et, Lr, ierr)
         if (ierr /= 0) then
            return
         end if
         FILENAME=trim(s% log_directory) // '/' // 'LINA_period_growth.data'
         open(15,file=trim(FILENAME),status='unknown')
         do I=1,NMODES
            s% rsp_LINA_periods(i) = PERS(I)
            s% rsp_LINA_growth_rates(i) = ETO(I)
            write(*,'(I3,2X,99e16.5)') I-1, &
               PERS(I)/86400.d0,ETO(I)
            write(15,'(I3,2X,99e16.5)') I-1, &
               PERS(I)/86400.d0,ETO(I)
         enddo
         close(15)
         s% RSP_have_set_velocities = .true.
      else      
         PERS(1:NMODES) = 0d0
         VEL0(1:NMODES) = 0d0
         do I=1,NZN
            do j=1,nmodes
               VEL(I,J) = 0d0
            end do
         enddo         
      end if

 11   continue
  
5568  format(1X,1P,5E15.6)

 444  format(F6.3,tr2,f8.2,tr2,f7.2,tr2,d9.3) 
      if (s% RSP_trace_RSP_build_model) then
         write(*,*) '*** done creating initial model ***'
         write(*,*)
      end if
      ! recall that w is actually Et = w**2 at this point
      call set_build_vars(s,M,DM,DM_BAR,R,Vol,T,w,Lr,Lc)

      IWORK=0
      TEFF = TE
      RSTA = R(NZN)
      do I=1,NZN
         kk = NZN+1-i
         s% L(kk)=L
         s% L_start(kk)=0.0d0  !SMOLEC!
         s% v(kk)=0d0
      end do
      s% L_center=L
      if(ALFA.eq.0.d0) EFL0=0.d0          
      s% rsp_period=s% RSP_default_PERIODLIN
      if (is_bad(s% rsp_period)) then
         write(*,1) 'rsp_period', s% rsp_period
         stop 'rsp_build read_model'
      end if
      amix1 = s% RSP_fraction_1st_overtone
      amix2 = s% RSP_fraction_2nd_overtone
      if((AMIX1+AMIX2).gt.1.d0) write(*,*) 'AMIX DO NOT ADD UP RIGHT' 
      if (.not. s% use_RSP_new_start_scheme) then      
         PERIODLIN=PERS(s% RSP_mode_for_setting_PERIODLIN+1)         
         s% rsp_period=PERIODLIN            
         s% v_center = 0d0
         do I=1,NZN
            s% v(NZN+1-i)=1.0d5*s% RSP_kick_vsurf_km_per_sec* &
              ((1.0d0-AMIX1-AMIX2)*VEL(I,1)+AMIX1*VEL(I,2)+AMIX2*VEL(I,3))
         enddo      
      end if
      
      end subroutine do_rsp_build


      subroutine STAH(s,MX,L,TE,H0,dmN0,TH0,NZT,NZN,ierr)
   !     INTEGRATE STATIC ENVELOPE
   !     CONVECTIVE FLUX INCLUDED (WUCHTERL & FEUCHTINGER 1998)
   !     TURBULENT PRESSERE AND OVERSHOOTING NEGLECTED 
   !     (ALPHA_P = ALPHA_T = 0)
   !     DIFFUSION APPROXIMATION FOR RADIATIVE TRANSFER
   !     HYDROGEN ZONE DEPTH (NZT, IN ZONES), AND TEMPERATURE(TH0) FIXED
   !     ARGUMENTS.. M(GM), L(CGS), TE,TH0(DEG K)

         type (star_info), pointer :: s
         real(dp), intent(in) :: MX,L,TE,TH0
         real(dp), intent(out) :: H0,dmN0
         integer, intent(in) :: NZT
         integer, intent(inout) :: NZN
         integer, intent(out) :: ierr
      
         real(dp) :: dmN,dm_0,H,Psurf,DDT
         real(dp) :: OPVV,OPTT,POM
         real(dp) :: GPF
         real(dp) :: DTN,DTLAST,RS00,RIN
         real(dp) :: F2,F1,D,HH,TT,RAMA,dmL
         real(dp) :: T4_1,RM,T4_0,WE,TNL,dmNL
         real(dp) :: FACQ,HH1,HH2
         integer :: N,N1,N2,I,ITIN,dmN_cnt,NCHANG,IG,H_cnt
         real(dp) :: HP_0,HP_1,IGR_0,IGR_1,PII,w_0
         real(dp) :: Lr_0,Lc_0,SVEL_0,HSTART,tau_sum,TH0_tol,TIN_tol, &
            dmN_too_large, dmN_too_small, H_too_large, H_too_small
         logical :: adjusting_dmN, in_photosphere, in_outer_env, &
            have_dmN_too_large, have_dmN_too_small, &
            have_H_too_large, have_H_too_small, have_T
         
         include 'formats'
         ierr = 0
         IG  = 0
         HH1 = 1.005d0!1.005
         HH2 = 1.005d0!1.005
         NCHANG = 0
         HSTART = 1.01d0 ! s% RSP_hstart
         FACQ = 1.005d0 ! 1.005
         ITIN = 1
         tau_sum = 0
         adjusting_dmN = .true.
         in_photosphere = .true.
         in_outer_env = .true.
         have_dmN_too_large = .false.
         have_dmN_too_small = .false.
         have_H_too_large = .false.
         have_H_too_small = .false.
         
         TH0_tol = s% RSP_T_anchor_tolerance
         TIN_tol = s% RSP_T_inner_tolerance
      
         call ZNVAR(s,H,dmN,L,TE,MX,ierr)
         if (ierr /= 0) return
         dmN_cnt = 1
         H_cnt = 1
         
         start_from_top_loop: do
            if (s% RSP_trace_RSP_build_model) write(*,*) 'call setup_outer_zone'
            call setup_outer_zone(ierr)
            if (ierr /= 0) return
            if (s% RSP_trace_RSP_build_model) write(*,*) 'call store_N'
            call store_N
            zone_loop: do
               R_1=pow(R_1**3-3.d0*V_0*dm_0/P4,1.d0/3.d0)
               N=N-1  
               if (s% RSP_trace_RSP_build_model) write(*,*) 'zone_loop', N, T_0, TIN
               if (N.eq.0 .or. T_0 >= TIN) then
                  if (s% RSP_trace_RSP_build_model) write(*,*) 'call next_H'
                  call next_H ! sets HH
                  if (N.eq.0 .and. abs(T(1)-TIN).lt.TIN*TIN_tol) then
                     s% M_center = M_0 - dm_0
                     s% star_mass = s% RSP_mass
                     s% mstar = s% star_mass*SUNM
                     s% xmstar = s% mstar - s% M_center
                     s% M_center = s% mstar - s% xmstar ! this is how it is set when read file
                     s% L_center = s% RSP_L*SUNL
                     s% R_center = pow(r(1)**3 - Vol(1)*dm(1)/P43, 1d0/3d0)
                     s% v_center = 0                     
                     if (s% RSP_trace_RSP_build_model) &
                        write(*,*) '   inner dm growth scale', HH
                     exit start_from_top_loop ! done
                  end if
                  if (H_cnt >= s% RSP_max_inner_scale_tries) then
                     write(*,*) 'failed to find inner dm scaling to satisify tolerance for T_inner'
                     write(*,*) 'you might try increasing RSP_T_inner_tolerance'
                     ierr = -1
                     return
                     !stop 1
                  end if
                  if (s% RSP_trace_RSP_build_model) write(*,*) 'call prepare_for_new_H'
                  call prepare_for_new_H
                  cycle zone_loop
               end if
               if (s% RSP_trace_RSP_build_model) write(*,*) 'call setup_next_zone'
               call setup_next_zone
               ! pick T to make Lr + Lc = L
               if (s% RSP_trace_RSP_build_model) write(*,*) 'call get_T'
               have_T = get_T(ierr)
               if (ierr /= 0) return
               if (.not. have_T) then
                  call failed    
                  DDT = dmN/1.d3
                  dmN = dmN-DDT
                  cycle start_from_top_loop            
               end if
               if (s% RSP_trace_RSP_build_model) write(*,*) 'call get_V'
               call get_V(ierr)
               if (ierr /= 0) return
               if((NZN-N+1.eq.NZT .or. T_0 >= TH0) &
                     .and. adjusting_dmN) then
                  if (s% RSP_trace_RSP_build_model) &
                     write(*,*) 'call next_dmN', dmN_cnt, NZN-N+1, T_0, TH0, abs(T_0-TH0), TH0_tol*TH0
                  if (dmN_cnt >= s% RSP_max_outer_dm_tries) then
                     write(*,*) 'failed to find outer dm to satisify tolerance for T_anchor'
                     write(*,*) 'you might try increasing RSP_T_anchor_tolerance'
                     ierr = -1
                     return
                     !stop 1
                  end if
                  call next_dmN
                  if(NZN-N+1.eq.NZT.and.abs(T_0-TH0).lt.TH0_tol*TH0) then
                     adjusting_dmN = .false.
                     if (s% RSP_trace_RSP_build_model) &
                        write(*,*) '           outer dm/Msun', dmN/SUNM
                  end if
                  ! one last repeat with final dmN
                  !if (s% RSP_trace_RSP_build_model) write(*,*) 'cycle start_from_top_loop', dmN_cnt, dmN/SUNM
                  cycle start_from_top_loop            
               end if
               if (s% RSP_trace_RSP_build_model) write(*,*) 'call store_N'
               call store_N
            end do zone_loop  
         end do start_from_top_loop
         
         s% R_center=R_1; H0=H; dmN0=dmN     
         if(N.ne.0) call change_NZN
      
      contains
      
      subroutine report_location_of_photosphere
         real(dp) :: tau, dtau
         integer :: i
         include 'formats'
         tau = 0d0 ! surface currently at tau=0
         do i=NZN,1,-1
            dtau = dm(i)*K(i)/(4d0*pi*r(i)**2)
            tau = tau + dtau
            if (tau >= 2d0/3d0) then
               write(*,2) 'cells from surface tau=0 to tau=2/3', &
                  NZN+1-i, tau-dtau, 2d0/3d0, tau, T(i), TE
               return
            end if
         end do
      end subroutine report_location_of_photosphere
      
      subroutine setup_outer_zone(ierr)
         use rsp_eval_eos_and_kap, only: get_surf_P_T_kap
         integer, intent(out) :: ierr
         real(dp) :: tau_surf, kap_guess, T_surf, kap_surf, Teff_atm
         include 'formats'
         dm_0=dmN*FSUB
         M_0=MX
         dm_bar_0=(dm_0/2.d0)
         if(.not.RSP_eddi) then !     EXACT GREY RELATION
            WE=TE**4
            T4_0=WE*sqrt(3.d0)/4.d0               !0.4330127018d0 
            T_0= pow(sqrt(3.d0)/4.d0,0.25d0)*TE !0.811194802d0*TE
         else !     EDDINGTON APPROXIMATION
            WE=TE**4
            T4_0=WE*0.5d0 ! T4_0=WE*1.0d0/2.d0
            T_0=pow(0.5d0, 0.25d0)*TE ! T_0= pow(1.0d0/2.d0,0.25d0)*TE
         endif      
         RM=sqrt(L/(P4*SIG*WE))
         R_1=RM
         if (s% RSP_use_atm_grey_with_kap_for_Psurf) then
            tau_surf = s% RSP_tau_surf_for_atm_grey_with_kap
            kap_guess = 1d-2
            call get_surf_P_T_kap(s, &
               MX, RM, L, tau_surf, kap_guess, &
               T_surf, Psurf, kap_surf, Teff_atm, ierr)
            !write(*,*) 'Psurf from atm, Prad_surf', Psurf, crad*T_0*T_0*T_0*T_0/3d0
            if (ierr /= 0) then
               write(*,*) 'failed in get_surf_P_T_kap'
               return
            end if
         else if (s% RSP_use_Prad_for_Psurf) then
            Psurf = crad*T_0*T_0*T_0*T_0/3d0
         else
            Psurf = 0d0
         end if
         Psurf_from_atm = Psurf
         P_0=Psurf+G*M_0*dm_bar_0/(P4*R_1**4)
         call EOP(s,-1,T_0,P_0,V_0, &
                  E_0,CP_0,QQ_0,SVEL_0,OP_0,ierr)
         if (ierr /= 0) return
         w_0=0.d0
         Lc_0=0.d0
         Lr_0=L
         N=NZN
      end subroutine setup_outer_zone
      
      subroutine get_V(ierr)
         integer, intent(out) :: ierr
         ierr = 0
         T_0=sqrt(sqrt(T4_0))
         call EOP(s,N,T_0,P_0,V_0, &
                  E_0,CP_0,QQ_0,SVEL_0,OP_0,ierr)
         if (ierr /= 0) return
         if(N.ne.NZN)then
            call CFLUX(HP_0,IGR_0,Lc_0,w_0,GPF,N)
            if(Lc_0.ge.L) then 
               write(*,*) 'trouble!',I
               stop
            endif
         endif
      end subroutine get_V

      subroutine next_dmN
         real(dp), parameter :: search_factor = 2d0
         dmN_cnt = dmN_cnt+1
         dmNL=dmN
         if (T_0 < TH0) then ! dmN is too small
            dmN_too_small = dmN
            have_dmN_too_small = .true.
            if (.not. have_dmN_too_large) then
               dmN = dmN * search_factor
               DDT = dmN - dmNL
               return
            end if
         else ! T_0 > TH0, dmN is too large
            dmN_too_large = dmN
            have_dmN_too_large = .true.
            if (.not. have_dmN_too_small) then
               dmN = dmN / search_factor
               DDT = dmN - dmNL
               return
            end if
         end if
         ! search using bounds
         ! just bisect since for dmN too large, stop short of target cell
         dmN = 0.5d0*(dmN_too_large + dmN_too_small)
         DDT = dmN - dmNL
      end subroutine next_dmN

      subroutine next_H ! same scheme as next_dmN.  bound and bisect.
         real(dp), parameter :: search_factor = 1.05d0
         real(dp) :: HH_prev
         HH_prev = HH
         H_cnt = H_cnt+1
         !write(*,*) 'next_H have_H_too_large', have_H_too_large, H_too_large
         !write(*,*) 'next_H have_H_too_small', have_H_too_small, H_too_small
         if (T_0 < TIN) then ! H is too small
            H_too_small = H
            have_H_too_small = .true.
            if (.not. have_H_too_large) then
               HH = H * search_factor
               !write(*,*) 'too small next_H HH, HH_prev', HH, HH_prev
               return
            end if
         else ! T_0 > TIN, H is too large
            H_too_large = H
            have_H_too_large = .true.
            if (.not. have_H_too_small) then
               HH = H / search_factor
               !write(*,*) 'too large next_H HH, HH_prev', HH, HH_prev
               return
            end if
         end if
         ! search using bounds.  keep it simple. 
         ! just bisect since for H too large, stop short of target cell
         HH = 0.5d0*(H_too_large + H_too_small)
         !write(*,*) 'next_H HH, HH_prev', HH, HH_prev
         !if (abs(HH - HH_prev) < 1d-6*HH) stop 'next_H'
      end subroutine next_H
      
      subroutine store_N
         real(dp) :: dtau
         R(N)  = R_1                               
         P(N)  = P_0
         Vol(N)  = V_0
         E(N)  = E_0
         T(N)  = T_0
         K(N) = OP_0
         M(N) = M_0
         dm(N) = dm_0
         dm_bar(N) = dm_bar_0
         Hp_face(N) = HP_0
         Y_face(N)= IGR_0
         CPS(N) = CP_0
         QQS(N) = QQ_0
         Lc(N) = Lc_0
         Lr(N) = Lr_0
         w(N) = w_0
         dtau = dm(N)*K(N)/(P4*R(N)**2)
         if (in_photosphere .and. T(N) >= TE) then
            if (s% RSP_testing) &
               write(*,*) 'nz phot, tau_sum, T', &
                  NZN-N, tau_sum, tau_sum+dtau, T(N+1), TE, T(N)
            in_photosphere = .false.
         end if
         tau_sum = tau_sum + dtau
      end subroutine store_N
      
      subroutine setup_next_zone      
         P_1  = P_0
         V_1  = V_0
         OP_1 = OP_0
         T4_1  = T4_0
         M_0   = M_0-dm_0
         dmL = dm_0   !dmL IS dm IN A PREVIOUS ZONE
         HP_1  = HP_0
         IGR_1 = IGR_0
         CP_1  = CP_0
         QQ_1  = QQ_0
         T_1 = sqrt(sqrt(T4_1))
         E_1 = E_0
         !     RESET dm FOR THE INNNER ZONES
         if(N.eq.NZN-1) then
            dm_0=dmN
         else if (T_1 > TH0) then
            if (in_outer_env) then
               in_outer_env = .false.
               if (s% RSP_testing) write(*,*) 'nz outer', NZN-N, T_1, TH0
            end if
            dm_0=dm_0*H
         end if         
         dm_bar_0=(dm_0+dmL)/2.d0
         P_0=P_1+G*M_0*dm_bar_0/(P4*R_1**4)         
      end subroutine setup_next_zone
      
      real(dp) function eval_T_residual(ierr)
         integer, intent(out) :: ierr
         call EOP(s,0, &
            T_0,P_0,V_0,E_0,CP_0,QQ_0,SVEL_0,OP_0,ierr)     
         if (ierr /= 0) return 
         call CFLUX(HP_0,IGR_0,Lc_0,w_0,GPF,N)
         TT=4.d0*SIG*P4**2*R_1**4/(3.d0*dm_bar_0*L)
         T4_0 = T_0**4
         Lr_0=TT*(T4_0/OP_0-T4_1/OP_1)/ &
               (1.d0-dlog(OP_0/OP_1)/dlog(T4_0/T4_1))*L
         eval_T_residual = (Lr_0 + Lc_0)/L - 1d0               
      end function eval_T_residual
      
      real(dp) function get_T_residual(lnT, dfdx, lrpar, rpar, lipar, ipar, ierr)
         ! returns with ierr = 0 if was able to evaluate f and df/dx at x
         ! if df/dx not available, it is okay to set it to 0
         use const_def, only: dp
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(in) :: lnT
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr         
         ierr = 0
         dfdx = 0
         T_0 = exp(lnT)
         get_T_residual = eval_T_residual(ierr)
      end function get_T_residual
      
      logical function get_T_estimate()
         use num_lib, only: safe_root_with_brackets
         real(dp) :: Tmax, epsx, epsy, residual, lnT, &
            lnT_min, lnT_max, resid_T_min, resid_T_max, dfdx
         integer :: i, n, ierr
         integer, parameter :: lrpar=1, lipar=1, imax=50
         real(dp), target :: rpar_target(lrpar)
         integer, target :: ipar_target(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         include 'formats'
         Tmax = 0.99d0*(3d0*P_0/crad)**0.25d0 ! Prad must be < P_0
         lnT_min = log(T_1)
         lnT_max = log(Tmax)
         resid_T_min = get_T_residual(lnT_min, &
            dfdx, lrpar, rpar, lipar, ipar, ierr)
         resid_T_max = get_T_residual(lnT_max, &
            dfdx, lrpar, rpar, lipar, ipar, ierr)
         epsx = 1d-5
         epsy = 1d-5
         ierr = 0
         lnT = safe_root_with_brackets( &
            get_T_residual, lnT_min, lnT_max, &
            resid_T_min, resid_T_max, &
            imax, epsx, epsy, &
            lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) then
            write(*,2) 'T_1', N, T_1
            write(*,2) 'Tmax', N, Tmax
            write(*,2) 'resid_T_min', N, resid_T_min
            write(*,2) 'resid_T_max', N, resid_T_max
            return
            stop 'get_T failed in root find'
         end if
         T_0 = exp(lnT)
         residual = get_T_residual(lnT, &
            dfdx, lrpar, rpar, lipar, ipar, ierr)
         get_T_estimate = .true.
      end function get_T_estimate
      
      logical function get_T(ierr)
         integer, intent(out) :: ierr
         real(dp) :: Prad
         ierr = 0
         if (get_T_estimate()) then
            ! use T_0 from get_T_estimate as initial guess
            T4_0 = T_0*T_0*T_0*T_0
         else
            T4_0=T4_1*1.2d0
            T_0=sqrt(sqrt(T4_0))
            if (s% RSP_testing) then
               do
                  Prad = Radiation_Pressure(T_0)
                  if (Prad < P_0) exit
                  write(*,*) '1 reduce T_0 in get_T', N
                  T_0 = 0.99d0*T_0
               end do
            end if
         end if

         Lc_loop1: do ! reduce T if Lc >= L
            call EOP(s,N, &
               T_0,P_0,V_0,E_0,CP_0,QQ_0,SVEL_0,OP_0,ierr)  
            if (ierr /= 0) return    
            call CFLUX(HP_0,IGR_0,Lc_0,w_0,GPF,N)
            if(Lc_0.lt.L) exit Lc_loop1
            T_0=T_0-10.d0*(Lc_0/L) 
            if (T_0 <= 0) then
               ierr = -1
               return
               stop 'T_0 <= 0 in Lc_loop1'
            end if
            T4_0=T_0**4
         end do Lc_loop1
         
         TT=4.d0*SIG*P4**2*R_1**4/(3.d0*dm_bar_0*L)
         Lr_0=TT*(T4_0/OP_0-T4_1/OP_1)/ &
               (1.d0-dlog(OP_0/OP_1)/dlog(T4_0/T4_1))*L
         F1=(Lr_0+Lc_0)/L-1.d0  
         
         D=T4_0/1.d3
         I=0
         T1_loop: do ! adjust T to make Lr + Lc = L
            I=I+1
            if(I.gt.10000) then
               get_T = .false.
               if (s% RSP_testing) then
                  write(*,*) 'failed get_T', N, T_0
                  stop 'get_T'
               end if
               return
            end if
            do while (abs(D/T4_0).gt.0.5d0) 
               D=(T4_0/2.d0)*(D/abs(D))
            end do
            T4_0 = T4_0-D
            T_0=sqrt(sqrt(T4_0))
            Lc_loop: do ! reduce T if Lc >= L
               call EOP(s,N, &
                   T_0,P_0,V_0,E_0,CP_0,QQ_0,SVEL_0,OP_0,ierr)
               if (ierr /= 0) return
               call CFLUX(HP_0,IGR_0,Lc_0,w_0,GPF,N)
               if (Lc_0.lt.L) exit Lc_loop
               T_0=T_0-10.d0*(Lc_0/L)
               if (T_0 <= 0) then
                  ierr = -1
                  return
                  stop 'T_0 <= 0 in Lc_loop'
               end if
               T4_0=T_0**4
            end do Lc_loop
            Lr_0=TT*(T4_0/OP_0-T4_1/OP_1)/ &
                  (1.d0-dlog(OP_0/OP_1)/dlog(T4_0/T4_1))*L
            F2=(Lr_0+Lc_0)/L-1.d0  
            if(ABS(F2).lt.PREC .or. F2.eq.F1) exit T1_loop
            D=F2*(T4_0-T4_0)/(F2-F1)
            F1=F2
            T4_0=T4_0
         end do T1_loop
         
         get_T = .true.
         if (s% RSP_testing) write(*,*) 'done get_T', N, T_0
      end function get_T
      
      subroutine failed
         write(*,*) 'NO CONVERGENCE IN STA INNER LOOP',I
         if((.not.adjusting_dmN) .OR. dmN_cnt.gt.1 .or. IG > 54) then
            write(*,*) 'zone ',N,'IGR= ',IGR_0 
            stop
         end if
         dmN=dmN/4.d0
         write(*,*) 'PHOENIX CONDITION'
         IG=IG+1
      end subroutine failed
      
      subroutine prepare_for_new_H
         ITIN = ITIN+1
         ! STRATOWE WARTOSCI DLA ITERACJI PONIZEJ TH0
         N1 = NZN-NZT+1
         R_1 = R(N1)
         V_0  = Vol(N1)
         P_0  = P(N1)
         T_0  = T(N1)
         T4_0  = T_0**4
         M_0  = M(N1)
         dm_0 = dm(N1)
         OP_0 = K(N1)
         N  = N1
         H  = HH
         HP_0 = Hp_face(N1)
         IGR_0 = Y_face(N1)
         QQ_0 = QQS(N1)
         CP_0 = CPS(N1)
         Lc_0 = Lc(N1)
         E_0 = E(N1)
      end subroutine prepare_for_new_H

      subroutine change_NZN
         NZN=NZN-N
         do I=1,NZN
            R(I)=R(I+N)
            P(I)=P(I+N)
            Vol(I)=Vol(I+N)
            E(I)=E(I+N)
            T(I)=T(I+N)
            K(I)=K(I+N)
            M(I)=M(I+N)
            dm(I)=dm(I+N)
            dm_bar(I)=dm_bar(I+N)
            Hp_face(I) = Hp_face(I+N)
            Y_face(I) = Y_face(I+N)
            CPS(I) = CPS(I+N)
            QQS(I) = QQS(I+N)
            Lc(I) = Lc(I+N)
            Lr(I) = Lr(I+N)
            w(I) = w(I+N)
         end do  
      end subroutine change_NZN

      end subroutine STAH


      subroutine ZNVAR(s,H,dmN,L,TE,M,ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      real(dp) :: H,dmN,L,TE,M,ha,Psurf
      real(dp) :: T0,R,TAU0,Pdm,CKP,P,PZ,V,OP,DU1,DU2,CP,QQ, &
         dtau, kap, alfa, G_M_dtau_div_R2, Prad, Pgas_0, Pgas_1, &
         dP_dV, dkap_dV, xx, residual, d_residual_dlnV, &
         dkap_dlnV, dlnV, dP_dlnV, lnV
      integer :: I
      H=HHFAC
      ierr = 0
      if(.not.RSP_eddi) then !     EXACT GREY RELATION
         T0= pow(sqrt(3.d0)/4.d0,0.25d0)*TE !0.811194802d0*TE
      else !     EDDINGTON APPROXIMATION
         T0= pow(0.5d0, 0.25d0)*TE ! T0= pow(1.0d0/2.d0,0.25d0)*TE 
      endif      
      if (s% RSP_use_Prad_for_Psurf) then
         Psurf = crad*T0*T0*T0*T0/3d0
      else
         Psurf = 0d0
      end if
      R=sqrt(L/(4.d0*PI*SIG))/TE**2
      TAU0=0.003d0

      ! dtau is optical depth we want to middle of 1st cell; ~ 0.003
      ! dtau = dm*kap/(4*pi*R^2)
      ! mass of 1st cell will be 2*dm, so dm is mass above center of cell.
      ! R is derived from L and Teff.
      ! T0 is derived from Teff. it is temperature at center of 1st cell.
      ! Psurf is pressure just outside 1st cell; either 0 or, better, Prad(T0).
      ! need to find P, pressure at center of 1st cell
      ! P = Psurf + G*M*dm/(4*pi*R^4)
      ! substitute for dm to get
      ! P = Psurf + G*M*dtau/(R^2*kap)
      ! kap depends on P, so need to solve this implicit equation iteratively.
      ! once find P and kap, can evaluate dm = 4*pi*R^2*dtau/kap
      
      ! note that since T0 is fixed, we can only change P by changing Pgas.
      ! for most cases this is not a problem since Prad << Pgas.
      ! but for massive blue stars, we can have Prad >> Pgas.
      ! to handle both cases well, we need to iterate on Pgas instead of P
      ! that will let us make sure we don't create guesses that imply Pgas < 0.
      
      dtau = TAU0 ! s% RSP_outer_dtau_target
      ! make rough initial guess for opacity based on T0
      if (T0 < 4700d0) then
         kap = 1d-3
      else if (T0 > 5100d0) then
         kap = 1d-2*(1d0 + X)
      else
         alfa = (T0 - 4700)/(5100 - 4700)
         kap = alfa*1d-2*(1d0 + X) + (1d0-alfa)*1d-3
      end if
      G_M_dtau_div_R2 = G*M*dtau/R**2
      Prad = crad*T0*T0*T0*T0/3d0
      Pgas_0 = G_M_dtau_div_R2/kap
      P = Pgas_0 + Prad ! initial guess for P
      call EOP(s,0,T0,P,V,xx,xx,xx,xx,xx,ierr) ! initial V
      if (ierr /= 0) return
      !write(*,*) 'init T0,P,V', T0,P,V
      do I=1,25
         call mesa_eos_kap(s,-4, &
            T0,V,P,dP_dV,xx,xx,xx,xx,xx,xx,xx,xx,xx,xx,kap,dkap_dV,xx,ierr)
         if (ierr /= 0) return
         residual = P - (Prad + G_M_dtau_div_R2/kap)
         if (abs(residual) < 1d-6*P) exit ! done
         dP_dlnV = dP_dV*V
         lnV = log(V)
         dkap_dlnV = dkap_dV*V
         d_residual_dlnV = dP_dlnV + G_M_dtau_div_R2*dkap_dlnV/kap**2
         dlnV = - residual/d_residual_dlnV
         V = exp(lnV + dlnV)
         !write(*,*) 'T0, new V, P, kap, residual, dP_dV', i, T0, V, P, kap, residual, dP_dV
      end do
      
      !write(*,*) 'V, P, kap, residual', i, V, P, kap, residual
      
      dmN = 4*pi*R**2*dtau/kap
      if (s% RSP_testing) write(*,*) 'initial dmN', dmN/SUNM
      !stop 
      end subroutine ZNVAR


      subroutine CFLUX(HP_0,IGR_0,Lc_0,OMEGA_0,GPF,N)
      implicit none

      real(dp) :: POM,POM2,HP_0,IGR_0,OMEGA_0,PII,Lc_0
      real(dp) :: FF,GG,GPF,ENT
      real(dp) :: AA,BB,CC,DELTA
      integer :: N
!-
      if(ALFA.eq.0.d0) then
         OMEGA_0=0.d0
         Lc_0=0.d0
         return
      endif

!     PRESSURE SCALE HEIGHT
      HP_0=R_1**2/(G*M_0)*(P_0*V_0+P_1*V_1)/2.d0
      POM=P4*R_1**2*HP_0/dm_bar_0/(V_0+V_1)*2.d0

!     SUPERADIABATIC GRADIENT, Y
      IGR_0=POM*((QQ_0/CP_0+QQ_1/CP_1)/2.d0*(P_1-P_0) &
            -(dlog(T_1)-dlog(T_0)))!hyt!

      POM=sqrt(2.d0/3.d0)*0.5d0
      ENT=(E_0+P_0*V_0)/T_0+(E_1+P_1*V_1)/T_1
      FF=POM*ENT
      POM=ALFAS*ALFA
      POM2=0.5d0*(CP_0+CP_1)
      GG=POM*POM2*IGR_0
      GPF=GG/FF
   

!     BOTTOM BOUNDARY CONDITION FOR CONVECTION
      if(N.le.IBOTOM)then 
         Lc_0=0.d0
         OMEGA_0=0.d0
         return
      endif
      if(IGR_0.gt.0.d0)then
         if(.true. .or. GAMMAR.eq.0.d0)then   ! gammar breaks Cep 11.5M model  BP
           OMEGA_0=sqrt(ALFA/CEDE*FF*GPF* &
             (T_0*P_0*QQ_0/CP_0+T_1*P_1*QQ_1/CP_1)*0.5d0)
         else
            AA=CEDE/ALFA
            POM=(GAMMAR**2/ALFA**2)*4.d0*SIG
            BB=(POM/HP_0)* &
               (T_0**3*V_0**2/CP_0/OP_0 &
                 +T_1**3*V_1**2/CP_1/OP_1)*0.5d0
            CC=-(T_0*P_0*QQ_0/CP_0 &
                     +T_1*P_1*QQ_1/CP_1)*0.5d0*FF*GPF
            DELTA=BB**2-4.d0*AA*CC
            if(DELTA.le.0.d0) then
               write(*,*) 'CFLUX: Error! : Y>0, but no solution found'
               stop
            endif
            OMEGA_0=(-BB+sqrt(DELTA))/(2.d0*AA)
         endif
         PII=FF*OMEGA_0*GPF
         Lc_0=P4*R_1**2*(T_0/V_0+T_1/V_1)*0.5d0*PII*(ALFAC/ALFAS)
!                                                 (REPLACES ALFAS BY ALFAC)
      else
         OMEGA_0=0.d0
         Lc_0=0.d0
      endif
      end subroutine CFLUX

            
      end module rsp_build
