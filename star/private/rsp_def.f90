! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton & The MESA Team
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

      module rsp_def
      use const_def, only: dp, qp, pi, crad, clight, ln10
      use math_lib
      use utils_lib, only: is_bad
      use star_def, only: star_info

      implicit none
      
      integer, parameter :: MAX_NZN = 1501
      
      real(dp), parameter :: f_Edd_isotropic = 1d0/3d0, f_Edd_free_stream = 1d0

      real(dp) :: rsp_tau_factor, rsp_min_dr_div_cs, rsp_min_rad_diff_time, Psurf_from_atm
      integer :: i_min_dr_div_cs, i_min_rad_diff_time
      
      real(dp), pointer, dimension(:) :: &
         dVol_dr_00, dVol_dr_in, &
         d_egas_dVol, d_egas_dT, d_egas_dr_00, d_egas_dr_in, &
         d_Pg_dVol, d_Pg_dT, d_Pg_dr_00, d_Pg_dr_in, &   
         dK_dVol, dK_dT, dK_dr_00, dK_dr_in, &
         dQQ_dVol, dQQ_dT, dQQ_dr_00, dQQ_dr_in, &
         dCp_dVol, dCp_dT, dCp_dr_00, dCp_dr_in, &
         d_Pr_dVol, d_Pr_der, d_Pr_dr_00, d_Pr_dr_in, &   

         dHp_dr_out, dHp_dr_00, dHp_dr_in, &
         dHp_dVol_00, dHp_dVol_out, &
         dHp_dT_00, dHp_dT_out, &
         dHp_der_00, dHp_der_out, &
         
         dY_dr_in, dY_dr_00, dY_dr_out, &
         dY_dVol_00, dY_dVol_out, &
         dY_dT_00, dY_dT_out, &
         dY_der_00, dY_der_out, &
         
         dPII_dr_in, dPII_dr_00, dPII_dr_out, &
         dPII_dVol_00, dPII_dVol_out, &
         dPII_dT_00, dPII_dT_out, &
         dPII_der_00, dPII_der_out, &
         
         d_avQ_dr_00, d_avQ_dr_in, &
         d_avQ_dVol, d_avQ_dT, d_avQ_der, &
         
         dPt_dr_00, dPt_dr_in, dPt_dVol_00, dPt_dw_00, &
         
         dChi_dr_in2, dChi_dr_in, dChi_dr_00, dChi_dr_out, &  
         dChi_dVol_in, dChi_dVol_00, dChi_dVol_out, &
         dChi_dT_in, dChi_dT_00, dChi_dT_out, &
         dChi_der_in, dChi_der_00, dChi_der_out, &
         dChi_dw_00, &
                
         dEq_dr_out, dEq_dr_00, dEq_dr_in, dEq_dr_in2, &
         dEq_dVol_out, dEq_dVol_00, dEq_dVol_in, &
         dEq_dT_out, dEq_dT_00, dEq_dT_in, &
         dEq_der_out, dEq_der_00, dEq_der_in, &
         dEq_dw_00, &
         
         dC_dr_in2, dC_dr_in, dC_dr_00, dC_dr_out, &         
         dC_dVol_in, dC_dVol_00, dC_dVol_out, &
         dC_dT_in, dC_dT_00, dC_dT_out, &
         dC_der_in, dC_der_00, dC_der_out, &
         dC_dw_00, &

         photo_T, photo_r, photo_Vol, photo_w, photo_opacity, photo_QQ, &
         photo_Pgas, photo_Prad, photo_egas, photo_erad, photo_Cp, &
         photo_v, photo_M, photo_dm, photo_dm_bar, photo_csound

      ! arrays for LAPACK must be declared at compile time. legacy fortran issue.
      integer, parameter :: NV=5
      integer, parameter :: HD_DIAG=2*NV+1, LD_HD=4*NV+1, LD_ABB=6*NV+1, LPSZ=NV*MAX_NZN+1
      real(dp) DX(LPSZ), HR(LPSZ), HD(LD_HD, LPSZ), ABB(LD_ABB, LPSZ)
      integer IPVT(LPSZ)      
      
      integer, parameter :: LD_LLL = 4*MAX_NZN
      real(dp), dimension(LD_LLL, LD_LLL) :: LLL, VLx, VRx
      integer :: ISORTx(LD_LLL)
      real(dp) :: WORKx(4*LD_LLL), WRx(LD_LLL), WIx(LD_LLL)
      
      ! for rsp_eval_eos_and_kap
      real(dp), pointer :: xa(:)
      real(dp) :: X, Z, Y, abar, zbar, z53bar, XC, XN, XO, Xne

      ! these for for rsp.f90 period and work calculations
      real(dp) :: ETOT, EGRV, ETHE, EKIN, EDE_start, ECON, &
         TE, ELSTA, TEFF, E0, TT1, TE_start, T0, UN, ULL, &
         RMAX, LMAX, LMIN, EKMAX, EKMIN, EKMAXL, EKDEL, &
         RSTA, RMIN, PERIODL, PERIODLIN, &
         PDVWORK, FASE0
      real(dp), pointer, dimension(:) :: &
         PPP0, PPQ0, PPT0, PPC0, VV0, &
         WORK, WORKQ, WORKT, WORKC               
      integer :: INSIDE, IWORK, ID, NSTART, FIRST, &
         run_num_retries_prev_period, prev_cycle_run_num_steps, &
         run_num_iters_prev_period
      
      ! for maps
      logical :: writing_map, done_writing_map
      integer, parameter :: max_map_cols = 200
      integer :: map_ids(max_map_cols), num_map_cols
      character(256) :: map_col_names(max_map_cols)
      
      ! marsaglia and zaman random number generator. period is 2**43 with
      ! 900 million different sequences. the state of the generator (for restarts)
      integer, parameter :: rn_u_len=97
      integer :: rn_i97, rn_j97
      real(dp) :: rn_u(rn_u_len), rn_c, rn_cd, rn_cm

      ! from const
      real(dp) :: G, SIG, SUNL, SUNM, SUNR, CL, P43, P4      
      
      ! these are set from inlist
      real(dp) :: ALFA, ALFAP, ALFAM, ALFAT, ALFAS, ALFAC, CEDE, GAMMAR
      real(dp) :: THETA, THETAT, THETAQ, THETAU, THETAE, WTR, WTC, WTT, GAM
      real(dp) :: THETA1, THETAT1, THETAQ1, THETAU1, THETAE1, WTR1, WTC1, WTT1, GAM1
      real(dp) :: EFL0, CQ, ZSH, kapE_factor, kapP_factor
      integer :: NZN, IBOTOM
      
      integer :: ITOP ! below non convective region at surface
      

      contains
      
      
      subroutine init_def(s)
         use const_def, only: standard_cgrav, boltz_sigma, &
            Lsun, Msun, Rsun
         type (star_info), pointer :: s
         
         P4=4.d0*PI
         P43=P4/3.d0
         
         G=standard_cgrav
         SIG=boltz_sigma
         SUNL=Lsun
         SUNM=Msun
         SUNR=Rsun
         CL=4d0*(4d0*PI)**2*SIG/3d0

         ALFA = s% RSP_alfa
         ALFAP = s% RSP_alfap
         ALFAM = s% RSP_alfam
         ALFAT = s% RSP_alfat
         ALFAS = s% RSP_alfas
         ALFAC = s% RSP_alfac
         CEDE = s% RSP_alfad
         GAMMAR = s% RSP_gammar

         ALFAP = ALFAP*2.d0/3.d0
         ALFAS = ALFAS*(1.d0/2.d0)*sqrt(2.d0/3.d0)
         ALFAC = ALFAC*(1.d0/2.d0)*sqrt(2.d0/3.d0)
         CEDE  = CEDE*(8.d0/3.d0)*sqrt(2.d0/3.d0)
         GAMMAR = GAMMAR*2.d0*sqrt(3.d0)
         
         call turn_on_time_weighting(s)

         CQ = s% RSP_cq
         ZSH = s% RSP_zsh
         EFL0 = s% RSP_efl0
         kapE_factor = 1d0 ! s% RSP_kapE_factor
         kapP_factor = 1d0 ! s% RSP_kapP_factor
         
         if (ALFA == 0.d0) EFL0=0.d0
         
         writing_map = .false.
            
      end subroutine init_def
      
         
      subroutine turn_off_time_weighting(s)
         type (star_info), pointer :: s
         THETA = 1d0
         THETAT = 1d0
         THETAQ = 1d0
         THETAE = 1d0
         THETAU = 1d0
         WTR = 1d0
         WTC = 1d0
         WTT = 1d0
         GAM = 1d0
         GAM1=0d0
         WTR1=0d0
         WTC1=0d0
         WTT1=0d0
         THETA1 =0d0
         THETAT1=0d0
         THETAQ1=0d0
         THETAE1=0d0
         THETAU1=0d0
      end subroutine turn_off_time_weighting
      
      
      subroutine turn_on_time_weighting(s)
         type (star_info), pointer :: s
         THETA = s% RSP_theta
         THETAT = s% RSP_thetat
         THETAQ = s% RSP_thetaq
         THETAE = s% RSP_thetae
         THETAU = s% RSP_thetau
         WTR = s% RSP_wtr
         WTC = s% RSP_wtc
         WTT = s% RSP_wtt
         GAM = s% RSP_gam
         GAM1=1.d0-GAM
         WTR1=1.d0-WTR
         WTC1=1.d0-WTC
         WTT1=1.d0-WTT
         THETA1 =1.d0-THETA
         THETAT1=1.d0-THETAT
         THETAQ1=1.d0-THETAQ
         THETAE1=1.d0-THETAE
         THETAU1=1.d0-THETAU
      end subroutine turn_on_time_weighting
      
      
      subroutine init_allocate(s,nz)
         !use rsp_eddfac, only: eddfac_allocate
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         integer :: n
         NZN = nz
         if (NZN > MAX_NZN) then
            write(*,*) 'NZN > MAX_NZN', NZN, MAX_NZN
            stop 'rsp init_allocate'
         end if
         IBOTOM = NZN/s% RSP_nz_div_IBOTOM
         n = NZN + 1 ! room for ghost cell
         allocate(xa(s% species), &
            dVol_dr_00(n), dVol_dr_in(n), &
            d_egas_dVol(n), d_egas_dT(n), d_egas_dr_00(n), d_egas_dr_in(n), &
            d_Pg_dVol(n), d_Pg_dT(n), d_Pg_dr_00(n), d_Pg_dr_in(n), &   
            d_Pr_dVol(n), d_Pr_der(n), d_Pr_dr_00(n), d_Pr_dr_in(n), &   
            dK_dVol(n), dK_dT(n), dK_dr_00(n), dK_dr_in(n), &
            dQQ_dVol(n), dQQ_dT(n), dQQ_dr_00(n), dQQ_dr_in(n), &
            dCp_dVol(n), dCp_dT(n), dCp_dr_00(n), dCp_dr_in(n), &
            dHp_dr_out(n), dHp_dr_00(n), dHp_dr_in(n), &
            dHp_dVol_00(n), dHp_dVol_out(n), &
            dHp_dT_00(n), dHp_dT_out(n), dHp_der_00(n), dHp_der_out(n), &         
            dY_dr_in(n), dY_dr_00(n), dY_dr_out(n), &
            dY_dVol_00(n), dY_dVol_out(n), &
            dY_dT_00(n), dY_dT_out(n), dY_der_00(n), dY_der_out(n), &         
            dPII_dr_in(n), dPII_dr_00(n), dPII_dr_out(n), &
            dPII_dVol_00(n), dPII_dVol_out(n), &
            dPII_dT_00(n), dPII_dT_out(n), dPII_der_00(n), dPII_der_out(n), &         
            d_avQ_dr_00(n), d_avQ_dr_in(n), d_avQ_dVol(n), d_avQ_dT(n), d_avQ_der(n), &         
            dPt_dr_00(n), dPt_dr_in(n), dPt_dVol_00(n), dPt_dw_00(n), &         
            dChi_dr_in2(n), dChi_dr_in(n), dChi_dr_00(n), dChi_dr_out(n), &  
            dChi_dVol_in(n), dChi_dVol_00(n), dChi_dVol_out(n), &
            dChi_dT_in(n), dChi_dT_00(n), dChi_dT_out(n), &
            dChi_der_in(n), dChi_der_00(n), dChi_der_out(n), &
            dChi_dw_00(n), &                
            dEq_dr_out(n), dEq_dr_00(n), dEq_dr_in(n), dEq_dr_in2(n), &
            dEq_dVol_out(n), dEq_dVol_00(n), dEq_dVol_in(n), &
            dEq_dT_out(n), dEq_dT_00(n), dEq_dT_in(n), &
            dEq_der_out(n), dEq_der_00(n), dEq_der_in(n), &
            dEq_dw_00(n), &         
            dC_dr_in2(n), dC_dr_in(n), dC_dr_00(n), dC_dr_out(n), &         
            dC_dVol_in(n), dC_dVol_00(n), dC_dVol_out(n), &
            dC_dT_in(n), dC_dT_00(n), dC_dT_out(n), &
            dC_der_in(n), dC_der_00(n), dC_der_out(n), dC_dw_00(n), &
            photo_T(n), photo_r(n), photo_Vol(n), photo_w(n), photo_opacity(n), photo_QQ(n), &
            photo_Pgas(n), photo_Prad(n), photo_egas(n), photo_erad(n), photo_Cp(n), &
            photo_v(n), photo_M(n), photo_dm(n), photo_dm_bar(n), photo_csound(n), &
            PPP0(n), PPQ0(n), PPT0(n), PPC0(n), VV0(n), &
            WORK(n), WORKQ(n), WORKT(n), WORKC(n))
      end subroutine init_allocate
      
      
      subroutine init_free(s)
         type (star_info), pointer :: s
         deallocate(xa, &
            dVol_dr_00, dVol_dr_in, &
            d_egas_dVol, d_egas_dT, d_egas_dr_00, d_egas_dr_in, &
            d_Pg_dVol, d_Pg_dT, d_Pg_dr_00, d_Pg_dr_in, &   
            d_Pr_dVol, d_Pr_der, d_Pr_dr_00, d_Pr_dr_in, &   
            dK_dVol, dK_dT, dK_dr_00, dK_dr_in, &
            dQQ_dVol, dQQ_dT, dQQ_dr_00, dQQ_dr_in, &
            dCp_dVol, dCp_dT, dCp_dr_00, dCp_dr_in, &
            dHp_dr_out, dHp_dr_00, dHp_dr_in, &
            dHp_dVol_00, dHp_dVol_out, &
            dHp_dT_00, dHp_dT_out, &
            dHp_der_00, dHp_der_out, &         
            dY_dr_in, dY_dr_00, dY_dr_out, &
            dY_dVol_00, dY_dVol_out, &
            dY_dT_00, dY_dT_out, dY_der_00, dY_der_out, &         
            dPII_dr_in, dPII_dr_00, dPII_dr_out, &
            dPII_dVol_00, dPII_dVol_out, &
            dPII_dT_00, dPII_dT_out, dPII_der_00, dPII_der_out, &         
            d_avQ_dr_00, d_avQ_dr_in, d_avQ_dVol, d_avQ_dT, d_avQ_der, &         
            dPt_dr_00, dPt_dr_in, dPt_dVol_00, dPt_dw_00, &         
            dChi_dr_in2, dChi_dr_in, dChi_dr_00, dChi_dr_out, &  
            dChi_dVol_in, dChi_dVol_00, dChi_dVol_out, &
            dChi_dT_in, dChi_dT_00, dChi_dT_out, &
            dChi_der_in, dChi_der_00, dChi_der_out, &
            dChi_dw_00, &                
            dEq_dr_out, dEq_dr_00, dEq_dr_in, dEq_dr_in2, &
            dEq_dVol_out, dEq_dVol_00, dEq_dVol_in, &
            dEq_dT_out, dEq_dT_00, dEq_dT_in, &
            dEq_der_out, dEq_der_00, dEq_der_in, &
            dEq_dw_00, &         
            dC_dr_in2, dC_dr_in, dC_dr_00, dC_dr_out, &         
            dC_dVol_in, dC_dVol_00, dC_dVol_out, &
            dC_dT_in, dC_dT_00, dC_dT_out, &
            dC_der_in, dC_der_00, dC_der_out, dC_dw_00, &
            photo_T, photo_r, photo_Vol, photo_w, photo_csound, &
            photo_Pgas, photo_Prad, photo_egas, photo_erad, photo_Cp, photo_QQ, &
            photo_v, photo_M, photo_dm, photo_dm_bar, photo_opacity, &
            PPP0, PPQ0, PPT0, PPC0, VV0, &
            WORK, WORKQ, WORKT, WORKC)
      end subroutine init_free
      
      
      real(dp) function rsp_phase_time0()
         rsp_phase_time0 = TT1
      end function rsp_phase_time0
      
      
      real(dp) function rsp_WORK(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         rsp_WORK = WORK(k)
      end function rsp_WORK
      
      
      real(dp) function rsp_WORKQ(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         rsp_WORKQ = WORKQ(k)
      end function rsp_WORKQ
      
      
      real(dp) function rsp_WORKT(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         rsp_WORKT = WORKT(k)
      end function rsp_WORKT
      
      
      real(dp) function rsp_WORKC(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         rsp_WORKC = WORKC(k)
      end function rsp_WORKC
      
      
      subroutine rsp_photo_out(s, iounit)
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer :: n
         include 'formats'
         !call copy_from_xh_to_rsp(s, NZN) ! resynch before photo
         n = NZN + 1
         write(iounit) NZN
         write(iounit) xa(1:s% species), &
            X, Z, Y, abar, zbar, z53bar, XC, XN, XO, Xne, IBOTOM, &
            s% RSP_num_periods, s% RSP_dt, s% RSP_period, s% rsp_DeltaR, &
            s% rsp_DeltaMag, s% rsp_GRPDV, s% rsp_GREKM, s% rsp_GREKM_avg_abs, &
            rsp_tau_factor, rsp_min_dr_div_cs, rsp_min_rad_diff_time, &
            i_min_dr_div_cs, i_min_rad_diff_time, Psurf_from_atm, &
            s% Fr(1:n), s% Lc(1:n), s% Lt(1:n), s% Y_face(1:n), &
            s% Pt(1:n), s% Chi(1:n), s% COUPL(1:n), s% avQ(1:n), &
            s% T(1:n), s% r(1:n), s% Vol(1:n), s% RSP_w(1:n), &
            s% Pgas(1:n), s% Prad(1:n), s% csound(1:n), s% Cp(1:n), &
            s% egas(1:n), s% erad(1:n), s% opacity(1:n), s% QQ(1:n), &
            s% v(1:n), s% M(1:n), s% dm(1:n), s% dm_bar(1:n), &
            ETOT, EGRV, ETHE, EKIN, EDE_start, ECON, &
            TE, ELSTA, TEFF, E0, TT1, TE_start, T0, UN, ULL, &
            RMAX, LMAX, LMIN, EKMAX, EKMIN, EKMAXL, EKDEL, &
            RSTA, RMIN, PERIODL, PERIODLIN, &
            PDVWORK, FASE0, INSIDE, IWORK, ID, NSTART, FIRST, &
            s% rsp_LINA_periods(1:3), s% rsp_LINA_growth_rates(1:3), &
            run_num_retries_prev_period, prev_cycle_run_num_steps, &
            run_num_iters_prev_period, writing_map, &
            PPP0(1:n), PPQ0(1:n), PPT0(1:n), PPC0(1:n), VV0(1:n), &
            WORK(1:n), WORKQ(1:n), WORKT(1:n), WORKC(1:n)
      end subroutine rsp_photo_out
      
      
      subroutine rsp_photo_in(s, iounit, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer :: n, k, i
         include 'formats'
         call init_def(s) 
         ierr = 0
         read(iounit, iostat=ierr) NZN
         if (ierr /= 0) stop 'read failed in rsp_photo_in'
         s% nz = NZN
         call init_allocate(s,NZN)
         n = NZN + 1
         read(iounit, iostat=ierr) xa(1:s% species), &
            X, Z, Y, abar, zbar, z53bar, XC, XN, XO, Xne, IBOTOM, &
            s% RSP_num_periods, s% RSP_dt, s% RSP_period, s% rsp_DeltaR, &
            s% rsp_DeltaMag, s% rsp_GRPDV, s% rsp_GREKM, s% rsp_GREKM_avg_abs, &
            rsp_tau_factor, rsp_min_dr_div_cs, rsp_min_rad_diff_time, &
            i_min_dr_div_cs, i_min_rad_diff_time, Psurf_from_atm, &
            s% Fr(1:n), s% Lc(1:n), s% Lt(1:n), s% Y_face(1:n), &
            s% Pt(1:n), s% Chi(1:n), s% COUPL(1:n), s% avQ(1:n), &
            photo_T(1:n), photo_r(1:n), photo_Vol(1:n), photo_w(1:n), &
            photo_Pgas(1:n), photo_Prad(1:n), photo_csound(1:n), photo_Cp(1:n), &
            photo_egas(1:n), photo_erad(1:n), photo_opacity(1:n), photo_QQ(1:n), &
            photo_v(1:n), photo_M(1:n), photo_dm(1:n), photo_dm_bar(1:n), &
            ETOT, EGRV, ETHE, EKIN, EDE_start, ECON, &
            TE, ELSTA, TEFF, E0, TT1, TE_start, T0, UN, ULL, &
            RMAX, LMAX, LMIN, EKMAX, EKMIN, EKMAXL, EKDEL, &
            RSTA, RMIN, PERIODL, PERIODLIN, &
            PDVWORK, FASE0, INSIDE, IWORK, ID, NSTART, FIRST, &
            s% rsp_LINA_periods(1:3), s% rsp_LINA_growth_rates(1:3), &
            run_num_retries_prev_period, prev_cycle_run_num_steps, &
            run_num_iters_prev_period, writing_map, &
            PPP0(1:n), PPQ0(1:n), PPT0(1:n), PPC0(1:n), VV0(1:n), &
            WORK(1:n), WORKQ(1:n), WORKT(1:n), WORKC(1:n)
         if (ierr /= 0) stop 'read failed in rsp_photo_in'
         if (writing_map) then
            write(*,*) 'sorry, cannot use photo to resume writing map data file'
            writing_map = .false.
         end if
      end subroutine rsp_photo_in
      
      
      subroutine finish_after_build_model(s)
         type (star_info), pointer :: s
         integer :: k, ierr
         ! restore bit-for-bit erad = crad*T**4*Vol
         include 'formats'
         do k = 1, NZN
            s% erad(k) = crad*s% T(k)**4*s% Vol(k)
            s% Prad(k) = s% f_Edd(k)*s% erad(k)/s% Vol(k)
         end do
      end subroutine finish_after_build_model
      
      
      subroutine finish_rsp_photo_in(s)
         use star_utils, only: set_rmid
         type (star_info), pointer :: s
         integer :: k, ierr
         ! restore bit-for-bit same as before photo
         include 'formats'
         do k = 1, NZN
            s% T(k) = photo_T(k)
            s% Pgas(k) = photo_Pgas(k)
            s% Prad(k) = photo_Prad(k)
            s% P(k) = s% Pgas(k) + s% Prad(k)
            s% egas(k) = photo_egas(k)
            s% erad(k) = photo_erad(k)
            s% opacity(k) = photo_opacity(k)
            s% csound(k) = photo_csound(k)
            s% Cp(k) = photo_Cp(k)
            s% QQ(k) = photo_QQ(k)
            s% r(k) = photo_r(k)
            s% Vol(k) = photo_Vol(k)
            s% RSP_w(k) = photo_w(k)
            s% v(k) = photo_v(k)
            s% M(k) = photo_M(k)
            s% dm(k) = photo_dm(k)
            s% dm_bar(k) = photo_dm_bar(k)
         end do
         call set_rmid(s, 1, NZN, ierr)
         if (ierr /= 0) then
            write(*,*) 'finish_rsp_photo_in failed in set_rmid'
            stop 'finish_rsp_photo_in'
         end if
      end subroutine finish_rsp_photo_in
      
      
      subroutine set_build_vars(s, &
            m, dm, dm_bar, r, Vol, T, RSP_Et, Lr, Lc)
         type (star_info), pointer :: s
         real(dp), dimension(:), intent(in) :: &
            m, dm, dm_bar, r, Vol, T, RSP_Et, Lr, Lc
         integer :: k, i
         include 'formats'
         do i=1, NZN
            k = NZN+1 - i
            s% m(k) = m(i)
            s% dm(k) = dm(i)
            s% dm_bar(k) = dm_bar(i)
            s% r(k) = r(i)
            s% Vol(k) = Vol(i)
            s% T(k) = T(i)
            s% RSP_Et(k) = RSP_Et(i)
            s% RSP_w(k) = sqrt(s% RSP_Et(k))
            s% Fr(k) = Lr(i)/(4d0*pi*s% r(k)**2)
            s% erad(k) = crad*s% T(k)**4*s% Vol(k)
            s% L(k) = Lr(i) + Lc(i)
            if (is_bad(s% L(k))) then
               write(*,2) 'L Lr Lc', k, s% L(k), Lr(i), Lc(i)
               stop 'set_build_vars'
            end if
         end do
      end subroutine set_build_vars
      
      
      subroutine set_star_vars(s, ierr)
         use star_utils, only: normalize_dqs, set_qs, set_m_and_dm, set_dm_bar
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: sum_dm
         integer :: i, k
         include 'formats'
         ierr = 0
         sum_dm = 0d0         
         do k=1, NZN
            sum_dm = sum_dm + s% dm(k)
         end do
         do k=1, NZN
            s% dq(k) = s% dm(k)/sum_dm
         end do
         s% xmstar = sum_dm
         if (s% nz /= NZN) then
            write(*,2) 'NZN', NZN
            write(*,2) 's% nz', s% nz
            stop 'bad nz'
         end if
         if (.not. s% do_normalize_dqs_as_part_of_set_qs) then
            call normalize_dqs(s, s% nz, s% dq, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'normalize_dqs failed in rsp_def set_star_vars'
               return
            end if
         end if
         call set_qs(s, s% nz, s% q, s% dq, ierr)
         if (ierr /= 0) stop 'failed in set_qs'         
         s% m(1) = s% mstar
         do k=2,s% nz
            s% dm(k-1) = s% dq(k-1)*s% xmstar
            s% m(k) = s% m(k-1) - s% dm(k-1)
         end do
         k = s% nz
         s% dm(k) = s% m(k) - s% m_center         
         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)                  
         do k=1, NZN
            if (k==NZN) then
               s% Vol(k)=P43/s% dm(k)*(s% r(k)**3 - s% R_center**3)
               s% rmid(k) = 0.5d0*(s% r(k) + s% R_center)
            else
               s% Vol(k)=P43/s% dm(k)*(s% r(k)**3 - s% r(k+1)**3)
               s% rmid(k) = 0.5d0*(s% r(k) + s% r(k+1))
            end if
            if (is_bad(s% Vol(k)))then
               write(*, 2) 's% Vol(k)', k, s% Vol(k)
               stop 'set_star_vars'
            end if            
            s% rho(k) = 1d0/s% Vol(k)
            s% lnd(k) = log(s% rho(k))
            s% rho(k) = exp(s% lnd(k))
            s% Vol(k) = 1d0/s% rho(k)  
            s% xh(s% i_lnd, k) = s% lnd(k)            
            s% lnT(k) = log(s% T(k))
            s% xh(s% i_lnT, k) = s% lnT(k)            
            s% lnR(k) = log(s% r(k))
            s% xh(s% i_lnR, k) = s% lnR(k)            
            s% RSP_Et(k) = s% RSP_w(k)*s% RSP_w(k)
            s% xh(s% i_etrb_RSP, k) = s% RSP_Et(k)               
            s% xh(s% i_v, k) = s% v(k)            
         end do
      end subroutine set_star_vars
      
      
      subroutine copy_from_xh_to_rsp(s, nz_new) 
         ! do this when load a file and after remesh
         type (star_info), pointer :: s
         integer, intent(in) :: nz_new
         integer :: k
         real(qp) :: q1, q2, q3, q4
         include 'formats'
         if (nz_new > 0) NZN = nz_new
         do k=NZN,1,-1
            s% lnd(k) = s% xh(s% i_lnd,k)
            s% rho(k) = exp(s% lnd(k))
            s% lnT(k) = s% xh(s% i_lnT,k)
            s% T(k) = exp(s% lnT(k))
            s% lnR(k) = s% xh(s% i_lnR,k)
            s% r(k) = exp(s% lnR(k))
            s% RSP_Et(k) = s% xh(s% i_etrb_RSP,k)
            s% RSP_w(k) = sqrt(s% RSP_Et(k))
            s% Fr(k) = s% xh(s% i_Fr_RSP,k)
            s% v(k) = s% xh(s% i_v,k)
            if (k == NZN) then ! center
               s% Vol(k)=P43/s% dm(k)*(s% r(k)**3 - s% R_center**3)
               s% rmid(k) = 0.5d0*(s% r(k) + s% R_center)
               if (s% Vol(k) <= 0d0 .or. is_bad(s% Vol(k))) then
                  write(*,2) 'copy from xh to rsp s% Vol(k) r00 r_center dm', k, &
                     s% Vol(k), s% r(k), s% R_center, s% dm(k)
                  stop 'copy_from_xh_to_rsp'
               end if
            else
               q1 = P43/s% dm(k)
               q2 = s% r(k)
               q3 = s% r(k+1)
               q4 = q1*(q2**3 - q3**3)
               s% Vol(k) = dble(q4)
               s% rmid(k) = 0.5d0*(s% r(k) + s% r(k+1))
               if (s% Vol(k) <= 0d0 .or. is_bad(s% Vol(k))) then
                  write(*,2) 'copy from xh to rsp s% Vol(k) r00 rp1 dm', k, &
                     s% Vol(k), s% r(k), s% r(k+1), s% dm(k)
                  stop 'copy_from_xh_to_rsp'
               end if
            end if
            s% erad(k) = s% xh(s% i_erad_RSP,k)
            s% Prad(k) = f_Edd_isotropic*s% erad(k)/s% Vol(k)
         end do
      end subroutine copy_from_xh_to_rsp


      subroutine check_for_T_or_P_inversions(s,str)
         type (star_info), pointer :: s      
         character (len=*), intent(in) :: str
         integer :: k
         logical :: okay
         include 'formats'
         okay = .true.
         do k=2, s% nz
            if (s% T(k) <= s% T(k-1)) then
               write(*,3) trim(str) // ' T inversion', k, s% model_number, s% T(k), s% T(k-1)
               okay = .false.
            end if
         end do
         do k=2, s% nz
            if (s% P(k) <= s% P(k-1)) then
               write(*,3) trim(str) // ' P inversion', k, s% model_number, s% P(k), s% P(k-1), &
                  s% Pt(k), s% Pt(k-1), s% avQ(k), s% avQ(k-1), s% v(k+1), s% v(k), s% v(k-1)
               okay = .false.
            end if
         end do
         if (.not. okay) stop 'rsp_one_step check_for_T_or_P_inversions'
      end subroutine check_for_T_or_P_inversions


      subroutine check_R(s,str)
         type (star_info), pointer :: s      
         character (len=*), intent(in) :: str
         integer :: k
         real(dp) :: V
         include 'formats'
         do k=s% nz,1,-1
            if (k == s% nz) then
               if (s% r(k) <= s% r_center) then
                  write(*,3) trim(str) // ' bad r', k, s% model_number, s% r(k), s% r_center
                  stop 'rsp_one_step check_R'
               end if
            else
               if (s% r(k) <= s% r(k+1)) then
                  write(*,3) trim(str) // ' bad r', k, s% model_number, s% r(k), s% r(k+1)
                  stop 'rsp_one_step check_R'
               end if
            end if
         end do
      end subroutine check_R
      
      
      subroutine rsp_dump_for_debug(s)
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         write(*,*) 'rsp_dump_for_debug'
         write(*,2) 'R_center', s% model_number, s% R_center
         write(*,2) 'xmstar', s% model_number, s% xmstar
         write(*,2) 'M/Msun', s% model_number, s% star_mass
         write(*,2) 's% R_center', s% model_number, s% R_center
         write(*,4) 'species', s% model_number, s% species
         write(*,2) 'm_center', s% model_number, s% M_center
         write(*,2) 'mstar', s% model_number, s% mstar
         write(*,2) 'L_center', s% model_number, s% L_center
         write(*,2) 'X', s% model_number, X
         write(*,2) 'Z', s% model_number, Z
         write(*,2) 'Y', s% model_number, Y
         write(*,2) 'abar', s% model_number, abar
         write(*,2) 'zbar', s% model_number, zbar
         write(*,2) 'z53bar', s% model_number, z53bar
         write(*,2) 'XC', s% model_number, XC
         write(*,2) 'XN', s% model_number, XN
         write(*,2) 'XO', s% model_number, XO
         write(*,2) 'Xne', s% model_number, Xne
         write(*,2) 'dt dt_next', s% model_number, s% dt
         write(*,2) 's% rsp_period', s% model_number, s% rsp_period
         write(*,2) 'dt', s% model_number, s% dt
         do k=1,s% nz
            write(*,2) '% Vol(k)', k, s% Vol(k)
            write(*,2) 's% v(k)', k, s% v(k)
            write(*,2) 's% r(k)', k, s% r(k)
            write(*,2) 's% dm(k)', k, s% dm(k)
            write(*,2) 's% RSP_w(k)', k, s% RSP_w(k)
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,2) 's% erad(k)', k, s% erad(k)
            write(*,2) 's% Prad(k)', k, s% Prad(k)
            write(*,2) 's% Fr(k)', k, s% Fr(k)
            !write(*,2) '', k, 
         end do
         !stop 'rsp_dump_for_debug'      
      end subroutine rsp_dump_for_debug
      
      
      subroutine cleanup_for_LINA( &
            s, M, DM, DM_BAR, R, Vol, T, RSP_Et, P, ierr)
         use star_utils, only: normalize_dqs, set_qs, set_m_and_dm, set_dm_bar
         type (star_info), pointer :: s
         real(dp), intent(inout), dimension(:) :: &
            M, DM, DM_BAR, R, Vol, T, RSP_Et, P
         integer, intent(out) :: ierr
         
         integer :: I, k
         
         include 'formats'
         
         ! get
         do i=1,NZN
            k = NZN+1-i 
            s% m(k) = M(i)
            s% dq(k) = DM(i)/s% xmstar
            s% r(k) = R(i)
            s% Vol(k) = Vol(i)
            s% T(k) = T(i)
            s% RSP_w(k) = sqrt(RSP_Et(i))
            s% P(k) = P(i)
            s% Prad(k) = crad*s% T(k)**4/3d0
            s% Pgas(k) = s% P(k) - s% Prad(k)
         end do                    
         s% dq(s% nz) = (s% m(NZN) - s% M_center)/s% xmstar
         
         ! fix
         if (.not. s% do_normalize_dqs_as_part_of_set_qs) then
            call normalize_dqs(s, NZN, s% dq, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'normalize_dqs failed in rsp_def cleanup_for_LINA'
               return
            end if
         end if
         call set_qs(s, NZN, s% q, s% dq, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_qs'
            stop 'build_rsp_model'
         end if
         call set_m_and_dm(s)
         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)
         s% Y_face(1:s% nz) = 0d0
         s% grada(1:s% nz) = 0d0
         s% rmid_start(1:s% nz) = -1d0
         s% Fr(1:s% nz) = 0d0
         s% Lc(1:s% nz) = 0d0
         s% Lt(1:s% nz) = 0d0
         s% csound(1:s% nz) = 0d0
         call copy_results(s)
         
         ! put back
         do i=1,NZN
            k = NZN+1-i 
            M(i) = s% m(k)
            DM(i) = s% dm(k)
            DM_BAR(i) = s% dm_bar(k)
            R(i) = s% r(k)
            Vol(i) = s% Vol(k)
            T(i) = s% T(k)
            RSP_Et(i) = s% RSP_w(k)**2
         end do                    
      
      end subroutine cleanup_for_LINA
      
      
      subroutine copy_results(s)
         use star_utils, only: set_rmid
         use const_def, only: convective_mixing, no_mixing, qp
         type (star_info), pointer :: s
         integer :: i, k, ierr
         real(dp) :: RSP_efl0_2
         real(qp) :: q1, q2, q3, q4
      
         RSP_efl0_2 = EFL0**2
         do i=1, NZN
         
            k = NZN+1 - i
            s% xh(s% i_v,k) = s% v(k)
            s% xh(s% i_erad_RSP,k) = s% erad(k)
            s% xh(s% i_Fr_RSP,k) = s% Fr(k)
            
            ! some tweaks needed for bit-for-bit with photos
            
            ! sqrt(w**2) /= original w, so need to redo
            s% RSP_Et(k) = s% RSP_w(k)**2
            s% xh(s% i_etrb_RSP,k) = s% RSP_Et(k)               
            s% RSP_w(k) = sqrt(s% xh(s% i_etrb_RSP,k))
            
            ! exp(log(r)) /= original r, so need to redo
            s% lnR(k) = log(s% r(k))
            s% xh(s% i_lnR,k) = s% lnR(k)
            s% r(k) = exp(s% xh(s% i_lnR,k))
            
            ! exp(log(T)) /= original T, so need to redo
            s% lnT(k) = log(s% T(k))
            s% xh(s% i_lnT,k) = s% lnT(k)
            s% T(k) = exp(s% xh(s% i_lnT,k))
            
            s% P(k) = s% Pgas(k) + s% Prad(k)
            if (k > 1) s% gradT(k) = &
               s% Y_face(k) + 0.5d0*(s% grada(k-1) + s% grada(k))
            
         end do
         s% gradT(1) = s% gradT(2)

         call set_rmid(s, 1, NZN, ierr)
         if (ierr /= 0) then
            write(*,*) 'copy_results failed in set_rmid'
            stop 'copy_results'
         end if
         
         do i=1, NZN         
            k = NZN+1 - i            
            ! revise Vol and rho using revised r
            if (i==1) then
               s% Vol(k)=P43/s% dm(k)*(s% r(k)**3 - s% R_center**3)
            else
               q1 = P43/s% dm(k)
               q2 = s% r(k)
               q3 = s% r(k+1)
               q4 = q1*(q2**3 - q3**3)
               s% Vol(k) = dble(q4)
            end if
            s% rho(k) = 1d0/s% Vol(k)            
            ! exp(log(rho)) /= original rho, so need to redo
            s% lnd(k) = log(s% rho(k))    
            s% xh(s% i_lnd,k) = s% lnd(k)
            s% rho(k) = exp(s% xh(s% i_lnd,k))
            s% L(k) = 4d0*pi*s% r(k)**2*s% Fr(k) + s% Lc(k) + s% Lt(k)
            if (s% RSP_w(k) > 1d4) then ! arbitrary cut
               s% mixing_type(k) = convective_mixing
            else
               s% mixing_type(k) = no_mixing
            end if      
         end do
         
         ! set some things for mesa output reporting
         i = 1
         s% rho_face(i) = s% rho(i)
         s% P_face_ad(i)%val = s% P(i)
         s% csound_face(i) = s% csound(i)
         do i = 2,NZN
            s% rho_face(i) = 0.5d0*(s% rho(i) + s% rho(i-1))
            s% P_face_ad(i)%val = 0.5d0*(s% P(i) + s% P(i-1))
            s% csound_face(i) = 0.5d0*(s% csound(i) + s% csound(i-1))
         end do
         
         ! these are necessary to make files consistent with photos.
         s% R_center = pow(s% r(NZN)**3 - s% Vol(NZN)*s% dm(NZN)/P43, 1d0/3d0)
         s% M_center = s% mstar - s% xmstar

      end subroutine copy_results


      end module rsp_def
      
