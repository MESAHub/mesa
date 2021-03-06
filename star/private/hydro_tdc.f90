! ***********************************************************************
!
!   Copyright (C) 2010-2020  Bill Paxton & The MESA Team
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

      module hydro_tdc

      use star_private_def
      use const_def
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use accurate_sum_auto_diff_star_order1
      use star_utils, only: em1, e00, ep1

      implicit none

      private
      public :: do1_tdc_L_eqn, do1_turbulent_energy_eqn, compute_Eq_cell, &
         compute_Uq_face, set_etrb_start_vars, reset_etrb_using_L
      
      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, &
         x_ALFAS = (1.d0/2.d0)*sqrt(2.d0/3.d0), &
         x_ALFAC = (1.d0/2.d0)*sqrt(2.d0/3.d0), &
         x_CEDE  = (8.d0/3.d0)*sqrt(2.d0/3.d0), &
         x_GAMMAR = 2.d0*sqrt(3.d0)

!         RSP       TDC
!         ALFA  =>  TDC_alfa
!         ALFAP =>  TDC_alfap*x_ALFAP
!         ALFAM =>  TDC_alfam
!         ALFAT =>  TDC_alfat
!         ALFAS =>  x_ALFAS
!         ALFAC =>  x_ALFAC
!         CEDE  =>  x_CEDE
!         GAMMAR => TDC_alfar*x_GAMMAR

      contains
      

      subroutine do1_tdc_L_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         type(auto_diff_real_star_order1) :: L_expected, L_actual, resid
         real(dp) :: scale, residual, e_avg, L_start_max
         logical :: test_partials
         include 'formats'

         test_partials = (k == s% solver_test_partials_k)
         !test_partials = .false.

         ierr = 0
         L_expected = compute_L_face(s, k, ierr)
         if (ierr /= 0) return        
         L_actual = wrap_L_00(s, k)  
         L_start_max = maxval(s% L_start(1:s% nz))
         scale = 1d0/L_start_max
         if (is_bad(scale)) then
            write(*,2) 'do1_tdc_L_eqn scale', k, scale
            stop 'do1_tdc_L_eqn'
         end if
         resid = (L_expected - L_actual)*scale         
         residual = resid%val
         s% equ(s% i_equL, k) = residual         
         if (test_partials) then
            s% solver_test_partials_val = residual 
            if (.false.) then
               write(*,2) 'L_expected%val', k, L_expected%val
               write(*,2) 'L_actual%val', k, L_actual%val
               write(*,2) 'diff', k, L_expected%val - L_actual%val
               write(*,2) 'diff*scale', k, (L_expected%val - L_actual%val)*scale
               write(*,2) 's% Lc(k)', k, s% Lc(k)
               write(*,2) 's% w(k)', k, s% w(k)
               write(*,2) 's% w(k-1)', k-1, s% w(k-1)
            end if
         end if
         
         if (skip_partials) return
         call save_eqn_residual_info(s, k, nvar, s% i_equL, resid, 'do1_tdc_L_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_w
            s% solver_test_partials_dval_dx = resid%d1Array(i_w_00)
            write(*,*) 'do1_tdc_L_eqn', s% solver_test_partials_var, k, s% nz
         end if      
      end subroutine do1_tdc_L_eqn
      
      
      function compute_Hp_cell(s, k, ierr) result(Hp_cell) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell
         type(auto_diff_real_star_order1) :: r_mid, r_00, r_p1, &
            P_00, d_00, P_m1, d_m1, alt_Hp_cell, alfa
         real(dp) :: cgrav_00, cgrav_p1, cgrav_mid, m_00, m_p1, m_mid
         include 'formats'
         ierr = 0
         r_00 = wrap_opt_time_center_r_00(s, k)
         cgrav_00 = s% cgrav(k)
         m_00 = s% m(k)
         d_00 = wrap_d_00(s, k)
         P_00 = wrap_P_00(s, k)
         r_p1 = wrap_opt_time_center_r_p1(s, k)
         if (k < s% nz) then
            cgrav_p1 = s% cgrav(k+1)
            m_p1 = s% m(k+1)
         else
            cgrav_p1 = s% cgrav(k)
            m_p1 = s% m_center
         end if
         cgrav_mid = 0.5d0*(cgrav_00 + cgrav_p1)
         m_mid = 0.5d0*(m_00 + m_p1)
         r_mid = 0.5d0*(r_00 + r_p1)
         Hp_cell = pow2(r_mid)*P_00 / (d_00*cgrav_mid*m_mid)
         if (s% alt_scale_height_flag) then
            ! consider sound speed*hydro time scale as an alternative scale height
            ! (this comes from Eggleton's code.)
            alt_Hp_cell = sqrt(P_00/cgrav_mid)/d_00
            if (alt_Hp_cell%val < Hp_cell%val) then ! blend
               alfa = pow2(alt_Hp_cell/Hp_cell) ! 0 <= alfa%val < 1
               Hp_cell = alfa*Hp_cell + (1d0 - alfa)*alt_Hp_cell
            end if
         end if
      end function compute_Hp_cell
      
      
      function compute_Hp_face(s, k, ierr) result(Hp_face) ! cm
         ! no time_centering for Hp
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face
         type(auto_diff_real_star_order1) :: &
            r_00, P_00, d_00, P_m1, d_m1, P_div_rho, &
            d_face, P_face, alt_Hp_face, alfa
         integer :: j
         include 'formats'
         ierr = 0
         r_00 = wrap_opt_time_center_r_00(s, k)
         d_00 = wrap_d_00(s, k)
         P_00 = wrap_p_00(s, k)
         if (k > 1) then
            d_m1 = wrap_d_m1(s, k)
            P_m1 = wrap_P_m1(s, k)
            P_div_rho = 0.5d0*(P_00/d_00 + P_m1/d_m1)
            Hp_face = pow2(r_00)*P_div_rho/(s% cgrav(k)*s% m(k))
            if (s% alt_scale_height_flag) then
               ! consider sound speed*hydro time scale as an alternative scale height
               ! (this comes from Eggleton's code.)
               d_face = 0.5d0*(d_00 + d_m1)
               P_face = 0.5d0*(P_00 + P_m1)
               alt_Hp_face = sqrt(P_face/s% cgrav(k))/d_face
               if (alt_Hp_face%val < Hp_face%val) then ! blend
                  alfa = pow2(alt_Hp_face/Hp_face) ! 0 <= alfa%val < 1
                  Hp_face = alfa*Hp_face + (1d0 - alfa)*alt_Hp_face
               end if
            end if
         else ! surface
            P_div_rho = P_00/d_00
            Hp_face = pow2(r_00)*P_div_rho/(s% cgrav(k)*s% m(k))
         end if
      end function compute_Hp_face
      

      subroutine do1_turbulent_energy_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: set_energy_eqn_scal, save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         real(dp) :: dt, dm, dt_div_dm, scal, residual
         integer :: i_detrb_dt, i_v
         type(auto_diff_real_star_order1) :: resid_ad, &
            d_turbulent_energy_ad, PtdV_ad, dt_dLt_dm_ad, dt_C_ad, dt_Eq_ad
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         logical :: non_turbulent_cell, test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         call init
         
         non_turbulent_cell = (s% TDC_alfa == 0d0 .or. & ! TDC_alfa == 0d0 means purely radiative
             k == 1 .or. k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent)
         
         if (non_turbulent_cell) then
             
            resid_ad = wrap_w_00(s,k) - s% TDC_w_min_for_damping
            
         else
         
            call setup_d_turbulent_energy(ierr); if (ierr /= 0) return ! erg g^-1 = cm^2 s^-2
            call setup_PtdV_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_C_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_Eq_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)
         
            ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
            esum_ad = d_turbulent_energy_ad + PtdV_ad + dt_dLt_dm_ad - dt_C_ad - dt_Eq_ad ! erg g^-1
         
            resid_ad = esum_ad ! convert back to auto_diff_real_star_order1
            resid_ad = resid_ad*scal/s% dt ! to make residual unitless, must cancel out the dt in scal
         
         end if

         residual = resid_ad%val
         s% equ(i_detrb_dt, k) = residual

         if (is_bad(residual)) then
!$omp critical (hydro_equ_turbulent_crit1)
            write(*,2) 'turbulent energy eqn residual', k, residual
            write(*,2) 'det', k, d_turbulent_energy_ad%val
            write(*,2) 'PtdV', k, PtdV_ad%val
            write(*,2) 'dt_dLt_dm', k, dt_dLt_dm_ad%val
            write(*,2) 'dt_C', k, dt_C_ad%val
            write(*,2) 'dt_Eq', k, dt_Eq_ad%val
            stop 'get1_turbulent_energy_eqn'
!$omp end critical (hydro_equ_turbulent_crit1)
         end if

         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         
         if (skip_partials) return
         call save_eqn_residual_info(s, k, nvar, i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var
         end if      

         contains
         
         subroutine init
            i_detrb_dt = s% i_detrb_dt
            i_v = s% i_v
            dt = s% dt
            dm = s% dm(k)
            dt_div_dm = dt/dm
         end subroutine init
         
         subroutine setup_d_turbulent_energy(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: w_00
            ierr = 0
            w_00 = wrap_w_00(s,k) ! wrap_dxh_w = w_00 - s% w_start(k)
            d_turbulent_energy_ad = wrap_dxh_w(s,k)*(w_00 + s% w_start(k))
         end subroutine setup_d_turbulent_energy
         
         ! PtdV_ad = Pt_ad*dV_ad
         subroutine setup_PtdV_ad(ierr) ! erg g^-1
            use star_utils, only: calc_Pt_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Pt_ad, dV_ad, d_00
            call calc_Pt_ad_tw(s, k, Pt_ad, ierr)
            if (ierr /= 0) return
            d_00 = wrap_d_00(s,k)
            dV_ad = 1d0/d_00 - 1d0/s% rho_start(k)
            PtdV_ad = Pt_ad*dV_ad ! erg cm^-3 cm^-3 g^-1 = erg g^-1
         end subroutine setup_PtdV_ad

         subroutine setup_dt_dLt_dm_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr            
            type(auto_diff_real_star_order1) :: Lt_00, Lt_p1, dLt_ad
            real(dp) :: Lt_00_start, Lt_p1_start
            include 'formats'
            ierr = 0
            Lt_00 = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            if (k < s% nz) then
               Lt_p1 = compute_Lt(s, k+1, ierr)
               if (ierr /= 0) return
               Lt_p1 = shift_p1(Lt_p1)
            else
               Lt_p1 = 0d0
               Lt_p1_start = 0d0
            end if
            if (s% using_velocity_time_centering .and. &
                  s% include_L_in_velocity_time_centering) then
               Lt_00_start = s% Lt_start(k)
               if (k < s% nz) then
                  Lt_p1_start = s% Lt_start(k+1)
               else
                  Lt_p1_start = 0d0
               end if
               dLt_ad = 0.5d0*(Lt_00 + Lt_00_start) - 0.5d0*(Lt_p1 + Lt_p1_start)
            else
               dLt_ad = Lt_00 - Lt_p1
            end if
            dt_dLt_dm_ad = dt*dLt_ad/dm
         end subroutine setup_dt_dLt_dm_ad
         
         subroutine setup_dt_C_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: C
            C = compute_C(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_C_ad = dt*C
         end subroutine setup_dt_C_ad
                  
         subroutine setup_dt_Eq_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Eq_cell
            Eq_cell = compute_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_Eq_ad = dt*Eq_cell
         end subroutine setup_dt_Eq_ad
      
      end subroutine do1_turbulent_energy_eqn

      
      function compute_Y_face(s, k, ierr) result(Y_face) ! superadiabatic gradient [unitless]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face
         type(auto_diff_real_star_order1) :: Hp_face, Y1, Y2, QQ_div_Cp_face, &
            r_00, d_00, P_00, Cp_00, T_00, chiT_00, chiRho_00, QQ_00, lnT_00, &
            r_m1, d_m1, P_m1, Cp_m1, T_m1, chiT_m1, chiRho_m1, QQ_m1, lnT_m1
         real(dp) :: dm_bar
         include 'formats'
         ierr = 0
         if (k == 1 .or. s% TDC_alfa == 0d0) then
            Y_face = 0d0
            s% Y_face(k) = 0d0
            return
         end if
         dm_bar = s% dm_bar(k)
         Hp_face = compute_Hp_face(s,k,ierr)
         if (ierr /= 0) return
         r_00 = wrap_opt_time_center_r_00(s, k)
         d_00 = wrap_d_00(s, k)
         P_00 = wrap_P_00(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         T_00 = wrap_T_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
         lnT_00 = wrap_lnT_00(s,k)
         
         r_m1 = wrap_opt_time_center_r_m1(s, k)
         d_m1 = wrap_d_m1(s, k)
         P_m1 = wrap_P_m1(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         T_m1 = wrap_T_m1(s, k)
         chiT_m1 = wrap_chiT_m1(s, k)
         chiRho_m1 = wrap_chiRho_m1(s, k)
         QQ_m1 = chiT_m1/(d_m1*T_m1*chiRho_m1)
         lnT_m1 = wrap_lnT_m1(s,k)
         QQ_div_Cp_face = 0.5d0*(QQ_00/Cp_00 + QQ_m1/Cp_m1)
         ! QQ units (g cm^-3 K)^-1 = g^-1 cm^3 K^-1
         ! Cp units erg g^-1 K^-1 = g cm^2 s^-2 g^-1 K^-1 = cm^2 s^-2 K^-1
         ! QQ/Cp units = (g^-1 cm^3 K^-1)/(cm^2 s^-2 K^-1)
         !  = g^-1 cm^3 K^-1 cm^-2 s^2 K
         !  = g^-1 cm s^2
         ! P units = erg cm^-3 = g cm^2 s^-2 cm^-3 = g cm^-1 s^-2
         ! QQ/Cp*P is unitless.
         
         Y1 = QQ_div_Cp_face*(P_m1 - P_00) - (lnT_m1 - lnT_00)
         ! Y1 unitless
         
         Y2 = 4d0*pi*pow2(r_00)*Hp_face*2d0/(1/d_00 + 1/d_m1)/dm_bar
         ! units = cm^2 cm / (cm^3 g^-1) / g
         !       = cm^2 cm cm^-3 g g^-1 = unitless
         
         Y_face = Y1*Y2 ! unitless
         s% Y_face(k) = Y_face%val

      end function compute_Y_face
      
      
      function compute_PII_face(s, k, ierr) result(PII_face) ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Cp_00, Cp_m1, Cp_face, Y_face
         real(dp) :: ALFAS, ALFA
         include 'formats'
         ierr = 0
         if (k == 1 .or. k == s% nz .or. s% TDC_alfa == 0d0) then
            PII_face = 0d0
            s% PII(k) = 0d0
            return
         end if
         Y_face = compute_Y_face(s, k, ierr)
         if (ierr /= 0) return
         Cp_00 = wrap_Cp_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         Cp_face = 0.5d0*(Cp_00 + Cp_m1) ! ergs g^-1 K^-1
         ALFAS = x_ALFAS
         ALFA = s% TDC_alfa
         PII_face = ALFAS*ALFA*Cp_face*Y_face
         s% PII(k) = PII_face%val
      end function compute_PII_face
      
      
      function compute_d_v_div_r(s, k, ierr) result(d_v_div_r) ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         r_00 = wrap_opt_time_center_r_00(s,k)
         v_00 = wrap_opt_time_center_v_00(s,k)
         r_p1 = wrap_opt_time_center_r_p1(s,k)
         v_p1 = wrap_opt_time_center_v_p1(s,k)
         d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1
      end function compute_d_v_div_r
      
      
      function compute_Chi_cell(s, k, ierr) result(Chi_cell) ! eddy viscosity energy (Kuhfuss 1986) [erg]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_rho2, r6_cell, d_v_div_r, Hp_cell, w_00, d_00, r_00, r_p1
         real(dp) :: f, ALFAM_ALFA
         include 'formats'
         ierr = 0
         ALFAM_ALFA = s% TDC_alfam*s% TDC_alfa
         Hp_cell = compute_Hp_cell(s, k, ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r(s, k, ierr)
         if (ierr /= 0) return
         w_00 = wrap_w_00(s,k)
         d_00 = wrap_d_00(s,k)
         f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
         w_rho2 = w_00*pow2(d_00)
         r_00 = wrap_opt_time_center_r_00(s,k)
         r_p1 = wrap_opt_time_center_r_p1(s,k)
         r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
         Chi_cell = f*w_rho2*r6_cell*d_v_div_r*Hp_cell
         ! units = g^-1 cm s^-1 g^2 cm^-6 cm^6 s^-1 cm
         !       = g cm^2 s^-2
         !       = erg
         s% Chi(k) = Chi_cell%val
         
         !s% xtra1_array(k) = f
         !s% xtra2_array(k) = safe_log10(w_rho2%val)
         !s% xtra3_array(k) = safe_log10(r6_cell%val)
         !s% xtra4_array(k) = d_v_div_r%val
         !s% xtra5_array(k) = Hp_cell%val

      end function compute_Chi_cell

      
      function compute_Eq_cell(s, k, ierr) result(Eq_cell) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: d_v_div_r, Chi_cell
         include 'formats'
         ierr = 0
         Chi_cell = compute_Chi_cell(s,k,ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r(s, k, ierr)
         if (ierr /= 0) return
         Eq_cell = 4d0*pi*Chi_cell*d_v_div_r/s% dm(k) ! erg s^-1 g^-1
         !s% xtra4_array(k) = 4d0*pi*d_v_div_r%val/s% dm(k)
         !s% xtra5_array(k) = Chi_cell%val
         s% Eq(k) = Eq_cell%val
      end function compute_Eq_cell


      function compute_Uq_face(s, k, ierr) result(Uq_face) ! cm s^-2, acceleration
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Uq_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_00, Chi_out, r_00
         include 'formats'
         ierr = 0         
         Chi_00 = compute_Chi_cell(s,k,ierr)
         if (k > 1) then
            Chi_out = compute_Chi_cell(s,k-1,ierr)
            Chi_out = shift_m1(Chi_out)
            if (ierr /= 0) return
         else
            Chi_out = 0d0
         end if
         r_00 = wrap_opt_time_center_r_00(s,k)
         Uq_face = 4d0*pi*(Chi_out - Chi_00)/(s% dm_bar(k)*r_00)   
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
         s% Uq(k) = Uq_face%val
      end function compute_Uq_face


      function compute_Source(s, k, ierr) result(Source) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, P_00, Cp_00, chiT_00, chiRho_00, QQ_00, &
            Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell
         include 'formats'
         ierr = 0
         w_00 = wrap_w_00(s, k)
         T_00 = wrap_T_00(s, k)                  
         d_00 = wrap_d_00(s, k)         
         P_00 = wrap_P_00(s, k)         
         Cp_00 = wrap_Cp_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            
         Hp_face_00 = compute_Hp_face(s,k,ierr)
         if (ierr /= 0) return
         Hp_face_p1 = compute_Hp_face(s,k+1,ierr)
         if (ierr /= 0) return
         Hp_face_p1 = shift_p1(Hp_face_p1)

         PII_face_00 = compute_PII_face(s, k, ierr)
         if (ierr /= 0) return
         PII_face_p1 = compute_PII_face(s, k+1, ierr)
         if (ierr /= 0) return
         PII_face_p1 = shift_p1(PII_face_p1)

         PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         Source = PII_div_Hp_cell*w_00*T_00*P_00*QQ_00/Cp_00  
         
         ! PII units same as Cp = erg g^-1 K^-1
         ! P*QQ/Cp is unitless (see Y_face)
         ! Source units = (erg g^-1 K^-1) cm^-1 cm s^-1 K
         !     = erg g^-1 s^-1

         s% SOURCE(k) = Source%val

         if (is_bad(Source%val)) then
            !$omp critical (hydro_equ_turbulent_crit3)
            write(*,2) 'Source', k, Source%val
            write(*,2) 'Hp_face_00', k, Hp_face_00%val
            write(*,2) 'Hp_face_p1', k, Hp_face_p1%val
            write(*,2) 'PII_face_00', k, PII_face_00%val
            write(*,2) 'PII_face_p1', k, PII_face_p1%val
            write(*,2) 'Cp_00', k, Cp_00%val
            write(*,2) 'QQ_00', k, QQ_00%val
            write(*,2) 'P_00', k, P_00%val
            write(*,2) 'T_00', k, T_00%val
            write(*,2) 'w_00', k, w_00%val
            stop 'compute_Source'
            !$omp end critical (hydro_equ_turbulent_crit3)
         end if

      end function compute_Source


      function compute_D(s, k, ierr) result(D) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell, w_00, dw3
         include 'formats'
         real(dp) :: alpha
         ierr = 0
         alpha = s% TDC_alfa
         Hp_cell = compute_Hp_cell(s,k,ierr)
         if (ierr /= 0) return
         w_00 = wrap_w_00(s,k)
         dw3 = pow3(w_00 - s% TDC_w_min_for_damping) ! keep nonzero partial even when w = 0
         D = (x_CEDE/alpha)*dw3/Hp_cell
         ! units cm^3 s^-3 cm^-1 = cm^2 s^-3 = erg g^-1 s^-1
         s% DAMP(k) = D%val
      end function compute_D


      function compute_Dr(s, k, ierr) result(Dr) ! erg g^-1 s^-1 = cm^2 s^-3
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Cp_00, kap_00, Hp_cell, POM2
         real(dp) :: gammar, alpha, POM
         include 'formats'
         ierr = 0
         alpha = s% TDC_alfa
         gammar = s% TDC_alfar*x_GAMMAR
         if (gammar == 0d0) then
            Dr = 0d0
            s% DAMPR(k) = 0d0
            return
         end if
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Hp_cell = compute_Hp_cell(s,k,ierr)
         if (ierr /= 0) return
         POM = 4d0*boltz_sigma*(gammar/alpha)**2 ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00) 
            ! K^3 / ((g cm^-3)^2 (erg g^-1 K^-1) (cm^2 g^-1))
            ! K^3 / (cm^-4 erg K^-1) = K^4 cm^4 erg^-1
         Dr = POM*POM2*pow2(w_00)/pow2(Hp_cell)
         ! (erg cm^-2 K^-4 s^-1) (K^4 cm^4 erg^-1) cm^2 s^-2 cm^-2
         ! cm^2 s^-3 = erg g^-1 s^-1
         s% DAMPR(k) = Dr%val
      end function compute_Dr


      function compute_C(s, k, ierr) result(C) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: C
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source, D, Dr
         Source = compute_Source(s, k, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         s% COUPL(k) = C%val
      end function compute_C


      function compute_L_face(s, k, ierr) result(L_face) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr, Lc, Lt
         call compute_L_terms(s, k, L_face, Lr, Lc, Lt, ierr)
         return
      end function compute_L_face


      function compute_Lr(s, k, ierr) result(Lr) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_00, T400, Erad, T_m1, T4m1, &
            kap_00, kap_m1, kap_face, diff_T4_div_kap, BW, BK
         real(dp) :: alfa
         include 'formats'
         ierr = 0
         r_00 = wrap_r_00(s,k) ! not time centered for luminosity
         area = 4d0*pi*pow2(r_00)
         T_00 = wrap_T_00(s,k)
         T400 = pow4(T_00)
         if (k == 1) then ! Lr(1) proportional to Erad in cell(1)
            Erad = crad * T400
            Lr = s% TDC_Lsurf_factor * area * clight * Erad
            s% Lr(k) = Lr%val
            return
         end if
         T_m1 = wrap_T_m1(s,k)
         T4m1 = pow4(T_m1)            
         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         kap_00 = wrap_kap_00(s,k)
         kap_m1 = wrap_kap_m1(s,k)
         kap_face = alfa*kap_00 + (1d0 - alfa)*kap_m1
         diff_T4_div_kap = (T4m1 - T400)/kap_face
         if (s% TDC_use_Stellingwerf_Lr) then ! RSP style
            BW = log(T4m1/T400)
            if (abs(BW%val) > 1d-20) then
               BK = log(kap_m1/kap_00)
               if (abs(1d0 - BK%val/BW%val) > 1d-15 .and. abs(BW%val - BK%val) > 1d-15) then
                  diff_T4_div_kap = (T4m1/kap_m1 - T400/kap_00)/(1d0 - BK/BW)
               end if
            end if
         end if
         Lr = -crad*clight/3d0*diff_T4_div_kap*pow2(area)/s% dm_bar(k)       
         ! units (erg cm^-3 K^-4) (cm s^-1) (K^4 cm^-2 g cm^4) g^-1 = erg s^-1  
         s% Lr(k) = Lr%val
      end function compute_Lr


      function compute_Lc(s, k, ierr) result(Lc) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lc_w_face_factor
         Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc


      function compute_Lc_terms(s, k, Lc_w_face_factor, ierr) result(Lc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc, Lc_w_face_factor
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area, &
            T_m1, T_00, d_m1, d_00, w_m1, w_00, T_rho_face, PII_face, w_face
         real(dp) :: ALFAC, ALFAS
         include 'formats'
         ierr = 0
         if (k == 1) then
            Lc = 0d0
            Lc_w_face_factor = 1
            return
         end if
         r_00 = wrap_r_00(s, k) ! not time centered for luminosity
         area = 4d0*pi*pow2(r_00)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)         
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         w_m1 = wrap_w_m1(s, k)
         w_00 = wrap_w_00(s, k)
         T_rho_face = 0.5d0*(T_m1*d_m1 + T_00*d_00)
         PII_face = compute_PII_face(s, k, ierr)
         w_face = 0.5d0*(w_m1 + w_00)
         ALFAC = x_ALFAC
         ALFAS = x_ALFAS
         Lc_w_face_factor = area*(ALFAC/ALFAS)*T_rho_face*PII_face
         ! units = cm^2 K g cm^-3 ergs g^-1 K^-1 = ergs cm^-1
         Lc = w_face*Lc_w_face_factor
         ! units = cm s^-1 ergs cm^-1 = ergs s^-1
      end function compute_Lc_terms


      function compute_Lt(s, k, ierr) result(Lt) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lt
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area2, d_m1, d_00, rho2_face, Hp_face, w_m1, w_00
         real(dp) :: alpha, alpha_t
         include 'formats'
         ierr = 0
         if (k == 1) then
            Lt = 0d0
            return
         end if
         alpha = s% TDC_alfa
         alpha_t = s% TDC_alfat
         r_00 = wrap_r_00(s,k) ! not time centered for luminosity         
         area2 = (4d0*pi)**2*pow4(r_00)
         d_m1 = wrap_d_m1(s,k)
         d_00 = wrap_d_00(s,k)
         rho2_face = 0.5d0*(pow2(d_00) + pow2(d_m1))
         w_m1 = wrap_w_m1(s,k)
         w_00 = wrap_w_00(s,k)
         Hp_face = compute_Hp_face(s,k,ierr)
         if (ierr /= 0) return
         Lt = -2d0/3d0*alpha*alpha_t * area2 * Hp_face * rho2_face * &
            (pow3(w_m1) - pow3(w_00))/s% dm_bar(k)         
         ! units = cm^4 cm g^2 cm^-6 cm^3 s^-3 g^-1 = g cm^2 s^-3 = erg s^-1
         s% Lt(k) = Lt%val
      end function compute_Lt


      subroutine compute_L_terms(s, k, L, Lr, Lc, Lt, ierr)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L, Lr, Lc, Lt
         integer, intent(out) :: ierr         
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            L = 0d0
            L%val = s% L_center
            Lr = 0d0
            Lc = 0d0
            Lt = 0d0
            return
         end if
         Lr = compute_Lr(s, k, ierr)
         if (ierr /= 0) return
         if (k == 1) then            
            Lc = 0d0
            Lt = 0d0
         else
            Lc = compute_Lc(s, k, ierr)
            if (ierr /= 0) return
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
         end if
         L = Lr + Lc + Lt
         if (abs(Lt%val)/max(1d-99,abs(L%val)) > 1d-2) then
            s% mixing_type(k) = overshoot_mixing
         else if (abs(Lc%val)/max(1d-99,abs(L%val)) > 1d-2) then
            s% mixing_type(k) = convective_mixing
         else
            s% mixing_type(k) = no_mixing
         end if
      end subroutine compute_L_terms


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, op_err
         type(auto_diff_real_star_order1) :: Lt
         include 'formats'
         ierr = 0
         do k=1,s%nz
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_start(k) = Lt%val  
            s% w_start(k) = s% w(k)
         end do         
      end subroutine set_etrb_start_vars
      
      
      subroutine reset_etrb_using_L(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr   
         integer :: k, nz
         real(dp) :: Lc_val, w_00
         type(auto_diff_real_star_order1) :: &
            Lc_w_face_factor, L, Lr, Lc, Lt
         real(dp), allocatable :: w_face(:)
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         if (s% TDC_alfa == 0d0) return ! no convection
         nz = s% nz
         allocate(w_face(nz))
         w_face(1) = 0d0
         do k=2, nz
            Lr = compute_Lr(s, k, ierr)
            if (ierr /= 0) stop 'failed in compute_Lr'
            Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
            if (ierr /= 0) stop 'failed in compute_Lc_terms'
            Lc_val = s% L(k) - Lr%val ! assume Lt = 0 for this
            if (abs(Lc_w_face_factor%val) < 1d-20) then
               w_face(k) = 0d0
            else
               w_face(k) = Lc_val/Lc_w_face_factor%val
            end if
         end do
         do k=1, nz
            if (k < nz) then
               w_00 = max(0.5d0*(w_face(k) + w_face(k+1)), 0d0)
            else
               w_00 = max(w_face(k), 0d0)
            end if
            s% w(k) = w_00
            s% xh(s%i_w,k) = s% w(k)
            call compute_L_terms(s, k, L, Lr, Lc, Lt, ierr) ! redo with new w(k)
            if (ierr /= 0) stop 'failed in compute_L reset_wturb_using_L'
         end do
         if (dbg) stop 'reset_etrb_using_L'
      end subroutine reset_etrb_using_L


      end module hydro_tdc

