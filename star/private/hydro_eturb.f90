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

      module hydro_eturb

      use star_private_def
      use const_def
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use star_utils, only: em1, e00, ep1

      implicit none

      private
      public :: do1_turbulent_energy_eqn, do1_Eturb_L_eqn, &
         set_eturb_start_vars, reset_eturb_using_L, calc_Eq_18, calc_Uq_18 
      

      contains
      

      subroutine do1_Eturb_L_eqn( &
            s, k, xscale, equ, L_start_max, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         real(dp), intent(in) :: L_start_max
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         real(dp), dimension(nvar) :: d_dm1, d_d00, d_dp1      
         include 'formats'
         call get1_Eturb_L_eqn( &
            s, k, xscale, equ, L_start_max, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_Eturb_L_eqn', k
            return
         end if         
         if (skip_partials) return         
         call store_partials(s, k, xscale, s% i_equL, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do1_Eturb_L_eqn

      
      subroutine get1_Eturb_L_eqn( &  
            s, k, xscale, equ, L_start_max, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         use star_utils, only: unpack_res18_partials
         use accurate_sum_auto_diff_18var_order1
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         real(dp), intent(in) :: L_start_max
         logical, intent(in) :: skip_partials
         real(dp), dimension(nvar), intent(out) :: d_dm1, d_d00, d_dp1      
         integer, intent(out) :: ierr
         type(auto_diff_real_18var_order1) :: L, Lr, Lc, Lt, L_actual, res18
         real(dp) :: scale, residual, e_avg
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         call compute_L(s, k, L, Lr, Lc, Lt, ierr)         
         if (ierr /= 0) return        
         L_actual = wrap_L_00(s, k)  
         scale = 1d0/L_start_max
         res18 = (L - L_actual)*scale         
         residual = res18%val
         equ(s% i_equL, k) = residual
         
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         if (skip_partials) return
         call unpack_res18_partials(s, k, nvar, xscale, s% i_equL, &
            res18, d_dm1, d_d00, d_dp1)

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)
            write(*,*) 'get1_Eturb_L_eqn', s% solver_test_partials_var
         end if      
      end subroutine get1_Eturb_L_eqn
      

      subroutine do1_turbulent_energy_eqn( &
            s, k, xscale, equ, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         real(dp), dimension(nvar) :: d_dm1, d_d00, d_dp1      
         include 'formats'
         call get1_turbulent_energy_eqn( &
            s, k, xscale, equ, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_turbulent_energy_eqn', k
            return
         end if         
         if (skip_partials) return         
         call store_partials(s, k, xscale, s% i_deturb_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do1_turbulent_energy_eqn

      
      subroutine get1_turbulent_energy_eqn( &  
            s, k, xscale, equ, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         use star_utils, only: calc_Pt_18_tw
         use accurate_sum_auto_diff_18var_order1
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         real(dp), dimension(nvar), intent(out) :: d_dm1, d_d00, d_dp1      
         integer, intent(out) :: ierr
         
         real(dp) :: dt, dm, dt_div_dm, residual, scale
         integer :: i_deturb_dt, i_eturb, i_lnd, i_lnT, i_lnR, i_v
         type(auto_diff_real_18var_order1) :: resid_18, &
            dEturb_18, PtdV_18, dt_dLt_dm_18, dt_C_18, dt_Eq_18
         type(accurate_auto_diff_real_18var_order1) :: esum_18
         real(dp) :: dm_m1, dm_00, dm_p1, m_00, cgrav_00
         type(auto_diff_real_18var_order1) :: Source, D, Dr, g_00, g_p1, area_00, area_p1, g_cell_00, h_00
         type(auto_diff_real_18var_order1) :: d_m1, d_00, d_p1, P_m1, P_00, P_p1
         type(auto_diff_real_18var_order1) :: eturb_m1, eturb_00, eturb_p1, T_m1, T_00, kap_m1, kap_00
         type(auto_diff_real_18var_order1) :: r_m1, r_00, r_p1, v_m1, v_00, v_p1
         type(auto_diff_real_18var_order1) :: ChiRho_00, ChiT_00, Cp_00, s_m1, s_00, s_p1
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         call init
         
         if (s% Eturb_alfa == 0d0) then
            equ(i_deturb_dt, k) = s% Eturb(k) - min_Eturb
            if (.not. skip_partials) d_d00(i_eturb) = 1d0
            return
         end if
         
         call setup_intermediates
         call setup_dEturb_18(ierr); if (ierr /= 0) return         
         call setup_PtdV_18(ierr); if (ierr /= 0) return         
         call setup_dt_dLt_dm_18(ierr); if (ierr /= 0) return         
         call setup_dt_C_18(ierr); if (ierr /= 0) return         
         call setup_dt_Eq_18(ierr); if (ierr /= 0) return    
         
         ! sum terms in esum_18 using accurate_auto_diff_real_18var_order1
         esum_18 = dEturb_18 + PtdV_18 + dt_dLt_dm_18 - dt_C_18 - dt_Eq_18
         
         resid_18 = esum_18 ! convert back to auto_diff_real_18var_order1
         scale = 1d0/s% energy_start(k)
         resid_18 = scale*resid_18
         residual = resid_18%val
         equ(i_deturb_dt, k) = residual

         if (is_bad(residual)) then
!$omp critical (hydro_equ_turbulent_crit1)
            write(*,2) 'turbulent energy eqn residual', k, residual
            write(*,2) 'dEturb', k, dEturb_18%val
            write(*,2) 'PtdV', k, PtdV_18%val
            write(*,2) 'dt_dLt_dm', k, dt_dLt_dm_18%val
            write(*,2) 'dt_C', k, dt_C_18%val
            write(*,2) 'dt_Eq', k, dt_Eq_18%val
            stop 'get1_turbulent_energy_eqn'
!$omp end critical (hydro_equ_turbulent_crit1)
         end if
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         if (skip_partials) return
         call unpack_res18(resid_18)

         if (test_partials) then
            s% solver_test_partials_var = i_lnT
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)
            write(*,*) 'get1_turbulent_energy_eqn', s% solver_test_partials_var, &
               s% Eturb(k), s% Eturb_start(k)
         end if      

         contains
         
         subroutine init
            i_deturb_dt = s% i_deturb_dt
            i_eturb = s% i_eturb
            i_lnd = s% i_lnd
            i_lnT = s% i_lnT
            i_lnR = s% i_lnR
            i_v = s% i_v
            dt = s% dt
            dm = s% dm(k)
            dt_div_dm = dt/dm
            d_dm1 = 0d0; d_d00 = 0d0; d_dp1 = 0d0
         end subroutine init
         
         subroutine setup_intermediates         
            if (k > 1) then
               dm_m1 = s%dm(k-1)
            else
               dm_m1 = 0d0
            end if
            dm_00 = s%dm(k)
            d_m1 = wrap_d_m1(s, k)
            d_00 = wrap_d_00(s, k)
            d_p1 = wrap_d_p1(s, k) ! Gradients vanish at the center, so d(nz+1) == d(nz).
            P_m1 = wrap_P_m1(s, k)
            P_00 = wrap_P_00(s, k)
            P_p1 = wrap_P_p1(s, k)
            s_m1 = wrap_s_m1(s, k)
            s_00 = wrap_s_00(s, k)
            s_p1 = wrap_s_p1(s, k) ! Gradients vanish at the center, so s(nz+1) == s(nz).
            eturb_m1 = wrap_eturb_m1(s, k)
            eturb_00 = wrap_eturb_00(s, k)
            eturb_p1 = wrap_eturb_p1(s, k)
            T_m1 = wrap_T_m1(s, k)
            T_00 = wrap_T_00(s, k)
            kap_m1 = wrap_kap_m1(s, k)
            kap_00 = wrap_kap_00(s, k)
            ChiRho_00 = wrap_ChiRho_00(s, k)
            ChiT_00 = wrap_ChiT_00(s, k)
            Cp_00 = wrap_Cp_00(s, k)
            r_m1 = wrap_r_m1(s, k)
            r_00 = wrap_r_00(s, k)
            r_p1 = wrap_r_p1(s, k) ! Set by wrap routine to r_center when k == nz.
            v_m1 = wrap_v_m1(s, k)
            v_00 = wrap_v_00(s, k)
            v_p1 = wrap_v_p1(s, k) ! Set by wrap routine to zero when k == nz.
            ! Compute areas on faces (used by L_turb)
            area_00 = 4d0 * pi * pow2(r_00)
            area_p1 = 4d0 * pi * pow2(r_p1)
            ! Compute gravity on faces (used by L_turb)
            m_00 = s% m(k)
            cgrav_00 = s% cgrav(k)
            g_00 = m_00 * cgrav_00 / pow2(r_00)
            if (r_p1 > 0d0) then
               g_p1 = (m_00 - dm_00) * cgrav_00 /  pow2(r_p1)
            else
               g_p1 = 0d0
            end if            
            g_cell_00 = (m_00 - 0.5d0 * dm_00) * cgrav_00 / (0.5d0 * (r_00 + r_p1))
         end subroutine setup_intermediates
         
         subroutine setup_dEturb_18(ierr)
            integer, intent(out) :: ierr
            dEturb_18 = 0d0
            dEturb_18%val = s% dxh_eturb(k) ! Eturb = Eturb_start + dxh_eturb
            dEturb_18%d1Array(i_eturb_00) = 1d0
         end subroutine setup_dEturb_18
         
         ! PtdV_18 = Pt_18*dV_18
         subroutine setup_PtdV_18(ierr)
            use star_utils, only: calc_Pt_18_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: Pt_18, dV_18
            real(dp) :: DV, d_DV_dlnd
            call calc_Pt_18_tw(s, k, Pt_18, ierr)
            if (ierr /= 0) return
            dV_18 = 1d0/d_00 - 1d0/s% rho_start(k)
            PtdV_18 = Pt_18*dV_18
         end subroutine setup_PtdV_18

         subroutine setup_dt_dLt_dm_18(ierr)
            integer, intent(out) :: ierr            
            type(auto_diff_real_18var_order1) :: Lt_00, Lt_p1, dLt_18
            real(dp) :: Lt_00_start, Lt_p1_start
            include 'formats'
            ierr = 0
            Lt_00 = compute_Lt(s, k, &
               dm_m1, dm_00, g_00, area_00, P_m1, P_00, d_m1, d_00, eturb_m1, eturb_00, ierr)
            if (ierr /= 0) return
            if (k < s% nz) then
               Lt_p1 = compute_Lt(s, k+1, &
                  dm_00, dm_p1, g_p1, area_p1, P_00, P_p1, d_00, d_p1, eturb_00, eturb_p1, ierr)
               if (ierr /= 0) return
            else
               Lt_p1 = 0d0
               Lt_p1_start = 0d0
            end if
            if (s% use_Fraley_time_centering .and. s% include_L_in_Fraley_time_centering) then
               Lt_00_start = s% Lt_start(k)
               if (k < s% nz) then
                  Lt_p1_start = s% Lt_start(k+1)
               else
                  Lt_p1_start = 0d0
               end if
               dLt_18 = 0.5d0*(Lt_00 + Lt_00_start) - 0.5d0*(Lt_p1 + Lt_p1_start)
            else
               dLt_18 = Lt_00 - Lt_p1
            end if
            dt_dLt_dm_18 = dt*dLt_18/dm
         end subroutine setup_dt_dLt_dm_18
         
         subroutine setup_dt_C_18(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: C
            C = compute_C(s, k, ierr)
            if (ierr /= 0) return
            dt_C_18 = dt*C
         end subroutine setup_dt_C_18
                  
         subroutine setup_dt_Eq_18(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: h_00, epsilon_q
            ! Compute scale height on cell k (used by epsilon_q)
            h_00 = compute_h_00(s, k, P_00, g_cell_00, d_00, ierr)
            if (ierr /= 0) return
            ! Compute epsilon_q
            epsilon_q = compute_epsilon_q(s, k, &
               v_p1, v_00, r_p1, r_00, d_00, dm_00, eturb_00, h_00, ierr)
            dt_Eq_18 = dt*epsilon_q
         end subroutine setup_dt_Eq_18

         subroutine unpack_res18(res18)
            use star_utils, only: unpack_res18_partials
            type(auto_diff_real_18var_order1) :: res18            
            call unpack_res18_partials(s, k, nvar, xscale, i_deturb_dt, &
               res18, d_dm1, d_d00, d_dp1)
         end subroutine unpack_res18
      
      end subroutine get1_turbulent_energy_eqn
      
      
      function compute_C(s, k, ierr) result(C)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         ! Outputs
         type(auto_diff_real_18var_order1) :: C
         integer, intent(out) :: ierr
         real(dp) :: dm_m1, dm_00, dm_p1, m_00, cgrav
         type(auto_diff_real_18var_order1) :: Source, D, Dr, &
            eturb_00, ChiT_00, ChiRho_00, Cp_00, s_m1, s_00, s_p1, &
            d_m1, d_00, d_p1, T_00, r_00, r_p1, h_00, kap_00, &
            g_cell_00, P_00
         include 'formats'
         ierr = 0
         
         if (k > 1) then
            dm_m1 = s%dm(k-1)
            s_m1 = wrap_s_m1(s, k)
            d_m1 = wrap_d_m1(s, k)
         else
            dm_m1 = 0d0
            s_m1 = 0d0
            d_m1 = 0d0
         end if
         if (k < s% nz) then
            dm_p1 = s%dm(k+1)
            s_p1 = wrap_s_p1(s, k)
            d_p1 = wrap_d_p1(s, k)
         else
            dm_p1 = 0d0
            s_p1 = 0d0
            d_p1 = 0d0
         end if
         
         r_p1 = wrap_r_p1(s, k)

         dm_00 = s%dm(k)
         eturb_00 = wrap_eturb_00(s, k)
         ChiT_00 = wrap_ChiT_00(s, k)
         ChiRho_00 = wrap_ChiRho_00(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         s_00 = wrap_s_00(s, k)
         d_00 = wrap_d_00(s, k)
         T_00 = wrap_T_00(s, k)
         r_00 = wrap_r_00(s, k)
         kap_00 = wrap_kap_00(s, k)
         
         m_00 = s% m(k)
         cgrav = s% cgrav(k)
         P_00 = wrap_P_00(s, k)
         g_cell_00 = (m_00 - 0.5d0 * dm_00) * cgrav / (0.5d0 * (r_00 + r_p1))
         h_00 = compute_h_00(s, k, P_00, g_cell_00, d_00, ierr)
         if (ierr /= 0) return
         
         Source = compute_Source(s, k, &
            eturb_00, d_00, d_p1, T_00, P_00, h_00, r_00, r_p1, ChiT_00, ChiRho_00, Cp_00, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, h_00, eturb_00, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, T_00, d_00, Cp_00, kap_00, h_00, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         
         s% COUPL(k) = C%val

         if (is_bad(C%val)) then
!$omp critical (hydro_equ_turbulent_crit2)
            write(*,2) 'C', k, C%val
            write(*,2) 'Source', k, Source%val
            write(*,2) 'D', k, D%val
            write(*,2) 'Dr', k, Dr%val
            stop 'compute_C'
!$omp end critical (hydro_equ_turbulent_crit2)
         end if
         
      end function compute_C


      !! Calculates the rate of turbulent energy generation from turbulent viscosity in cell k.
      !!
      !! @param v_p1 Velocity of face k+1
      !! @param v_00 Velocity of face k
      !! @param r_p1 Radius of face k+1
      !! @param r_00 Radius of face k
      !! @param d_00 Density of face k
      !! @param dm_00 Mass of cell k
      !! @param w_00 Turbulent velocity in cell k
      !! @param h_00 Pressure scale height in cell k
      function compute_epsilon_q_adam(s, k, v_p1, v_00, &
                                 r_p1, r_00, d_00, dm_00, w_00, h_00, ierr) result(epsilon_q)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_00
         type(auto_diff_real_18var_order1), intent(in) :: v_p1, v_00, r_p1, r_00, d_00, w_00, h_00

         ! Intermediates
         type(auto_diff_real_18var_order1) :: dr_00, d_v_div_r_dr_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: epsilon_q
         integer, intent(out) :: ierr
         
         real(dp) :: alpha, alpha_m

         ierr = 0
         alpha = s% Eturb_alfa
         alpha_m = s% Eturb_alfam

         dr_00 = dm_00 / (4d0 * pi * d_00 * (pow2(r_00) + r_00*r_p1 + pow2(r_p1)))
         d_v_div_r_dr_00 = (v_p1 / r_p1 - v_00 / r_00) / dr_00
         epsilon_q = (4d0 / 3d0) * alpha * alpha_m * h_00 * w_00 * pow2(d_v_div_r_dr_00)

         s% Eq(k) = epsilon_q%val

      end function compute_epsilon_q_adam

      function compute_epsilon_q(s, k, & ! following RSP's definition
            v_p1, v_00, r_p1, r_00, d_00, dm_00, eturb_00, h_00, ierr) result(epsilon_q)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_00
         type(auto_diff_real_18var_order1), intent(in) :: v_p1, v_00, r_p1, r_00, d_00, eturb_00, h_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: epsilon_q
         integer, intent(out) :: ierr

         type(auto_diff_real_18var_order1) :: w_rho2, r6_cell, d_v_div_r, Chi
         real(dp) :: alpha, alpha_m
         include 'formats'

         ierr = 0
         alpha = s% Eturb_alfa
         alpha_m = s% Eturb_alfam

         w_rho2 = sqrt(eturb_00)*pow2(d_00)
         r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
         if (k < s% nz) then
            d_v_div_r = v_00/r_00 - v_p1/r_p1
         else
            d_v_div_r = v_00/r_00
         end if
         Chi = (16d0/3d0)*pi*alpha*alpha_m*w_rho2*r6_cell*h_00*d_v_div_r/dm_00         
         epsilon_q = pi4*Chi*d_v_div_r/dm_00

         s% Eq(k) = epsilon_q%val

      end function compute_epsilon_q

      !! Calculates the radiative dissipation of turbulent energy in cell k.
      !!
      !! @param T_00 Temperature in cell k
      !! @param d_00 Density in cell k
      !! @param Cp_00 Specific heat in cell k
      !! @param kap_00 Opacity in cell k
      !! @param h_00 Pressure scale height in cell k
      function compute_Dr_adam(s, k, T_00, d_00, cp_00, kap_00, h_00, ierr) result(Dr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: T_00, d_00, cp_00, kap_00, h_00
         type(auto_diff_real_18var_order1) :: Dr
         integer, intent(out) :: ierr
         
         real(dp) :: gammar, alpha
         ierr = 0
         gammar = s% Eturb_alfar
         alpha = s% Eturb_alfa
         Dr = (4d0 * boltz_sigma * pow2(gammar) / pow2(alpha)) * (pow3(T_00) / (pow2(d_00) + cp_00 * kap_00 * h_00))

         s% DAMPR(k) = Dr%val

      end function compute_Dr_adam

      function compute_Dr(s, k, & ! following RSP's definition
            T_00, d_00, cp_00, kap_00, h_00, ierr) result(Dr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: T_00, d_00, cp_00, kap_00, h_00
         type(auto_diff_real_18var_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_18var_order1) :: eturb_00
         real(dp) :: gammar, alpha
         include 'formats'
         ierr = 0
         alpha = s% Eturb_alfa
         if (alpha == 0d0) then
            Dr = 0d0
         else
            gammar = s% Eturb_alfar
            eturb_00 = wrap_eturb_00(s,k)
            Dr = (4d0 * boltz_sigma * pow2(gammar) / pow2(alpha)) * pow3(T_00) * eturb_00 / &
                  (pow2(d_00) * Cp_00 * kap_00 * pow2(h_00))
         end if
         s% DAMPR(k) = Dr%val
      end function compute_Dr

      !! Calculates turbulent dissipation rate of turbulent energy in cell k.
      !!
      !! @param h_00 Pressure scale height in cell k
      !! @param w_00 Turbulent velocity in cell k
      function compute_D_adam(s, k, h_00, w_00, ierr) result(D)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: h_00, w_00
         type(auto_diff_real_18var_order1) :: D
         integer, intent(out) :: ierr

         real(dp) :: alpha
         ierr = 0
         alpha = s% Eturb_alfa
         D = 0d0
         if (alpha /= 0d0) D = pow(w_00,1.5d0) / (alpha * h_00)
         s% DAMP(k) = D%val

      end function compute_D_adam

      function compute_D(s, k, h_00, eturb_00, ierr) result(D) ! following RSP's definition
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: h_00, eturb_00
         type(auto_diff_real_18var_order1) :: D
         integer, intent(out) :: ierr
         include 'formats'

         real(dp) :: alpha
         ierr = 0
         alpha = s% Eturb_alfa
         if (alpha == 0d0) then
            D = 0d0
         else
            D = (pow(eturb_00,1.5d0) - min_Eturb_pow_1_pt_5)/(alpha * h_00)
         end if
         s% DAMP(k) = D%val

      end function compute_D

      !! Calculates the source term for turbulent energy.
      !!
      !! From MESA V: S = alpha * w_turb * T * rho * Q * Y_sag / h
      !!
      !! Q here is the thermal expansion coefficient, given by
      !!
      !! Q = ChiT / (rho * T * ChiRho)
      !!
      !! So we can write
      !!
      !! S = alpha * w_turb * (ChiT / ChiRho) * Y_sag / h
      !!
      !! Y_sag = -(h / c_p) ds/dr
      !!
      !! So we can further simplify
      !!
      !! S = alpha * w_turb * (ChiT / ChiRho) * (-1/c_p) ds/dr
      !!
      !! Because s (entropy) is defined on cells ds/dr is naturally defined on faces.
      !! S (the source) is defined on cells though, so we have to average ds/dr onto cells.
      !! The way Eturb handles this is by computing Y*cp on faces (i.e. -h*ds/dr, using h averaged onto faces)
      !! then averaging that quantity divided by h_face onto cells. 
      !! This is equivalent to just averaging ds/dr between the two adjacent faces, so we do that.
      !!
      !! @param alpha Eturb alpha parameter
      !! @param dm_m1 Mass of cell k-1
      !! @param dm_00 Mass of cell k
      !! @param dm_p1 Mass of cell k+1
      !! @param w_00 Turbulent velocity in cell k
      !! @param ChiT_00 dlnP/dlnT|rho in cell k
      !! @param ChiRho_00 dlnP/dlnRho|T in cell k
      !! @param Cp_00 Specific heat in cell k
      !! @param s_m1 Entropy in cell k-1
      !! @param s_00 Entropy in cell k
      !! @param s_p1 Entropy in cell k+1
      !! @param d_m1 Density of cell k-1
      !! @param d_00 Density of cell k
      !! @param d_p1 Density of cell k+1
      !! @param r_00 Radius of face k
      !! @param r_p1 Radius of face k+1
      function compute_Source_adam(s, k, &
            dm_m1, dm_00, dm_p1, &
            w_00, ChiT_00, ChiRho_00, Cp_00, s_m1, s_00, s_p1, &
            d_m1, d_00, d_p1, r_00, r_p1, ierr) result(Source)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_m1, dm_00, dm_p1
         type(auto_diff_real_18var_order1), intent(in) :: w_00, ChiT_00, ChiRho_00, Cp_00
         type(auto_diff_real_18var_order1), intent(in) :: s_m1, s_00, s_p1, d_m1, d_00, d_p1, r_00, r_p1

         ! Intermediates
         type(auto_diff_real_18var_order1) :: ds_dr_00, ds_dr_p1, ds_dr, dr_bar_00, dr_bar_p1

         ! Outputs
         type(auto_diff_real_18var_order1) :: Source
         integer, intent(out) :: ierr
         
         real(dp) :: alpha
         ierr = 0
         alpha = s% Eturb_alfa
         
         dr_bar_00 = 0.5d0 * (dm_m1 / d_m1 + dm_00 / d_00) / (4d0 * pi * pow2(r_00))
         dr_bar_p1 = 0.5d0 * (dm_00 / d_00 + dm_p1 / d_p1) / (4d0 * pi * pow2(r_p1))

         ds_dr_00 = (s_00 - s_m1) / dr_bar_00
         ds_dr_p1 = (s_p1 - s_00) / dr_bar_p1
         ds_dr = 0.5d0 * (ds_dr_00 + ds_dr_p1)
         Source = alpha * w_00 * (ChiT_00 / ChiRho_00) * (-1d0 / Cp_00) * ds_dr

         s% SOURCE(k) = Source%val

      end function compute_Source_adam

      function compute_Source(s, k, & ! following RSP's definition
            eturb_00, d_00, d_p1, T_00, P_00, h_00, r_00, r_p1, &
            ChiT_00, ChiRho_00, Cp_00, ierr) result(Source)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: &
            eturb_00, d_00, d_p1, T_00, P_00, h_00, r_00, r_p1, &
            ChiT_00, ChiRho_00, Cp_00
         ! Outputs
         type(auto_diff_real_18var_order1) :: Source
         integer, intent(out) :: ierr
         
         type(auto_diff_real_18var_order1) :: &
            T_m1, T_p1, d_m1, P_m1, P_p1, chiT_m1, chiT_p1, chiRho_m1, chiRho_p1, &
            Cp_m1, Cp_p1, QQ_m1, QQ_00, QQ_p1, lnT_m1, lnT_00, lnT_p1, &
            QQ_div_Cp_face_00, QQ_div_Cp_face_p1, QQ_div_Cp, Y1_00, Y1_p1, Y_cell, PII
         real(dp) :: alpha
         include 'formats'
         ierr = 0
         alpha = s% Eturb_alfa
         if (alpha == 0d0) then
            Source = 0d0
         else
         
            T_m1 = wrap_T_m1(s, k)
            T_p1 = wrap_T_p1(s, k)                  
            d_m1 = wrap_d_m1(s, k)         
            P_m1 = wrap_P_m1(s, k)
            P_p1 = wrap_P_p1(s, k)
            chiT_m1 = wrap_chiT_m1(s, k)
            chiT_p1 = wrap_chiT_p1(s, k)         
            chiRho_m1 = wrap_chiRho_m1(s, k)
            chiRho_p1 = wrap_chiRho_p1(s, k)
            Cp_m1 = wrap_Cp_m1(s, k)
            Cp_p1 = wrap_Cp_p1(s, k)        
            QQ_m1 = chiT_m1/(d_m1*T_m1*chiRho_m1) ! thermal expansion coefficient
            QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            QQ_p1 = chiT_p1/(d_p1*T_p1*chiRho_p1)         
            lnT_m1 = 0d0
            if (k > 1) then
               lnT_m1%val = s% lnT(k-1)
               lnT_m1% d1Array(i_lnT_m1) = 1d0
            end if
            lnT_00 = 0d0
            lnT_00%val = s% lnT(k)
            lnT_00% d1Array(i_lnT_00) = 1d0
            lnT_p1 = 0d0
            if (k < s% nz) then
               lnT_p1%val = s% lnT(k+1)
               lnT_p1% d1Array(i_lnT_p1) = 1d0
            end if
            
            if (k > 1) then
               QQ_div_Cp_face_00 = 0.5d0*(QQ_m1/Cp_m1 + QQ_00/Cp_00)
               Y1_00 = QQ_div_Cp_face_00*(P_m1 - P_00) - (lnT_m1 - lnT_00)
            else
               QQ_div_Cp_face_00 = 0d0
               Y1_00 = 0d0
            end if
            
            if (k < s% nz) then
               QQ_div_Cp_face_p1 = 0.5d0*(QQ_00/Cp_00 + QQ_p1/Cp_p1)
               Y1_p1 = QQ_div_Cp_face_p1*(P_00 - P_p1) - (lnT_00 - lnT_p1)
            else
               QQ_div_Cp_face_p1 = 0d0
               Y1_p1 = 0d0
            end if
               
            QQ_div_Cp = 0.5d0*(QQ_div_Cp_face_00 + QQ_div_Cp_face_p1)

            Y_cell = 0.5d0*(Y1_00*pow2(r_00) + Y1_p1*pow2(r_p1))

            Source = alpha * Y_cell * d_00 * T_00 * P_00 * Cp_00 * QQ_div_Cp * sqrt(eturb_00) / s% dm(k)

         if (is_bad(Source%val)) then
!$omp critical (hydro_equ_turbulent_crit3)
            write(*,2) 'Source', k, Source%val
            write(*,2) 'Y_cell', k, Y_cell%val
            write(*,2) 'QQ_div_Cp', k, QQ_div_Cp%val
            write(*,2) 'QQ_div_Cp_face_00', k, QQ_div_Cp_face_00%val
            write(*,2) 'QQ_div_Cp_face_p1', k, QQ_div_Cp_face_p1%val
            write(*,2) 'd_00', k, d_00%val
            write(*,2) 'T_00', k, T_00%val
            write(*,2) 'P_00', k, P_00%val
            write(*,2) 'Cp_00', k, Cp_00%val
            write(*,2) 'eturb_00', k, eturb_00%val
            write(*,2) 'Y1_00', k, Y1_00%val
            write(*,2) 'Y1_p1', k, Y1_p1%val
            write(*,2) 'r_00', k, r_00%val
            write(*,2) 'r_p1', k, r_p1%val
            stop 'compute_Source'
!$omp end critical (hydro_equ_turbulent_crit3)
         end if
            
         end if

         s% SOURCE(k) = Source%val

      end function compute_Source

      !! Compute the turbulent viscosity for face k.
      !!
      !! Uq = (1/rho r^3) d/dr [(4/3) alpha_m alpha rho Hp wturb r^3 d/dr (v/r)]
      !!
      !! @param r_m1 Radius of face k-1
      !! @param r_00 Radius of face k
      !! @param r_p1 Radius of face k+1
      !! @param v_m1 Velocity of face k-1
      !! @param v_00 Velocity of face k
      !! @param v_p1 Velocity of face k+1
      !! @param m_m1 Mass coordinate of face k-1
      !! @param m_00 Mass coordinate of face k
      !! @param m_p1 Mass coordinate of face k+1
      !! @param cgrav_m1 Gravitational constant of face k-1
      !! @param cgrav_00 Gravitational constant of face k
      !! @param cgrav_p1 Gravitational constant of face k+1
      !! @param dm_m1 Mass of cell k-1
      !! @param dm_00 Mass of cell k
      !! @param d_m1 Density of cell k-1
      !! @param d_00 Density of cell k
      !! @param P_m1 Pressure of cell k-1
      !! @param P_00 Pressure of cell k
      !! @param w_m1 Turbulent velocty of cell k-1
      !! @param w_00 Turbulent velocity of cell k
      function compute_Uq_adam(s, k, r_m1, r_00, r_p1, v_m1, v_00, v_p1, m_m1, m_00, m_p1, &
                          cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00, d_m1, d_00, &
                          P_m1, P_00, w_m1, w_00, ierr) result(Uq)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: m_m1, m_00, m_p1, cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00
         type(auto_diff_real_18var_order1), intent(in) :: r_m1, r_00, r_p1, v_m1, v_00, v_p1
         type(auto_diff_real_18var_order1), intent(in) :: d_m1, d_00, P_m1, P_00, w_m1, w_00

         ! Intermediates
         type(auto_diff_real_18var_order1) :: dr_00, dr_m1, d_v_div_r_dr_00, d_v_div_r_dr_m1
         type(auto_diff_real_18var_order1) :: rmid_00, rmid_m1, mmid_00, mmid_m1, g_00, g_m1
         type(auto_diff_real_18var_order1) :: h_00, h_m1, flux_00, flux_m1, rho_inv_avg, area, dr_bar

         ! Outputs
         type(auto_diff_real_18var_order1) :: Uq
         
         real(dp) :: alpha, alpha_m, cgrav_mid_00, cgrav_mid_m1
         integer, intent(out) :: ierr
         ierr = 0
         
         alpha = s% Eturb_alfa
         alpha_m = s% Eturb_alfam

         dr_00 = dm_00 / (4d0 * pi * d_00 * (pow2(r_00) + r_00*r_p1 + pow2(r_p1)))
         dr_m1 = dm_m1 / (4d0 * pi * d_m1 * (pow2(r_m1) + r_m1*r_00 + pow2(r_00)))

         d_v_div_r_dr_00 = (v_p1 / r_p1 - v_00 / r_00) / dr_00
         d_v_div_r_dr_m1 = (v_00 / r_00 - v_00 / r_00) / dr_m1
         
         cgrav_mid_00 = 0.5d0 * (cgrav_p1 + cgrav_00)
         rmid_00 = 0.5d0 * (r_00 + r_p1)
         mmid_00 = 0.5d0 * (m_00 + m_p1)
         g_00 = cgrav_mid_00 * mmid_00 / pow2(rmid_00) ! gravity in cell k

         cgrav_mid_m1 = 0.5d0 * (cgrav_00 + cgrav_m1)
         rmid_m1 = 0.5d0 * (r_m1 + r_00)
         mmid_m1 = 0.5d0 * (m_m1 + m_00)
         g_m1 = cgrav_mid_m1 * mmid_m1 / pow2(rmid_m1) ! gravity in cell k-1

         h_00 = P_00 / (d_00 * g_00) ! Scale height in cell k
         h_m1 = P_m1 / (d_m1 * g_m1) ! Scale height in cell k-1

         flux_00 = (4d0 / 3d0) * alpha * alpha_m * d_00 * h_00 * w_00 * pow3(rmid_00) * d_v_div_r_dr_00 ! Flux through cell k
         flux_m1 = (4d0 / 3d0) * alpha * alpha_m * d_m1 * h_m1 * w_m1 * pow3(rmid_m1) * d_v_div_r_dr_m1 ! Flux through cell k-1

         rho_inv_avg = 0.5d0 * (1d0 / d_m1 + 1d0 / d_00)

         area = 4d0 * pi * pow2(r_00)
         dr_bar = 0.5d0 * (dm_m1 / d_m1 + dm_00 / d_00) / area
         Uq = rho_inv_avg * (flux_m1 - flux_00) / (dr_bar * pow3(r_00))

         s% Uq(k) = Uq%val

      end function compute_Uq_adam

      function compute_Uq(s, k, & ! following RSP's definition
            r_m1, r_00, r_p1, v_m1, v_00, v_p1, m_m1, m_00, m_p1, &
            cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00, d_m1, d_00, &
            P_m1, P_00, eturb_m1, eturb_00, ierr) result(Uq)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: m_m1, m_00, m_p1, cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00
         type(auto_diff_real_18var_order1), intent(in) :: r_m1, r_00, r_p1, v_m1, v_00, v_p1
         type(auto_diff_real_18var_order1), intent(in) :: d_m1, d_00, P_m1, P_00, eturb_m1, eturb_00
         ! Outputs
         type(auto_diff_real_18var_order1) :: Uq
         integer, intent(out) :: ierr
         real(dp) :: alpha, alpha_m
         type(auto_diff_real_18var_order1) :: w_m1, w_00, r6_00, &
            cgrav_mid_00, rmid_00, mmid_00, g_00, h_00, &
            w_rho2_00, r6_cell_00, d_v_div_r_00, Chi_00, &
            cgrav_mid_m1, rmid_m1, mmid_m1, g_m1, h_m1, &
            w_rho2_m1, r6_cell_m1, d_v_div_r_m1, Chi_m1
         include 'formats'
         ierr = 0
         
         alpha = s% Eturb_alfa
         alpha_m = s% Eturb_alfam
         if (alpha == 0d0 .or. alpha_m == 0d0 .or. k == 1) then
         
            Uq = 0d0
            
         else
         
            w_m1 = sqrt(eturb_m1)
            w_00 = sqrt(eturb_00)
         
            cgrav_mid_00 = 0.5d0 * (cgrav_p1 + cgrav_00)
            rmid_00 = 0.5d0 * (r_00 + r_p1)
            mmid_00 = 0.5d0 * (m_00 + m_p1)
            g_00 = cgrav_mid_00 * mmid_00 / pow2(rmid_00) ! gravity in cell k

            cgrav_mid_m1 = 0.5d0 * (cgrav_00 + cgrav_m1)
            rmid_m1 = 0.5d0 * (r_m1 + r_00)
            mmid_m1 = 0.5d0 * (m_m1 + m_00)
            g_m1 = cgrav_mid_m1 * mmid_m1 / pow2(rmid_m1) ! gravity in cell k-1

            h_00 = P_00 / (d_00 * g_00) ! Scale height in cell k
            h_m1 = P_m1 / (d_m1 * g_m1) ! Scale height in cell k-1
         
            r6_00 = pow6(r_00)
            w_rho2_00 = w_00*pow2(d_00)
            r6_cell_00 = 0.5d0*(r6_00 + pow6(r_p1))
            d_v_div_r_00 = v_00/r_00 - v_p1/r_p1
            Chi_00 = (16d0/3d0)*pi*alpha*alpha_m*w_rho2_00*r6_cell_00*h_00*d_v_div_r_00/dm_00         

            w_rho2_m1 = w_m1*pow2(d_m1)
            r6_cell_m1 = 0.5d0*(pow6(r_m1) + r6_00)
            d_v_div_r_m1 = v_m1/r_m1 - v_00/r_00
            Chi_m1 = (16d0/3d0)*pi*alpha*alpha_m*w_rho2_m1*r6_cell_m1*h_m1*d_v_div_r_m1/dm_m1         
         
            Uq = pi4*(Chi_m1 - Chi_00)/(s% dm_bar(k)*r_00)
         
         if (.false.) then
            
            write(*,2) 'eturb_00%val', k, eturb_00%val
            if (is_bad(eturb_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 eturb_00', k, eturb_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            write(*,2) 'w_00%val', k, w_00%val
            if (is_bad(w_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 w_00', k, w_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(d_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 d_00', k, d_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(w_rho2_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 w_rho2_00', k, w_rho2_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(r6_cell_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 r6_cell_00', k, r6_cell_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(h_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 h_00', k, h_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(d_v_div_r_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 d_v_div_r_00', k, d_v_div_r_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

            if (is_bad(Chi_m1%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 Chi_m1', k, Chi_m1%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if
            if (is_bad(Chi_00%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 Chi_00', k, Chi_00%d1Array(i_lnd_m1)
               stop 'compute_Uq'
            end if

         end if
            
         end if

         s% Uq(k) = Uq%val

      end function compute_Uq

      !! Computes the radiative luminosity at face k
      !!
      !! @param k Face index
      !! @param dm_bar Face mass
      !! @param r Face radius
      !! @param T_m1 Temperature of cell k-1
      !! @param T_00 Temperature of cell k
      !! @param kap_m1 Opacity of cell k-1
      !! @param kap_00 Opacity of cell k
      !! @param Lr Radiative Luminosity
      function compute_Lr(s, k, dm_bar, area, T_m1, T_00, kap_m1, kap_00, ierr) result(Lr)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_bar
         type(auto_diff_real_18var_order1), intent(in) :: area, T_m1, T_00, kap_m1, kap_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: Lr
         integer, intent(out) :: ierr

         ! Intermediates
         type(auto_diff_real_18var_order1) :: Erad, T400, T4m1, d_P_rad, kap_face
         real(dp) :: alfa
         
         include 'formats'
         
         ierr = 0

         if (k == 1) then ! Lr(1) proportional to Erad in cell(1)
         
            Erad = crad * pow4(T_00)
            Lr = s% Eturb_Lsurf_factor * area * clight * Erad
            
         else 
         
            T400 = pow4(T_00)
            T4m1 = pow4(T_m1)
            d_P_rad = crad/3d0*(T4m1 - T400)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            kap_face = alfa*kap_00 + (1d0 - alfa)*kap_m1
            Lr = -d_P_rad/dm_bar*clight*pow2(area)/kap_face
            
         end if
         
         s% Lr(k) = Lr%val

      end function compute_Lr

      !! Computes the convective luminosity at face k
      !!
      !! @param dm_m1 Mass of cell k-1
      !! @param dm_00 Mass of cell k
      !! @param r Face radius
      !! @param T_m1 Temperature of cell k-1
      !! @param T_00 Temperature of cell k
      !! @param d_m1 Density of cell k-1
      !! @param d_00 Density of cell k
      !! @param P_m1 Pressure of cell k-1
      !! @param P_00 Pressure of cell k
      !! @param s_m1 Entropy of cell k-1
      !! @param s_00 Entropy of cell k
      !! @param Lc Convective Luminosity
      function compute_Lc_adam(s, k, dm_m1, dm_00, g, area, &
                          T_m1, T_00, d_m1, d_00, P_m1, P_00, s_m1, s_00, w_m1, w_00, ierr) result(Lc)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_m1, dm_00
         type(auto_diff_real_18var_order1), intent(in) :: &
            g, area, T_m1, T_00, d_m1, d_00, P_m1, P_00, s_m1, s_00, w_m1, w_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: Lc
         integer, intent(out) :: ierr

         ! Intermediates
         real(dp) :: alpha
         type(auto_diff_real_18var_order1) :: dr_bar, w_turb_face, T_rho_face, h_face, ds_dr, Y_sag_over_cp_face

         ierr = 0
         alpha = s% Eturb_alfa

         h_face = 0.5d0 * (P_m1 / d_m1 + P_00 / d_00) / g
         dr_bar = 0.5d0 * (dm_m1 / d_m1 + dm_00 / d_00) / area

         w_turb_face = 0.5d0 * (w_m1 + w_00)
         T_rho_face = 0.5d0 * (T_00 * d_00 + T_m1 * d_m1)

         ds_dr = -(s_00 - s_m1) / dr_bar ! Minus because r(k) < r(k-1)
         Y_sag_over_cp_face = -h_face * ds_dr

         Lc = alpha * T_rho_face * w_turb_face * Y_sag_over_cp_face

         s% Lc(k) = Lc%val

      end function compute_Lc_adam

      function compute_Lc(s, k, ierr) result(Lc) ! copy the RSP version
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         ! Outputs
         type(auto_diff_real_18var_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_18var_order1) :: Lc_w_face_factor
         Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc

      function compute_Lc_terms(s, k, Lc_w_face_factor, ierr) result(Lc) ! copy the RSP version
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         ! Outputs
         type(auto_diff_real_18var_order1) :: Lc, Lc_w_face_factor
         integer, intent(out) :: ierr
         ! Intermediates
         real(dp) :: alpha, dm_bar, cgrav, m
         type(auto_diff_real_18var_order1) :: &
            T_m1, lnT_m1, d_m1, P_m1, w_m1, chiT_m1, chiRho_m1, Cp_m1, QQ_m1, &
            T_00, lnT_00, d_00, P_00, w_00, chiT_00, chiRho_00, Cp_00, QQ_00, r_00, &
            P_div_rho_face, Hp_face, QQ_div_Cp_face, Y1, avg_Vol, area, Y2, Y_face, &
            Cp_face, PII_face, w_face, T_rho_face
         include 'formats'

         ierr = 0

         ! Get reals
         alpha = s% Eturb_alfa
         if (alpha <= 0d0 .or. k == 1) then
            area = 0d0
            T_rho_face = 0d0
            PII_face = 0d0
            w_face = 0d0
            Lc = 0d0
            return
         end if
         
         dm_bar = s% dm_bar(k)
         cgrav = s% cgrav(k)
         m = s% m(k)

         ! Wrap auto_diff variables
         r_00 = wrap_r_00(s, k)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)         
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         P_m1 = wrap_P_m1(s, k)
         P_00 = wrap_P_00(s, k)
         w_m1 = sqrt(wrap_eturb_m1(s, k))
         w_00 = sqrt(wrap_eturb_00(s, k))
         chiT_m1 = wrap_chiT_m1(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_m1 = wrap_chiRho_m1(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         lnT_m1 = 0d0
         lnT_m1%val = s% lnT(k-1)
         lnT_m1% d1Array(i_lnT_m1) = 1d0
         lnT_00 = 0d0
         lnT_00%val = s% lnT(k)
         lnT_00% d1Array(i_lnT_00) = 1d0

         QQ_m1 = chiT_m1/(d_m1*T_m1*chiRho_m1) ! thermal expansion coefficient
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
         
         P_div_rho_face = 0.5d0*(P_00/d_00 + P_m1/d_m1)
         Hp_face = P_div_rho_face*pow2(r_00)/(cgrav*m)
         QQ_div_Cp_face = 0.5d0*(QQ_m1/Cp_m1 + QQ_00/Cp_00)
         Y1 = QQ_div_Cp_face*(P_m1 - P_00) - (lnT_m1 - lnT_00)
         avg_Vol = 0.5d0*(1d0/d_00 + 1d0/d_m1)
         area = 4d0*pi*pow2(r_00)
         Y2 = area*Hp_face/(avg_Vol*dm_bar)
         Y_face = Y1*Y2
         
         Cp_face = 0.5d0*(Cp_m1 + Cp_00)
         PII_face = alpha*Cp_face*Y_face
         w_face = 0.5d0*(w_m1 + w_00)
         T_rho_face = 0.5d0*(T_m1*d_m1 + T_00*d_00)
         Lc_w_face_factor = area*T_rho_face*PII_face
         Lc = Lc_w_face_factor*w_face

      end function compute_Lc_terms

      !! Computes the turbulent luminosity (i.e. transport of turbulent kinetic energy) at face k
      !!
      !! @param dm_m1 Mass of cell k-1
      !! @param dm_00 Mass of cell k
      !! @param r Face radius
      !! @param P_m1 Pressure of cell k-1
      !! @param P_00 Pressure of cell k
      !! @param d_m1 Density of cell k-1
      !! @param d_00 Density of cell k
      !! @param w_m1 Turbulent velocity of cell k-1
      !! @param w_00 Turbulent velocity of cell k
      function compute_Lt_adam(s, k, dm_m1, dm_00, g, area, &
                          P_m1, P_00, d_m1, d_00, w_m1, w_00, ierr) result(Lt_18)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_m1, dm_00
         type(auto_diff_real_18var_order1), intent(in) :: g, area, P_m1, P_00, d_m1, d_00, w_m1, w_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: Lt_18
         integer, intent(out) :: ierr
         
         type(auto_diff_real_18var_order1) :: rho_h_face, dr_bar, wturb_d_eturb_dr, Ft
         real(dp) :: alpha, alpha_t
         ierr = 0
         alpha = s% Eturb_alfa
         alpha_t = s% Eturb_alfat
         
         rho_h_face = 0.5d0 * (P_m1 + P_00) / g

         dr_bar = 0.5d0 * (dm_m1 / d_m1 + dm_00 / d_00) / area
         wturb_d_eturb_dr = -(2d0/3d0) * (pow(w_00,1.5d0) - pow(w_m1,1.5d0)) / dr_bar ! Minus because r(k) < r(k-1)
         Ft = -alpha * alpha_t * rho_h_face * wturb_d_eturb_dr

         Lt_18 = area * Ft

         s% Lt(k) = Lt_18%val

      end function compute_Lt_adam

      function compute_Lt(s, k, & ! copy the RSP version
            dm_m1, dm_00, g, area, &
            P_m1, P_00, d_m1, d_00, eturb_m1, eturb_00, ierr) result(Lt_18)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dm_m1, dm_00
         type(auto_diff_real_18var_order1), intent(in) :: &
            g, area, P_m1, P_00, d_m1, d_00, eturb_m1, eturb_00

         ! Outputs
         type(auto_diff_real_18var_order1) :: Lt_18
         integer, intent(out) :: ierr
         
         type(auto_diff_real_18var_order1) :: &
            r_00, P_div_rho_face, Hp_face, rho2_face
         real(dp) :: alpha, alpha_t
         include 'formats'
         ierr = 0

         alpha = s% Eturb_alfa
         alpha_t = s% Eturb_alfat
         if (alpha <= 0d0 .or. alpha_t <= 0d0 .or. k == 1) then
            Lt_18 = 0d0
         else
            r_00 = wrap_r_00(s,k)
            P_div_rho_face = 0.5d0*(P_00/d_00 + P_m1/d_m1)
            Hp_face = P_div_rho_face*pow2(r_00)/(s% cgrav(k)*s% m(k))
            rho2_face = 0.5d0*(pow2(d_00) + pow2(d_m1))
            Lt_18 = -2d0/3d0*alpha*alpha_t*pow2(area)*Hp_face*rho2_face*&
               (pow(eturb_m1,1.5d0) - pow(eturb_00,1.5d0))/s% dm_bar(k)            
         end if

         s% Lt(k) = Lt_18%val

      end function compute_Lt

      !! Compute the luminosity L = L_rad + L_conv + L_turb.
      !!
      !! @param s Star pointer
      !! @param k Face index
      subroutine compute_L(s, k, L, Lr, Lc, Lt, ierr)
         ! Inputs
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k

         ! Outputs
         type(auto_diff_real_18var_order1), intent(out) :: L, Lr, Lc, Lt
         integer, intent(out) :: ierr         

         ! Intermediates
         real(dp) :: cgrav, dm_bar, m, dm_m1, dm_00, unused
         type(auto_diff_real_18var_order1) :: r_00, T_m1, T_00, kap_m1, kap_00, d_m1, d_00
         type(auto_diff_real_18var_order1) :: P_m1, P_00, s_m1, s_00, eturb_m1, eturb_00, g, area
         type(auto_diff_real_18var_order1) :: g_cell_00, h_00, r_p1, v_00, v_p1
         
         include 'formats'
         ierr = 0
         
         ! Get reals
         dm_bar = s% dm_bar(k)

         ! Wrap auto_diff variables
         r_00 = wrap_r_00(s, k)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)
         kap_m1 = wrap_kap_m1(s, k)
         kap_00 = wrap_kap_00(s, k)

         area = 4d0 * pi * pow2(r_00)

         ! Radiative luminosity               
         Lr = compute_Lr(s, k, dm_bar, area, T_m1, T_00, kap_m1, kap_00, ierr)
         if (ierr /= 0) return

         ! Turbulent and convective luminosities.
         if (k == 1) then            
            Lc = 0d0
            Lt = 0d0
         else
            m = s% m(k)
            cgrav = s% cgrav(k)

            g = m * cgrav / pow2(r_00)

            r_p1 = wrap_r_p1(s, k)

            dm_m1 = s%dm(k-1)
            dm_00 = s%dm(k)

            d_m1 = wrap_d_m1(s, k)
            d_00 = wrap_d_00(s, k)
            P_m1 = wrap_P_m1(s, k)
            P_00 = wrap_P_00(s, k)
            s_m1 = wrap_s_m1(s, k)
            s_00 = wrap_s_00(s, k)

            eturb_m1 = wrap_eturb_m1(s, k)
            eturb_00 = wrap_eturb_00(s, k)

            v_00 = wrap_v_00(s, k)
            v_p1 = wrap_v_p1(s, k) ! Set by wrap routine to zero when k == nz.

            !Lc = compute_Lc_adam(s, k, dm_m1, dm_00, g, area, T_m1, T_00, d_m1, d_00, P_m1, P_00, s_m1, s_00, w_m1, w_00, ierr)
            Lc = compute_Lc(s, k, ierr)
            if (ierr /= 0) return
            
            Lt = compute_Lt(s, k, dm_m1, dm_00, g, area, P_m1, P_00, d_m1, d_00, eturb_m1, eturb_00, ierr)
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

      end subroutine compute_L


      subroutine calc_Eq_18(s, k, Eq_18, ierr) ! used by energy equation
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(out) :: Eq_18
         integer, intent(out) :: ierr
         real(dp) :: m_00, dm_00, cgrav
         type(auto_diff_real_18var_order1) :: &
            v_p1, v_00, r_p1, r_00, d_00, P_00, eturb_00, g_cell_00, h_00
         include 'formats'

         ierr = 0
         m_00 = s% m(k)
         dm_00 = s% dm(k)
         cgrav = s% cgrav(k)
         r_00 = wrap_r_00(s, k)
         r_p1 = wrap_r_p1(s, k) ! Set by wrap routine to r_center when k == nz.
         P_00 = wrap_P_00(s, k)
         d_00 = wrap_d_00(s, k)
         eturb_00 = wrap_eturb_00(s, k)
         v_00 = wrap_v_00(s, k)
         v_p1 = wrap_v_p1(s, k) ! Set by wrap routine to zero when k == nz.
         
         ! Compute scale height on cell k (used by epsilon_q)
         g_cell_00 = (m_00 - 0.5d0 * dm_00) * cgrav / (0.5d0 * (r_00 + r_p1))
         h_00 = compute_h_00(s, k, P_00, g_cell_00, d_00, ierr)
         if (ierr /= 0) return 
         
         Eq_18 = compute_epsilon_q(s, k, v_p1, v_00, &
                                 r_p1, r_00, d_00, dm_00, eturb_00, h_00, ierr)

      end subroutine calc_Eq_18
      

      subroutine calc_Uq_18(s, k, Uq_18, ierr) ! used by momentum equation
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(out) :: Uq_18
         integer, intent(out) :: ierr
         real(dp) :: m_m1, m_00, m_p1, cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00
         type(auto_diff_real_18var_order1) :: r_m1, r_00, r_p1, v_m1, v_00, v_p1
         type(auto_diff_real_18var_order1) :: d_m1, d_00, P_m1, P_00, eturb_m1, eturb_00
         include 'formats'

         ierr = 0
         
         if (k == 1) then
            Uq_18 = 0d0
            return
         end if
         
         cgrav_m1 = s% cgrav(k-1)
         m_m1 = s% m(k-1)
         dm_m1 = s% dm(k-1)
         cgrav_00 = s% cgrav(k)
         m_00 = s% m(k)
         dm_00 = s% dm(k)
         if (k == s% nz) then
            cgrav_p1 = cgrav_00
            m_p1 = s% m_center
         else
            cgrav_p1 = s% cgrav(k+1)
            m_p1 = s% m(k+1)
         end if

         r_m1 = wrap_r_m1(s, k)
         r_00 = wrap_r_00(s, k)
         r_p1 = wrap_r_p1(s, k) ! Set by wrap routine to r_center when k == nz.

         v_m1 = wrap_v_m1(s, k)
         v_00 = wrap_v_00(s, k)
         v_p1 = wrap_v_p1(s, k) ! Set by wrap routine to zero when k == nz.

         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)

         eturb_m1 = wrap_eturb_m1(s, k)
         eturb_00 = wrap_eturb_00(s, k)

         P_00 = wrap_P_00(s, k)
         P_m1 = wrap_P_m1(s, k)

         Uq_18 = compute_Uq(s, k, r_m1, r_00, r_p1, v_m1, v_00, v_p1, m_m1, m_00, m_p1, &
                          cgrav_m1, cgrav_00, cgrav_p1, dm_m1, dm_00, d_m1, d_00, &
                          P_m1, P_00, eturb_m1, eturb_00, ierr)

      end subroutine calc_Uq_18


      !! Calculates the pressure scale height in cell k
      !!
      !! @param P_00 Pressure in cell k
      !! @param g_00 Gravity in cell k
      !! @param d_00 Density in cell k
      function compute_h_00(s, k, P_00, g_cell_00, d_00, ierr) result(h_00)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: P_00, g_cell_00, d_00
         type(auto_diff_real_18var_order1) :: h_00
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         h_00 = P_00 / (d_00 * g_cell_00)
      end function compute_h_00


      subroutine set_eturb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, op_err
         logical :: time_center
         include 'formats'
         ierr = 0
         time_center = (s% use_Fraley_time_centering .and. s% include_L_in_Fraley_time_centering)
         do k=1,s%nz
            call set1_eturb_start_vars(k, op_err) 
            if (op_err /= 0) ierr = op_err  
         end do
         contains
         subroutine set1_eturb_start_vars(k, ierr)   
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: L, Lr, Lc, Lt
            include 'formats'
            ierr = 0               
            if (time_center) then
               call compute_L(s, k, L, Lr, Lc, Lt, ierr)
               if (ierr /= 0) return
               s% Lt_start(k) = Lt%val  
            else
               s% Lt_start(k) = 0d0  
            end if
            s% Eturb_start(k) = s% Eturb(k)
         end subroutine set1_eturb_start_vars
         
      end subroutine set_eturb_start_vars
      
      
      subroutine reset_eturb_using_L(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr   
         integer :: k, i_eturb, nz
         real(dp) :: alpha, dm_bar, Lc_val, w_00
         type(auto_diff_real_18var_order1) :: &
            r_00, area, T_m1, T_00, kap_m1, kap_00, &
            Lc_w_face_factor, L, Lr, Lc, Lt
         real(dp), allocatable :: w_face(:)
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         alpha = s% Eturb_alfa
         if (alpha == 0d0) return
         nz = s% nz
         allocate(w_face(nz))
         i_eturb = s% i_eturb
         w_face(1) = 0d0
         do k=2, nz
            dm_bar = s% dm_bar(k)
            r_00 = wrap_r_00(s,k)
            area = 4d0*pi*pow2(r_00)
            T_m1 = wrap_T_m1(s,k)
            T_00 = wrap_T_00(s,k)
            kap_m1 = wrap_kap_m1(s,k)
            kap_00 = wrap_kap_00(s,k)
            Lr = compute_Lr(s, k, dm_bar, area, T_m1, T_00, kap_m1, kap_00, ierr)
            if (ierr /= 0) stop 'failed in compute_Lr'
            Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
            if (ierr /= 0) stop 'failed in compute_Lc_terms'
            Lc_val = s% L(k) - Lr%val ! assume Lt = 0
            w_face(k) = Lc_val/Lc_w_face_factor%val
         end do
         do k=1, nz
            if (k < nz) then
               w_00 = 0.5d0*(w_face(k) + w_face(k+1))
            else
               w_00 = w_face(k)
            end if
            if (w_00 < 0d0) w_00 = 0d0
            s% xh(i_eturb,k) = max(min_eturb, pow2(w_00))
            s% Eturb(k) = s% xh(i_eturb,k)
            call compute_L(s, k, L, Lr, Lc, Lt, ierr)
            if (ierr /= 0) stop 'failed in compute_L reset_wturb_using_L'
            if (dbg) write(*,2) 'L_eturb/L, Lc_eturb/L, Lr_eturb/L, Eturb/energy, mlt_vc/cs', k, &
               L%val/s% L(k), Lc%val/s% L(k), Lr%val/s% L(k), &
               s% Eturb(k)/s% energy(k), s% mlt_vc(k)/s% csound(k)
         end do
         if (dbg) stop 'reset_eturb_using_L'
      end subroutine reset_eturb_using_L


      end module hydro_eturb

