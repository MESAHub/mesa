! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module hydro_rsp2

      use star_private_def
      use const_def, only: dp, boltz_sigma, pi, clight, crad, ln10
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use hydro_gradient_support, only: get_RSP2_alfa_beta_face_weights, &
         get_rsp2_face_eos, get_rsp2_Lrad_coeff, get_rsp2_thermal_gradient
      use accurate_sum_auto_diff_star_order1
      use star_utils
      use tdc_hydro, only: get_TDC_Hp_face, &
         get_TDC_mixing_length_face, get_TDC_mixing_length_cell

      implicit none

      private
      public :: do1_rsp2_L_eqn
      public :: do1_turbulent_energy_eqn
      public :: do1_rsp2_flux_eqn
      public :: compute_Source, compute_D, compute_Dr
      public :: compute_Source_div_w
      public :: compute_Eq_cell
      public :: compute_Uq_face, compute_Uq_dm_cell
      public :: set_RSP2_vars
      public :: rsp2_flux_residual, set_etrb_start_vars
      public :: RSP2_adjust_vars_before_call_solver
      public :: get_RSP2_alfa_beta_face_weights
      public :: do1_rsp2_moment_eqns, rsp2_moment_rhs
      public :: rsp2_zero_moments, init_rsp2_moments
      public :: rsp2_dormant_moments, rsp2_local_w_equation
      public :: remesh_rsp2_moments, interpolate_rsp2_face
      public :: rsp2_remesh_w_face

      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, &  ! Ptrb
         x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3, &  ! PII_face and Lc
         x_ALFAC = (1.d0/2.d0)*sqrt_2_div_3, &  ! Lc
         x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3, &  ! DAMP
         x_GAMMAR = 2.d0*sqrt(3.d0)  ! DAMPR

      ! Kuhfuss local MLT calibration, as in Braun et al. (2026), section 2.
      ! Pi = <v_r*s>; Phi = <s*s>. Using the full variance does not change its decay rate.
      ! RSP2_alfa_pi and RSP2_alfa_phi multiply these coefficients.
      real(dp), parameter :: &
         x_ALFAPI = 6d0*sqrt_2_div_3, &
         x_ALFAPHI = 4d0*sqrt_2_div_3

      contains

      subroutine set_RSP2_vars(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: x
         integer :: k, op_err
         include 'formats'
         ierr = 0
         op_err = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
         do k=1,s%nz
            op_err = 0
            x = get_TDC_Hp_face(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            s% Hp_face(k) = x%val
            s% scale_height_ad(k) = x
            s% scale_height(k) = x%val
            x = get_TDC_mixing_length_face(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            s% Lambda_ad(k) = x
            s% mlt_mixing_length(k) = x%val
            x = compute_RSP2_gradT(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            s% gradT_ad(k) = x
            s% gradT(k) = x%val
            s% mlt_gradT(k) = x%val
            s% Y_face_ad(k) = wrap_Y_00(s, k)
            s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k)
            op_err = 0
            x = compute_PII_face(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            !Pvsc           skip?
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 1', s% model_number
            return
         end if
         !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
         do k=1,s% nz
            op_err = 0
            x = compute_Eq_cell(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            if (s% u_flag) then
               x = compute_Uq_dm_cell(s, k, op_err)
            else if (s% v_flag) then
               x = compute_Uq_face(s, k, op_err)
            end if
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            op_err = 0
            x = compute_C(s, k, op_err)  ! COUPL
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_L_face(s, k, op_err)  ! Lr, Lt, Lc
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_Lrad_coeff(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
            else if (k > 1 .and. x%val > 0d0) then
               s% gradr_ad(k) = (wrap_L_00(s,k) - s% Lt_ad(k))/x
               s% gradr(k) = s% gradr_ad(k)%val
            end if
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 2', s% model_number
            return
         end if
         do k = 1, s% RSP2_num_outermost_cells_forced_nonturbulent
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0;
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
         do k = s% nz + 1 - int(s% nz/s% RSP2_nz_div_IBOTOM) , s% nz
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0;
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
      end subroutine set_RSP2_vars


      subroutine do1_rsp2_L_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) ::  &
            L_expected, L_actual,resid
         type(accurate_auto_diff_real_star_order1) :: L_sum
         real(dp) :: scale, residual
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (.not. s% RSP2_flag) then
            ierr = -1
            return
         end if

         ierr = 0
         if (k /= 1) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'do1_rsp2_L_eqn requires k == 1', k
            return
         end if
         L_sum = s% Lr_ad(k)
         L_sum = L_sum + s% Lc_ad(k)
         L_sum = L_sum + s% Lt_ad(k)
         L_expected = L_sum
         L_actual = wrap_L_00(s, k)
         scale = 1d0/max(abs(L_expected%val),abs(L_actual%val),1d0)
         if (is_bad(scale)) then
            write(*,2) 'do1_rsp2_L_eqn scale', k, scale
            call mesa_error(__FILE__,__LINE__,'do1_rsp2_L_eqn')
         end if
         resid = (L_expected - L_actual)*scale

         residual = resid%val
         s% equ(s% i_equL, k) = residual
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if

         call save_eqn_residual_info(s, k, nvar, s% i_equL, resid, 'do1_rsp2_L_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_rsp2_L_eqn', s% solver_test_partials_var
         end if
      end subroutine do1_rsp2_L_eqn


      function rsp2_flux_residual(s, k) result(resid)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: resid, L_expected, L_actual
         type(accurate_auto_diff_real_star_order1) :: L_sum
         real(dp) :: scale

         if (k == 1) then
            ! The surface luminosity has its own physical boundary condition.
            resid = wrap_Y_00(s, k)
            return
         end if
         L_sum = s% Lr_ad(k)
         L_sum = L_sum + s% Lc_ad(k)
         L_sum = L_sum + s% Lt_ad(k)
         L_expected = L_sum
         L_actual = wrap_L_00(s, k)
         ! Use the TDC luminosity scale, fixed during solver iterations.
         if (s% solver_iter == 0) then
            scale = max(1d0, abs(s% L(k)), 1d-3*maxval(abs(s% L(1:s% nz))))
         else
            scale = max(1d0, abs(s% L_start(k)), 1d-3*maxval(abs(s% L_start(1:s% nz))))
         end if
         resid = (L_expected - L_actual)/scale
      end function rsp2_flux_residual


      subroutine do1_rsp2_flux_eqn(s, k, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: resid

         ierr = 0
         resid = rsp2_flux_residual(s, k)
         s% equ(s% i_rsp2_flux, k) = resid%val
         call save_eqn_residual_info(s, k, nvar, s% i_rsp2_flux, &
            resid, 'do1_rsp2_flux_eqn', ierr)
      end subroutine do1_rsp2_flux_eqn


      logical function rsp2_zero_moments(s, k) result(zero_moments)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         zero_moments = k <= 1 .or. k >= s% nz .or. s% mixing_length_alpha == 0d0 .or. &
            k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
      end function rsp2_zero_moments


      logical function rsp2_dormant_moments(s, k) result(dormant)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         dormant = .true.
         if (rsp2_zero_moments(s,k)) return
         dormant = get_etrb(s,k) == 0d0 .and. get_etrb(s,k-1) == 0d0 .and. &
            s% Pi(k) == 0d0 .and. s% Phi(k) == 0d0
      end function rsp2_dormant_moments


      subroutine interpolate_rsp2_face( &
            s, i_var, n_old, x_old, n_new, x_new, face_old, face_new, work, ierr)
         use interp_1d_lib, only: interpolate_vector_pm
         type(star_info), pointer :: s
         integer, intent(in) :: i_var, n_old, n_new
         real(dp), intent(in) :: x_old(:), x_new(:), face_old(:)
         real(dp), intent(inout) :: face_new(:)
         real(dp), pointer :: work(:)
         integer, intent(out) :: ierr
         integer :: j, k

         ierr = 0
         if (i_var == s% i_Phi) then
            if (any(is_bad(face_old(1:n_old))) .or. minval(face_old(1:n_old)) < 0d0) then
               s% retry_message = 'invalid RSP2 Phi before remesh'
               ierr = -1
               return
            end if
         end if
         call interpolate_vector_pm( &
            n_old, x_old, n_new, x_new, face_old, face_new, work, 'RSP2 face remesh', ierr)
         if (ierr /= 0) return
         j = 1
         do k=1,n_new
            do while (j < n_old-1)
               if (x_new(k) < x_old(j+1)) exit
               j = j+1
            end do
            ! Preserve old faces exactly and protect Phi from cubic evaluation roundoff.
            if (x_new(k) == x_old(j)) then
               face_new(k) = face_old(j)
            else if (x_new(k) == x_old(j+1)) then
               face_new(k) = face_old(j+1)
            else if (i_var == s% i_Phi) then
               face_new(k) = min(max(face_new(k), min(face_old(j),face_old(j+1))), &
                  max(face_old(j),face_old(j+1)))
            end if
         end do
      end subroutine interpolate_rsp2_face


      real(dp) function rsp2_remesh_w_face(s, k, nz, dq, w) result(w_face)
         type(star_info), pointer :: s
         integer, intent(in) :: k, nz
         real(dp), intent(in) :: dq(:), w(:)
         real(dp) :: alfa

         w_face = 0d0
         if (k <= 1 .or. k >= nz) return
         alfa = 0.5d0
         if (s% RSP2_use_mass_interp_face_values) alfa = dq(k-1)/(dq(k-1) + dq(k))
         w_face = sqrt(alfa*pow2(w(k)) + (1d0-alfa)*pow2(w(k-1)))
      end function rsp2_remesh_w_face


      subroutine remesh_rsp2_moments(s, nz, dq, xh, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), intent(in) :: dq(:)
         real(dp), intent(inout) :: xh(:,:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: alfa, beta, etrb_face

         ierr = 0
         ! Use the new mesh and remapped energy before rebuilding derived variables.
         do k=1,nz
            if (k <= 1 .or. k >= nz .or. s% mixing_length_alpha == 0d0 .or. &
                  k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
                  k > nz - int(nz/s% RSP2_nz_div_IBOTOM)) then
               xh(s% i_Pi,k) = 0d0
               xh(s% i_Phi,k) = 0d0
               cycle
            end if
            alfa = 0.5d0
            if (s% RSP2_use_mass_interp_face_values) alfa = dq(k-1)/(dq(k-1) + dq(k))
            beta = 1d0 - alfa
            etrb_face = alfa*pow2(xh(s% i_w,k)) + beta*pow2(xh(s% i_w,k-1))
            ! Zero kinetic energy cannot carry entropy flux. Retain the entropy variance.
            if (etrb_face == 0d0) xh(s% i_Pi,k) = 0d0
         end do
         if (minval(xh(s% i_Phi,1:nz)) < 0d0) then
            ierr = -1
            s% retry_message = 'negative RSP2 Phi after remesh'
         end if
      end subroutine remesh_rsp2_moments


      function rsp2_buoyancy_face(s, k, ierr) result(buoyancy)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: buoyancy, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face

         ierr = 0
         buoyancy = 0d0
         if (rsp2_zero_moments(s,k)) return
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         ! -(1/rho)*dP/dr times the isobaric expansion per unit entropy.
         buoyancy = (ChiT_face/(ChiRho_face*Cp_face))* &
            (4d0*pi*pow2(wrap_r_00(s,k))/(0.5d0*(s% dm(k-1) + s% dm(k))))* &
            (wrap_Peos_00(s,k) - wrap_Peos_m1(s,k))
      end function rsp2_buoyancy_face


      function rsp2_moment_source(s, k, divide_by_w, ierr) result(source)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: divide_by_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: source, buoyancy, Pi_face, etrb_face, w_face, w_cell
         real(dp) :: alfa, beta
         integer :: j

         ierr = 0
         source = 0d0
         w_cell = wrap_w_00(s,k)
         do j=k,min(k+1,s% nz)
            if (rsp2_zero_moments(s,j)) cycle
            call get_RSP2_alfa_beta_face_weights(s,j,alfa,beta)
            etrb_face = alfa*pow2(wrap_w_00(s,j)) + beta*pow2(wrap_w_m1(s,j))
            Pi_face = wrap_Pi_00(s,j)
            if (etrb_face%val == 0d0) then
               if (Pi_face%val == 0d0) cycle
               ierr = -1
               s% retry_message = 'RSP2 entropy flux without turbulent energy'
               return
            end if
            w_face = sqrt(etrb_face)
            buoyancy = rsp2_buoyancy_face(s,j,ierr)
            if (ierr /= 0) return
            if (j > k) then
               w_face = shift_p1(w_face)
               buoyancy = shift_p1(buoyancy)
               Pi_face = shift_p1(Pi_face)
            end if
            if (divide_by_w) then
               source = source + 0.5d0*buoyancy*div_by_w(Pi_face,w_face)
            else
               source = source + 0.5d0*div_by_w(w_cell,w_face)*buoyancy*Pi_face
            end if
         end do
      end function rsp2_moment_source


      subroutine rsp2_moment_rhs(s, k, Pi_rhs, Phi_rhs, ierr, use_time_centering)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Pi_rhs, Phi_rhs
         integer, intent(out) :: ierr
         logical, intent(in), optional :: use_time_centering
         type(auto_diff_real_star_order1) :: Pi_face, Phi_face, &
            etrb_face, w_face, Cp_face, T_face, rho_face, kap_face, Lambda_face, &
            entropy_gradient, inverse_rad_time, buoyancy, inverse_dr, strain, &
            P_face, ChiRho_face, ChiT_face, grad_ad, gradL
         real(dp) :: alfa, beta

         ierr = 0
         Pi_rhs = 0d0
         Phi_rhs = 0d0
         if (rsp2_zero_moments(s,k)) return
         call get_RSP2_alfa_beta_face_weights(s,k,alfa,beta)
         etrb_face = alfa*pow2(wrap_w_00(s,k)) + beta*pow2(wrap_w_m1(s,k))
         w_face = 0d0
         if (etrb_face%val > 0d0) w_face = sqrt(etrb_face)
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         Lambda_face = get_TDC_mixing_length_face(s,k,ierr)
         if (ierr /= 0) return
         if (Lambda_face%val <= 0d0) then
            ierr = -1
            return
         end if
         inverse_dr = 4d0*pi*pow2(wrap_r_00(s,k))*rho_face/(0.5d0*(s% dm(k-1) + s% dm(k)))
         call get_rsp2_thermal_gradient(s, k, grad_ad, gradL, entropy_gradient, ierr, use_time_centering)
         if (ierr /= 0) return
         buoyancy = rsp2_buoyancy_face(s,k,ierr)
         if (ierr /= 0) return
         inverse_rad_time = 4d0*boltz_sigma*pow2(s% RSP2_alfar*x_GAMMAR)*pow3(T_face)/ &
            (Cp_face*kap_face*pow2(rho_face)*pow2(Lambda_face))
         if (s% u_flag) then
            strain = inverse_dr*(wrap_u_m1(s,k) - wrap_u_00(s,k))
         else
            ! Interpolate the two cell velocity gradients to their common face.
            strain = alfa*(wrap_v_00(s,k) - wrap_v_p1(s,k))/(wrap_r_00(s,k) - wrap_r_p1(s,k)) + &
               beta*(wrap_v_m1(s,k) - wrap_v_00(s,k))/(wrap_r_m1(s,k) - wrap_r_00(s,k))
         end if
         Pi_face = wrap_Pi_00(s,k)
         Phi_face = wrap_Phi_00(s,k)
         Pi_rhs = (2d0/3d0)*etrb_face*entropy_gradient + buoyancy*Phi_face - &
            (s% RSP2_alfa_pi*x_ALFAPI*w_face/Lambda_face + inverse_rad_time + strain)*Pi_face
         Phi_rhs = 2d0*entropy_gradient*Pi_face - &
            (s% RSP2_alfa_phi*x_ALFAPHI*w_face/Lambda_face + 2d0*inverse_rad_time)*Phi_face
      end subroutine rsp2_moment_rhs


      subroutine do1_rsp2_moment_eqns(s, k, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         integer :: i_flux, i_variance
         type(auto_diff_real_star_order1) :: Pi_rhs, Phi_rhs, Pi_resid, Phi_resid
         real(dp) :: Pi_scale, Phi_scale

         ierr = 0
         if (.not. s% RSP2_3equation_flag) return
         i_flux = s% i_Pi
         i_variance = s% i_Phi
         Pi_scale = s% Pi_scale(k)
         Phi_scale = s% Phi_scale(k)
         Pi_resid = wrap_Pi_00(s,k)
         Phi_resid = wrap_Phi_00(s,k)
         if (.not. rsp2_zero_moments(s,k) .and. &
               .not. (rsp2_dormant_moments(s,k) .and. &
                  s% xh_start(i_flux,k) == 0d0 .and. s% xh_start(i_variance,k) == 0d0)) then
            call rsp2_moment_rhs(s,k,Pi_rhs,Phi_rhs,ierr)
            if (ierr /= 0) return
            ! Use the same implicit source time level as COUPL in the w equation.
            Pi_resid = Pi_resid - s% xh_start(i_flux,k) - s% dt*Pi_rhs
            Phi_resid = Phi_resid - s% xh_start(i_variance,k) - s% dt*Phi_rhs
         end if
         Pi_resid = Pi_resid/Pi_scale
         Phi_resid = Phi_resid/Phi_scale
         s% equ(i_flux,k) = Pi_resid%val
         s% equ(i_variance,k) = Phi_resid%val
         call save_eqn_residual_info(s,k,nvar,i_flux,Pi_resid,'dPi_dt',ierr)
         if (ierr /= 0) return
         call save_eqn_residual_info(s,k,nvar,i_variance,Phi_resid,'dPhi_dt',ierr)
      end subroutine do1_rsp2_moment_eqns


      subroutine init_rsp2_moments(s, Lc_old, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: Lc_old(:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: alfa, beta, etrb_face
         type(auto_diff_real_star_order1) :: &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, &
            gradL, entropy_gradient

         ierr = 0
         do k=1,s% nz
            s% Pi(k) = 0d0
            s% Phi(k) = 0d0
            if (.not. rsp2_zero_moments(s,k)) then
               call get_rsp2_thermal_gradient(s,k,grad_ad,gradL,entropy_gradient,ierr)
               if (ierr /= 0) return
               ! Do not infer finite entropy fluctuations from a quiet stable layer.
               if (entropy_gradient%val > 0d0 .and. Lc_old(k) > 0d0) then
                  call get_RSP2_alfa_beta_face_weights(s,k,alfa,beta)
                  etrb_face = alfa*pow2(s% w(k)) + beta*pow2(s% w(k-1))
                  call get_rsp2_face_eos( &
                     s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
                  if (ierr /= 0) return
                  ! Preserve the old heat flux; seed a fully correlated parcel variance.
                  s% Pi(k) = Lc_old(k)/(4d0*pi*pow2(s% r(k))*rho_face%val*T_face%val)
                  if (etrb_face == 0d0 .and. s% Pi(k) /= 0d0) then
                     ierr = -1
                     s% retry_message = 'cannot initialize RSP2 moments without turbulent energy'
                     return
                  end if
                  if (etrb_face > 0d0) s% Phi(k) = &
                     1.5d0*pow2(s% Pi(k))/etrb_face
               end if
            end if
            s% xh(s% i_Pi,k) = s% Pi(k)
            s% xh(s% i_Phi,k) = s% Phi(k)
         end do
      end subroutine init_rsp2_moments


      logical function rsp2_local_w_equation(s, k) result(local)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         ! Use the energy represented by the accepted state, including underflow to zero.
         local = s% RSP2_source_seed == 0d0 .and. s% RSP2_alfat == 0d0 .and. &
            get_etrb_start(s,k) == 0d0 .and. (.not. s% u_flag .or. s% RSP2_alfam == 0d0)
      end function rsp2_local_w_equation


      function div_by_w(x, w) result(quotient)
         type(auto_diff_real_star_order1), intent(in) :: x, w
         type(auto_diff_real_star_order1) :: quotient

         ! Avoid 1/w**2, which can overflow when the quotient derivatives are finite.
         quotient%val = x%val/w%val
         quotient%d1Array = (x%d1Array - quotient%val*w%d1Array)/w%val
      end function div_by_w


      subroutine do1_turbulent_energy_eqn(s, k, nvar, ierr)
         use star_utils, only: calc_Ptrb_ad_tw, set_energy_eqn_scal, save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            d_turbulent_energy_ad, Ptrb_dV_ad, dt_C_ad, dt_Eq_ad, &
            source_div_w_ad, D_div_w_ad, Dr_div_w_ad, Eq_div_w_ad, C_div_w_ad, &
            Ptrb_dV_div_w_ad
         type(auto_diff_real_star_order1) :: w_00, Ptrb_div_etrb, dV_ad
         type(auto_diff_real_star_order1) :: tst, resid_ad, dt_dLt_dm_ad
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         logical :: non_turbulent_cell, positive_branch, test_partials
         real(dp) :: residual, scal, P_theta, Ptrb_start
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         w_00 = wrap_w_00(s,k)

         non_turbulent_cell = &
            s% mixing_length_alpha == 0d0 .or. &
            k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
         if (.not. s% RSP2_flag) then
            resid_ad = w_00 - s% w_start(k)  ! just hold w constant when not using RSP2
         else if (non_turbulent_cell) then
            resid_ad = w_00/s% csound(k)  ! make w = 0
         else
            call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)
            positive_branch = rsp2_local_w_equation(s,k)
            if (positive_branch .or. (s% RSP2_3equation_flag .and. &
                  w_00%val > 0d0 .and. w_00%val < s% csound(k))) then
               source_div_w_ad = compute_Source_div_w(s, k, ierr)
               if (ierr /= 0) return
               if (.not. s% RSP2_3equation_flag .and. s% RSP2_source_seed /= 0d0) &
                  source_div_w_ad = source_div_w_ad + div_by_w(s% RSP2_source_seed*source_div_w_ad,w_00)
               Eq_div_w_ad = compute_Eq_div_w_cell(s, k, ierr)
               if (ierr /= 0) return

               ! Cancel known w factors before differentiating the divided row.
               call calc_Ptrb_ad_tw(s, k, Ptrb_dV_ad, Ptrb_div_etrb, ierr)
               if (ierr /= 0) return
               dV_ad = 1d0/wrap_d_00(s,k) - 1d0/s% rho_start(k)
               P_theta = 1d0
               if (s% using_velocity_time_centering .and. &
                     s% include_P_in_velocity_time_centering) &
                  P_theta = s% P_theta_for_velocity_time_centering
               Ptrb_dV_div_w_ad = P_theta*Ptrb_div_etrb*w_00*dV_ad
               D_div_w_ad = compute_D_div_w(s, k, ierr)
               if (ierr /= 0) return
               Dr_div_w_ad = compute_Dr_div_w(s, k, ierr)
               if (ierr /= 0) return
               C_div_w_ad = source_div_w_ad - D_div_w_ad - Dr_div_w_ad

               esum_ad = w_00
               esum_ad = esum_ad + Ptrb_dV_div_w_ad
               esum_ad = esum_ad - s% dt*C_div_w_ad
               esum_ad = esum_ad - s% dt*Eq_div_w_ad
               if (.not. positive_branch) then
                  Ptrb_start = s% RSP2_alfap*(2d0/3d0)*s% rho_start(k)*get_etrb_start(s,k)
                  esum_ad = esum_ad + div_by_w( &
                     (1d0 - P_theta)*Ptrb_start*dV_ad - get_etrb_start(s,k), w_00)
                  call setup_dt_dLt_dm_ad(ierr, .true.); if (ierr /= 0) return
                  esum_ad = esum_ad + dt_dLt_dm_ad
               end if
               resid_ad = esum_ad
               if (positive_branch .and. w_00%val <= resid_ad%val) resid_ad = w_00
               resid_ad = resid_ad*s% csound(k)*scal/s% dt
            else
               call setup_d_turbulent_energy(ierr); if (ierr /= 0) return  ! erg g^-1 = cm^2 s^-2
               call setup_Ptrb_dV_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_C_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_Eq_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
               esum_ad = d_turbulent_energy_ad
               esum_ad = esum_ad + Ptrb_dV_ad
               esum_ad = esum_ad + dt_dLt_dm_ad
               esum_ad = esum_ad - dt_C_ad
               esum_ad = esum_ad - dt_Eq_ad  ! erg g^-1
               resid_ad = esum_ad
               if (w_00%val == 0d0 .and. resid_ad%val == 0d0 .and. get_etrb_start(s,k) == 0d0) then
                  resid_ad = w_00*s% csound(k)
               end if
               resid_ad = resid_ad*scal/s%dt
            end if

         end if

         residual = resid_ad%val
         s% equ(s% i_detrb_dt, k) = residual

         if (test_partials) then
            tst = residual
            s% solver_test_partials_val = tst%val
            if (s% solver_iter == 12) &
               write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k), tst%val
         end if

         call save_eqn_residual_info(s, k, nvar, s% i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)     ! xi0 good , xi1 partial 0, xi2 good.  Af horrible.'
            write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k)/ln10, tst%val
         end if

         contains

         subroutine setup_d_turbulent_energy(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            ierr = 0
            d_turbulent_energy_ad = wrap_etrb_00(s,k) - get_etrb_start(s,k)
         end subroutine setup_d_turbulent_energy

         ! Ptrb_dV_ad = Ptrb_ad*dV_ad
         subroutine setup_Ptrb_dV_ad(ierr)  ! erg g^-1
            use star_utils, only: calc_Ptrb_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Ptrb_ad, PT0, dV_ad, d_00
            call calc_Ptrb_ad_tw(s, k, Ptrb_ad, PT0, ierr)
            if (ierr /= 0) return
            d_00 = wrap_d_00(s,k)
            dV_ad = 1d0/d_00 - 1d0/s% rho_start(k)
            Ptrb_dV_ad = Ptrb_ad*dV_ad  ! erg cm^-3 cm^-3 g^-1 = erg g^-1
         end subroutine setup_Ptrb_dV_ad

         subroutine setup_dt_dLt_dm_ad(ierr, divide_by_w)
            integer, intent(out) :: ierr
            logical, intent(in), optional :: divide_by_w
            type(auto_diff_real_star_order1) :: Lt_00, Lt_p1
            real(dp) :: L_theta, Lt_start_p1
            include 'formats'
            ierr = 0
            if (s% using_velocity_time_centering .and. &
                     s% include_L_in_velocity_time_centering) then
               L_theta = s% L_theta_for_velocity_time_centering
            else
               L_theta = 1d0
            end if
            if (present(divide_by_w)) then
               if (divide_by_w) then
                  Lt_00 = compute_Lt(s, k, ierr, k)
                  if (ierr /= 0) return
                  Lt_p1 = compute_Lt(s, k+1, ierr, k)
                  if (ierr /= 0) return
                  if (k < s% nz) Lt_p1 = shift_p1(Lt_p1)
                  dt_dLt_dm_ad = L_theta*(Lt_00 - Lt_p1)*s% dt/s% dm(k)
                  if (L_theta /= 1d0) then
                     Lt_start_p1 = 0d0
                     if (k < s% nz) Lt_start_p1 = s% Lt_start(k+1)
                     Lt_00 = (1d0 - L_theta)*(s% Lt_start(k) - Lt_start_p1)
                     dt_dLt_dm_ad = dt_dLt_dm_ad + div_by_w(Lt_00,w_00)*s% dt/s% dm(k)
                  end if
                  return
               end if
            end if
            Lt_00 = L_theta*s% Lt_ad(k) + (1d0 - L_theta)*s% Lt_start(k)
            if (k == s% nz) then
               Lt_p1 = 0d0
            else
               Lt_p1 = L_theta*shift_p1(s% Lt_ad(k+1)) + (1d0 - L_theta)*s% Lt_start(k+1)
               if (ierr /= 0) return
            end if
            dt_dLt_dm_ad = (Lt_00 - Lt_p1)*s%dt/s%dm(k)
         end subroutine setup_dt_dLt_dm_ad

         subroutine setup_dt_C_ad(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: C
            ierr = 0
            C = s% COUPL_ad(k)  ! compute_C(s, k, ierr) ! erg g^-1 s^-1
            dt_C_ad = s%dt*C
         end subroutine setup_dt_C_ad

         subroutine setup_dt_Eq_ad(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Eq_cell
            ierr = 0
            Eq_cell = s% Eq_ad(k)  ! compute_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
            dt_Eq_ad = s%dt*Eq_cell
         end subroutine setup_dt_Eq_ad

      end subroutine do1_turbulent_energy_eqn





      function compute_PII_face(s, k, ierr) result(PII_face)  ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face, w_face
         real(dp) :: alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            PII_face = 0d0
            return
         end if
         if (k == 1 .or. k == s% nz .or. s% mixing_length_alpha == 0d0 .or. &
               k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
               k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            PII_face = 0d0
            s% PII(k) = 0d0
            s% PII_ad(k) = 0d0
            return
         end if
         Y_face = s% Y_face_ad(k)
         if (ierr /= 0) return
         if (s% RSP2_3equation_flag) then
            ! PII remains a diagnostic; the moment equations use <v_r*s> directly.
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            w_face = alfa*wrap_w_00(s,k) + beta*wrap_w_m1(s,k)
            PII_face = 0d0
            if (w_face%val > 0d0) PII_face = div_by_w(wrap_Pi_00(s,k),w_face)
         else
            PII_face = compute_PII_from_Y(s, k, Y_face, ierr)
            if (ierr /= 0) return
         end if

         s% PII(k) = PII_face%val
         s% PII_ad(k) = PII_face
         if (k == -2 .and. s% PII(k) < 0d0) then
            write(*,2) 's% PII(k)', k, s% PII(k)
            write(*,2) 'Y_face', k, Y_face%val
            call mesa_error(__FILE__,__LINE__,'compute_PII_face')
         end if
      end function compute_PII_face


      function compute_PII_from_Y(s, k, Y_face, ierr) result(PII_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: Y_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: PII_face
         type(auto_diff_real_star_order1) :: &
            Cp_00, Cp_m1, Cp_face, Lambda_face, Hp_face
         real(dp) :: alfa, beta

         ierr = 0
         Cp_00 = wrap_Cp_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         Cp_face = alfa*Cp_00 + beta*Cp_m1  ! ergs g^-1 K^-1
         Lambda_face = get_TDC_mixing_length_face(s, k, ierr)
         if (ierr /= 0) return
         Hp_face = get_TDC_Hp_face(s, k, ierr)
         if (ierr /= 0) return
         PII_face = x_ALFAS*(Lambda_face/Hp_face)*Cp_face*Y_face
      end function compute_PII_from_Y


      function compute_Source(s, k, ierr) result(Source)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: w_00, Source_div_w
         include 'formats'
         ierr = 0
         if (s% RSP2_3equation_flag) then
            Source = rsp2_moment_source(s,k,.false.,ierr)
            if (ierr /= 0) return
            s% SOURCE(k) = Source%val
            return
         end if
         w_00 = wrap_w_00(s, k)
         Source_div_w = compute_Source_div_w(s, k, ierr)
         if (ierr /= 0) return
         Source = (w_00 + s% RSP2_source_seed)*Source_div_w

         if (k==-109) then
            write(*,3) 'RSP2 Source w source_div_w', k, s% solver_iter, &
               Source%val, w_00%val, Source_div_w%val
         end if
         s% SOURCE(k) = Source%val

      end function compute_Source


      function compute_Source_div_w(s, k, ierr) result(Source_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            T_00, d_00, Peos_00, Cp_00, chiT_00, chiRho_00, QQ_00, &
            Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell, &
            P_QQ_div_Cp
         ierr = 0
         if (s% RSP2_3equation_flag) then
            Source_div_w = rsp2_moment_source(s,k,.true.,ierr)
            return
         end if
         T_00 = wrap_T_00(s, k)
         d_00 = wrap_d_00(s, k)
         Peos_00 = wrap_Peos_00(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)

         Hp_face_00 = get_TDC_Hp_face(s, k, ierr)
         if (ierr /= 0) return
         PII_face_00 = s% PII_ad(k)

         if (k == s% nz) then
            PII_div_Hp_cell = PII_face_00/Hp_face_00
         else
            Hp_face_p1 = shift_p1(get_TDC_Hp_face(s, k+1, ierr))
            if (ierr /= 0) return
            PII_face_p1 = shift_p1(s% PII_ad(k+1))
            PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         end if

         ! Peos_00*QQ_00/Cp_00 = grad_ad if all perfect.
         !grad_ad_00 = wrap_grad_ad_00(s, k)
         P_QQ_div_Cp = Peos_00*QQ_00/Cp_00  ! use this to be same as RSP
         Source_div_w = PII_div_Hp_cell*T_00*P_QQ_div_Cp
      end function compute_Source_div_w


      function compute_D(s, k, ierr) result(D)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: D_div_w, w_00
         ierr = 0
         w_00 = wrap_w_00(s,k)
         D_div_w = compute_D_div_w(s, k, ierr)
         if (ierr /= 0) return
         D = w_00*D_div_w
         s% DAMP(k) = D%val
      end function compute_D


      function compute_D_div_w(s, k, ierr) result(D_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lambda_cell, w_00
         ierr = 0
         if (s% mixing_length_alpha == 0d0) then
            D_div_w = 0d0
         else
            Lambda_cell = get_TDC_mixing_length_cell(s, k, ierr)
            if (ierr /= 0) return
            w_00 = wrap_w_00(s,k)
            D_div_w = (s% RSP2_alfad*x_CEDE)* &
               pow2(w_00)/Lambda_cell
         end if
      end function compute_D_div_w


      function compute_Dr(s, k, ierr) result(Dr)  ! erg g^-1 s^-1 = cm^2 s^-3
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Dr_div_w, w_00
         ierr = 0
         w_00 = wrap_w_00(s,k)
         Dr_div_w = compute_Dr_div_w(s, k, ierr)
         if (ierr /= 0) return
         Dr = w_00*Dr_div_w
         s% DAMPR(k) = Dr%val
      end function compute_Dr


      function compute_Dr_div_w(s, k, ierr) result(Dr_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Cp_00, kap_00, Lambda_cell, POM2
         real(dp) :: gammar, alpha, POM
         ierr = 0
         Dr_div_w = 0d0
         if (s% RSP2_3equation_flag) return
         alpha = s% mixing_length_alpha
         gammar = s% RSP2_alfar*x_GAMMAR
         if (gammar == 0d0 .or. alpha == 0d0) then
            Dr_div_w = 0d0
            return
         end if
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Lambda_cell = get_TDC_mixing_length_cell(s, k, ierr)
         if (ierr /= 0) return
         POM = 4d0*boltz_sigma*pow2(gammar)  ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00)
         Dr_div_w = w_00*POM*POM2/pow2(Lambda_cell)
      end function compute_Dr_div_w


      function compute_d_v_div_r(s, k, use_time_centering) result(d_v_div_r)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: use_time_centering
         type(auto_diff_real_star_order1) :: d_v_div_r, v_00, v_p1, r_00, r_p1

         v_00 = wrap_v_00(s, k)
         v_p1 = wrap_v_p1(s, k)
         r_00 = wrap_r_00(s, k)
         r_p1 = wrap_r_p1(s, k)
         if (use_time_centering) then
            if (s% using_velocity_time_centering .or. &
                  .not. s% use_P_d_1_div_rho_form_of_work) then
               v_00 = 0.5d0*(v_00 + s% v_start(k))
               if (k < s% nz) v_p1 = 0.5d0*(v_p1 + s% v_start(k+1))
            end if
            r_00 = wrap_opt_time_center_r_00(s, k)
            r_p1 = wrap_opt_time_center_r_p1(s, k)
         end if
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1
      end function compute_d_v_div_r


      function compute_Chi_div_w_cell(s, k, ierr) result(Chi_div_w)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_div_w, Lambda_cell, r6_cell

         ierr = 0
         Chi_div_w = 0d0
         if (s% mixing_length_alpha == 0d0 .or. s% RSP2_alfam == 0d0 .or. &
               k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
               k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) return
         Lambda_cell = get_TDC_mixing_length_cell(s, k, ierr)
         if (ierr /= 0) return
         r6_cell = 0.5d0*(pow6(wrap_r_00(s,k)) + pow6(wrap_r_p1(s,k)))
         Chi_div_w = (16d0/3d0)*pi*s% RSP2_alfam*pow2(wrap_d_00(s,k))* &
            r6_cell*Lambda_cell*compute_d_v_div_r(s, k, .false.)/s% dm(k)
      end function compute_Chi_div_w_cell


      function compute_Chi_cell(s, k, ierr) result(Chi_cell)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_cell

         Chi_cell = compute_Chi_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         Chi_cell = Chi_cell*wrap_w_00(s, k)
      end function compute_Chi_cell


      function compute_Eq_div_w_cell(s, k, ierr) result(Eq_div_w)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_div_w, Chi_div_w, Eq_p1

         ierr = 0
         Eq_div_w = 0d0
         if (s% mixing_length_alpha == 0d0 .or. s% RSP2_alfam == 0d0) return
         if (s% u_flag) then
            Eq_div_w = 0.5d0*compute_Eq_face(s, k, ierr, k)
            if (ierr /= 0) return
            Eq_p1 = compute_Eq_face(s, k+1, ierr, k)
            if (ierr /= 0) return
            if (k < s% nz) Eq_p1 = shift_p1(Eq_p1)
            Eq_div_w = Eq_div_w + 0.5d0*Eq_p1
            return
         end if
         Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         Eq_div_w = 4d0*pi*Chi_div_w*compute_d_v_div_r(s, k, .true.)/s% dm(k)
      end function compute_Eq_div_w_cell


      function compute_Eq_cell(s, k, ierr) result(Eq_cell)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_cell, Eq_p1

         ierr = 0
         Eq_cell = 0d0
         if (s% mixing_length_alpha == 0d0 .or. s% RSP2_alfam == 0d0) then
            Eq_cell = 0d0
         else if (s% u_flag) then
            ! Each cell receives half of each adjacent face's specific heating.
            Eq_cell = 0.5d0*compute_Eq_face(s, k, ierr)
            if (ierr /= 0) return
            Eq_p1 = compute_Eq_face(s, k+1, ierr)
            if (ierr /= 0) return
            if (k < s% nz) Eq_p1 = shift_p1(Eq_p1)
            Eq_cell = Eq_cell + 0.5d0*Eq_p1
         else if (s% v_flag) then
            Eq_cell = compute_Eq_div_w_cell(s, k, ierr)
            if (ierr /= 0) return
            Eq_cell = Eq_cell*wrap_w_00(s, k)
         end if
         s% Eq(k) = Eq_cell%val
         s% Eq_ad(k) = Eq_cell
      end function compute_Eq_cell


      function compute_Uq_face(s, k, ierr) result(Uq_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Uq_face, Chi_00, Chi_m1, r_00
         real(dp) :: dm_face

         ierr = 0
         Uq_face = 0d0
         Chi_00 = 0d0
         if (s% mixing_length_alpha /= 0d0 .and. s% RSP2_alfam /= 0d0) then
            Chi_00 = compute_Chi_cell(s, k, ierr)
            if (ierr /= 0) return
            Chi_m1 = 0d0
            if (k > 1) then
               Chi_m1 = shift_m1(compute_Chi_cell(s, k-1, ierr))
               if (ierr /= 0) return
            end if
            r_00 = wrap_opt_time_center_r_00(s, k)
            if (s% use_mass_corrections) then
               dm_face = 0.5d0*s% dm(k)*s% mass_correction(k)
               if (k > 1) dm_face = dm_face + &
                  0.5d0*s% dm(k-1)*s% mass_correction(k-1)
            else
               dm_face = 0.5d0*s% dm(k)
               if (k > 1) dm_face = dm_face + 0.5d0*s% dm(k-1)
            end if
            ! Mask stresses, retaining their force at the edge of a turbulent region.
            Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*dm_face)
         end if
         s% Chi(k) = Chi_00%val
         s% Chi_ad(k) = Chi_00
         s% Uq(k) = Uq_face%val
      end function compute_Uq_face


      function compute_d_u_div_r_face(s, k, use_time_centering) result(d_u_div_r)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: use_time_centering
         type(auto_diff_real_star_order1) :: d_u_div_r, u_out, u_in
         real(dp) :: r_out, r_in

         d_u_div_r = 0d0
         if (k == 1) return
         if (k > s% nz) then
            if (s% R_center <= 0d0) return
            ! Boundary results use the innermost cell's AD slots.
            u_out = wrap_u_00(s, s% nz)
            u_in = s% v_center
            r_out = s% rmid_start(s% nz)
            r_in = s% R_center
            if (use_time_centering .and. (s% using_velocity_time_centering .or. &
                  .not. s% use_P_d_1_div_rho_form_of_work)) &
               u_out = 0.5d0*(u_out + s% u_start(s% nz))
         else
            u_out = wrap_u_m1(s, k)
            u_in = wrap_u_00(s, k)
            r_out = s% rmid_start(k-1)
            r_in = s% rmid_start(k)
            if (use_time_centering .and. (s% using_velocity_time_centering .or. &
                  .not. s% use_P_d_1_div_rho_form_of_work)) then
               u_out = 0.5d0*(u_out + s% u_start(k-1))
               u_in = 0.5d0*(u_in + s% u_start(k))
            end if
         end if
         d_u_div_r = u_out/r_out - u_in/r_in
      end function compute_d_u_div_r_face


      function compute_Chi_face(s, k, ierr, k_div_w) result(Chi_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         integer, intent(in), optional :: k_div_w  ! divide by this cell's w
         type(auto_diff_real_star_order1) :: Chi_face, rho_face, r_face, Lambda_face, w_face
         real(dp) :: alfa, beta, dm_face

         ierr = 0
         Chi_face = 0d0
         if (k == 1 .or. s% mixing_length_alpha == 0d0 .or. s% RSP2_alfam == 0d0) return
         if (k > s% nz) then
            if (.not. s% TDC_include_inner_boundary_eddy_viscosity .or. &
                  s% R_center <= 0d0 .or. int(s% nz/s% RSP2_nz_div_IBOTOM) > 0) return
            rho_face = wrap_d_00(s, s% nz)
            r_face = s% R_center
            Lambda_face = get_TDC_mixing_length_face(s, s% nz, ierr)
            w_face = wrap_w_00(s, s% nz)
            if (present(k_div_w)) w_face = 1d0
            dm_face = 0.5d0*s% dm(s% nz)
         else
            if (k <= s% RSP2_num_outermost_cells_forced_nonturbulent + 1 .or. &
                  k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) return
            rho_face = get_rho_face(s, k)
            r_face = wrap_r_00(s, k)
            Lambda_face = get_TDC_mixing_length_face(s, k, ierr)
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            w_face = alfa*wrap_w_00(s,k) + beta*wrap_w_m1(s,k)
            if (present(k_div_w)) then
               ! Cancel the local contribution before forming w_face/w_cell.
               if (k_div_w == k) then
                  w_face = alfa + div_by_w(beta*wrap_w_m1(s,k),wrap_w_00(s,k))
               else
                  w_face = div_by_w(alfa*wrap_w_00(s,k),wrap_w_m1(s,k)) + beta
               end if
            end if
            dm_face = 0.5d0*(s% dm(k) + s% dm(k-1))
         end if
         if (ierr /= 0) return
         Chi_face = (16d0/3d0)*pi*s% RSP2_alfam*pow2(rho_face)*pow6(r_face)* &
            Lambda_face*w_face*compute_d_u_div_r_face(s, k, .false.)/dm_face
      end function compute_Chi_face


      function compute_Eq_face(s, k, ierr, k_div_w) result(Eq_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         integer, intent(in), optional :: k_div_w
         type(auto_diff_real_star_order1) :: Eq_face, Chi_face
         real(dp) :: dm_face

         ierr = 0
         Eq_face = 0d0
         if (k == 1) return
         Chi_face = compute_Chi_face(s, k, ierr, k_div_w)
         if (ierr /= 0) return
         if (k > s% nz) then
            dm_face = 0.5d0*s% dm(s% nz)
         else
            dm_face = 0.5d0*(s% dm(k) + s% dm(k-1))
         end if
         Eq_face = 4d0*pi*Chi_face*compute_d_u_div_r_face(s, k, .true.)/dm_face
      end function compute_Eq_face


      function compute_Uq_dm_cell(s, k, ierr) result(Uq_dm_cell)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Uq_dm_cell, Chi_00, Chi_p1

         ierr = 0
         Uq_dm_cell = 0d0
         Chi_00 = 0d0
         if (s% mixing_length_alpha /= 0d0 .and. s% RSP2_alfam /= 0d0) then
            Chi_00 = compute_Chi_face(s, k, ierr)
            if (ierr /= 0) return
            Chi_p1 = compute_Chi_face(s, k+1, ierr)
            if (ierr /= 0) return
            if (k < s% nz) Chi_p1 = shift_p1(Chi_p1)
            ! Return force; the Riemann momentum row divides by dm.
            Uq_dm_cell = 4d0*pi*(Chi_00 - Chi_p1)/s% rmid_start(k)
         end if
         s% Chi(k) = Chi_00%val
         s% Chi_ad(k) = Chi_00
         s% Uq(k) = Uq_dm_cell%val/s% dm(k)
      end function compute_Uq_dm_cell


      function compute_C(s, k, ierr) result(C)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: C
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source, D, Dr
         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            if (k >= 1 .and. k <= s% nz) then
               s% SOURCE(k) = 0d0
               s% DAMP(k) = 0d0
               s% DAMPR(k) = 0d0
               s% COUPL(k) = 0d0
               s% COUPL_ad(k) = 0d0
            end if
            C = 0d0
            return
         end if
         Source = compute_Source(s, k, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         s% COUPL(k) = C%val
         s% COUPL_ad(k) = C
      end function compute_C


      function compute_L_face(s, k, ierr) result(L_face)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr, Lc, Lt
         call compute_L_terms(s, k, L_face, Lr, Lc, Lt, ierr)
      end function compute_L_face


      subroutine compute_L_terms(s, k, L, Lr, Lc, Lt, ierr)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L, Lr, Lc, Lt
         type(accurate_auto_diff_real_star_order1) :: L_sum
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
         L_sum = Lr
         L_sum = L_sum + Lc
         L_sum = L_sum + Lt
         L = L_sum
         s% Lr_ad(k) = Lr
         s% Lc_ad(k) = Lc
         s% Lt_ad(k) = Lt
      end subroutine compute_L_terms


      function compute_Lr(s, k, ierr) result(Lr)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_00, T400, Erad, Lrad_coeff, gradT
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lr = s% L_center
         else
            r_00 = wrap_r_00(s,k)  ! not time centered
            area = 4d0*pi*pow2(r_00)
            T_00 = wrap_T_00(s,k)
            T400 = pow4(T_00)
            if (k == 1) then
               if (s% RSP2_use_L_eqn_at_surface) then
                  Erad = crad*T400
                  Lr = s% RSP2_Lsurf_factor*area*clight*Erad
               else
                  Lr = wrap_L_00(s,k)
               end if
            else
               Lrad_coeff = compute_Lrad_coeff(s, k, ierr)
               if (ierr /= 0) return
               gradT = s% gradT_ad(k)
               Lr = Lrad_coeff*gradT
            end if
         end if
         s% Lr(k) = Lr%val
      end function compute_Lr


      function compute_Lrad_coeff(s, k, ierr) result(Lrad_coeff)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lrad_coeff
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_face, rho_face, kap_face, Hp_face, krad
         real(dp) :: alfa, beta
         include 'formats'
         ierr = 0

         if (s% RSP2_3equation_flag) then
            Lrad_coeff = get_rsp2_Lrad_coeff(s, k, ierr)
            return
         end if

         r_00 = wrap_r_00(s,k)
         area = 4d0*pi*pow2(r_00)
         Hp_face = get_TDC_Hp_face(s, k, ierr)
         if (ierr /= 0) return
         if (k == 1) then
            T_face = wrap_T_00(s,k)
            rho_face = wrap_d_00(s,k)
            kap_face = wrap_kap_00(s,k)
         else
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            T_face = alfa*wrap_T_00(s,k) + beta*wrap_T_m1(s,k)
            rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
            kap_face = alfa*wrap_kap_00(s,k) + beta*wrap_kap_m1(s,k)
         end if
         krad = 4d0*crad*clight*pow3(T_face)/(3d0*kap_face*rho_face)
         Lrad_coeff = area*krad*T_face/Hp_face
      end function compute_Lrad_coeff


      function compute_Lc(s, k, ierr) result(Lc)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lc_div_w_face
         Lc = compute_Lc_terms(s, k, Lc_div_w_face, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc


      function compute_Lc_terms(s, k, Lc_div_w_face, ierr) result(Lc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc, Lc_div_w_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area, &
            T_m1, T_00, d_m1, d_00, w_m1, w_00, T_rho_face, PII_face, w_face, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face
         real(dp) :: ALFAC, ALFAS, alfa, beta
         include 'formats'
         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Lc = 0d0
            Lc_div_w_face = 1
            return
         end if
         r_00 = wrap_r_00(s, k)
         area = 4d0*pi*pow2(r_00)
         if (s% RSP2_3equation_flag) then
            Lc = 0d0
            Lc_div_w_face = 0d0
            if (rsp2_zero_moments(s,k)) return
            call get_rsp2_face_eos( &
               s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
            if (ierr /= 0) return
            Lc = area*rho_face*T_face*wrap_Pi_00(s,k)
            return
         end if
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         w_m1 = wrap_w_m1(s, k)
         w_00 = wrap_w_00(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         PII_face = s% PII_ad(k)  ! compute_PII_face(s, k, ierr)
         w_face = alfa*w_00 + beta*w_m1
         ALFAC = x_ALFAC
         ALFAS = x_ALFAS
         Lc_div_w_face = area*(ALFAC/ALFAS)*T_rho_face*PII_face
         ! units = cm^2 K g cm^-3 ergs g^-1 K^-1 = ergs cm^-1
         Lc = w_face*Lc_div_w_face
         ! units = cm s^-1 ergs cm^-1 = ergs s^-1
         if (k == -458) then
            write(*,2) 'Lc%val', k, Lc%val
            write(*,2) 'w_face%val', k, w_face%val
            write(*,2) 'Lc_div_w_face', k, Lc_div_w_face%val
            write(*,2) 'PII_face%val', k, PII_face%val
            write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k,
            !write(*,2) '', k,
            call mesa_error(__FILE__,__LINE__,'compute_Lc_terms')
         end if
      end function compute_Lc_terms


      function compute_RSP2_gradT(s, k, ierr) result(gradT)
         use turb_support, only: get_TDC_dynamical_gradL
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: gradT, grad_ad, gradL, dynamical_gradL, entropy_gradient
         real(dp) :: alfa, beta

         ierr = 0
         gradT = 0d0
         if (s% RSP2_3equation_flag) then
            call get_rsp2_thermal_gradient(s, k, grad_ad, gradL, entropy_gradient, ierr)
            if (ierr /= 0) return
         else
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            grad_ad = alfa*wrap_grad_ad_00(s,k) + beta*wrap_grad_ad_m1(s,k)
            gradL = grad_ad
            if (s% use_Ledoux_criterion) &
               gradL = gradL + s% gradL_composition_term(k)
            if (k > 1 .and. s% TDC_use_dynamical_gradL) then
               call get_TDC_dynamical_gradL(s, k, gradL, dynamical_gradL, ierr)
               if (ierr /= 0) return
               gradL = dynamical_gradL
            end if
         end if
         s% grada_face_ad(k) = grad_ad
         s% grada_face(k) = grad_ad%val
         s% gradL_ad(k) = gradL
         s% gradL(k) = gradL%val

         gradT = gradL + wrap_Y_00(s, k)
      end function compute_RSP2_gradT


      function compute_Lt(s, k, ierr, k_div_w) result(Lt)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lt
         integer, intent(out) :: ierr
         integer, intent(in), optional :: k_div_w  ! divide by this cell's w
         type(auto_diff_real_star_order1) :: r_00, area2, d_m1, d_00, &
            rho2_face, Lambda_face, w_m1, w_00, w_face, etrb_m1, etrb_00
         real(dp) :: alpha_t, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lt = 0d0
            return
         end if
         alpha_t = s% RSP2_alfat
         if (k == 1 .or. s% mixing_length_alpha == 0d0 .or. alpha_t == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Lt = 0d0
            s% Lt(k) = 0d0
            return
         end if
         r_00 = wrap_r_00(s,k)
         area2 = pow2(4d0*pi*pow2(r_00))
         d_m1 = wrap_d_m1(s,k)
         d_00 = wrap_d_00(s,k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         rho2_face = alfa*pow2(d_00) + beta*pow2(d_m1)
         w_m1 = wrap_w_m1(s,k)
         w_00 = wrap_w_00(s,k)
         w_face = alfa*w_00 + beta*w_m1
         if (present(k_div_w)) then
            if (k_div_w == k) then
               w_face = alfa + div_by_w(beta*w_m1,w_00)
            else
               w_face = div_by_w(alfa*w_00,w_m1) + beta
            end if
         end if
         etrb_m1 = pow2(w_m1)
         etrb_00 = pow2(w_00)
         Lambda_face = get_TDC_mixing_length_face(s, k, ierr)
         if (ierr /= 0) return
         Lt = - alpha_t * area2 * rho2_face * Lambda_face * w_face * (etrb_m1 - etrb_00) / s% dm_bar(k)
         ! units = (cm^4) (g^2 cm^-6) (cm) (cm s^-1) (ergs g^-1) g^-1 = erg s^-1
         if (.not. present(k_div_w)) s% Lt(k) = Lt%val
      end function compute_Lt


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: Lt, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face
         real(dp) :: L_scale, max_L, flux_reference, velocity_reference, area_rho
         include 'formats'
         ierr = 0
         max_L = maxval(abs(s% L(1:s% nz)))
         do k=1,s%nz
            s% Y_face_start(k) = s% Y_face(k)
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_ad(k) = Lt
            s% Lt_start(k) = Lt%val
            s% w_start(k) = s% w(k)
            if (s% RSP2_3equation_flag) then
               call get_rsp2_face_eos( &
                  s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
               if (ierr /= 0) return
               L_scale = max(1d0, abs(s% L(k)), 1d-3*max_L)
               area_rho = 4d0*pi*pow2(s% r(k))*rho_face%val
               flux_reference = L_scale/(area_rho*T_face%val)
               velocity_reference = pow(L_scale/area_rho,1d0/3d0)
               s% Pi_scale(k) = max(abs(s% Pi(k)), flux_reference)
               s% Phi_scale(k) = max(s% Phi(k), &
                  1.5d0*pow2(flux_reference/velocity_reference))
            end if
         end do
      end subroutine set_etrb_start_vars


      subroutine RSP2_adjust_vars_before_call_solver(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, pass, k_lo, k_hi, k_first, k_last, k_step
         real(dp) :: velocity_guess, source_coeff, linear_coeff, available_energy, &
            discr, soln, w_initial, L_theta
         type(auto_diff_real_star_order1) :: source, damping, rad_damping, Eq_cell, &
            Lt_00, Lt_p1, Ptrb, Ptrb_div_etrb, dV, rhs, buoyancy

         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. s% dt <= 0d0) return
         k_lo = s% RSP2_num_outermost_cells_forced_nonturbulent + 1
         k_hi = s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
         if (k_lo > k_hi) return

         ! A remaining entropy variance can generate flux at negligible kinetic energy.
         if (s% RSP2_3equation_flag) then
            do k=max(2,k_lo),min(k_hi,s% nz-1)
               if (rsp2_zero_moments(s,k) .or. s% Phi(k) <= 0d0) cycle
               buoyancy = rsp2_buoyancy_face(s,k,ierr)
               if (ierr /= 0) return
               velocity_guess = s% dt*abs(buoyancy%val)*sqrt(s% Phi(k))
               if (velocity_guess == 0d0) cycle
               ! Include w below roundoff relative to the buoyancy impulse.
               if (max(s% w(k),s% w(k-1)) > epsilon(1d0)*velocity_guess) cycle
               s% w(k) = velocity_guess
               if (k-1 >= k_lo) s% w(k-1) = velocity_guess
               s% Pi(k) = s% Pi(k) + s% dt*buoyancy%val*s% Phi(k)
            end do
         end if

         ! Include neighboring w derivatives when shear heats a quiet region.
         ! A diagonal test alone misses collective production with u velocity.
         do k=k_lo,k_hi
            if (s% w(k) /= 0d0) cycle
            source = compute_Source_div_w(s,k,ierr)
            if (ierr /= 0) return
            Eq_cell = compute_Eq_cell(s,k,ierr)
            if (ierr /= 0) return
            source_coeff = source%val + Eq_cell%d1Array(i_w_m1) + &
               Eq_cell%d1Array(i_w_00) + Eq_cell%d1Array(i_w_p1)
            if (source_coeff <= 0d0) cycle
            s% w(k) = s% dt*source_coeff
         end do

         L_theta = 1d0
         if (s% using_velocity_time_centering .and. s% include_L_in_velocity_time_centering) &
            L_theta = s% L_theta_for_velocity_time_centering
         ! Sweep in both directions so imported energy also seeds empty cells.
         ! These are trial values only; accepted state and time weights are unchanged.
         do pass=1,2
            k_first = k_lo
            k_last = k_hi
            k_step = 1
            if (pass == 2) then
               k_first = k_hi
               k_last = k_lo
               k_step = -1
            end if
            do k=k_first,k_last,k_step
               source = compute_Source(s,k,ierr)
               if (ierr /= 0) return
               damping = compute_D(s,k,ierr)
               if (ierr /= 0) return
               rad_damping = compute_Dr(s,k,ierr)
               if (ierr /= 0) return
               Eq_cell = compute_Eq_cell(s,k,ierr)
               if (ierr /= 0) return
               Lt_00 = compute_Lt(s,k,ierr)
               if (ierr /= 0) return
               Lt_p1 = compute_Lt(s,k+1,ierr)
               if (ierr /= 0) return
               if (k < s% nz) Lt_p1 = shift_p1(Lt_p1)
               Lt_00 = L_theta*Lt_00 + (1d0-L_theta)*s% Lt_start(k)
               if (k < s% nz) Lt_p1 = L_theta*Lt_p1 + (1d0-L_theta)*s% Lt_start(k+1)
               call calc_Ptrb_ad_tw(s,k,Ptrb,Ptrb_div_etrb,ierr)
               if (ierr /= 0) return
               dV = 1d0/wrap_d_00(s,k) - 1d0/s% rho_start(k)
               rhs = s% dt*(source - damping - rad_damping + Eq_cell + (Lt_p1-Lt_00)/s% dm(k)) - Ptrb*dV
               linear_coeff = rhs%d1Array(i_w_00)
               available_energy = pow2(s% w_start(k)) + rhs%val - linear_coeff*s% w(k)
               if (available_energy < 0d0 .or. is_bad(available_energy) .or. is_bad(linear_coeff)) cycle
               discr = sqrt(pow2(linear_coeff) + 4d0*available_energy)
               if (is_bad(discr)) cycle
               soln = 0d0
               if (linear_coeff >= 0d0) then
                  soln = 0.5d0*(linear_coeff + discr)
               else
                  soln = 2d0*available_energy/(discr - linear_coeff)
               end if
               if (soln <= s% w(k) .or. is_bad(soln)) cycle
               w_initial = s% w(k)
               s% w(k) = soln
               if (s% RSP2_report_adjust_w) &
                  write(*,'(a,i7,2es16.7)') 'RSP2 initial w ', k, w_initial, soln
            end do
         end do
      end subroutine RSP2_adjust_vars_before_call_solver

      end module hydro_rsp2
