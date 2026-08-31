! ***********************************************************************
!
!   Copyright (C) 2010-2025  Ebraheem Farag & The MESA Team
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

module tdc_hydro

   use star_private_def
   use const_def, only: dp, boltz_sigma, pi, clight, crad, ln10
   use utils_lib, only: is_bad
   use auto_diff
   use auto_diff_support
   use star_utils
   use reconstructed_face_support, only: &
      get_reconstructed_scale_height_ad, get_reconstructed_hse_scale_height_ad

   implicit none

   private
   public :: &
      compute_tdc_Uq_face, compute_tdc_Eq_cell, compute_tdc_Eq_div_w_face, &
      compute_tdc_Eq_div_w_inner_boundary, &
      get_TDC_alfa_beta_face_weights, get_TDC_mixing_length_face, &
      set_viscosity_vars_TDC, compute_tdc_Uq_dm_cell

contains

   ! This routine is called to initialize eq and uq for TDC.
   subroutine set_viscosity_vars_TDC(s, ierr)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: x
      integer :: k, op_err
      include 'formats'
      ierr = 0
      op_err = 0

      if (.not. (s%v_flag .or. s%u_flag)) then ! set values 0 if not using v_flag or u_flag.
         do k = 1, s%nz
            s%Eq(k) = 0d0; s%Eq_ad(k) = 0d0
            s%Chi(k) = 0d0; s%Chi_ad(k) = 0d0
            s%Uq(k) = 0d0
         end do
         return
      end if

      !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
      do k = 1, s%nz
         if (s%use_face_reconstruction) then
            x = get_TDC_Hp_face(s, k, op_err)
            if (op_err /= 0) then
               !$OMP ATOMIC WRITE
               ierr = op_err
            else if (s%Hp_face(k) <= 0d0) then
               s%Hp_face(k) = x%val
            end if
         else if (s%Hp_face(k) <= 0d0) then
            ! mlt_vars has not updated scale_height yet.
            s%Hp_face(k) = get_scale_height_face_val(s, k)
         end if
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) then
         if (s%report_ierr) write (*, 2) 'failed in set_viscosity_vars_TDC loop 1', s%model_number
         return
      end if
      !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
      do k = 1, s%nz
         if (s% v_flag) then
            x = compute_tdc_Eq_cell(s, k, op_err)
         else
            x = compute_tdc_Eq_div_w_face(s, k, op_err)
         end if
         if (op_err /= 0) then
            !$OMP ATOMIC WRITE
            ierr = op_err
            cycle
         end if
         if (s% v_flag) then
            x = compute_tdc_Uq_face(s, k, op_err)
         else if (s% u_flag) then
            x = compute_tdc_Uq_dm_cell(s, k, op_err)
         end if
         if (op_err /= 0) then
            !$OMP ATOMIC WRITE
            ierr = op_err
         end if
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) then
         if (s%report_ierr) write (*, 2) 'failed in set_viscosity_vars_TDC loop 2', s%model_number
         return
      end if
   end subroutine set_viscosity_vars_TDC

   subroutine get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(out) :: alfa, beta
      ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
      if (k == 1) call mesa_error(__FILE__, __LINE__, 'bad k==1 for get_TDC_alfa_beta_face_weights')
      if (s%TDC_hydro_use_mass_interp_face_values) then
         alfa = s%dq(k - 1)/(s%dq(k - 1) + s%dq(k))
         beta = 1d0 - alfa
      else
         alfa = 0.5d0
         beta = 0.5d0
      end if
   end subroutine get_TDC_alfa_beta_face_weights


   function get_TDC_Hp_face(s, k, ierr) result(Hp_face)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Hp_face

      ierr = 0
      if (s%use_face_reconstruction) then
         call get_reconstructed_scale_height_ad(s, k, Hp_face, ierr)
      else
         Hp_face = get_scale_height_face(s, k)
      end if
   end function get_TDC_Hp_face


   function get_TDC_mixing_length_face(s, k, ierr) result(Lambda_face)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Lambda_face, Hp_face

      ierr = 0
      if (s% harmonic_dissipation_length_beta > 0d0) then
         call get_reconstructed_hse_scale_height_ad(s, k, Hp_face, ierr)
         if (ierr /= 0) return
      else
         Hp_face = get_TDC_Hp_face(s, k, ierr)
         if (ierr /= 0) return
      end if
      Lambda_face = get_mlt_mixing_length( &
         s, Hp_face, wrap_r_00(s,k), s%mixing_length_alpha)
   end function get_TDC_mixing_length_face


   function get_TDC_mixing_length_cell(s, k, ierr) result(Lambda_cell)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Lambda0, Lambda1, Lambda_cell

      ierr = 0
      Lambda0 = get_TDC_mixing_length_face(s, k, ierr)
      if (ierr /= 0) return
      Lambda1 = 0d0
      if (k < s%nz) then
         Lambda1 = shift_p1(get_TDC_mixing_length_face(s, k+1, ierr))
         if (ierr /= 0) return
      end if
      Lambda_cell = 0.5d0*(Lambda0 + Lambda1)
   end function get_TDC_mixing_length_cell


   function compute_Chi_div_w_cell(s, k, ierr) result(Chi_div_w)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_div_w
      type(auto_diff_real_star_order1) :: &
         rho2, r6_cell, d_v_div_r, Lambda_cell, d_00, r_00, r_p1
      real(dp) :: f, ALFAM

      ierr = 0
      Chi_div_w = 0d0

      if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
         ALFAM = s%TDC_alpha_M
      else
         ALFAM = 0d0
      end if

      if (ALFAM == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) return

      Lambda_cell = get_TDC_mixing_length_cell(s, k, ierr)
      if (ierr /= 0) return
      d_v_div_r = compute_d_v_div_r(s, k, ierr)
      if (ierr /= 0) return

      d_00 = wrap_d_00(s, k)
      f = (16d0/3d0)*pi*ALFAM/s%dm(k)
      rho2 = pow2(d_00)
      r_00 = wrap_r_00(s, k)
      r_p1 = wrap_r_p1(s, k)
      r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
      Chi_div_w = f*rho2*r6_cell*d_v_div_r*Lambda_cell
   end function compute_Chi_div_w_cell


   ! Cell stress used by v_flag momentum and energy.
   function compute_Chi_cell(s, k, ierr) result(Chi_cell)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_cell, Chi_div_w, w_cell

      ierr = 0
      Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
      if (ierr /= 0) return

      if (s% okay_to_set_mlt_vc .and. &
            s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
         if (k < s%nz) then
            w_cell = 0.5d0*(s%mlt_vc_old(k) + s%mlt_vc_old(k+1))/sqrt_2_div_3
         else
            w_cell = 0.5d0*s%mlt_vc_old(k)/sqrt_2_div_3
         end if
      else
         if (k < s%nz) then
            w_cell = 0.5d0*(s%mlt_vc_ad(k) + &
               shift_p1(s%mlt_vc_ad(k+1)))/sqrt_2_div_3
         else
            w_cell = 0.5d0*s%mlt_vc_ad(k)/sqrt_2_div_3
         end if
      end if
      Chi_cell = Chi_div_w*w_cell
   end function compute_Chi_cell

   ! Face stress used by u_flag momentum and energy.
   function compute_Chi_div_w_face(s, k, ierr) result(Chi_div_w)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_div_w
      type(auto_diff_real_star_order1) :: &
         rho2, r6_face, d_v_div_r, Lambda_face, d_00, r_00
      real(dp) :: f, ALFAM, dmbar

      ierr = 0
      Chi_div_w = 0d0

      if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
         ALFAM = s%TDC_alpha_M
      else
         ALFAM = 0d0
      end if

      if (ALFAM /= 0d0 .and. &
            k > s%TDC_num_outermost_cells_forced_nonturbulent + 1 .and. &
            k <= s%nz - s%TDC_num_innermost_cells_forced_nonturbulent) then
         Lambda_face = get_TDC_mixing_length_face(s, k, ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r_face(s, k, ierr)
         if (ierr /= 0) return

         dmbar = 0.5d0*(s%dm(k) + s%dm(k-1))
         d_00 = get_rho_face(s, k)
         f = (16d0/3d0)*pi*ALFAM/dmbar
         rho2 = pow2(d_00)
         r_00 = wrap_r_00(s, k)
         r6_face = pow6(r_00)
         Chi_div_w = f*rho2*r6_face*d_v_div_r*Lambda_face
      end if
   end function compute_Chi_div_w_face


   function compute_d_u_div_r_inner_boundary( &
         s, use_time_centering) result(d_u_div_r)
      type(star_info), pointer :: s
      logical, intent(in) :: use_time_centering
      type(auto_diff_real_star_order1) :: d_u_div_r, u_inner
      integer :: nz

      nz = s%nz
      if (use_time_centering) then
         if (s% using_velocity_time_centering .or. &
               .not. s% use_P_d_1_div_rho_form_of_work) then
            ! Match the exact finite-step kinetic-energy work, which uses
            ! (u + u_start)/2 even when optional velocity time centering is off.
            u_inner = 0.5d0*(wrap_u_00(s, nz) + s%u_start(nz))
         else
            u_inner = wrap_u_00(s, nz)
         end if
      else
         u_inner = wrap_u_00(s, nz)
      end if
      d_u_div_r = u_inner/s%rmid_start(nz) - s%v_center/s%R_center
   end function compute_d_u_div_r_inner_boundary


   function compute_Chi_div_w_inner_boundary(s, ierr) result(Chi_div_w)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         Chi_div_w, Lambda_boundary, d_u_div_r, rho_inner
      real(dp) :: ALFAM, dm_boundary
      integer :: nz

      ierr = 0
      Chi_div_w = 0d0
      if (.not. s%TDC_include_inner_boundary_eddy_viscosity .or. &
            s%R_center <= 0d0 .or. .not. s%u_flag .or. &
            s%TDC_num_innermost_cells_forced_nonturbulent > 0d0) return

      if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
         ALFAM = s%TDC_alpha_M
      else
         ALFAM = 0d0
      end if
      if (ALFAM == 0d0 .or. s%mixing_length_alpha == 0d0) return

      nz = s%nz
      ! Extrapolate the innermost turbulent state to the excised boundary.
      Lambda_boundary = get_TDC_mixing_length_face(s, nz, ierr)
      if (ierr /= 0) return
      d_u_div_r = compute_d_u_div_r_inner_boundary(s, .false.)
      rho_inner = wrap_d_00(s, nz)
      ! The boundary is one half-cell in mass from the innermost cell center.
      dm_boundary = 0.5d0*s%dm(nz)
      Chi_div_w = (16d0/3d0)*pi*ALFAM/dm_boundary*pow2(rho_inner)* &
         pow6(s%R_center)*d_u_div_r*Lambda_boundary
   end function compute_Chi_div_w_inner_boundary


   function compute_tdc_Eq_div_w_inner_boundary(s, ierr) result(Eq_div_w)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         Eq_div_w, Chi_div_w, d_u_div_r
      real(dp) :: dm_boundary

      ierr = 0
      Eq_div_w = 0d0
      if (.not. s%TDC_include_inner_boundary_eddy_viscosity .or. &
            s%R_center <= 0d0 .or. .not. s%u_flag .or. &
            s%TDC_num_innermost_cells_forced_nonturbulent > 0d0) return

      Chi_div_w = compute_Chi_div_w_inner_boundary(s, ierr)
      if (ierr /= 0) return
      d_u_div_r = compute_d_u_div_r_inner_boundary(s, .true.)
      dm_boundary = 0.5d0*s%dm(s%nz)
      Eq_div_w = 4d0*pi*Chi_div_w*d_u_div_r/dm_boundary
   end function compute_tdc_Eq_div_w_inner_boundary


   function compute_tdc_Eq_cell(s, k, ierr) result(Eq_cell)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         Eq_cell, Chi_div_w, d_v_div_r, w_cell

      ierr = 0
      Eq_cell = 0d0
      if (s%mixing_length_alpha /= 0d0) then
         Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r_opt_time_center(s, k, ierr)
         if (ierr /= 0) return
         if (k < s%nz) then
            w_cell = 0.5d0*(s%mlt_vc_ad(k) + &
               shift_p1(s%mlt_vc_ad(k+1)))/sqrt_2_div_3
         else
            w_cell = 0.5d0*s%mlt_vc_ad(k)/sqrt_2_div_3
         end if
         Eq_cell = 4d0*pi*Chi_div_w*w_cell*d_v_div_r/s%dm(k)
      end if
      s%Eq(k) = Eq_cell%val
      s%Eq_ad(k) = Eq_cell
   end function compute_tdc_Eq_cell


   function compute_tdc_Eq_div_w_cell(s, k, ierr) result(Eq_div_w)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Eq_div_w, Chi_div_w, d_v_div_r

      ierr = 0
      Eq_div_w = 0d0
      if (s%mixing_length_alpha == 0d0) return

      Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
      if (ierr /= 0) return
      d_v_div_r = compute_d_v_div_r_opt_time_center(s, k, ierr)
      if (ierr /= 0) return
      Eq_div_w = 4d0*pi*Chi_div_w*d_v_div_r/s%dm(k)
   end function compute_tdc_Eq_div_w_cell


   function compute_tdc_Eq_div_w_face(s, k, ierr) result(Eq_div_w)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Eq_div_w
      type(auto_diff_real_star_order1) :: &
         Eq_div_w_m1, d_v_div_r, Chi_div_w, w_face, w_momentum
      real(dp) :: dmbar

      ierr = 0
      Eq_div_w = 0d0
      Chi_div_w = 0d0
      if (s% v_flag) then
         Eq_div_w = compute_tdc_Eq_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         if (k > 1) then
            Eq_div_w_m1 = shift_m1(compute_tdc_Eq_div_w_cell(s, k-1, ierr))
            if (ierr /= 0) return
            Eq_div_w = (s%dm(k)*Eq_div_w + &
               s%dm(k-1)*Eq_div_w_m1)/(s%dm(k) + s%dm(k-1))
         end if
      else if (s% u_flag) then
         if (s%mixing_length_alpha /= 0d0 .and. &
               k <= s%nz - s%TDC_num_innermost_cells_forced_nonturbulent) then
            Chi_div_w = compute_Chi_div_w_face(s, k, ierr)
            if (ierr /= 0) return
            d_v_div_r = compute_d_v_div_r_opt_time_center_face(s, k, ierr)
            if (ierr /= 0) return
            if (k > 1) then
               dmbar = 0.5d0*(s%dm(k) + s%dm(k-1))
            else
               dmbar = 0.5d0*s%dm(k)
            end if
            Eq_div_w = 4d0*pi*Chi_div_w*d_v_div_r/dmbar
         end if

         w_face = s%mlt_vc_ad(k)/sqrt_2_div_3
         if (s% okay_to_set_mlt_vc .and. &
               s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
            w_momentum = s%mlt_vc_old(k)/sqrt_2_div_3
         else
            w_momentum = w_face
         end if
         s%Chi(k) = Chi_div_w%val*w_momentum%val
         s%Chi_ad(k) = Chi_div_w*w_momentum
         s%Eq(k) = Eq_div_w%val*w_face%val
         s%Eq_ad(k) = Eq_div_w*w_face
      end if
   end function compute_tdc_Eq_div_w_face

   ! for v_flag only. face centered Uq for hydro_momentum
   function compute_tdc_Uq_face(s, k, ierr) result(Uq_face) !(v_flag only)  ! cm s^-2, acceleration
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Uq_face
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_00
      real(dp) :: dm_face
      include 'formats'
      ierr = 0
      Chi_00 = 0d0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
         Uq_face = 0d0
      else
         r_00 = wrap_opt_time_center_r_00(s, k)

         Chi_00 = compute_Chi_cell(s, k, ierr)
         if (ierr /= 0) return

         if (k > 1) then
            Chi_m1 = shift_m1(compute_Chi_cell(s, k-1, ierr))
            if (ierr /= 0) return
         else
            Chi_m1 = 0d0
         end if

         if (s%use_mass_corrections) then
            dm_face = 0.5d0*s%dm(k)*s%mass_correction(k)
            if (k > 1) dm_face = dm_face + &
               0.5d0*s%dm(k-1)*s%mass_correction(k-1)
         else
            dm_face = 0.5d0*s%dm(k)
            if (k > 1) dm_face = dm_face + 0.5d0*s%dm(k-1)
         end if
         ! Match the dual-face mass used by the momentum and kinetic-energy equations.
         Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*dm_face)

         if (k == -56) then
            write (*, 3) 'TDC Uq chi_m1 chi_00 r', k, s%solver_iter, &
               Uq_face%val, Chi_m1%val, Chi_00%val, r_00%val
         end if

      end if
      s%Chi(k) = Chi_00%val
      s%Chi_ad(k) = Chi_00
      ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
      s%Uq(k) = Uq_face%val
   end function compute_tdc_Uq_face

   ! for u_flag only. cell centered Uq as source for Reimann flux.
   function compute_tdc_Uq_dm_cell(s, k, ierr) result(Uq_cell)  ! cm s^-2, acceleration
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_00, Chi_p1, w_00, w_p1, r_cell, Uq_cell
      include 'formats'
      ierr = 0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
         Uq_cell = 0d0
      else
         ! Use the same lagged cell radius as the u/r strain.
         r_cell = s% rmid_start(k)

         if (s% okay_to_set_mlt_vc .and. &
            s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
            w_00 = s% mlt_vc_old(k)/sqrt_2_div_3
         else
            w_00 = s% mlt_vc_ad(k)/sqrt_2_div_3
         end if

         Chi_00 = compute_Chi_div_w_face(s, k, ierr) * w_00
         if (ierr /= 0) return

         if (k < s% nz) then
            if (s% okay_to_set_mlt_vc .and. &
               s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
               w_p1 = s% mlt_vc_old(k+1)/sqrt_2_div_3
            else
               w_p1 = shift_p1(s% mlt_vc_ad(k+1))/sqrt_2_div_3
            end if

            Chi_p1 = shift_p1(compute_Chi_div_w_face(s, k+1, ierr))*w_p1
            if (ierr /= 0) return
         else
            Chi_p1 = compute_Chi_div_w_inner_boundary(s, ierr)
            if (ierr /= 0) return
            if (s% okay_to_set_mlt_vc .and. &
                  s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
               w_p1 = s%mlt_vc_old(k)/sqrt_2_div_3
            else
               w_p1 = s%mlt_vc_ad(k)/sqrt_2_div_3
            end if
            Chi_p1 = Chi_p1*w_p1
         end if

         ! Return force; Riemann momentum divides by dm.
         Uq_cell = 4d0*pi*(Chi_00 - Chi_p1)/r_cell
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2 [g], acceleration*mass = Force

         if (k == -56) then
            write (*, 3) 'TDC Uq chi_m1 chi_00 r', k, s%solver_iter, &
               Uq_cell%val, Chi_p1%val, Chi_00%val, r_cell%val
         end if

      end if
      s%Uq(k) = Uq_cell%val/s%dm(k)
   end function compute_tdc_Uq_dm_cell


! all the forms of d(v/r)/dr, below
   function compute_d_v_div_r(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.
      v_00 = wrap_v_00(s, k)
      v_p1 = wrap_v_p1(s, k)
      r_00 = wrap_r_00(s, k)
      r_p1 = wrap_r_p1(s, k)
      if (r_p1%val == 0d0) r_p1 = 1d0
      d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_p1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_p1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r

   function compute_d_v_div_r_opt_time_center(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
      include 'formats'
      ierr = 0
      v_00 = wrap_v_00(s, k)
      v_p1 = wrap_v_p1(s, k)
      if (s% using_velocity_time_centering .or. &
            .not. s% use_P_d_1_div_rho_form_of_work) then
         ! Match the exact finite-step kinetic-energy work, which uses
         ! (v + v_start)/2 even when optional velocity time centering is off.
         v_00 = 0.5d0*(v_00 + s%v_start(k))
         if (k < s%nz) v_p1 = 0.5d0*(v_p1 + s%v_start(k+1))
      end if
      r_00 = wrap_opt_time_center_r_00(s, k)
      r_p1 = wrap_opt_time_center_r_p1(s, k)
      if (r_p1%val == 0d0) r_p1 = 1d0
      d_v_div_r = v_00/r_00 - v_p1/r_p1  ! units s^-1
   end function compute_d_v_div_r_opt_time_center

   function compute_d_v_div_r_face(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_m1, r_00, r_m1
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.

     if (s% v_flag) then
         v_00 = 0.5d0*(wrap_v_00(s, k) + wrap_v_p1(s, k))
         v_m1 = 0.5d0*(wrap_v_00(s, k) + wrap_v_m1(s, k))
     else if(s% u_flag) then
         v_00 = wrap_u_00(s,k)
         v_m1 = wrap_u_m1(s,k)
      end if

      if (s% v_flag) then
         r_00 = 0.5d0*(wrap_r_00(s, k) + wrap_r_p1(s, k))
         r_m1 = 0.5d0*(wrap_r_00(s, k) + wrap_r_m1(s, k))
      else if(s% u_flag) then
         ! Lag the cell-centered radius to retain tridiagonality.
         r_00 = s% rmid_start(k)
         if (k > 1) then
            r_m1 = s% rmid_start(k-1)
         else
            r_m1 = 1d0
         end if
      end if

      if (r_00%val == 0d0) r_00 = 1d0
      if (r_m1%val == 0d0) r_m1 = 1d0
      d_v_div_r = v_m1/r_m1 - v_00/r_00 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_m1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_m1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r_face

   function compute_d_v_div_r_opt_time_center_face(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_m1, r_00, r_m1
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.

     if (s% v_flag) then
         v_00 = 0.5d0 *(wrap_opt_time_center_v_00(s, k) + wrap_opt_time_center_v_p1(s, k))
         v_m1 = 0.5d0*(wrap_opt_time_center_v_00(s, k) + wrap_opt_time_center_v_m1(s, k))
     else if(s% u_flag) then
         v_00 = wrap_u_00(s,k)
         v_m1 = wrap_u_m1(s,k)
         if (s% using_velocity_time_centering .or. &
               .not. s% use_P_d_1_div_rho_form_of_work) then
            ! Match the exact finite-step kinetic-energy work, which uses
            ! (u + u_start)/2 even when optional velocity time centering is off.
            v_00 = 0.5d0*(v_00 + s%u_start(k))
            if (k > 1) v_m1 = 0.5d0*(v_m1 + s%u_start(k-1))
         end if
      end if

      if (s% v_flag) then
         r_00 = 0.5d0*(wrap_opt_time_center_r_00(s, k) + wrap_opt_time_center_r_p1(s, k))
         r_m1 = 0.5d0*(wrap_opt_time_center_r_00(s, k) + wrap_opt_time_center_r_m1(s, k))
      else if(s% u_flag) then
         ! Lag the cell-centered radius to retain tridiagonality.
         r_00 = s% rmid_start(k)
         if (k > 1) then
            r_m1 = s% rmid_start(k-1)
         else
            r_m1 = 1d0
         end if
      end if

      if (r_00%val == 0d0) r_00 = 1d0
      if (r_m1%val == 0d0) r_m1 = 1d0
      d_v_div_r = v_m1/r_m1 - v_00/r_00 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_m1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_m1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r_opt_time_center_face

end module tdc_hydro
