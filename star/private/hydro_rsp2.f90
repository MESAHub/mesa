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
      public :: compute_Source, compute_D, compute_Dr, compute_Lt_center
      public :: compute_Source_div_w
      public :: compute_Eq_cell
      public :: compute_Uq_face, compute_Uq_dm_cell
      public :: set_RSP2_vars
      public :: rsp2_flux_residual, set_etrb_start_vars
      public :: RSP2_adjust_vars_before_call_solver
      public :: get_RSP2_alfa_beta_face_weights
      public :: do1_rsp2_moment_eqns, rsp2_moment_rhs
      public :: init_rsp2_moments, check_rsp2_moments
      public :: rsp2_dormant_moments, rsp2_local_w_equation
      public :: remap_rsp2, interpolate_rsp2_face

      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, &  ! Ptrb
         x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3, &  ! PII_face and Lc
         x_ALFAC = (1.d0/2.d0)*sqrt_2_div_3, &  ! Lc
         x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3, &  ! DAMP
         x_GAMMAR = 2.d0*sqrt(3.d0)  ! DAMPR

      ! Phi is the full entropy variance. Pi decays at half the sum of variance losses.
      ! The unit Phi multiplier retains the Kuhfuss local variance decay coefficient.
      real(dp), parameter :: x_ALFAPHI = 4d0*sqrt_2_div_3

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


      function rsp2_flux_residual(s, k, perturb_convective_flux) result(resid)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in), optional :: perturb_convective_flux
         type(auto_diff_real_star_order1) :: resid, L_expected, L_actual, Lc, Lt
         type(accurate_auto_diff_real_star_order1) :: L_sum
         real(dp) :: scale

         if (k == 1) then
            ! The surface luminosity has its own physical boundary condition.
            resid = wrap_Y_00(s, k)
            return
         end if
         Lc = s% Lc_ad(k)
         Lt = s% Lt_ad(k)
         if (present(perturb_convective_flux)) then
            if (.not. perturb_convective_flux) then
               Lc = Lc%val
               Lt = Lt%val
            end if
         end if
         L_sum = s% Lr_ad(k)
         L_sum = L_sum + Lc
         L_sum = L_sum + Lt
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


      logical function rsp2_dormant_moments(s, k) result(dormant)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer :: j
         dormant = .true.
         if (rsp2_zero_w(s,k)) return
         dormant = get_etrb(s,k) == 0d0 .and. s% Pi(k) == 0d0 .and. s% Phi(k) == 0d0
         if (.not. dormant .or. s% RSP2_alfat == 0d0) return
         ! An empty face can receive moments from either active neighbor.
         do j=max(1,k-1),min(s% nz,k+1)
            if (rsp2_zero_w(s,j)) cycle
            if (s% w(j) /= 0d0 .or. s% Pi(j) /= 0d0 .or. s% Phi(j) /= 0d0) dormant = .false.
         end do
      end function rsp2_dormant_moments


      logical function rsp2_moments_valid(w, Pi_value, Phi_value) result(valid)
         real(dp), intent(in) :: w, Pi_value, Phi_value
         real(dp) :: limit
         valid = .false.
         if (is_bad(w) .or. is_bad(Pi_value) .or. is_bad(Phi_value)) return
         if (w < 0d0 .or. Phi_value < 0d0) return
         limit = sqrt_2_div_3*w*sqrt(Phi_value)
         valid = abs(Pi_value) <= (1d0+64d0*epsilon(1d0))*limit
      end function rsp2_moments_valid


      subroutine check_rsp2_moments(s,ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         if (.not. s% RSP2_3equation_flag) return
         do k=1,s% nz
            if (rsp2_moments_valid(s% w(k),s% Pi(k),s% Phi(k))) cycle
            ierr = -1
            write(s% retry_message,'(a,i0)') 'invalid RSP3 moments at face ', k
            return
         end do
      end subroutine check_rsp2_moments


      subroutine interpolate_rsp2_face( &
            n_old, x_old, n_new, x_new, face_old, face_new, work, ierr)
         use interp_1d_lib, only: interpolate_vector_pm
         integer, intent(in) :: n_old, n_new
         real(dp), intent(in) :: x_old(:), x_new(:), face_old(:)
         real(dp), intent(inout) :: face_new(:)
         real(dp), pointer :: work(:)
         integer, intent(out) :: ierr
         integer :: j, k

         ierr = 0
         call interpolate_vector_pm( &
            n_old, x_old, n_new, x_new, face_old, face_new, work, 'RSP2 face remesh', ierr)
         if (ierr /= 0) return
         j = 1
         do k=1,n_new
            do while (j < n_old-1)
               if (x_new(k) < x_old(j+1)) exit
               j = j+1
            end do
            ! Preserve old faces exactly.
            if (x_new(k) == x_old(j)) then
               face_new(k) = face_old(j)
            else if (x_new(k) == x_old(j+1)) then
               face_new(k) = face_old(j+1)
            end if
         end do
      end subroutine interpolate_rsp2_face


      subroutine remap_rsp2(s,nz_old,dq_old,xh_old,nz,dq,xh,ierr)
         use const_def, only: qp
         type(star_info), pointer :: s
         integer, intent(in) :: nz_old, nz
         real(dp), intent(in) :: dq_old(:), xh_old(:,:), dq(:)
         real(dp), intent(inout) :: xh(:,:)
         integer, intent(out) :: ierr
         real(qp) :: old_edge(nz_old+2), new_edge(nz+2), &
            old_value(3,nz_old+1), new_integral(3,nz+1), overlap, energy0, energy1
         real(dp) :: dm_face
         integer :: k, j, recipient, first, last

         ierr = 0
         if (any(is_bad(dq_old(1:nz_old))) .or. any(is_bad(dq(1:nz))) .or. &
               any(dq_old(1:nz_old) <= 0d0) .or. any(dq(1:nz) <= 0d0)) then
            ierr = -1
            s% retry_message = 'invalid mass in RSP2 remap'
            return
         end if
         old_value = 0.0_qp
         old_value(1,1:nz_old) = real(xh_old(s% i_w,1:nz_old),qp)**2
         if (s% RSP2_3equation_flag) then
            old_value(2,1:nz_old) = xh_old(s% i_Pi,1:nz_old)
            old_value(3,1:nz_old) = xh_old(s% i_Phi,1:nz_old)
         end if
         if (any(is_bad(real(old_value,dp))) .or. any(old_value(3,:) < 0.0_qp) .or. &
               any(xh_old(s% i_w,1:nz_old) < 0d0)) then
            ierr = -1
            s% retry_message = 'invalid moment in RSP2 remap'
            return
         end if
         if (s% RSP2_3equation_flag) then
            do k=1,nz_old
               if (rsp2_moments_valid(xh_old(s% i_w,k),xh_old(s% i_Pi,k),xh_old(s% i_Phi,k))) cycle
               ierr = -1
               s% retry_message = 'invalid covariance in RSP3 remap'
               return
            end do
         end if
         ! Quadruple coordinates retain overlap widths in very small zones.
         call face_edges(nz_old,dq_old,old_edge)
         call face_edges(nz,dq,new_edge)
         if (abs(new_edge(nz+2)-old_edge(nz_old+2)) > &
               128.0_qp*epsilon(1d0)*old_edge(nz_old+2)) then
            ierr = -1
            s% retry_message = 'different total mass in RSP2 remap'
            return
         end if
         new_edge = new_edge*(old_edge(nz_old+2)/new_edge(nz+2))
         new_integral = 0.0_qp
         j = 1
         do k=1,nz+1
            do while (j <= nz_old+1)
               overlap = min(new_edge(k+1),old_edge(j+1)) - max(new_edge(k),old_edge(j))
               if (overlap > 0.0_qp) new_integral(:,k) = new_integral(:,k) + overlap*old_value(:,j)
               if (old_edge(j+1) >= new_edge(k+1)) exit
               j = j+1
            end do
         end do
         first = max(2,s% RSP2_num_outermost_cells_forced_nonturbulent+2)
         last = nz-int(nz/s% RSP2_nz_div_IBOTOM)
         if (s% mixing_length_alpha == 0d0) last = 0
         if (first > last) then
            if (any(new_integral /= 0.0_qp)) then
               ierr = -1
               s% retry_message = 'no active face for RSP2 remap'
               return
            end if
         else
            do k=1,nz+1
               if (k >= first .and. k <= last) cycle
               ! Boundary volumes cannot store moments. Transfer their integrals together.
               recipient = min(last,max(first,k))
               new_integral(:,recipient) = new_integral(:,recipient) + new_integral(:,k)
               new_integral(:,k) = 0.0_qp
            end do
         end if
         energy0 = sum((old_edge(2:nz_old+2)-old_edge(1:nz_old+1))*old_value(1,:))
         energy1 = 0.0_qp
         do k=1,nz
            dm_face = 0.5d0*dq(k)
            if (k > 1) dm_face = dm_face + 0.5d0*dq(k-1)
            xh(s% i_w,k) = real(sqrt(new_integral(1,k)/dm_face),dp)
            energy1 = energy1 + real(dm_face,qp)*real(xh(s% i_w,k),qp)**2
            if (s% RSP2_3equation_flag) then
               xh(s% i_Pi,k) = real(new_integral(2,k)/dm_face,dp)
               xh(s% i_Phi,k) = real(new_integral(3,k)/dm_face,dp)
               if (.not. rsp2_moments_valid(xh(s% i_w,k),xh(s% i_Pi,k),xh(s% i_Phi,k))) then
                  ierr = -1
                  s% retry_message = 'invalid remapped RSP3 covariance'
                  return
               end if
            end if
         end do
         s% mesh_adjust_Eturb_conservation = real(abs(energy1-energy0)/max(energy0,tiny(energy0)),dp)

         contains

         subroutine face_edges(n,dq,edge)
            integer, intent(in) :: n
            real(dp), intent(in) :: dq(:)
            real(qp), intent(out) :: edge(:)
            real(qp) :: mass
            integer :: i
            edge(1) = 0.0_qp
            mass = 0.0_qp
            do i=1,n
               edge(i+1) = mass + 0.5_qp*real(dq(i),qp)
               mass = mass + real(dq(i),qp)
            end do
            edge(n+2) = mass
         end subroutine face_edges

      end subroutine remap_rsp2


      function rsp2_buoyancy_face(s, k, ierr) result(buoyancy)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: buoyancy, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face

         ierr = 0
         buoyancy = 0d0
         if (rsp2_zero_w(s,k)) return
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
         type(auto_diff_real_star_order1) :: source, buoyancy, Pi_face

         ierr = 0
         source = 0d0
         if (rsp2_zero_w(s,k)) return
         buoyancy = rsp2_buoyancy_face(s,k,ierr)
         if (ierr /= 0) return
         Pi_face = wrap_Pi_00(s,k)
         if (divide_by_w) then
            if (s% w(k) <= 0d0) then
               if (s% Pi(k) /= 0d0) ierr = -1
               return
            end if
            source = buoyancy*div_by_w(Pi_face,wrap_w_00(s,k))
         else
            source = buoyancy*Pi_face
         end if
      end function rsp2_moment_source


      subroutine rsp2_moment_coefficients( &
            s, k, buoyancy, entropy_gradient, Lambda_face, inverse_rad_time, ierr, use_time_centering)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: &
            buoyancy, entropy_gradient, Lambda_face, inverse_rad_time
         integer, intent(out) :: ierr
         logical, intent(in), optional :: use_time_centering
         type(auto_diff_real_star_order1) :: &
            Cp_face, T_face, rho_face, kap_face, P_face, ChiRho_face, ChiT_face, grad_ad, gradL

         ierr = 0
         buoyancy = 0d0
         entropy_gradient = 0d0
         Lambda_face = 1d0
         inverse_rad_time = 0d0
         if (rsp2_zero_w(s,k)) return
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         Lambda_face = get_TDC_mixing_length_face(s,k,ierr)
         if (ierr /= 0) return
         if (Lambda_face%val <= 0d0) then
            ierr = -1
            return
         end if
         call get_rsp2_thermal_gradient(s, k, grad_ad, gradL, entropy_gradient, ierr, use_time_centering)
         if (ierr /= 0) return
         buoyancy = rsp2_buoyancy_face(s,k,ierr)
         if (ierr /= 0) return
         inverse_rad_time = 4d0*boltz_sigma*pow2(s% RSP2_alfar*x_GAMMAR)*pow3(T_face)/ &
            (Cp_face*kap_face*pow2(rho_face)*pow2(Lambda_face))
      end subroutine rsp2_moment_coefficients


      subroutine rsp2_moment_rhs(s, k, Pi_rhs, Phi_rhs, ierr, use_time_centering)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Pi_rhs, Phi_rhs
         integer, intent(out) :: ierr
         logical, intent(in), optional :: use_time_centering
         type(auto_diff_real_star_order1) :: Pi_face, Phi_face, w_face, &
            buoyancy, entropy_gradient, Lambda_face, inverse_rad_time, Lt, &
            Pi_outer, Phi_outer, Pi_inner, Phi_inner
         real(dp) :: C_Pi, C_Phi, dm_face

         ierr = 0
         Pi_rhs = 0d0
         Phi_rhs = 0d0
         if (rsp2_zero_w(s,k)) return
         call rsp2_moment_coefficients( &
            s,k,buoyancy,entropy_gradient,Lambda_face,inverse_rad_time,ierr,use_time_centering)
         if (ierr /= 0) return
         C_Phi = s% RSP2_alfa_phi*x_ALFAPHI
         C_Pi = 0.5d0*(s% RSP2_alfad*x_CEDE + C_Phi)*s% RSP2_alfa_pi
         w_face = wrap_w_00(s,k)
         Pi_face = wrap_Pi_00(s,k)
         Phi_face = wrap_Phi_00(s,k)
         ! Pair buoyancy production and decay with the two isotropic variances.
         Pi_rhs = (2d0/3d0)*pow2(w_face)*entropy_gradient + (buoyancy/3d0)*Phi_face - &
            (C_Pi*w_face/Lambda_face + inverse_rad_time)*Pi_face
         Phi_rhs = 2d0*entropy_gradient*Pi_face - &
            (C_Phi*w_face/Lambda_face + 2d0*inverse_rad_time)*Phi_face
         if (s% RSP2_alfat == 0d0) return
         Lt = compute_Lt_center(s,k-1,ierr,Pi_outer,Phi_outer)
         if (ierr /= 0) return
         Lt = compute_Lt_center(s,k,ierr,Pi_inner,Phi_inner)
         if (ierr /= 0) return
         dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
         Pi_rhs = Pi_rhs - (shift_m1(Pi_outer) - Pi_inner)/dm_face
         Phi_rhs = Phi_rhs - (shift_m1(Phi_outer) - Phi_inner)/dm_face
      end subroutine rsp2_moment_rhs


      subroutine do1_rsp2_moment_eqns(s, k, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         integer :: i_flux, i_variance
         type(auto_diff_real_star_order1) :: Pi_rhs, Phi_rhs, Pi_resid, Phi_resid, &
            work, work_new, work_start
         real(dp) :: Pi_scale, Phi_scale

         ierr = 0
         if (.not. s% RSP2_3equation_flag) return
         i_flux = s% i_Pi
         i_variance = s% i_Phi
         Pi_scale = s% Pi_scale(k)
         Phi_scale = s% Phi_scale(k)
         Pi_resid = wrap_Pi_00(s,k)
         Phi_resid = wrap_Phi_00(s,k)
         if (.not. rsp2_zero_w(s,k) .and. &
               .not. (rsp2_dormant_moments(s,k) .and. &
                  s% xh_start(i_flux,k) == 0d0 .and. s% xh_start(i_variance,k) == 0d0)) then
            call rsp2_moment_rhs(s,k,Pi_rhs,Phi_rhs,ierr)
            if (ierr /= 0) return
            ! Use the same implicit source time level as COUPL in the w equation.
            work = calc_Ptrb_work_face(s,k,work_new=work_new,work_start=work_start)
            if (work_start%val >= 1d0) then
               ierr = -1
               s% retry_message = 'RSP3 pressure work exceeds old turbulent energy'
               return
            end if
            Pi_resid = Pi_resid - s% xh_start(i_flux,k) - s% dt*Pi_rhs + &
               0.5d0*work_new*wrap_Pi_00(s,k) + &
               work_start/(1d0 + sqrt(1d0-work_start))*s% xh_start(i_flux,k)
            Phi_resid = Phi_resid - s% xh_start(i_variance,k) - s% dt*Phi_rhs
            Pi_scale = Pi_scale + s% dt*s% Pi_rate_scale(k)
            Phi_scale = Phi_scale + s% dt*s% Phi_rate_scale(k)
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
         real(dp) :: etrb_face
         type(auto_diff_real_star_order1) :: &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, &
            gradL, entropy_gradient

         ierr = 0
         do k=1,s% nz
            s% Pi(k) = 0d0
            s% Phi(k) = 0d0
            if (.not. rsp2_zero_w(s,k)) then
               call get_rsp2_thermal_gradient(s,k,grad_ad,gradL,entropy_gradient,ierr)
               if (ierr /= 0) return
               ! Do not infer finite entropy fluctuations from a quiet stable layer.
               if (entropy_gradient%val > 0d0 .and. Lc_old(k) > 0d0) then
                  etrb_face = pow2(s% w(k))
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
         local = .not. s% RSP2_3equation_flag .and. s% RSP2_source_seed == 0d0 .and. &
            s% RSP2_alfat == 0d0 .and. get_etrb_start(s,k) == 0d0
      end function rsp2_local_w_equation


      function div_by_w(x, w) result(quotient)
         type(auto_diff_real_star_order1), intent(in) :: x, w
         type(auto_diff_real_star_order1) :: quotient

         ! Avoid 1/w**2, which can overflow when the quotient derivatives are finite.
         quotient%val = x%val/w%val
         quotient%d1Array = (x%d1Array - quotient%val*w%d1Array)/w%val
      end function div_by_w


      subroutine do1_turbulent_energy_eqn(s,k,nvar,ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: resid, w, source, damping, rad_damping, Eq_face, dLt_dm, work, old_energy
         type(accurate_auto_diff_real_star_order1) :: esum
         real(dp) :: scal, scal_outer, dm_face, w_start_over_w
         logical :: positive_branch, divided

         ierr = 0
         w = wrap_w_00(s,k)
         if (.not. s% RSP2_flag) then
            resid = w - s% w_start(k)
         else if (rsp2_zero_w(s,k)) then
            resid = w/s% csound(k)
         else
            call set_energy_eqn_scal(s,k,scal,ierr)
            if (ierr /= 0) return
            call set_energy_eqn_scal(s,k-1,scal_outer,ierr)
            if (ierr /= 0) return
            dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
            scal = 2d0*dm_face/(s% dm(k)/scal + s% dm(k-1)/scal_outer)
            positive_branch = rsp2_local_w_equation(s,k)
            divided = positive_branch .or. (s% RSP2_3equation_flag .and. &
               w%val > 0d0 .and. w%val < s% csound(k))
            if (divided) then
               source = compute_Source_div_w(s,k,ierr)
               if (ierr /= 0) return
               damping = compute_D_div_w(s,k,ierr)
               if (ierr /= 0) return
               rad_damping = compute_Dr_div_w(s,k,ierr)
               if (ierr /= 0) return
               Eq_face = compute_Eq_div_w_face(s,k,ierr)
               if (ierr /= 0) return
               work = calc_Ptrb_work_face(s,k,.true.)
               esum = w + work - s% dt*(source - damping - rad_damping + Eq_face)
               if (.not. positive_branch) then
                  ! Form the old energy quotient without squaring a small velocity.
                  w_start_over_w = s% w_start(k)/w%val
                  old_energy%val = s% w_start(k)*w_start_over_w
                  old_energy%d1Array = -pow2(w_start_over_w)*w%d1Array
                  esum = esum - old_energy
                  dLt_dm = rsp2_dLt_dm_face(s,k,ierr)
                  if (ierr /= 0) return
                  esum = esum + s% dt*div_by_w(dLt_dm,w)
               end if
               resid = esum
               if (positive_branch .and. w%val <= resid%val) resid = w
               resid = resid*s% csound(k)*scal/s% dt
            else
               source = compute_Source(s,k,ierr)
               if (ierr /= 0) return
               damping = compute_D(s,k,ierr)
               if (ierr /= 0) return
               rad_damping = compute_Dr(s,k,ierr)
               if (ierr /= 0) return
               Eq_face = compute_Eq_face(s,k,ierr)
               if (ierr /= 0) return
               dLt_dm = rsp2_dLt_dm_face(s,k,ierr)
               if (ierr /= 0) return
               work = calc_Ptrb_work_face(s,k)
               esum = (w-s% w_start(k))*(w+s% w_start(k))
               esum = esum + work + s% dt*(dLt_dm - source + damping + rad_damping - Eq_face)
               resid = esum
               if (w%val == 0d0 .and. resid%val == 0d0 .and. get_etrb_start(s,k) == 0d0) &
                  resid = w*s% csound(k)
               resid = resid*scal/s% dt
            end if
         end if
         s% equ(s% i_detrb_dt,k) = resid%val
         call save_eqn_residual_info(s,k,nvar,s% i_detrb_dt,resid,'detrb_dt',ierr)
      end subroutine do1_turbulent_energy_eqn


      function compute_PII_face(s, k, ierr) result(PII_face)  ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face, w_face
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            PII_face = 0d0
            return
         end if
         if (rsp2_zero_w(s,k)) then
            PII_face = 0d0
            s% PII(k) = 0d0
            s% PII_ad(k) = 0d0
            return
         end if
         Y_face = s% Y_face_ad(k)
         if (ierr /= 0) return
         if (s% RSP2_3equation_flag) then
            ! PII remains a diagnostic; the moment equations use <v_r*s> directly.
            w_face = wrap_w_00(s,k)
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


      function compute_Source_div_w(s, k, ierr) result(Source_div_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source_div_w, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, Hp_face

         ierr = 0
         Source_div_w = 0d0
         if (rsp2_zero_w(s,k)) return
         if (s% RSP2_3equation_flag) then
            Source_div_w = rsp2_moment_source(s,k,.true.,ierr)
            return
         end if
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         Hp_face = get_TDC_Hp_face(s,k,ierr)
         if (ierr /= 0) return
         Source_div_w = P_face*ChiT_face*s% PII_ad(k)/(rho_face*ChiRho_face*Cp_face*Hp_face)
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
         type(auto_diff_real_star_order1) :: Lambda_face, w_00
         ierr = 0
         if (rsp2_zero_w(s,k)) then
            D_div_w = 0d0
         else
            Lambda_face = get_TDC_mixing_length_face(s, k, ierr)
            if (ierr /= 0) return
            w_00 = wrap_w_00(s,k)
            D_div_w = (s% RSP2_alfad*x_CEDE)* &
               pow2(w_00)/Lambda_face
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


      function compute_Dr_div_w(s, k, ierr) result(Dr_div_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Dr_div_w, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, Lambda_face

         ierr = 0
         Dr_div_w = 0d0
         if (s% RSP2_3equation_flag .or. s% RSP2_alfar == 0d0 .or. rsp2_zero_w(s,k)) return
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         Lambda_face = get_TDC_mixing_length_face(s,k,ierr)
         if (ierr /= 0) return
         Dr_div_w = 4d0*boltz_sigma*pow2(s% RSP2_alfar*x_GAMMAR)*pow3(T_face)*wrap_w_00(s,k)/ &
            (Cp_face*kap_face*pow2(rho_face)*pow2(Lambda_face))
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
            r6_cell*Lambda_cell*compute_d_v_div_r(s, k, s% RSP2_3equation_flag)/s% dm(k)
      end function compute_Chi_div_w_cell


      function compute_Chi_cell(s, k, ierr) result(Chi_cell)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_cell
         Chi_cell = compute_Chi_div_w_cell(s,k,ierr)
         if (ierr /= 0) return
         Chi_cell = Chi_cell*0.5d0*(wrap_w_00(s,k) + wrap_w_p1(s,k))
      end function compute_Chi_cell


      function compute_Eq_div_w_cell(s,k,ierr) result(Eq_div_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_div_w, Chi_div_w
         Chi_div_w = compute_Chi_div_w_cell(s,k,ierr)
         if (ierr /= 0) return
         Eq_div_w = 4d0*pi*Chi_div_w*compute_d_v_div_r(s,k,.true.)/s% dm(k)
      end function compute_Eq_div_w_cell


      function compute_Eq_cell(s,k,ierr,Eq_inner) result(Eq_cell)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1), intent(out), optional :: Eq_inner
         type(auto_diff_real_star_order1) :: Eq_cell, Eq_p1

         Eq_p1 = 0d0
         if (s% v_flag) then
            ! Use the cell stress heating, as in TDC.
            Eq_cell = compute_Eq_div_w_cell(s,k,ierr)
            if (ierr /= 0) return
            Eq_cell = Eq_cell*0.5d0*(wrap_w_00(s,k) + wrap_w_p1(s,k))
         else
            Eq_cell = 0.5d0*compute_Eq_face(s,k,ierr)
            if (ierr /= 0) return
            if (k < s% nz) then
               Eq_p1 = 0.5d0*compute_Eq_face(s,k+1,ierr)
               if (ierr /= 0) return
               Eq_cell = Eq_cell + shift_p1(Eq_p1)
            end if
         end if
         ! Preserve the unshifted part for the energy row's k+2 derivatives.
         if (present(Eq_inner)) Eq_inner = Eq_p1
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


      function compute_Chi_div_w_face(s,k,ierr) result(Chi_div_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_div_w, rho_face, r_face, Lambda_face
         real(dp) :: dm_face

         ierr = 0
         Chi_div_w = 0d0
         if (s% RSP2_alfam == 0d0 .or. s% mixing_length_alpha == 0d0 .or. k <= 1) return
         if (k > s% nz) then
            if (.not. s% TDC_include_inner_boundary_eddy_viscosity .or. &
                  s% R_center <= 0d0 .or. int(s% nz/s% RSP2_nz_div_IBOTOM) > 0) return
            rho_face = wrap_d_00(s,s% nz)
            r_face = s% R_center
            Lambda_face = get_TDC_mixing_length_face(s,s% nz,ierr)
            dm_face = 0.5d0*s% dm(s% nz)
         else
            if (rsp2_zero_w(s,k)) return
            rho_face = get_rho_face(s,k)
            r_face = wrap_r_00(s,k)
            Lambda_face = get_TDC_mixing_length_face(s,k,ierr)
            dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
         end if
         if (ierr /= 0) return
         Chi_div_w = (16d0/3d0)*pi*s% RSP2_alfam*pow2(rho_face)*pow6(r_face)* &
            Lambda_face*compute_d_u_div_r_face(s,k,s% RSP2_3equation_flag)/dm_face
      end function compute_Chi_div_w_face

      function compute_Chi_face(s,k,ierr) result(Chi_face)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_face
         Chi_face = compute_Chi_div_w_face(s,k,ierr)
         if (ierr /= 0) return
         ! The envelope wall uses the adjacent interior viscosity.
         Chi_face = Chi_face*wrap_w_00(s,min(k,s% nz))
      end function compute_Chi_face


      function compute_Eq_div_w_face(s,k,ierr) result(Eq_div_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_div_w, Eq_outer, Chi_div_w
         real(dp) :: dm_face

         ierr = 0
         Eq_div_w = 0d0
         if (rsp2_zero_w(s,k) .or. s% RSP2_alfam == 0d0) return
         dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
         if (s% v_flag) then
            Eq_div_w = compute_Eq_div_w_cell(s,k,ierr)
            if (ierr /= 0) return
            Eq_outer = compute_Eq_div_w_cell(s,k-1,ierr)
            if (ierr /= 0) return
            Eq_div_w = (s% dm(k)*Eq_div_w + s% dm(k-1)*shift_m1(Eq_outer))/(2d0*dm_face)
         else if (s% u_flag) then
            Chi_div_w = compute_Chi_div_w_face(s,k,ierr)
            if (ierr /= 0) return
            Eq_div_w = 4d0*pi*Chi_div_w*compute_d_u_div_r_face(s,k,.true.)/dm_face
            if (k == s% nz) then
               Chi_div_w = compute_Chi_div_w_face(s,k+1,ierr)
               if (ierr /= 0) return
               ! Put wall dissipation into the last evolved face volume.
               Eq_div_w = Eq_div_w + 4d0*pi*Chi_div_w*compute_d_u_div_r_face(s,k+1,.true.)/dm_face
            end if
         end if
      end function compute_Eq_div_w_face

      function compute_Eq_face(s,k,ierr) result(Eq_face)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_face
         Eq_face = compute_Eq_div_w_face(s,k,ierr)
         if (ierr /= 0) return
         if (k >= 1 .and. k <= s% nz) Eq_face = Eq_face*wrap_w_00(s,k)
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
         if (rsp2_zero_w(s,k)) then
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
            s% Lc(k) = 0d0
            s% Lt(k) = 0d0
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
            T_m1, T_00, d_m1, d_00, w_00, T_rho_face, PII_face, w_face, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face
         real(dp) :: ALFAC, ALFAS, alfa, beta
         include 'formats'
         ierr = 0
         if (rsp2_zero_w(s,k)) then
            Lc = 0d0
            Lc_div_w_face = 1
            return
         end if
         r_00 = wrap_r_00(s, k)
         area = 4d0*pi*pow2(r_00)
         if (s% RSP2_3equation_flag) then
            Lc = 0d0
            Lc_div_w_face = 0d0
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
         w_00 = wrap_w_00(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         PII_face = s% PII_ad(k)  ! compute_PII_face(s, k, ierr)
         w_face = w_00
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


      function compute_Lt(s, k, ierr) result(Lt)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lt, Lt_outer, Lt_inner

         ierr = 0
         Lt = 0d0
         if (k >= 1 .and. k <= s% nz) s% Lt(k) = 0d0
         if (rsp2_zero_w(s,k)) return
         Lt_outer = compute_Lt_center(s,k-1,ierr)
         if (ierr /= 0) return
         Lt_outer = shift_m1(Lt_outer)
         Lt_inner = compute_Lt_center(s,k,ierr)
         if (ierr /= 0) return
         Lt = (s% dm(k)*Lt_outer + s% dm(k-1)*Lt_inner)/(s% dm(k-1) + s% dm(k))
         s% Lt(k) = Lt%val
      end function compute_Lt

      function compute_Lt_center(s,k,ierr,Pi_flux,Phi_flux,transport_coeff) result(Lt)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1), intent(out), optional :: Pi_flux, Phi_flux, transport_coeff
         type(auto_diff_real_star_order1) :: Lt, Lambda_cell, r_cell, w_outer, w_inner, coeff, w_cell

         ierr = 0
         Lt = 0d0
         if (present(Pi_flux)) Pi_flux = 0d0
         if (present(Phi_flux)) Phi_flux = 0d0
         if (present(transport_coeff)) transport_coeff = 0d0
         if (s% RSP2_alfat == 0d0 .or. rsp2_zero_w(s,k) .or. rsp2_zero_w(s,k+1)) return
         Lambda_cell = get_TDC_mixing_length_cell(s,k,ierr)
         if (ierr /= 0) return
         r_cell = pow(0.5d0*(pow3(wrap_r_00(s,k)) + pow3(wrap_r_p1(s,k))),1d0/3d0)
         w_outer = wrap_w_00(s,k)
         w_inner = wrap_w_p1(s,k)
         coeff = s% RSP2_alfat*pow2(4d0*pi*pow2(r_cell))*pow2(wrap_d_00(s,k))*Lambda_cell/s% dm(k)
         w_cell = 0.5d0*(w_outer + w_inner)
         Lt = -coeff*w_cell*(pow2(w_outer) - pow2(w_inner))
         if (present(transport_coeff)) transport_coeff = coeff
         if (.not. s% RSP2_3equation_flag) return
         if (present(Pi_flux)) Pi_flux = &
            -coeff*w_cell*(wrap_Pi_00(s,k) - shift_p1(wrap_Pi_00(s,k+1)))
         if (present(Phi_flux)) Phi_flux = &
            -coeff*w_cell*(wrap_Phi_00(s,k) - shift_p1(wrap_Phi_00(s,k+1)))
      end function compute_Lt_center

      function rsp2_dLt_dm_face(s,k,ierr) result(dLt_dm)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: dLt_dm, Lt_outer, Lt_inner
         real(dp) :: theta

         ierr = 0
         dLt_dm = 0d0
         if (rsp2_zero_w(s,k)) return
         Lt_outer = compute_Lt_center(s,k-1,ierr)
         if (ierr /= 0) return
         Lt_outer = shift_m1(Lt_outer)
         Lt_inner = compute_Lt_center(s,k,ierr)
         if (ierr /= 0) return
         theta = 1d0
         if (.not. s% RSP2_3equation_flag .and. &
               s% using_velocity_time_centering .and. s% include_L_in_velocity_time_centering) &
            theta = s% L_theta_for_velocity_time_centering
         dLt_dm = (theta*(Lt_outer - Lt_inner) + &
            (1d0-theta)*(s% Lt_center_start(k-1) - s% Lt_center_start(k)))/ &
            (0.5d0*(s% dm(k-1) + s% dm(k)))
      end function rsp2_dLt_dm_face


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: Lt, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, &
            buoyancy, entropy_gradient, Lambda_face, inverse_rad_time, transport_coeff
         real(dp) :: L_scale, max_L, flux_reference, velocity_reference, area_rho, &
            C_Pi, C_Phi, coeff, dm_face, energy_reference
         real(dp) :: w_reference(s% nz)
         include 'formats'
         ierr = 0
         if (s% RSP2_3equation_flag) then
            call check_rsp2_moments(s,ierr)
            if (ierr /= 0) return
         end if
         max_L = maxval(abs(s% L(1:s% nz)))
         w_reference = 0d0
         s% Pi_rate_scale(1:s% nz) = 0d0
         s% Phi_rate_scale(1:s% nz) = 0d0
         do k=1,s%nz
            s% Y_face_start(k) = s% Y_face(k)
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_ad(k) = Lt
            s% Lt_start(k) = Lt%val
            Lt = compute_Lt_center(s,k,ierr)
            if (ierr /= 0) return
            s% Lt_center_start(k) = Lt%val
            s% w_start(k) = s% w(k)
            if (s% RSP2_3equation_flag) then
               call get_rsp2_face_eos( &
                  s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
               if (ierr /= 0) return
               L_scale = max(1d0, abs(s% L(k)), 1d-3*max_L)
               area_rho = 4d0*pi*pow2(s% r(k))*rho_face%val
               flux_reference = L_scale/(area_rho*T_face%val)
               velocity_reference = pow(L_scale/area_rho,1d0/3d0)
               w_reference(k) = max(s% w(k),velocity_reference)
               s% Pi_scale(k) = max(abs(s% Pi(k)), flux_reference)
               s% Phi_scale(k) = max(s% Phi(k), &
                  1.5d0*pow2(flux_reference/velocity_reference))
            end if
         end do
         if (.not. s% RSP2_3equation_flag) return
         C_Phi = s% RSP2_alfa_phi*x_ALFAPHI
         C_Pi = 0.5d0*(s% RSP2_alfad*x_CEDE + C_Phi)*s% RSP2_alfa_pi
         do k=2,s% nz
            if (rsp2_zero_w(s,k)) cycle
            call rsp2_moment_coefficients( &
               s,k,buoyancy,entropy_gradient,Lambda_face,inverse_rad_time,ierr)
            if (ierr /= 0) return
            energy_reference = pow2(w_reference(k))
            s% Pi_rate_scale(k) = (2d0/3d0)*abs(entropy_gradient%val)*energy_reference + &
               abs(buoyancy%val)*s% Phi_scale(k)/3d0 + &
               (C_Pi*w_reference(k)/Lambda_face%val + inverse_rad_time%val)*s% Pi_scale(k)
            s% Phi_rate_scale(k) = 2d0*abs(entropy_gradient%val)*s% Pi_scale(k) + &
               (C_Phi*w_reference(k)/Lambda_face%val + 2d0*inverse_rad_time%val)*s% Phi_scale(k)
         end do
         do k=2,s% nz-1
            Lt = compute_Lt_center(s,k,ierr,transport_coeff=transport_coeff)
            if (ierr /= 0) return
            coeff = transport_coeff%val*0.5d0*(w_reference(k)+w_reference(k+1))
            if (coeff == 0d0) cycle
            dm_face = 0.5d0*(s% dm(k-1)+s% dm(k))
            s% Pi_rate_scale(k) = s% Pi_rate_scale(k) + coeff*(s% Pi_scale(k)+s% Pi_scale(k+1))/dm_face
            s% Phi_rate_scale(k) = s% Phi_rate_scale(k) + coeff*(s% Phi_scale(k)+s% Phi_scale(k+1))/dm_face
            dm_face = 0.5d0*(s% dm(k)+s% dm(k+1))
            s% Pi_rate_scale(k+1) = s% Pi_rate_scale(k+1) + coeff*(s% Pi_scale(k)+s% Pi_scale(k+1))/dm_face
            s% Phi_rate_scale(k+1) = s% Phi_rate_scale(k+1) + coeff*(s% Phi_scale(k)+s% Phi_scale(k+1))/dm_face
         end do
      end subroutine set_etrb_start_vars


      subroutine rsp2_predict_moments(s,k_lo,k_hi,ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k_lo, k_hi
         integer, intent(out) :: ierr
         integer :: n, i, j, l, k, iter, row, col, n_eq, bracket_iter
         integer, allocatable :: ipiv(:)
         real(dp), allocatable :: moments(:,:), next(:,:), scale(:,:), band(:,:), rhs(:), &
            buoyancy(:), entropy_gradient(:), Lambda(:), rad_time(:), Eq_div_w(:), &
            work_new(:), work_start(:), transport_coeff(:), conductance(:), dm_face(:)
         type(auto_diff_real_star_order1) :: b_ad, h_ad, Lambda_ad, rad_ad, &
            Eq_ad, work_ad, work_new_ad, work_start_ad, Lt_ad, coeff_ad
         real(dp) :: C_D, C_Pi, C_Phi, w, inverse_dt, entropy_per_velocity, &
            energy_decay, Phi_decay, iteration_rate, outer_rate, inner_rate, &
            block(3,3), rhs_face(3), row_scale, change, moment_norm, &
            lower_w, upper_w, trial_w, entropy_ratio, linear_rate, discr
         logical :: supported

         ierr = 0
         ! Outside the predictor's paired decay assumptions, keep the current guess.
         if (s% RSP2_alfa_pi < 1d0 .or. &
               (s% RSP2_alfad == 0d0 .and. s% RSP2_alfa_phi == 0d0)) return
         n = k_hi-k_lo+1
         n_eq = 3*n
         allocate(moments(3,n),next(3,n),scale(3,n),band(10,n_eq),rhs(n_eq),ipiv(n_eq), &
            buoyancy(n),entropy_gradient(n),Lambda(n),rad_time(n),Eq_div_w(n), &
            work_new(n),work_start(n),transport_coeff(0:n),conductance(0:n),dm_face(n),stat=ierr)
         if (ierr /= 0) return
         C_D = s% RSP2_alfad*x_CEDE
         C_Phi = s% RSP2_alfa_phi*x_ALFAPHI
         C_Pi = 0.5d0*(C_D+C_Phi)*s% RSP2_alfa_pi
         inverse_dt = 1d0/s% dt
         transport_coeff = 0d0
         do i=1,n
            k = k_lo+i-1
            call rsp2_moment_coefficients(s,k,b_ad,h_ad,Lambda_ad,rad_ad,ierr)
            if (ierr /= 0) return
            buoyancy(i) = b_ad%val
            entropy_gradient(i) = h_ad%val
            Lambda(i) = Lambda_ad%val
            rad_time(i) = rad_ad%val
            Eq_ad = compute_Eq_div_w_face(s,k,ierr)
            if (ierr /= 0) return
            Eq_div_w(i) = Eq_ad%val
            work_ad = calc_Ptrb_work_face(s,k,work_new=work_new_ad,work_start=work_start_ad)
            work_new(i) = work_new_ad%val
            work_start(i) = work_start_ad%val
            if (work_start(i) >= 1d0 .or. Eq_div_w(i) < 0d0) then
               ierr = -1
               s% retry_message = 'invalid RSP3 pressure or viscous work'
               return
            end if
            dm_face(i) = 0.5d0*(s% dm(k-1)+s% dm(k))
            if (.not. rsp2_moments_valid(s% w_start(k), &
                  s% xh_start(s% i_Pi,k),s% xh_start(s% i_Phi,k))) then
               ierr = -1
               s% retry_message = 'invalid RSP3 starting moments'
               return
            end if
            moments(:,i) = [pow2(s% w_start(k)),s% xh_start(s% i_Pi,k),s% xh_start(s% i_Phi,k)]
            if (rsp2_moments_valid(s% w(k),s% Pi(k),s% Phi(k))) &
               moments(:,i) = [pow2(s% w(k)),s% Pi(k),s% Phi(k)]
            scale(:,i) = [max(pow2(s% w(k)),1.5d0*pow2(s% Pi_scale(k)/sqrt(s% Phi_scale(k)))), &
               s% Pi_scale(k),s% Phi_scale(k)]
            if (i == n) cycle
            Lt_ad = compute_Lt_center(s,k,ierr,transport_coeff=coeff_ad)
            if (ierr /= 0) return
            transport_coeff(i) = coeff_ad%val
         end do
         do i=1,n
            if (moments(1,i) == 0d0 .and. Eq_div_w(i) > 0d0) then
               linear_rate = inverse_dt*(1d0+work_new(i))
               w = s% dt*Eq_div_w(i)
               if (C_D > 0d0) then
                  discr = sqrt(pow2(linear_rate)+4d0*C_D*Eq_div_w(i)/Lambda(i))
                  if (linear_rate >= 0d0) then
                     w = 2d0*Eq_div_w(i)/(linear_rate+discr)
                  else
                     w = 0.5d0*(discr-linear_rate)*Lambda(i)/C_D
                  end if
               else if (linear_rate > 0d0) then
                  w = Eq_div_w(i)/linear_rate
               end if
               moments(1,i) = pow2(w)
            end if
            supported = any(moments(:,i) /= 0d0)
            if (i > 1) supported = supported .or. &
               (transport_coeff(i-1) > 0d0 .and. any(moments(:,i-1) /= 0d0))
            if (i < n) supported = supported .or. &
               (transport_coeff(i) > 0d0 .and. any(moments(:,i+1) /= 0d0))
            if (.not. supported .or. buoyancy(i)*entropy_gradient(i) <= 0d0) cycle
            lower_w = sqrt(moments(1,i))
            if (positive_local_resolvent(lower_w,i)) cycle
            upper_w = max(2d0*lower_w,Lambda(i)* &
               max(inverse_dt,abs(work_new(i))*inverse_dt,rad_time(i), &
                  sqrt(abs(buoyancy(i)*entropy_gradient(i))))/max(C_D,C_Phi))
            do bracket_iter=1,100
               if (positive_local_resolvent(upper_w,i)) exit
               upper_w = 2d0*upper_w
            end do
            if (bracket_iter > 100 .or. is_bad(upper_w)) then
               ierr = -1
               s% retry_message = 'no positive RSP3 moment guess'
               return
            end if
            do bracket_iter=1,40
               trial_w = 0.5d0*(lower_w+upper_w)
               if (positive_local_resolvent(trial_w,i)) then
                  upper_w = trial_w
               else
                  lower_w = trial_w
               end if
            end do
            ! Select the positive nonlinear solution without changing the old state.
            w = 2d0*upper_w
            entropy_ratio = sqrt(3d0*abs(entropy_gradient(i)/buoyancy(i)))
            moments(:,i) = [w*w,sign(1d0,entropy_gradient(i))*(2d0/3d0)*w*w*entropy_ratio, &
               (2d0/3d0)*pow2(w*entropy_ratio)]
         end do
         ! One unit conversion gives a common covariance norm across the mesh.
         entropy_per_velocity = sqrt(maxval(scale(3,:)))/sqrt((2d0/3d0)*maxval(scale(1,:)))
         do iter=1,200
            conductance = 0d0
            do i=1,n-1
               conductance(i) = transport_coeff(i)*0.5d0*(sqrt(moments(1,i))+sqrt(moments(1,i+1)))
            end do
            band = 0d0
            do i=1,n
               k = k_lo+i-1
               w = sqrt(moments(1,i))
               energy_decay = C_D*w/Lambda(i) + work_new(i)*inverse_dt
               Phi_decay = C_Phi*w/Lambda(i) + 2d0*rad_time(i)
               ! Keep this bound independent of the iterate's turnover damping.
               ! Otherwise a weak growing solution can collapse onto the wrong root.
               iteration_rate = max(0d0, &
                  abs(buoyancy(i)*entropy_per_velocity/3d0 + entropy_gradient(i)/entropy_per_velocity) - &
                  min(work_new(i)*inverse_dt,2d0*rad_time(i)))
               outer_rate = conductance(i-1)/dm_face(i)
               inner_rate = conductance(i)/dm_face(i)
               block = 0d0
               block(1,1) = inverse_dt + energy_decay
               block(1,2) = -buoyancy(i)
               block(2,1) = -(2d0/3d0)*entropy_gradient(i)
               block(2,2) = inverse_dt + 0.5d0*work_new(i)*inverse_dt + C_Pi*w/Lambda(i) + rad_time(i)
               block(2,3) = -buoyancy(i)/3d0
               block(3,2) = -2d0*entropy_gradient(i)
               block(3,3) = inverse_dt + Phi_decay
               do j=1,3
                  block(j,j) = block(j,j) + iteration_rate + outer_rate + inner_rate
               end do
               rhs_face(1) = (1d0-work_start(i))*pow2(s% w_start(k))*inverse_dt + Eq_div_w(i)*w
               rhs_face(2) = sqrt(1d0-work_start(i))*s% xh_start(s% i_Pi,k)*inverse_dt
               rhs_face(3) = s% xh_start(s% i_Phi,k)*inverse_dt
               ! The iteration term cancels at a fixed point of the original three rows.
               rhs_face = rhs_face + iteration_rate*moments(:,i)
               do j=1,3
                  row = 3*(i-1)+j
                  row_scale = maxval(abs(block(j,:)*scale(:,i)))
                  if (i > 1) row_scale = max(row_scale,outer_rate*scale(j,i-1))
                  if (i < n) row_scale = max(row_scale,inner_rate*scale(j,i+1))
                  if (row_scale <= 0d0 .or. is_bad(row_scale)) then
                     ierr = -1
                     s% retry_message = 'invalid RSP3 predictor scale'
                     return
                  end if
                  rhs(row) = rhs_face(j)/row_scale
                  do l=1,3
                     col = 3*(i-1)+l
                     band(7+row-col,col) = block(j,l)*scale(l,i)/row_scale
                  end do
                  if (i > 1) band(10,row-3) = -outer_rate*scale(j,i-1)/row_scale
                  if (i < n) band(4,row+3) = -inner_rate*scale(j,i+1)/row_scale
               end do
            end do
            ! The three-moment stencil has scalar lower and upper bandwidth three.
            call DGBTRF(n_eq,n_eq,3,3,band,10,ipiv,ierr)
            if (ierr /= 0) then
               s% retry_message = 'RSP3 predictor factorization failed'
               return
            end if
            call DGBTRS('N',n_eq,3,3,1,band,10,ipiv,rhs,n_eq,ierr)
            if (ierr /= 0) then
               s% retry_message = 'RSP3 predictor solve failed'
               return
            end if
            next = reshape(rhs,[3,n])*scale
            do i=1,n
               if (is_bad(next(1,i)) .or. next(1,i) < 0d0) then
                  ierr = -1
               else if (.not. rsp2_moments_valid(sqrt(next(1,i)),next(2,i),next(3,i))) then
                  ierr = -1
               end if
               if (ierr /= 0) then
                  s% retry_message = 'invalid covariance in RSP3 predictor'
                  return
               end if
            end do
            change = 0d0
            do i=1,n
               moment_norm = max(maxval(abs(next(:,i))/scale(:,i)), &
                  maxval(abs(moments(:,i))/scale(:,i)),tiny(1d0))
               change = max(change,maxval(abs(next(:,i)-moments(:,i))/scale(:,i))/moment_norm)
            end do
            moments = next
            if (change <= 1d-8) exit
         end do
         s% w(k_lo:k_hi) = sqrt(moments(1,:))
         s% Pi(k_lo:k_hi) = moments(2,:)
         s% Phi(k_lo:k_hi) = moments(3,:)
         if (s% RSP2_report_adjust_w) &
            write(*,'(a,i6,es16.7)') 'RSP3 moment predictor ', min(iter,200), change

         contains

         logical function positive_local_resolvent(w,i) result(positive)
            real(dp), intent(in) :: w
            integer, intent(in) :: i
            real(dp) :: energy_rate, Pi_rate, Phi_rate, production
            energy_rate = inverse_dt*(1d0+work_new(i)) + C_D*w/Lambda(i)
            Pi_rate = inverse_dt*(1d0+0.5d0*work_new(i)) + C_Pi*w/Lambda(i) + rad_time(i)
            Phi_rate = inverse_dt + C_Phi*w/Lambda(i) + 2d0*rad_time(i)
            positive = .false.
            if (energy_rate <= 0d0 .or. Phi_rate <= 0d0) return
            production = (2d0/3d0)*buoyancy(i)*entropy_gradient(i)
            positive = Pi_rate > production/energy_rate + production/Phi_rate
         end function positive_local_resolvent
      end subroutine rsp2_predict_moments


      subroutine RSP2_adjust_vars_before_call_solver(s,ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, pass, k_lo, k_hi, k_first, k_last, k_step
         real(dp) :: source_coeff, linear_coeff, available_energy, discr, soln, w_initial
         type(auto_diff_real_star_order1) :: source, damping, rad_damping, Eq_face, dLt_dm, work, rhs

         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. s% dt <= 0d0) return
         k_lo = max(2,s% RSP2_num_outermost_cells_forced_nonturbulent+2)
         k_hi = s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
         if (k_lo > k_hi) return
         if (s% RSP2_3equation_flag) then
            call rsp2_predict_moments(s,k_lo,k_hi,ierr)
            ! The first solver evaluation uses the current cached fluxes.
            if (ierr == 0) call set_RSP2_vars(s,ierr)
            return
         end if
         do k=k_lo,k_hi
            if (s% w(k) /= 0d0) cycle
            source = compute_Source_div_w(s,k,ierr)
            if (ierr /= 0) return
            Eq_face = compute_Eq_div_w_face(s,k,ierr)
            if (ierr /= 0) return
            source_coeff = source%val + Eq_face%val
            if (source_coeff > 0d0) s% w(k) = s% dt*source_coeff
         end do
         ! Seed imported energy in both directions without changing accepted state.
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
               Eq_face = compute_Eq_face(s,k,ierr)
               if (ierr /= 0) return
               dLt_dm = rsp2_dLt_dm_face(s,k,ierr)
               if (ierr /= 0) return
               work = calc_Ptrb_work_face(s,k)
               rhs = s% dt*(source - damping - rad_damping + Eq_face - dLt_dm) - work
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
         ! Refresh fluxes after changing w, as on the three equation path.
         call set_RSP2_vars(s,ierr)
      end subroutine RSP2_adjust_vars_before_call_solver
      end module hydro_rsp2
