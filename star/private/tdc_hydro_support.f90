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

module tdc_hydro_support

   use star_private_def
   use const_def, only: dp, ln10, pi, lsun, rsun
   use utils_lib, only: is_bad
   use auto_diff
   use auto_diff_support
   use star_utils

   implicit none

   private
   public :: remesh_for_TDC_pulsations

contains

   subroutine remesh_for_TDC_pulsations(s, ierr)
      ! uses these controls
      !  TDC_hydro_nz = 190
      !  TDC_hydro_nz_outer = 40
      !  TDC_hydro_nz_inner = 20
      !  TDC_hydro_nz_T_gradient = 20
      !  TDC_hydro_T_anchor = 11d3
      !  TDC_hydro_dq_1_factor = 2d0
      use interp_1d_def, only: pm_work_size
      use interp_1d_lib, only: interpolate_vector_pm
      use hydro_vars, only: set_cgrav
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      integer :: k, j, nz_old, nz
      real(dp) :: xm_anchor, P_surf, T_surf, old_L1, old_r1, old_J, old_abs_J
      real(dp), allocatable, dimension(:) :: &
         xm_old, xm, xm_mid_old, xm_mid, v_old, v_new
      real(dp), pointer :: work1(:)  ! =(nz_old+1, pm_work_size)
      include 'formats'
      ierr = 0
      nz_old = s%nz
      nz = s%TDC_hydro_nz
      call validate_controls2
      if (ierr /= 0) return
      call setvars2(ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in setvars')
      old_J = 0d0
      old_abs_J = 0d0
      if (s%rotation_flag) then
         old_J = dot_product(s%dm_bar(1:nz_old), s%j_rot(1:nz_old))
         old_abs_J = dot_product(s%dm_bar(1:nz_old), abs(s%j_rot(1:nz_old)))
      end if
      old_L1 = s%L(1)
      old_r1 = s%r(1)
      call set_phot_info(s)  ! sets Teff
      call get_PT_surf2(P_surf, T_surf, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in get_PT_surf')
      allocate ( &
         xm_old(nz_old + 1), xm_mid_old(nz_old), v_old(nz_old + 1), &
         xm(nz + 1), xm_mid(nz), v_new(nz + 1), work1((nz_old + 1)*pm_work_size))
      call set_xm_old2
      call find_xm_anchor2
      if (ierr /= 0) then
         deallocate(work1)
         return
      end if
      call set_xm_new2
      if (ierr /= 0) then
         deallocate(work1)
         return
      end if
      call interpolate1_face_val2(s%i_lnR, log(max(1d0, s%r_center)))
      call check_new_lnR2
      call interpolate1_face_val2(s%i_lum, s%L_center)
      if (s%i_v /= 0) call interpolate1_face_val2(s%i_v, s%v_center)
      if (s%have_mlt_vc) then
         ! Preserve the lagged TDC state on the new face grid.
         v_old(1:nz_old) = s%mlt_vc(1:nz_old)
         if (s%r_center == 0d0) then
            v_old(nz_old + 1) = 0d0
         else
            ! A truncated envelope has no center symmetry condition.
            v_old(nz_old + 1) = s%mlt_vc(nz_old)
         end if
         call interpolate_vector_pm( &
            nz_old + 1, xm_old, nz + 1, xm, v_old, v_new, work1, &
            'remesh_for_TDC mlt_vc', ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC remesh mlt_vc interpolation failed')
         s%mlt_vc(1:nz) = v_new(1:nz)
      end if
      call set_new_lnd2
      call interpolate1_cell_val2(s%i_lnT)
      if (s%i_u /= 0) call interpolate1_cell_val2(s%i_u)
      do j = 1, s%species
         call remap1_xa2(j)
      end do
      if (s%rotation_flag) call remap_rotation2(old_J, old_abs_J)
      s%nz = nz
      call update_composition_info2
      call set_cgrav(s, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in set_cgrav')
      call revise_lnT_for_QHSE2(P_surf, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in revise_lnT_for_QHSE')
      if (s%rotation_flag) call set_rotation_seed2
      deallocate (work1)
      write (*, 1) 'new old L_surf/Lsun', s%xh(s%i_lum, 1)/Lsun, old_L1/Lsun
      write (*, 1) 'new old R_surf/Rsun', exp(s%xh(s%i_lnR, 1))/Rsun, old_r1/Rsun
      write (*, '(A)')

   contains

      subroutine validate_controls2
         integer :: nz_base
         include 'formats'
         nz_base = nz - s%TDC_hydro_nz_T_gradient
         if (s%RSP_flag .or. s%RSP2_flag) then
            write(*,'(A)') 'TDC remesh cannot be applied after enabling RSP or RSP2'
            ierr = -1
         else if (nz > nz_old) then
            write(*,3) 'TDC remesh cannot increase the number of zones', nz, nz_old
            ierr = -1
         else if (s%TDC_hydro_nz_T_gradient < 0) then
            write(*,1) 'TDC_hydro_nz_T_gradient must be nonnegative', &
               real(s%TDC_hydro_nz_T_gradient, dp)
            ierr = -1
         else if (s%TDC_hydro_nz_inner < 0) then
            write(*,1) 'TDC_hydro_nz_inner must be nonnegative', &
               real(s%TDC_hydro_nz_inner, dp)
            ierr = -1
         else if (s%TDC_hydro_nz_outer < 2) then
            write(*,1) 'TDC_hydro_nz_outer must be at least 2', &
               real(s%TDC_hydro_nz_outer, dp)
            ierr = -1
         else if (s%remesh_for_TDC_pulsations_log_core_zoning .and. &
               nz_base - s%TDC_hydro_nz_outer < 2) then
            write(*,3) 'TDC remesh needs at least two interior zones', &
               nz_base, s%TDC_hydro_nz_outer
            ierr = -1
         else if (.not. s%remesh_for_TDC_pulsations_log_core_zoning .and. &
               nz_base - s%TDC_hydro_nz_outer - s%TDC_hydro_nz_inner < 2) then
            write(*,3) 'TDC remesh needs at least two middle zones', &
               nz_base, s%TDC_hydro_nz_outer, s%TDC_hydro_nz_inner
            ierr = -1
         else if (.not. s%remesh_for_TDC_pulsations_log_core_zoning .and. &
               s%TDC_hydro_nz_inner > 0 .and. s%max_center_cell_dq <= 0d0) then
            write(*,1) 'max_center_cell_dq must be positive for inner TDC zoning', &
               s%max_center_cell_dq
            ierr = -1
         else if (s%TDC_hydro_dq_1_factor <= 0d0) then
            write(*,1) 'TDC_hydro_dq_1_factor must be positive', &
               s%TDC_hydro_dq_1_factor
            ierr = -1
         else if (s%TDC_hydro_T_anchor <= 0d0) then
            write(*,1) 'TDC_hydro_T_anchor must be positive', s%TDC_hydro_T_anchor
            ierr = -1
         end if
      end subroutine validate_controls2

      subroutine setvars2(ierr)
         use hydro_vars, only: unpack_xh, set_hydro_vars
         integer, intent(out) :: ierr
         logical, parameter :: &
            skip_basic_vars = .false., &
            skip_micro_vars = .false., &
            skip_m_grav_and_grav = .false., &
            skip_net = .true., &
            skip_neu = .true., &
            skip_kap = .false., &
            skip_grads = .true., &
            skip_rotation = .true., &
            skip_brunt = .true., &
            skip_other_cgrav = .true., &
            skip_mixing_info = .true., &
            skip_set_cz_bdy_mass = .true., &
            skip_mlt = .true., &
            skip_eos = .false.
         ierr = 0
         call unpack_xh(s, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in unpack_xh')
         call set_hydro_vars( &
            s, 1, nz_old, skip_basic_vars, &
            skip_micro_vars, skip_m_grav_and_grav, skip_eos, skip_net, skip_neu, &
            skip_kap, skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'remesh_for_TDC failed in set_hydro_vars')
      end subroutine setvars2

      subroutine get_PT_surf2(P_surf, T_surf, ierr)
         use atm_support, only: get_atm_PT
         real(dp), intent(out) :: P_surf, T_surf
         integer, intent(out) :: ierr
         real(dp) :: &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         logical, parameter :: skip_partials = .true.
         ierr = 0
         call set_phot_info(s)  ! sets s% Teff
         Teff = s%Teff
         call get_atm_PT( &  ! this uses s% opacity(1)
            s, s%tau_factor*s%tau_base, s%L(1), s%r(1), s%m(1), s%cgrav(1), skip_partials, &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'get_P_surf failed in get_atm_PT')
         P_surf = exp(lnP_surf)
         T_surf = exp(lnT_surf)
      end subroutine get_PT_surf2

      subroutine set_xm_old2
         xm_old(1) = 0d0
         do k = 2, nz_old
            xm_old(k) = xm_old(k - 1) + s%dm(k - 1)
         end do
         xm_old(nz_old + 1) = s%xmstar
         do k = 1, nz_old
            xm_mid_old(k) = xm_old(k) + 0.5d0*s%dm(k)
         end do
      end subroutine set_xm_old2

      subroutine find_xm_anchor2
         real(dp) :: lnT_anchor, xmm1, xm00, lnTm1, lnT00
         include 'formats'
         lnT_anchor = log(s%TDC_hydro_T_anchor)
         if (lnT_anchor <= s%xh(s%i_lnT, 1)) then
            write (*, 1) 'T_anchor < T_surf', s%TDC_hydro_T_anchor, exp(s%xh(s%i_lnT, 1))
            call mesa_error(__FILE__, __LINE__, 'find_xm_anchor')
         end if
         xm_anchor = xm_old(nz_old)
         do k = 2, nz_old
            if (s%xh(s%i_lnT, k) >= lnT_anchor) then
               xmm1 = xm_mid_old(k - 1)
               xm00 = xm_mid_old(k)
               lnTm1 = s%xh(s%i_lnT, k - 1)
               lnT00 = s%xh(s%i_lnT, k)
               xm_anchor = xmm1 + &
                           (xm00 - xmm1)*(lnT_anchor - lnTm1)/(lnT00 - lnTm1)
               if (is_bad(xm_anchor) .or. xm_anchor <= 0d0) then
                  write (*, 2) 'bad xm_anchor', k, xm_anchor, xmm1, xm00, lnTm1, lnT00, lnT_anchor, s%lnT(1)
                  call mesa_error(__FILE__, __LINE__, 'find_xm_anchor')
               end if
               return
            end if
         end do
         write(*,1) 'T_anchor exceeds the model temperature', s%TDC_hydro_T_anchor
         ierr = -1
      end subroutine find_xm_anchor2

      subroutine set_xm_new2  ! sets xm, dm, m, dq, q
         integer :: nz_outer, nz_inner, nz_T_gradient, nz_base, n_middle, k
         real(dp) :: dq_1_factor, dxm_outer, lnx, dlnx, base_dm, &
            rem_mass, center_dm, peak_dm, correction, ratio, H, H_inner
         real(dp) :: H_low, H_high, H_mid, f_low, f_high, f_mid
         integer :: iter
         include 'formats'
         nz_outer = s%TDC_hydro_nz_outer
         nz_inner = s%TDC_hydro_nz_inner
         nz_T_gradient = s%TDC_hydro_nz_T_gradient
         nz_base = nz - nz_T_gradient
         dq_1_factor = s%TDC_hydro_dq_1_factor
         ratio = max(2.5d0, s%mesh_max_allowed_ratio)
         dxm_outer = xm_anchor/(nz_outer - 1d0 + dq_1_factor)
         xm(1) = 0d0
         xm(2) = dxm_outer*dq_1_factor
         s%dm(1) = xm(2)
         do k = 3, nz_outer + 1
            xm(k) = xm(k - 1) + dxm_outer
            s%dm(k - 1) = dxm_outer
         end do

         if (.not. s%remesh_for_TDC_pulsations_log_core_zoning) then
            ! do rsp style core zoning with a power law on dq
            n_middle = nz_base - nz_outer - nz_inner
            rem_mass = s%xmstar - xm(nz_outer + 1)
            base_dm = dxm_outer

            if (nz_inner == 0) then
               ! Original single inward-increasing power law.
               H_low = 1.001d0
               H_high = 1.40d0
               f_low = base_dm*(1d0 - pow(H_low, real(n_middle, dp)))/(1d0 - H_low) - rem_mass
               f_high = base_dm*(1d0 - pow(H_high, real(n_middle, dp)))/(1d0 - H_high) - rem_mass
               if (f_low*f_high > 0d0) then
                  write(*,2) 'failed to bracket TDC core zoning ramp', &
                     n_middle, f_low, f_high
                  ierr = -1
                  return
               end if
               do iter = 1, 1000
                  H_mid = 0.5d0*(H_low + H_high)
                  f_mid = base_dm*(1d0 - pow(H_mid, real(n_middle, dp)))/(1d0 - H_mid) - rem_mass
                  if (abs(f_mid) < 1d-12*rem_mass) exit
                  if (f_low*f_mid <= 0d0) then
                     H_high = H_mid
                     f_high = f_mid
                  else
                     H_low = H_mid
                     f_low = f_mid
                  end if
               end do
               H = H_mid

               s%dm(nz_outer + 1) = base_dm
               do k = nz_outer + 2, nz_base - 1
                  s%dm(k) = H*s%dm(k - 1)
               end do
               s%dm(nz_base) = s%xmstar - sum(s%dm(1:nz_base - 1))

            else
               ! Add an inward-decreasing ramp that ends at max_center_cell_dq.
               H_low = 1d0
               H_high = ratio
               center_dm = s%max_center_cell_dq*s%xmstar
               H_low = max(H_low, &
                  pow(center_dm/base_dm, 1d0/real(n_middle - 1, dp)))
               if (H_low > H_high) then
                  write(*,2) 'TDC inner zoning requires a larger mesh ratio', &
                     H_low, H_high
                  ierr = -1
                  return
               end if
               f_low = double_ramp_mass2( &
                  H_low, base_dm, center_dm, n_middle, nz_inner) - rem_mass
               f_high = double_ramp_mass2( &
                  H_high, base_dm, center_dm, n_middle, nz_inner) - rem_mass
               if (f_low*f_high > 0d0) then
                  write(*,2) 'failed to bracket TDC core zoning ramp', &
                     n_middle, f_low, f_high
                  ierr = -1
                  return
               end if
               do iter = 1, 1000
                  H_mid = 0.5d0*(H_low + H_high)
                  f_mid = double_ramp_mass2( &
                     H_mid, base_dm, center_dm, n_middle, nz_inner) - rem_mass
                  if (abs(f_mid) < 1d-12*rem_mass) exit
                  if (f_low*f_mid <= 0d0) then
                     H_high = H_mid
                     f_high = f_mid
                  else
                     H_low = H_mid
                     f_low = f_mid
                  end if
               end do
               H = H_mid

               do k = nz_outer + 1, nz_outer + n_middle
                  s%dm(k) = base_dm*pow(H, real(k - nz_outer - 1, dp))
               end do

               peak_dm = s%dm(nz_outer + n_middle)
               H_inner = pow(peak_dm/center_dm, 1d0/real(nz_inner, dp))
               do k = nz_outer + n_middle + 1, nz_base
                  s%dm(k) = center_dm*pow(H_inner, real(nz_base - k, dp))
               end do

               correction = s%xmstar - sum(s%dm(1:nz_base))
               s%dm(nz_outer + n_middle) = &
                  s%dm(nz_outer + n_middle) + correction
            end if

         else ! use log zoning inward from anchor to core.
            lnx = log(xm(nz_outer + 1))
            if (is_bad(lnx)) then
               write (*, 2) 'bad lnx', nz_outer + 1, lnx, xm(nz_outer + 1)
               call mesa_error(__FILE__, __LINE__, 'set_xm_new')
            end if
            dlnx = (log(s%xmstar) - lnx)/(nz_base - nz_outer)
            do k = nz_outer + 2, nz_base
               lnx = lnx + dlnx
               xm(k) = exp(lnx)
               s%dm(k - 1) = xm(k) - xm(k - 1)
            end do
            s%dm(nz_base) = s%xmstar - xm(nz_base)

            ! enforce the last boundary at total mass
            xm(nz_base + 1) = s%xmstar

            ! recompute cell masses
            do k = nz_outer + 1, nz_base
               s%dm(k) = xm(k + 1) - xm(k)
            end do

         end if

         xm(1) = 0d0
         do k = 1, nz_base
            xm(k + 1) = xm(k) + s%dm(k)
         end do
         xm(nz_base + 1) = s%xmstar

         if (nz_T_gradient > 0) then
            call add_T_gradient_zones2(nz_outer, nz_base, nz_T_gradient)
            if (ierr /= 0) return
         end if

         do k = 1, nz
            if (s%dm(k) <= 0d0 .or. is_bad(s%dm(k))) then
               write(*,2) 'bad cell mass in TDC remesh', k, s%dm(k)
               ierr = -1
               return
            end if
         end do

         if ((.not. s%remesh_for_TDC_pulsations_log_core_zoning .and. nz_inner > 0) .or. &
               nz_T_gradient > 0) then
            do k = 2, nz
               if (s%dm(k) > ratio*(1d0 + 1d-10)*s%dm(k - 1) .or. &
                     s%dm(k - 1) > ratio*(1d0 + 1d-10)*s%dm(k)) then
                  write(*,2) 'bad adjacent cell mass ratio in TDC mesh', &
                     k, s%dm(k)/s%dm(k - 1), ratio
                  ierr = -1
                  return
               end if
            end do
         end if

         xm(1) = 0d0
         do k = 1, nz
            xm(k + 1) = xm(k) + s%dm(k)
         end do
         xm(nz + 1) = s%xmstar

         do k = 1, nz - 1
            xm_mid(k) = 0.5d0*(xm(k) + xm(k + 1))
         end do
         xm_mid(nz) = 0.5d0*(xm(nz) + s%xmstar)
         s%m(1) = s%mstar
         s%q(1) = 1d0
         s%dq(1) = s%dm(1)/s%xmstar
         do k = 2, nz
            s%m(k) = s%m(k - 1) - s%dm(k - 1)
            s%dq(k) = s%dm(k)/s%xmstar
            s%q(k) = s%q(k - 1) - s%dq(k - 1)
         end do
         call set_dm_bar(s, nz, s%dm, s%dm_bar)
      end subroutine set_xm_new2

      real(dp) function geometric_sum2(H, n) result(sum_H)
         real(dp), intent(in) :: H
         integer, intent(in) :: n

         if (abs(H - 1d0) < 1d-12) then
            sum_H = real(n, dp)
         else
            sum_H = (pow(H, real(n, dp)) - 1d0)/(H - 1d0)
         end if
      end function geometric_sum2

      real(dp) function double_ramp_mass2( &
            H, base_dm, center_dm, n_middle, n_inner) result(total_mass)
         real(dp), intent(in) :: H, base_dm, center_dm
         integer, intent(in) :: n_middle, n_inner
         real(dp) :: H_inner, peak_dm

         peak_dm = base_dm*pow(H, real(n_middle - 1, dp))
         H_inner = pow(peak_dm/center_dm, 1d0/real(n_inner, dp))
         total_mass = base_dm*geometric_sum2(H, n_middle) + &
            center_dm*geometric_sum2(H_inner, n_inner)
      end function double_ramp_mass2

      subroutine add_T_gradient_zones2(nz_outer, nz_base, nz_T_gradient)
         integer, intent(in) :: nz_outer, nz_base, nz_T_gradient
         integer :: i, i_base, i_old, j, k_new, max_points, npts, nz_core_base
         real(dp) :: dlnT_total, eps_xm, fraction, next_base, next_old, &
            next_xm, target
         real(dp), allocatable :: base_xm(:), lnT_monitor(:), &
            monitor(:), monitor_xm(:)
         include 'formats'

         nz_core_base = nz_base - nz_outer
         max_points = nz_old + nz_core_base + 2
         allocate(base_xm(nz_base + 1), lnT_monitor(max_points), &
            monitor(max_points), monitor_xm(max_points))
         base_xm = xm(1:nz_base + 1)
         eps_xm = 16d0*epsilon(1d0)*max(abs(s%xmstar), 1d0)

         ! Include both old temperature points and base-grid boundaries in the monitor.
         npts = 1
         monitor_xm(npts) = xm_anchor
         i_base = nz_outer + 2
         i_old = 1
         do while (i_old <= nz_old .and. xm_mid_old(i_old) <= xm_anchor + eps_xm)
            i_old = i_old + 1
         end do

         do
            next_base = huge(1d0)
            if (i_base <= nz_base) next_base = base_xm(i_base)
            next_old = huge(1d0)
            if (i_old <= nz_old) next_old = xm_mid_old(i_old)
            next_xm = min(next_base, next_old)
            if (next_xm >= s%xmstar - eps_xm) exit

            if (next_xm > monitor_xm(npts) + eps_xm) then
               npts = npts + 1
               monitor_xm(npts) = next_xm
            end if
            if (next_base <= next_xm + eps_xm) i_base = i_base + 1
            if (next_old <= next_xm + eps_xm) i_old = i_old + 1
         end do
         npts = npts + 1
         monitor_xm(npts) = s%xmstar

         do i = 1, npts
            lnT_monitor(i) = lnT_at_xm2(monitor_xm(i))
         end do
         monitor(1) = 0d0
         do i = 2, npts
            monitor(i) = monitor(i - 1) + &
               abs(lnT_monitor(i) - lnT_monitor(i - 1))
         end do
         dlnT_total = monitor(npts)

         if (dlnT_total > 0d0) then
            do i = 1, npts
               monitor(i) = base_mesh_coordinate2( &
                  monitor_xm(i), base_xm, nz_outer, nz_base) + &
                  nz_T_gradient*monitor(i)/dlnT_total
            end do
         else
            do i = 1, npts
               monitor(i) = base_mesh_coordinate2( &
                  monitor_xm(i), base_xm, nz_outer, nz_base)* &
                  real(nz_core_base + nz_T_gradient, dp)/real(nz_core_base, dp)
            end do
         end if
         monitor(1) = 0d0
         monitor(npts) = real(nz_core_base + nz_T_gradient, dp)

         xm(1:nz_outer + 1) = base_xm(1:nz_outer + 1)
         j = 2
         do k_new = nz_outer + 2, nz
            target = real(k_new - nz_outer - 1, dp)
            do while (j < npts .and. monitor(j) < target)
               j = j + 1
            end do
            if (monitor(j) <= monitor(j - 1)) then
               write(*,2) 'bad TDC temperature-gradient monitor', &
                  j, monitor(j - 1), monitor(j)
               ierr = -1
               deallocate(base_xm, lnT_monitor, monitor, monitor_xm)
               return
            end if
            fraction = (target - monitor(j - 1))/(monitor(j) - monitor(j - 1))
            xm(k_new) = monitor_xm(j - 1) + &
               fraction*(monitor_xm(j) - monitor_xm(j - 1))
         end do
         xm(nz + 1) = s%xmstar
         do i = nz_outer + 1, nz
            s%dm(i) = xm(i + 1) - xm(i)
         end do

         deallocate(base_xm, lnT_monitor, monitor, monitor_xm)
      end subroutine add_T_gradient_zones2

      real(dp) function base_mesh_coordinate2( &
            x, base_xm, nz_outer, nz_base) result(coordinate)
         real(dp), intent(in) :: x, base_xm(:)
         integer, intent(in) :: nz_outer, nz_base
         integer :: i

         if (x <= base_xm(nz_outer + 1)) then
            coordinate = 0d0
            return
         end if
         do i = nz_outer + 1, nz_base
            if (x <= base_xm(i + 1)) then
               coordinate = real(i - nz_outer - 1, dp) + &
                  (x - base_xm(i))/(base_xm(i + 1) - base_xm(i))
               return
            end if
         end do
         coordinate = real(nz_base - nz_outer, dp)
      end function base_mesh_coordinate2

      real(dp) function lnT_at_xm2(x) result(lnT_at_xm)
         real(dp), intent(in) :: x
         integer :: i
         real(dp) :: fraction

         if (nz_old == 1) then
            lnT_at_xm = s%xh(s%i_lnT, 1)
            return
         end if
         if (x <= xm_mid_old(1)) then
            i = 2
         else
            do i = 2, nz_old
               if (x <= xm_mid_old(i)) exit
            end do
            if (i > nz_old) i = nz_old
         end if
         fraction = (x - xm_mid_old(i - 1))/ &
            (xm_mid_old(i) - xm_mid_old(i - 1))
         lnT_at_xm = s%xh(s%i_lnT, i - 1) + fraction* &
            (s%xh(s%i_lnT, i) - s%xh(s%i_lnT, i - 1))
      end function lnT_at_xm2

      subroutine interpolate1_face_val2(i, cntr_val)
         integer, intent(in) :: i
         real(dp), intent(in) :: cntr_val
         do k = 1, nz_old
            v_old(k) = s%xh(i, k)
         end do
         v_old(nz_old + 1) = cntr_val
         call interpolate_vector_pm( &
            nz_old + 1, xm_old, nz + 1, xm, v_old, v_new, work1, 'remesh_for_TDC', ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC remesh face interpolation failed')
         do k = 1, nz
            s%xh(i, k) = v_new(k)
         end do
      end subroutine interpolate1_face_val2

      subroutine check_new_lnR2
         include 'formats'
         do k = 1, nz
            s%lnR(k) = s%xh(s%i_lnR, k)
            s%r(k) = exp(s%lnR(k))
         end do
         do k = 1, nz - 1
            if (s%r(k) <= s%r(k + 1)) then
               write (*, 2) 'bad r', k, s%r(k), s%r(k + 1)
               call mesa_error(__FILE__, __LINE__, 'check_new_lnR remesh TDC')
            end if
         end do
         if (s%r(nz) <= s%r_center) then
            write (*, 2) 'bad r center', nz, s%r(nz), s%r_center
            call mesa_error(__FILE__, __LINE__, 'check_new_lnR remesh TDC')
         end if
      end subroutine check_new_lnR2

      subroutine set_new_lnd2
         real(dp) :: vol, r300, r3p1
         include 'formats'
         do k = 1, nz
            r300 = pow3(s%r(k))
            if (k < nz) then
               r3p1 = pow3(s%r(k + 1))
            else
               r3p1 = pow3(s%r_center)
            end if
            vol = (4d0*pi/3d0)*(r300 - r3p1)
            s%rho(k) = s%dm(k)/vol
            s%lnd(k) = log(s%rho(k))
            s%xh(s%i_lnd, k) = s%lnd(k)
            if (is_bad(s%lnd(k))) then
               write (*, 2) 'bad lnd vol dm r300 r3p1', k, s%lnd(k), vol, s%dm(k), r300, r3p1
               call mesa_error(__FILE__, __LINE__, 'remesh for TDC')
            end if
         end do
      end subroutine set_new_lnd2

      subroutine interpolate1_cell_val2(i)
         integer, intent(in) :: i
         do k = 1, nz_old
            v_old(k) = s%xh(i, k)
         end do
         call interpolate_vector_pm( &
            nz_old, xm_mid_old, nz, xm_mid, v_old, v_new, work1, 'remesh_for_TDC', ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC remesh cell interpolation failed')
         do k = 1, nz
            s%xh(i, k) = v_new(k)
         end do
      end subroutine interpolate1_cell_val2

      subroutine remap1_xa2(j)
         integer, intent(in) :: j
         integer :: k_old, k_scan, k_new
         real(dp) :: overlap, species_mass

         v_old(1:nz_old) = s%xa(j, 1:nz_old)
         k_old = 1
         do k_new = 1, nz
            do while (k_old < nz_old .and. xm_old(k_old + 1) <= xm(k_new))
               k_old = k_old + 1
            end do
            species_mass = 0d0
            k_scan = k_old
            do while (k_scan <= nz_old .and. xm_old(k_scan) < xm(k_new + 1))
               ! A mass-overlap average conserves the total mass of each species.
               overlap = min(xm(k_new + 1), xm_old(k_scan + 1)) - &
                  max(xm(k_new), xm_old(k_scan))
               if (overlap > 0d0) species_mass = species_mass + overlap*v_old(k_scan)
               k_scan = k_scan + 1
            end do
            s%xa(j, k_new) = species_mass/s%dm(k_new)
         end do
      end subroutine remap1_xa2

      subroutine remap_rotation2(old_J, old_abs_J)
         real(dp), intent(in) :: old_J, old_abs_J
         real(dp) :: new_J
         include 'formats'

         v_old(1:nz_old) = s%j_rot(1:nz_old)
         call remap1_face_average2
         s%j_rot(1:nz) = v_new(1:nz)
         if (s%i_j_rot /= 0) s%xh(s%i_j_rot,1:nz) = s%j_rot(1:nz)

         v_old(1:nz_old) = s%omega(1:nz_old)
         call remap1_face_average2
         s%omega(1:nz) = v_new(1:nz)

         new_J = dot_product(s%dm_bar(1:nz), s%j_rot(1:nz))
         if (abs(new_J - old_J) > 1d-12*max(old_abs_J, tiny(1d0))) then
            write(*,1) 'relative angular momentum error in TDC remesh', &
               (new_J - old_J)/max(abs(old_J), old_abs_J, tiny(1d0))
            call mesa_error(__FILE__, __LINE__, 'TDC remesh failed to conserve angular momentum')
         end if
         s%total_angular_momentum = new_J
         s%total_abs_angular_momentum = &
            dot_product(s%dm_bar(1:nz), abs(s%j_rot(1:nz)))
      end subroutine remap_rotation2

      subroutine remap1_face_average2
         integer :: k_old, k_scan, k_new
         real(dp) :: new_outer, new_inner, old_outer, old_inner, &
            overlap, integral

         k_old = 1
         do k_new = 1, nz
            if (k_new == 1) then
               new_outer = 0d0
            else
               new_outer = xm_mid(k_new - 1)
            end if
            if (k_new == nz) then
               new_inner = s%xmstar
            else
               new_inner = xm_mid(k_new)
            end if

            do while (k_old < nz_old .and. xm_mid_old(k_old) <= new_outer)
               k_old = k_old + 1
            end do
            integral = 0d0
            k_scan = k_old
            do while (k_scan <= nz_old)
               if (k_scan == 1) then
                  old_outer = 0d0
               else
                  old_outer = xm_mid_old(k_scan - 1)
               end if
               if (k_scan == nz_old) then
                  old_inner = s%xmstar
               else
                  old_inner = xm_mid_old(k_scan)
               end if
               if (old_outer >= new_inner) exit
               overlap = min(new_inner, old_inner) - max(new_outer, old_outer)
               if (overlap > 0d0) integral = integral + overlap*v_old(k_scan)
               k_scan = k_scan + 1
            end do
            v_new(k_new) = integral/(new_inner - new_outer)
         end do
      end subroutine remap1_face_average2

      subroutine update_composition_info2
         use chem_lib, only: basic_composition_info
         real(dp) :: sumx

         do k = 1, nz
            call basic_composition_info( &
               s%species, s%chem_id, s%xa(1:s%species,k), &
               s%X(k), s%Y(k), s%Z(k), s%abar(k), s%zbar(k), s%z2bar(k), &
               s%z53bar(k), s%ye(k), s%mass_correction(k), sumx)
         end do
      end subroutine update_composition_info2

      subroutine set_rotation_seed2
         use hydro_rotation, only: w_div_w_roche_jrot

         do k = 1, nz
            s%w_div_w_crit_roche(k) = &
               w_div_w_roche_jrot(s%r(k), s%m(k), s%j_rot(k), s%cgrav(k), &
                  s%w_div_wcrit_max, s%w_div_wcrit_max2, s%w_div_wc_flag)
         end do
         if (s%i_w_div_wc /= 0) &
            s%xh(s%i_w_div_wc,1:nz) = s%w_div_w_crit_roche(1:nz)
      end subroutine set_rotation_seed2

      subroutine revise_lnT_for_QHSE2(P_surf, ierr)
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results
         use chem_def, only: chem_isos
         use eos_support, only: solve_eos_given_DP
         use eos_def, only: i_eta, i_lnfree_e
         use kap_def, only: num_kap_fracs
         use kap_support, only: get_kap
         real(dp), intent(in) :: P_surf
         integer, intent(out) :: ierr
         real(dp) :: logRho, logP, logT_guess, &
                     logT_tol, logP_tol, logT, P_m1, P_00, dm_face, &
                     kap_fracs(num_kap_fracs), kap, dlnkap_dlnRho, dlnkap_dlnT, &
                     old_kap, new_P_surf, new_T_surf
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: dres_dxa(num_eos_d_dxa_results, s%species)
         include 'formats'
         ierr = 0
         P_m1 = P_surf
         do k = 1, nz
            s%lnT(k) = s%xh(s%i_lnT, k)
            s%lnR(k) = s%xh(s%i_lnR, k)
            s%r(k) = exp(s%lnR(k))
         end do
         do k = 1, nz
            if (k < nz) then
               dm_face = s%dm_bar(k)
            else
               dm_face = 0.5d0*(s%dm(k - 1) + s%dm(k))
            end if
            P_00 = P_m1 + s%cgrav(k)*s%m(k)*dm_face/(4d0*pi*pow4(s%r(k)))
            logP = log10(P_00)  ! value for QHSE
            s%lnPeos(k) = logP*ln10
            s%Peos(k) = P_00
            logRho = s%lnd(k)/ln10
            logT_guess = s%lnT(k)/ln10
            logT_tol = 1d-11
            logP_tol = 1d-11
            call solve_eos_given_DP( &
               s, k, s%xa(:, k), &
               logRho, logP, logT_guess, logT_tol, logP_tol, &
               logT, res, d_dlnd, d_dlnT, dres_dxa, ierr)
            if (ierr /= 0) then
               write (*, 2) 'solve_eos_given_DP failed', k
               write (*, '(A)')
               write (*, 1) 'sum(xa)', sum(s%xa(:, k))
               do j = 1, s%species
                  write (*, 4) 'xa(j,k) '//trim(chem_isos%name(s%chem_id(j))), j, j + s%nvar_hydro, k, s%xa(j, k)
               end do
               write (*, 1) 'logRho', logRho
               write (*, 1) 'logP', logP
               write (*, 1) 'logT_guess', logT_guess
               write (*, 1) 'logT_tol', logT_tol
               write (*, 1) 'logP_tol', logP_tol
               write (*, '(A)')
               call mesa_error(__FILE__, __LINE__, 'revise_lnT_for_QHSE')
            end if
            s%lnT(k) = logT*ln10
            s%xh(s%i_lnT, k) = s%lnT(k)
            P_m1 = P_00

            if (k == 1) then  ! get opacity and recheck surf BCs
               call get_kap( &
                  s, k, s%zbar(k), s%xa(:, k), logRho, logT, &
                  res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
                  res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
                  kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, &
                  ierr)
               if (ierr /= 0) then
                  write (*, 2) 'get_kap failed', k
                  call mesa_error(__FILE__, __LINE__, 'revise_lnT_for_QHSE')
               end if
               old_kap = s%opacity(1)
               s%opacity(1) = kap  ! for use by atm surf PT
               call get_PT_surf2(new_P_surf, new_T_surf, ierr)
               if (ierr /= 0) then
                  write (*, 2) 'get_PT_surf failed', k
                  call mesa_error(__FILE__, __LINE__, 'revise_lnT_for_QHSE')
               end if
               write (*, 1) 'new old T_surf', new_T_surf, T_surf
               write (*, 1) 'new old P_surf', new_P_surf, P_surf
               write (*, 1) 'new old kap(1)', kap, old_kap
            end if

         end do
      end subroutine revise_lnT_for_QHSE2

   end subroutine remesh_for_TDC_pulsations

end module tdc_hydro_support
