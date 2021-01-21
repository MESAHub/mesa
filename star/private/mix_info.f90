! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
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


      module mix_info

      use const_def
      use num_lib
      use utils_lib
      use star_private_def


      implicit none

      private
      public :: set_mixing_info, set_RTI_mixing_info, do_smoothing_by_mass, &
         update_rotation_mixing_info, set_dPdr_dRhodr_info, get_convection_sigmas, &
         set_dxdt_mix, set_cz_bdy_mass


      contains


      subroutine set_mixing_info(s, skip_set_cz_bdy_mass, ierr)
         ! set convection variables cdc and conv_vel starting from local MLT results.
         ! overshooting can also be added.    and rotation mixing.
         use rates_def, only: i_rate
         use chem_def, only: ipp, icno, i3alf, ih1, ihe4, ic12
         use star_utils, only: start_time, update_time
         use overshoot, only: add_overshooting
         use predictive_mix, only: add_predictive_mixing
         type (star_info), pointer :: s
         logical, intent(in) :: skip_set_cz_bdy_mass
         integer, intent(out) :: ierr

         integer :: nz, i, k, max_conv_bdy, max_mix_bdy, k_dbg, k_Tmax, i_h1, i_he4, i_c12
         real(dp) :: c, rho_face, f, Tmax, conv_vel, min_conv_vel_for_convective_mixing_type, &
            region_bottom_q, region_top_q
         real(dp), pointer, dimension(:) :: eps_h, eps_he, eps_z, cdc_factor

         logical :: rsp_or_eturb, dbg

         integer(8) :: time0
         real(dp) :: total

         include 'formats'

         ierr = 0
         nz = s% nz

         dbg = .false.
         !dbg = .true.
         k_dbg = -1152
         
         min_conv_vel_for_convective_mixing_type = 1d0 ! make this a control parameter
         
         rsp_or_eturb = s% RSP_flag .or. s% Eturb_flag

         if (dbg) write(*, *) 'set_mixing_info'
         if (s% doing_timing) call start_time(s, time0, total)
         
         if (s% RTI_flag) then
            if (dbg) write(*,*) 'call set_RTI_mixing_info'
            call set_RTI_mixing_info(s, ierr)
            if (failed('set_RTI_mixing_info')) return
            if (dbg) write(*,*) 'call set_dPdr_dRhodr_info'
            call set_dPdr_dRhodr_info(s, ierr)
            if (failed('set_dPdr_dRhodr_info')) return
         end if

         nullify(eps_h, eps_he, eps_z, cdc_factor)

         max_conv_bdy = 10 ! will automatically be increased if necessary
         max_mix_bdy = 10 ! will automatically be increased if necessary

         s% num_conv_boundaries = 0
         if (.not. associated(s% conv_bdy_loc)) allocate(s% conv_bdy_loc(max_conv_bdy))
         if (.not. associated(s% conv_bdy_q)) allocate(s% conv_bdy_q(max_conv_bdy))
         if (.not. associated(s% top_conv_bdy)) allocate(s% top_conv_bdy(max_conv_bdy))
         if (.not. associated(s% burn_h_conv_region)) allocate(s% burn_h_conv_region(max_conv_bdy))
         if (.not. associated(s% burn_he_conv_region)) allocate(s% burn_he_conv_region(max_conv_bdy))
         if (.not. associated(s% burn_z_conv_region)) allocate(s% burn_z_conv_region(max_conv_bdy))

         s% num_mix_boundaries = 0
         if (.not. associated(s% mix_bdy_loc)) allocate(s% mix_bdy_loc(max_mix_bdy))
         if (.not. associated(s% mix_bdy_q)) allocate(s% mix_bdy_q(max_mix_bdy))
         if (.not. associated(s% top_mix_bdy)) allocate(s% top_mix_bdy(max_mix_bdy))
         if (.not. associated(s% burn_h_mix_region)) allocate(s% burn_h_mix_region(max_mix_bdy))
         if (.not. associated(s% burn_he_mix_region)) allocate(s% burn_he_mix_region(max_mix_bdy))
         if (.not. associated(s% burn_z_mix_region)) allocate(s% burn_z_mix_region(max_mix_bdy))

         call do_alloc(ierr)
         if (ierr /= 0) return
         
         if (.not. s% RSP_flag) then
            do k = 1, nz
               s% mixing_type(k) = s% mlt_mixing_type(k)
            end do
         end if
         
         cdc_factor(1) = 1d0
         do k = 2, nz
            rho_face = (s% dq(k-1)*s% rho(k) + s% dq(k)*s% rho(k-1))/&
                           (s% dq(k-1) + s% dq(k))
            f = 4d0*pi*s% r(k)*s% r(k)*rho_face
            cdc_factor(k) = f*f
         end do
         
         if (s% RSP_flag) then
            do k = 1, nz
               s% mixing_type(k) = no_mixing
               s% D_mix(k) = 0d0
               s% cdc(k) = 0d0
               s% conv_vel(k) = 0d0
            end do
         else if (s% conv_vel_flag .or. s% Eturb_flag) then
            do k = 1, nz
               if (s% Eturb_flag) then
                  s% conv_vel(k) = sqrt2*sqrt(s% Eturb(k))
                  if (s% conv_vel(k) >= min_conv_vel_for_convective_mixing_type) then
                     s% mixing_type(k) = convective_mixing
                  else
                     s% mixing_type(k) = no_mixing
                  end if
               end if
               s% D_mix(k) = s% conv_vel(k)*s% mlt_mixing_length(k)/3d0
               if (s% conv_vel_ignore_thermohaline) then
                  s% D_mix(k) = s% D_mix(k) + s% mlt_D_thrm(k)
               end if
               if (s% conv_vel_ignore_semiconvection) then
                  s% D_mix(k) = s% D_mix(k) + s% mlt_D_semi(k)
               end if
               s% cdc(k) = cdc_factor(k)*s% D_mix(k)
            end do
         else
            do k = 1, nz
               s% D_mix(k) = s% mlt_D(k)
               s% cdc(k) = s% mlt_cdc(k)
               s% conv_vel(k) = s% mlt_vc(k)
            end do
         end if
         
         call check('after get mlt_D')
         
         if (dbg) write(*,3) 'after copy mlt results', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         
         if (s% remove_mixing_glitches) then

            if (dbg) write(*, *) 'remove_mixing_glitches'

            if (dbg) write(*,3) 'call remove_tiny_mixing', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call remove_tiny_mixing(s, ierr)
            if (failed('remove_tiny_mixing')) return

            if (dbg) write(*,3) 'call remove_mixing_singletons', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call remove_mixing_singletons(s, ierr)
            if (failed('remove_mixing_singletons')) return

            if (dbg) write(*,3) 'call close_convection_gaps', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call close_convection_gaps(s, ierr)
            if (failed('close_convection_gaps')) return

            if (dbg) write(*,3) 'call close_thermohaline_gaps', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call close_thermohaline_gaps(s, ierr)
            if (failed('close_thermohaline_gaps')) return

            if (dbg) write(*,3) 'call remove_thermohaline_dropouts', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call remove_thermohaline_dropouts(s, ierr)
            if (failed('remove_thermohaline_dropouts')) return

            if (dbg) write(*,3) 'call close_semiconvection_gaps', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call close_semiconvection_gaps(s, ierr)
            if (failed('close_semiconvection_gaps')) return

            if (dbg) write(*,3) 'call remove_embedded_semiconvection', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call remove_embedded_semiconvection(s, ierr)
              if (failed('remove_embedded_semiconvection')) return

         end if
         
         call check('after get remove_mixing_glitches')

         if (dbg) write(*,3) 'call do_mix_envelope', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call do_mix_envelope(s)

         do k=1,s% nz
            eps_h(k) = s% eps_nuc_categories(ipp,k) + &
                       s% eps_nuc_categories(icno,k)
            eps_he(k) = s% eps_nuc_categories(i3alf,k)
            eps_z(k) = s% eps_nuc(k) - (eps_h(k) + eps_he(k))
         end do

         if (dbg) write(*,3) 'call set_mlt_cz_boundary_info', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call set_mlt_cz_boundary_info(s, ierr)
         if (failed('set_mlt_cz_boundary_info')) return

         if (dbg) write(*,3) 'call locate_convection_boundaries', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call locate_convection_boundaries( &
            s, nz, eps_h, eps_he, eps_z, s% mstar, &
            s% q, s% cdc, ierr)
         if (failed('locate_convection_boundaries')) return
         
         if (.not. rsp_or_eturb) then
            if (dbg) write(*,3) 'call add_predictive_mixing', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call add_predictive_mixing(s, ierr)
            if (failed('add_predictive_mixing')) return
         end if
         
         call check('after add_predictive_mixing')

         ! NB: re-call locate_convection_boundries to take into
         ! account changes from add_predictive_mixing

         if (dbg) write(*,3) 'call locate_convection_boundaries', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call locate_convection_boundaries( &
            s, nz, eps_h, eps_he, eps_z, s% mstar, &
            s% q, s% cdc, ierr)
         if (failed('locate_convection_boundaries')) return

         if (dbg) write(*,*) 'call locate_mixing_boundaries'
         ! need to call this before add_overshooting
         call locate_mixing_boundaries(s, eps_h, eps_he, eps_z, ierr)
         if (failed('locate_mixing_boundaries')) return

         if (.not. rsp_or_eturb) then
            if (dbg) write(*,3) 'call add_overshooting', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call add_overshooting(s, ierr)
            if (failed('add_overshooting')) return
         end if
         
         call check('after add_overshooting')

         if (dbg) write(*,3) 'call add_RTI_turbulence', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call add_RTI_turbulence(s, ierr)
         if (failed('add_RTI_turbulence')) return

         if (dbg) write(*,3) 'call s% other_after_set_mixing_info', &
            k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
         call s% other_after_set_mixing_info(s% id, ierr)
         if (failed('other_after_set_mixing_info')) return

         if (.not. skip_set_cz_bdy_mass) then
            if (dbg) write(*,3) 'call set_cz_bdy_mass', &
               k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)
            call set_cz_bdy_mass(s, ierr)
            if (failed('set_cz_bdy_mass')) return
         end if

         if (s% set_min_D_mix .and. s% ye(nz) >= s% min_center_Ye_for_min_D_mix) then
            do k=1,nz
               if (s% D_mix(k) >= s% min_D_mix) cycle
               if (s% m(k) > s% mass_upper_limit_for_min_D_mix*Msun) cycle
               if (s% m(k) < s% mass_lower_limit_for_min_D_mix*Msun) cycle
               s% D_mix(k) = s% min_D_mix
               s% mixing_type(k) = minimum_mixing
            end do
         end if

         if (s% set_min_D_mix_below_Tmax) then
            Tmax = -1
            k_Tmax = -1
            do k=1,nz
               if (s% T(k) > Tmax) then
                  Tmax = s% T(k)
                  k_Tmax = k
               end if
            end do
            do k=k_Tmax+1,nz
               if (s% D_mix(k) >= s% min_D_mix_below_Tmax) cycle
               s% D_mix(k) = s% min_D_mix_below_Tmax
               s% mixing_type(k) = minimum_mixing
            end do
         end if

         if (s% set_min_D_mix_in_H_He) then
            do k=1,nz
               if (s% X(k) + s% Y(k) < 0.5d0) exit
               if (s% D_mix(k) >= s% min_D_mix_in_H_He) cycle
               s% D_mix(k) = s% min_D_mix_in_H_He
               s% mixing_type(k) = minimum_mixing
            end do
         end if

         if (s% use_other_D_mix) then
            if (dbg) write(*,*) 'call other_D_mix'
            call s% other_D_mix(s% id, ierr)
            if (failed('other_D_mix')) return
         end if

         do k=1,nz
            s% D_mix_non_rotation(k) = s% D_mix(k)
         end do
         
         call check('before rotation_flag')

         if (s% rotation_flag) then

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do k=1,nz
                  write(*,3) 'before update_rotation_mixing_info D_mix', &
                     s% model_number, k, s% D_mix(k)
               end do
            end if

            if (dbg) write(*,*) 'call update_rotation_mixing_info'
            call update_rotation_mixing_info(s,ierr)
            if (failed('update_rotation_mixing_info')) return

            do k = 2, nz
               if (s% D_mix(k) < 1d-10) s% D_mix(k) = 0d0
               s% cdc(k) = s% D_mix(k)*cdc_factor(k)
               if (s% D_mix(k) /= 0 .and. s% mixing_type(k) == no_mixing) then
                  s% mixing_type(k) = rotation_mixing
               end if
            end do
            s% cdc(1) = s% cdc(2)

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do k=1,nz
                  write(*,3) 'after do rotation mixing D_mix', &
                     s% model_number, k, s% D_mix(k)
               end do
            end if

         end if
         
         call check('after update_rotation_mixing_info')
         
         if (.not. s% conv_vel_flag) then

            region_bottom_q = s% D_mix_zero_region_bottom_q
            region_top_q = s% D_mix_zero_region_top_q
            
            if (s% dq_D_mix_zero_at_H_He_crossover > 0d0) then
               i_h1 = s% net_iso(ih1)
               i_he4 = s% net_iso(ihe4)
               if (i_h1 > 0 .and. i_he4 > 0) then
                  do k=2,nz
                     if (s% xa(i_h1,k-1) > s% xa(i_he4,k-1) .and. &
                         s% xa(i_h1,k) <= s% xa(i_he4,k)) then ! crossover
                        region_bottom_q = &
                           s% q(k) - 0.5d0*s% dq_D_mix_zero_at_H_He_crossover
                        region_top_q = &
                           s% q(k) + 0.5d0*s% dq_D_mix_zero_at_H_He_crossover
                        exit
                     end if
                  end do
               end if
            end if
            
            if (region_bottom_q < region_top_q) then
               do k=2,nz
                  if (s% q(k) >= region_bottom_q .and. s% q(k) <= region_top_q) then
                     s% D_mix(k) = 0d0
                     s% mixing_type(k) = no_mixing
                  end if
               end do
            end if
            
            region_bottom_q = 1d99
            region_top_q = -1d99
            if (s% dq_D_mix_zero_at_H_C_crossover > 0d0) then
               i_h1 = s% net_iso(ih1)
               i_c12 = s% net_iso(ic12)
               if (i_h1 > 0 .and. i_c12 > 0) then
                  do k=2,nz
                     if (s% xa(i_h1,k-1) > s% xa(i_c12,k-1) .and. &
                         s% xa(i_h1,k) <= s% xa(i_c12,k)) then ! crossover
                        region_bottom_q = &
                           s% q(k) - 0.5d0*s% dq_D_mix_zero_at_H_C_crossover
                        region_top_q = &
                           s% q(k) + 0.5d0*s% dq_D_mix_zero_at_H_C_crossover
                        exit
                     end if
                  end do
               end if
            end if
            
            if (region_bottom_q < region_top_q) then
               do k=2,nz
                  if (s% q(k) >= region_bottom_q .and. s% q(k) <= region_top_q) then
                     s% D_mix(k) = 0d0
                     s% mixing_type(k) = no_mixing
                  end if
               end do
            end if
         
            ! as last thing, update conv_vel from D_mix and mixing length.
            do k=2,nz
               if (s% alpha_mlt(k)*s% scale_height(k) > 0) then
                  s% conv_vel(k) = &
                     3*s% D_mix(k)/(s% alpha_mlt(k)*s% scale_height(k))
               else
                  s% conv_vel(k) = 0
               end if
            end do
            
         end if

         ! set these just for plotting.  not used.
         s% mixing_type(1) = s% mixing_type(2)
         s% D_mix(1) = s% D_mix(2)
         if (.not. s% conv_vel_flag) then
            s% conv_vel(1) = 0d0
         end if

         call check('final')
         if (failed('set_mixing_info')) return

         if (dbg) write(*,3) 'done mixing', k_dbg, s% mixing_type(k_dbg), s% D_mix(k_dbg)

         call dealloc

         if (s% doing_timing) &
            call update_time(s, time0, total, s% time_set_mixing_info)

         s% have_mixing_info = .true.

         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            if (s% report_ierr .or. dbg) &
               write(*,*) 'set_mixing_info failed in call to ' // trim(str)
            failed = .true.
            call dealloc
         end function failed


         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               eps_h, nz, nz_alloc_extra, 'mix_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               eps_he, nz, nz_alloc_extra, 'mix_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               eps_z, nz, nz_alloc_extra, 'mix_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               cdc_factor, nz, nz_alloc_extra, 'mix_info', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

         subroutine check(str)
            character(len=*) :: str
            integer :: k
            include 'formats'
            do k = 1, s% nz
               if (is_bad_num(s% D_mix(k))) then
                  ierr = -1
                  if (s% report_ierr) then
                     write(*,3) trim(str) // ' mixing_type, D_mix', k, s% mixing_type(k), s% D_mix(k)
                     if (s% rotation_flag) then
                        if (is_bad_num(s% D_mix_non_rotation(k))) &
                           write(*,2) 's% D_mix_non_rotation(k)', k, s% D_mix_non_rotation(k)
                        if (is_bad_num(s% D_visc(k))) write(*,2) 's% D_visc(k)', k, s% D_visc(k)
                        if (is_bad_num(s% D_DSI(k))) write(*,2) 's% D_DSI(k)', k, s% D_DSI(k)
                        if (is_bad_num(s% D_SH(k))) write(*,2) 's% D_SH(k)', k, s% D_SH(k)
                        if (is_bad_num(s% D_SSI(k))) write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
                        if (is_bad_num(s% D_ES(k))) write(*,2) 's% D_ES(k)', k, s% D_ES(k)
                        if (is_bad_num(s% D_GSF(k))) write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
                        if (is_bad_num(s% D_ST(k))) write(*,2) 's% D_ST(k)', k, s% D_ST(k)
                     end if
                  end if
                  if (s% stop_for_bad_nums) stop 'set mixing info'
               end if
            end do
         end subroutine check

      end subroutine set_mixing_info


      subroutine set_mlt_cz_boundary_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, mt, mt1, mt2, nz
         real(dp) :: dg0, dg1

         include 'formats'

         ! NOTE: this routine must be called BEFORE overshooting is done.

         ! convective zone boundary is where gradL = gradr
         ! semi-convective zone boundary is where grada_face = gradr

         ierr = 0
         nz = s% nz
         s% cz_bdy_dq(1:nz) = 0d0
         
         if (s% rsp_flag .or. s% Eturb_flag) return ! don't have MLT info

         do k = 2, nz
            mt1 = s% mixing_type(k-1)
            mt2 = s% mixing_type(k)
            if (mt1 == mt2) cycle
            if (mt2 == convective_mixing .or. mt1 == convective_mixing) then
               dg0 = s% gradL(k-1) - s% gradr(k-1)
               dg1 = s% gradL(k) - s% gradr(k)
            else if (mt2 == semiconvective_mixing .or. mt1 == semiconvective_mixing) then
               dg0 = s% grada_face(k-1) - s% gradr(k-1)
               dg1 = s% grada_face(k) - s% gradr(k)
            else
               cycle
            end if
            if (dg0*dg1 >= 0) cycle
            s% cz_bdy_dq(k-1) = find0(0d0,dg0,s% dq(k-1),dg1)
            if (s% cz_bdy_dq(k-1) < 0d0 .or. s% cz_bdy_dq(k-1) > s% dq(k-1)) then
               write(*,2) 'bad cz_bdy_dq', k-1, s% cz_bdy_dq(k-1), s% dq(k-1)
               stop 'set_mlt_cz_boundary_info'
               ierr = -1
               return
            end if
         end do

      end subroutine set_mlt_cz_boundary_info


      subroutine set_cz_bdy_mass(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, j, nz
         logical, parameter :: dbg = .false.

         include 'formats'
         ierr = 0
         nz = s% nz

         s% cz_top_mass(:) = s% mstar
         s% cz_bot_mass(:) = s% mstar

         s% n_conv_regions = 0
         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         if (in_convective_region) then
            s% n_conv_regions = 1
            s% cz_bot_mass(1) = s% M_center
         end if

         if (dbg) write(*,*) 'initial in_convective_region', in_convective_region

         do k=nz-1, 2, -1
            if (in_convective_region) then
               if (s% mixing_type(k) /= convective_mixing) then ! top of convective region
                  s% cz_top_mass(s% n_conv_regions) = &
                     s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
                  in_convective_region = .false.
               end if
            else
               if (s% mixing_type(k) == convective_mixing) then ! bottom of convective region
                  if (s% n_conv_regions < max_num_mixing_regions) then
                     s% n_conv_regions = s% n_conv_regions + 1
                     s% cz_bot_mass(s% n_conv_regions) = &
                        s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
                  end if
                  in_convective_region = .true.
               end if
            end if
         end do
         if (in_convective_region) then
            s% cz_top_mass(s% n_conv_regions) = s% mstar
         end if
         
         s% have_new_cz_bdy_info = .true.

         if (dbg) then
            write(*,*)
            write(*,2) 'set_mixing_info s% n_conv_regions', s% n_conv_regions
            do j = 1, s% n_conv_regions
               write(*,2) 'conv region', j, s% cz_bot_mass(j)/Msun, s% cz_top_mass(j)/Msun
            end do
            write(*,*)
         end if

      end subroutine set_cz_bdy_mass


      subroutine locate_convection_boundaries( &
            s, nz, eps_h, eps_he, eps_z, mstar, q, cdc, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), dimension(:), intent(in) :: eps_h, eps_he, eps_z
         real(dp), intent(in) :: mstar
         real(dp), pointer, dimension(:) :: q, cdc
         integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, k_bot, i, j, iounit, max_conv_bdy
         real(dp) :: dgrad00, dgradp1, turnover_time, &
            bot_Hp, bot_r, top_Hp, top_r, dr

         logical :: dbg
         logical, parameter :: write_debug = .false.

         include 'formats'

         ierr = 0
         dbg = .false.

         if (write_debug) then
            open(newunit=iounit, file=trim('debug.data'), action='write', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'open debug.data failed'
               return
            end if
            write(*,*) 'write debug.data'
            write(iounit,*) 'nz', nz
            write(iounit,1) 'mstar', mstar
            do k=1,nz
               write(iounit,2) 'q', k, q(k)
               write(iounit,2) 'cdc', k, cdc(k)
               write(iounit,2) 'eps_h', k, eps_h(k)
               write(iounit,2) 'eps_he', k, eps_he(k)
               write(iounit,2) 'eps_z', k, eps_z(k)
            end do
         end if

         max_conv_bdy = size(s% conv_bdy_q, dim=1)
         s% conv_bdy_q(:) = 0
         s% conv_bdy_loc(:) = 0
         s% top_conv_bdy(:) = .false.
         s% burn_h_conv_region(:) = .false.
         s% burn_he_conv_region(:) = .false.
         s% burn_z_conv_region(:) = .false.
         bot_Hp = 0; bot_r = 0; top_Hp = 0; top_r = 0; dr = 0

         s% num_conv_boundaries = 0
         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         k_bot = nz
         turnover_time = 0

         do k=nz-1, 2, -1
            if (in_convective_region) then
               if (s% mixing_type(k) /= convective_mixing) then
                  call end_of_convective_region
               else
                  if(s% conv_vel(k).ne. 0d0) turnover_time = turnover_time + (s% rmid(k-1) - s% rmid(k))/s% conv_vel(k)
               end if
            else ! in non-convective region
               if (s% mixing_type(k) == convective_mixing) then ! start of a convective region
                  if (s% num_conv_boundaries == max_conv_bdy) then
                     call realloc(ierr)
                     if (ierr /= 0) then
                        return
                     end if
                  end if
                  s% num_conv_boundaries = s% num_conv_boundaries+1
                  i = s% num_conv_boundaries
                  k_bot = k+1
                  if (k == 1) then
                     s% conv_bdy_q(i) = 1
                  else ! bottom of region is between k+1 and k
                     s% conv_bdy_q(i) = s% q(k) - s% cz_bdy_dq(k)
                  end if
                  s% top_conv_bdy(i) = .false.
                  s% conv_bdy_loc(i) = k_bot
                     ! bottom of region is between k_bot and k_bot-1
                  in_convective_region = .true.
                  bot_r = s% r(k_bot)
                  !bot_Hp = scale_height(s,k_bot,.true.)
                  bot_Hp = s% scale_height(k_bot)
                  turnover_time = 0
               end if
            end if
         end do

         if (in_convective_region) then
            k = 1 ! end at top
            call end_of_convective_region
         end if

         if (write_debug) then
            write(iounit,*) 's% num_conv_boundaries', s% num_conv_boundaries
            do j=1,s% num_conv_boundaries
               write(iounit,2) 's% conv_bdy_q', j, s% conv_bdy_q(j)
               write(iounit,*) 's% top_conv_bdy', j, s% top_conv_bdy(j)
               write(iounit,*) 's% burn_h_conv_region', j, s% burn_h_conv_region(j)
               write(iounit,*) 's% burn_he_conv_region', j, s% burn_he_conv_region(j)
               write(iounit,*) 's% burn_z_conv_region', j, s% burn_z_conv_region(j)
               write(iounit,*) 's% conv_bdy_loc', j, s% conv_bdy_loc(j)
            end do
            close(iounit)
         end if

         if (dbg) then
            write(*,*) 's% num_conv_boundaries', s% num_conv_boundaries
            do j=1,s% num_conv_boundaries
               write(*,*) 's% conv_bdy_q', j, s% conv_bdy_q(j)
               write(*,*) 's% top_conv_bdy', j, s% top_conv_bdy(j)
               write(*,*) 's% burn_h_conv_region', j, s% burn_h_conv_region(j)
               write(*,*) 's% burn_he_conv_region', j, s% burn_he_conv_region(j)
               write(*,*) 's% burn_z_conv_region', j, s% burn_z_conv_region(j)
               write(*,*) 's% conv_bdy_loc', j, s% conv_bdy_loc(j)
               write(*,*) 'mixing type', s% mixing_type(s% conv_bdy_loc(j)-3:s% conv_bdy_loc(j)+3)
            end do
            write(*,*)
            write(*,3) 'mixing_type', 1152, s% mixing_type(1152)
            stop 'locate_convection_boundaries'
         end if


         contains


         subroutine end_of_convective_region()
            integer :: max_logT_loc, kk, op_err, mix_type
            real(dp) :: max_logT, max_X, max_Y, Hp, max_eps
            logical :: end_dbg

            9 format(a40, 3i7, 99(1pd26.16))
            include 'formats'

            in_convective_region = .false.

            end_dbg = .false.

            top_r = s% r(k)
            top_Hp = s% scale_height(k)
            dr = top_r - bot_r
            Hp = (bot_Hp + top_Hp)/2

            if (dr/Hp < s% prune_bad_cz_min_Hp_height .and. s% prune_bad_cz_min_Hp_height > 0) then
               max_eps = maxval(eps_h(k:k_bot) + eps_he(k:k_bot) + eps_z(k:k_bot))
               if (end_dbg) write(*,3) 'max_eps', k, k_bot, max_eps, &
                  exp10(s% prune_bad_cz_min_log_eps_nuc)
               if (max_eps < exp10(s% prune_bad_cz_min_log_eps_nuc) &
                     .and. all(s% mixing_type(k+1:k_bot-1) /= thermohaline_mixing)) then
                  do kk = k, k_bot ! this includes the radiative points at boundaries
                     call set_use_gradr(s,kk)
                     s% cdc(kk) = 0
                     s% D_mix(kk) = 0
                     if (.not. s% conv_vel_flag) s% conv_vel(kk) = 0
                     s% mixing_type(kk) = no_mixing
                  end do
                  if (s% num_conv_boundaries > 0) &
                     s% num_conv_boundaries = s% num_conv_boundaries-1
                  return
               end if
            end if

            if (s% num_conv_boundaries == max_conv_bdy) then
               call realloc(ierr)
               if (ierr /= 0) return
            end if

            s% num_conv_boundaries = s% num_conv_boundaries+1
            i = s% num_conv_boundaries

            ! check for burning in region
            max_logT = -1d99
            max_X = 0d0
            max_Y = 0d0
            max_logT_loc = 0
            do kk=k,min(nz,k_bot+1)
               if (s% lnT(kk)/ln10 > max_logT) then
                  max_logT = s% lnT(kk)/ln10
                  max_logT_loc = kk
               end if
               if (s% X(kk) > max_X) then
                  max_X = s% X(kk)
               end if
               if (s% Y(kk) > max_Y) then
                  max_Y = s% Y(kk)
               end if
            end do
            if (max_logT > s% burn_z_mix_region_logT &
                  .and. max_Y < s% max_Y_for_burn_z_mix_region) then            
               s% burn_z_conv_region(i) = .true.
               if (i > 1) s% burn_z_conv_region(i-1) = .true.
               !write(*,*) 'burn z mix region', i, &
               !   s% burn_z_mix_region_logT, max_logT, s% m(max_logT_loc)/Msun
            else if (max_logT > s% burn_he_mix_region_logT &
                  .and. max_X < s% max_X_for_burn_he_mix_region) then
               s% burn_he_conv_region(i) = .true.
               if (i > 1) s% burn_he_conv_region(i-1) = .true.
               !write(*,*) 'burn he mix region', i, &
               !   s% burn_he_mix_region_logT, max_logT, s% m(max_logT_loc)/Msun
            else if (max_logT > s% burn_h_mix_region_logT) then
               s% burn_h_conv_region(i) = .true.
               if (i > 1) s% burn_h_conv_region(i-1) = .true.
               !write(*,*) 'burn h mix region', i, &
               !   s% burn_h_mix_region_logT, max_logT, s% m(max_logT_loc)/Msun
            end if

            if (k == 1) then
               s% conv_bdy_q(i) = 1
            else
               ! top of region is between k+1 and k
               s% conv_bdy_q(i) = s% q(k) - s% cz_bdy_dq(k)
            end if
            s% top_conv_bdy(i) = .true.
            s% conv_bdy_loc(i) = k

         end subroutine end_of_convective_region


         subroutine realloc(ierr)
            use utils_lib
            integer, intent(out) :: ierr

            integer :: sz

            sz = size(s% conv_bdy_q, dim=1)

            ierr = 0
            max_conv_bdy = 2*(10+max_conv_bdy)

            call realloc_double(s% conv_bdy_q,max_conv_bdy,ierr)
            if (ierr /= 0) return

            call realloc_integer(s% conv_bdy_loc,max_conv_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% top_conv_bdy,max_conv_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_h_conv_region,max_conv_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_he_conv_region,max_conv_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_z_conv_region,max_conv_bdy,ierr)
            if (ierr /= 0) return

            s% conv_bdy_q(sz+1:max_conv_bdy) = 0
            s% conv_bdy_loc(sz+1:max_conv_bdy) = 0
            s% top_conv_bdy(sz+1:max_conv_bdy) = .false.
            s% burn_h_conv_region(sz+1:max_conv_bdy) = .false.
            s% burn_he_conv_region(sz+1:max_conv_bdy) = .false.
            s% burn_z_conv_region(sz+1:max_conv_bdy) = .false.

         end subroutine realloc

      end subroutine locate_convection_boundaries


      subroutine set_use_gradr(s,k)
         use mlt_info
         type (star_info), pointer :: s
         integer, intent(in) :: k
         call switch_to_no_mixing(s,k)
         call switch_to_radiative(s,k)
      end subroutine set_use_gradr


      subroutine locate_mixing_boundaries(s, eps_h, eps_he, eps_z, ierr)
         type (star_info), pointer :: s
         real(dp), dimension(:), intent(in) :: eps_h, eps_he, eps_z
         integer, intent(out) :: ierr

         logical :: in_mixing_region
         integer :: k, k_bot, i, j, iounit, max_mix_bdy, nz

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         nz = s% nz

         max_mix_bdy = size(s% mix_bdy_q, dim=1)
         s% mix_bdy_q(:) = 0
         s% mix_bdy_loc(:) = 0
         s% top_mix_bdy(:) = .false.
         s% burn_h_mix_region(:) = .false.
         s% burn_he_mix_region(:) = .false.
         s% burn_z_mix_region(:) = .false.

         s% num_mix_boundaries = 0
         s% num_mix_regions = 0
         in_mixing_region = (s% mixing_type(nz) /= no_mixing)
         k_bot = nz
         do k=nz-1, 2, -1
            if (in_mixing_region) then
               if (s% mixing_type(k) == no_mixing) call end_of_mixing_region
            else ! in non-mixing region
               if (s% mixing_type(k) /= no_mixing) then ! start of a mixing region
                  if (s% num_mix_boundaries == max_mix_bdy) then
                     call realloc(ierr)
                     if (ierr /= 0) return
                  end if
                  s% num_mix_boundaries = s% num_mix_boundaries+1
                  i = s% num_mix_boundaries
                  k_bot = k+1
                  if (k == 1) then
                     s% mix_bdy_q(i) = 1
                  else ! bottom of region is between k+1 and k
                     s% mix_bdy_q(i) = s% q(k) - s% cz_bdy_dq(k)
                  end if
                  s% top_mix_bdy(i) = .false.
                  s% mix_bdy_loc(i) = k_bot
                  in_mixing_region = .true.
               end if
            end if
         end do

         if (in_mixing_region) then
            k = 1 ! end at top
            call end_of_mixing_region
         end if

            
         !do i=1,s% num_conv_boundaries
         !   write(*,*) 'locate_mixing_boundaries region burn_h he z', i, &
         !      s% burn_h_mix_region(i), s% burn_he_conv_region(i), s% burn_z_conv_region(i)
         !end do

         contains


         subroutine end_of_mixing_region()
            integer :: kk, max_logT_loc
            real(dp) :: max_logT, max_X, max_Y

            9 format(a40, 3i7, 99(1pd26.16))
            include 'formats'

            in_mixing_region = .false.

            if (s% num_mix_boundaries == max_mix_bdy) then
               call realloc(ierr)
               if (ierr /= 0) return
            end if

            s% num_mix_regions = s% num_mix_regions+1
            s% num_mix_boundaries = s% num_mix_boundaries+1
            i = s% num_mix_boundaries

            ! check for burning in region
            max_logT = -1d99
            max_X = 0d0
            max_Y = 0d0
            max_logT_loc = 0
            do kk=k,min(nz,k_bot+1)
               if (s% lnT(kk)/ln10 > max_logT) then
                  max_logT = s% lnT(kk)/ln10
                  max_logT_loc = kk
               end if
               if (s% X(kk) > max_X) then
                  max_X = s% X(kk)
               end if
               if (s% Y(kk) > max_Y) then
                  max_Y = s% Y(kk)
               end if
            end do
            if (max_logT > s% burn_z_mix_region_logT &
                  .and. max_Y < s% max_Y_for_burn_z_mix_region) then            
               s% burn_z_mix_region(i) = .true.
               if (i > 1) s% burn_z_mix_region(i-1) = .true.
               !write(*,*) 'burn z mix region', i
            else if (max_logT > s% burn_he_mix_region_logT &
                  .and. max_X < s% max_X_for_burn_he_mix_region) then
               s% burn_he_mix_region(i) = .true.
               if (i > 1) s% burn_he_mix_region(i-1) = .true.
               !write(*,*) 'burn he mix region', i
            else if (max_logT > s% burn_h_mix_region_logT) then
               s% burn_h_mix_region(i) = .true.
               if (i > 1) s% burn_h_mix_region(i-1) = .true.
               !write(*,*) 'burn h mix region', i
            end if

            if (k == 1) then
               s% mix_bdy_q(i) = 1
            else
               ! top of region is between k+1 and k
               s% mix_bdy_q(i) = s% q(k) - s% cz_bdy_dq(k)
            end if
            s% top_mix_bdy(i) = .true.
            s% mix_bdy_loc(i) = k

         end subroutine end_of_mixing_region


         subroutine realloc(ierr)
            use utils_lib
            integer, intent(out) :: ierr

            integer :: sz

            sz = size(s% mix_bdy_q, dim=1)

            ierr = 0
            max_mix_bdy = 2*(10+max_mix_bdy)

            call realloc_double(s% mix_bdy_q,max_mix_bdy,ierr)
            if (ierr /= 0) return

            call realloc_integer(s% mix_bdy_loc,max_mix_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% top_mix_bdy,max_mix_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_h_mix_region,max_mix_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_he_mix_region,max_mix_bdy,ierr)
            if (ierr /= 0) return

            call realloc_logical(s% burn_z_mix_region,max_mix_bdy,ierr)
            if (ierr /= 0) return

            s% mix_bdy_q(sz+1:max_mix_bdy) = 0
            s% mix_bdy_loc(sz+1:max_mix_bdy) = 0
            s% top_mix_bdy(sz+1:max_mix_bdy) = .false.
            s% burn_h_mix_region(sz+1:max_mix_bdy) = .false.
            s% burn_he_mix_region(sz+1:max_mix_bdy) = .false.
            s% burn_z_mix_region(sz+1:max_mix_bdy) = .false.

         end subroutine realloc


      end subroutine locate_mixing_boundaries


      subroutine remove_tiny_mixing(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz
         logical, parameter :: dbg = .false.
         real(dp) :: tiny

         include 'formats'

         if (dbg) write(*,*) 'remove_tiny_mixing'

         ierr = 0
         nz = s% nz

         tiny = s% clip_D_limit
         do k=1,nz
            if (s% D_mix(k) < tiny) then
               s% cdc(k) = 0
               s% D_mix(k) = 0
               if (.not. s% conv_vel_flag) then
                  s% conv_vel(k) = 0
               end if
               s% mixing_type(k) = no_mixing
            end if
         end do

      end subroutine remove_tiny_mixing


      ! remove single point mixing or non-mixing regions
      ! NOTE: does not remove thermohaline singletons
      subroutine remove_mixing_singletons(s, ierr)
         !use star_utils, only: scale_height
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz
         logical, parameter :: dbg = .false.
         real(dp) :: lambda

         include 'formats'

         if (dbg) write(*,*) 'remove_mixing_singletons'

         ierr = 0
         nz = s% nz

         do k=2,nz-1
            if (s% cdc(k) == 0) then
               if (s% cdc(k-1) /= 0 .and. s% cdc(k+1) /= 0) then
                  s% cdc(k) = (s% cdc(k-1) + s% cdc(k+1))/2
                  s% D_mix(k) = s% cdc(k)/pow2(4*pi*s% r(k)*s% r(k)*s% rho(k))
                  lambda = s% alpha_mlt(k)* &
                     (s% scale_height(k-1) + s% scale_height(k+1))/2
                  if (.not. s% conv_vel_flag) then
                     s% conv_vel(k) = 3*s% D_mix(k)/lambda
                  end if
                  s% mixing_type(k) = max(s% mixing_type(k-1), s% mixing_type(k+1))
               if (dbg) write(*,3) 'remove radiative singleton', k, nz
               end if
            else if (s% okay_to_remove_mixing_singleton) then
               if (s% cdc(k-1) == 0 .and. s% cdc(k+1) == 0) then
                  call set_use_gradr(s,k)
                  s% cdc(k) = 0
                  s% D_mix(k) = 0
                  if (.not. s% conv_vel_flag) then
                     s% conv_vel(k) = 0
                  end if
                  s% mixing_type(k) = no_mixing
                  if (dbg) write(*,3) 'remove mixing singleton', k, nz
               end if
            end if
         end do

         if (s% cdc(1) == 0) then
            if (s% cdc(2) /= 0) then
               s% cdc(1) = s% cdc(2)
               s% D_mix(1) = s% D_mix(2)
               if (.not. s% conv_vel_flag) then
                  s% conv_vel(1) = s% conv_vel(2)
               end if
               s% mixing_type(1) = s% mixing_type(2)
               if (dbg) write(*,3) 'remove radiative singleton', 1, nz
            end if
         else
            if (s% cdc(2) == 0) then
               call set_use_gradr(s,1)
               s% cdc(1) = 0
               s% D_mix(1) = 0
               if (.not. s% conv_vel_flag) then
                  s% conv_vel(1) = 0
               end if
               s% mixing_type(1) = no_mixing
               if (dbg) write(*,2) 'remove mixing singleton', 1
            end if
         end if

         if (s% cdc(nz) == 0) then
            if (s% cdc(nz-1) /= 0) then
               s% cdc(nz) = s% cdc(nz-1)
               s% D_mix(nz) = s% D_mix(nz-1)
               if (.not. s% conv_vel_flag) s% conv_vel(nz) = s% conv_vel(nz-1)
               s% mixing_type(nz) = s% mixing_type(nz-1)
               if (dbg) write(*,2) 'remove radiative singleton: s% cdc(nz-1)', nz, s% cdc(nz-1)
            end if
         else
            if (s% cdc(nz-1) == 0) then
               call set_use_gradr(s,nz)
               s% cdc(nz) = 0
               s% D_mix(nz) = 0
               if(.not. s% conv_vel_flag) s% conv_vel(nz) = 0
               s% mixing_type(nz) = no_mixing
               if (dbg) write(*,2) 'remove mixing singleton: s% cdc(nz)', nz, s% cdc(nz)
            end if
         end if

      end subroutine remove_mixing_singletons


      subroutine close_convection_gaps(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call close_gaps(s, convective_mixing, s% min_convective_gap, ierr)
      end subroutine close_convection_gaps


      subroutine close_thermohaline_gaps(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call close_gaps(s, thermohaline_mixing, s% min_thermohaline_gap, ierr)
      end subroutine close_thermohaline_gaps


      subroutine close_semiconvection_gaps(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call close_gaps(s, semiconvective_mixing, s% min_semiconvection_gap, ierr)
      end subroutine close_semiconvection_gaps


      subroutine close_gaps(s, mix_type, min_gap, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: mix_type
         real(dp), intent(in) :: min_gap
         integer, intent(out) :: ierr

         integer :: k, kk, nz
         logical :: in_region, dbg
         real(dp) :: rtop, rbot, Hp
         integer :: ktop, kbot ! k's for gap
         include 'formats'

         dbg = .false.
         !dbg = (mix_type == convective_mixing)
         if (dbg) write(*,*) 'close_gaps convective_mixing'
         if (dbg) write(*,3) 'mixing_type', 1152, s% mixing_type(1152)
         ierr = 0
         if (min_gap < 0) return
         nz = s% nz
         in_region = (s% mixing_type(nz) == mix_type)
         rbot = 0
         kbot = nz
         do k=nz-1, 2, -1
            if (in_region) then
               if (s% mixing_type(k) /= mix_type) then ! end of region
                  kbot = k+1
                  rbot = s% r(kbot)
                  in_region = .false.
                  if (dbg) write(*,2) 'end of region', kbot, rbot
               end if
            else
               if (s% mixing_type(k) == mix_type) then ! start of region
                  ktop = k
                  rtop = s% r(ktop)
                  Hp = s% P(ktop)/(s% rho(ktop)*s% grav(ktop))
                  if (dbg) write(*,2) 'start of region', ktop, rtop
                  if (dbg) write(*,1) 'rtop - rbot < Hp*min_gap', (rtop - rbot) - Hp*min_gap, &
                     rtop - rbot, Hp*min_gap, Hp, min_gap, (rtop-rbot)/Hp
                  if (rtop - rbot < Hp*min_gap) then
                     if (kbot < nz) then
                        s% cdc(ktop+1:kbot-1) = (s% cdc(ktop) + s% cdc(kbot))/2
                        s% D_mix(ktop+1:kbot-1) = &
                           (s% D_mix(ktop) + s% D_mix(kbot))/2
                        if(.not. s% conv_vel_flag) s% conv_vel(ktop+1:kbot-1) = (s% conv_vel(ktop) + s% conv_vel(kbot))/2
                        s% mixing_type(ktop+1:kbot-1) = mix_type
                        if (dbg) write(*,3) 'close mixing gap', &
                              ktop+1, kbot-1, (rtop - rbot)/Hp, rtop - rbot, Hp
                     else
                        s% cdc(ktop+1:kbot) = s% cdc(ktop)
                        s% D_mix(ktop+1:kbot) = s% D_mix(ktop)
                        if(.not. s% conv_vel_flag) s% conv_vel(ktop+1:kbot) = s% conv_vel(ktop)
                        s% mixing_type(ktop+1:kbot) = mix_type
                        if (dbg) write(*,3) 'close mixing gap', &
                           ktop+1, kbot, (rtop - rbot)/Hp, rtop - rbot, Hp
                     end if
                  end if
                  in_region = .true.
               end if
            end if
         end do
         if (dbg) write(*,3) 'mixing_type', 1152, s% mixing_type(1152)
         if (dbg) write(*,*) 'done close_gaps'
         !if (dbg) stop 'done close_gaps'

      end subroutine close_gaps


      ! if find radiative region embedded in thermohaline,
      ! and max(gradL - grada) in region is < 1d-3
      ! and region height is < min_thermohaline_dropout
      ! then convert the region to thermohaline
      subroutine remove_thermohaline_dropouts(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz, j
         logical :: in_region
         real(dp) :: rtop, rbot, Hp, q_upper, q_lower, alfa, beta
         integer :: ktop, kbot ! k's for gap
         logical :: all_small
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         nz = s% nz
         rbot = s% r(nz)
         kbot = nz-1
         in_region = (s% mixing_type(kbot) == no_mixing)
         all_small = .false.
         do k=nz-2, 2, -1
            if (in_region) then
               if (s% mixing_type(k) == no_mixing) then ! check if okay
                  if (s% gradL(k) - s% grada_face(k) > s% max_dropout_gradL_sub_grada) &
                     all_small = .false.
               else ! end of radiative region
                  ktop = k+1
                  rtop = s% r(ktop)
                  Hp = s% P(ktop)/(s% rho(ktop)*s% grav(ktop))
                  q_upper = s% q(ktop-1)
                  q_lower = s% q(kbot+1)
                  if (rtop - rbot < Hp*s% min_thermohaline_dropout .and. &
                      s% mixing_type(ktop-1) == thermohaline_mixing .and. &
                      s% mixing_type(kbot+1) == thermohaline_mixing .and. &
                      q_upper - q_lower > 1d-20 .and. all_small) then
                     do j = ktop, kbot ! interpolate in q
                        alfa = (s% q(j) - q_lower)/(q_upper - q_lower)
                        beta = 1 - alfa
                        s% cdc(j) = alfa*s% cdc(ktop-1) + beta*s% cdc(kbot+1)
                        s% D_mix(j) = alfa*s% D_mix(ktop-1) + beta*s% D_mix(kbot+1)
                        if (.not. s% conv_vel_flag) then
                           s% conv_vel(j) = alfa*s% conv_vel(ktop-1) + beta*s% conv_vel(kbot+1)
                        end if
                        s% mixing_type(j) = thermohaline_mixing
                     end do
                  end if
                  in_region = .false.
               end if
            else
               if (s% mixing_type(k) == no_mixing) then ! start of region
                  kbot = k
                  rbot = s% r(kbot)
                  in_region = .true.
                  all_small = &
                     (s% gradL(k) - s% grada_face(k) <= s% max_dropout_gradL_sub_grada)
               end if
            end if
         end do

      end subroutine remove_thermohaline_dropouts


      subroutine remove_embedded_semiconvection(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz
         logical :: in_region
         integer :: kbot, ktop

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         if (.not. s% remove_embedded_semiconvection) return

         nz = s% nz

         in_region = check(nz)
         kbot = nz
         do k=nz-1, 2, -1
            if (in_region) then
               if (.not. check(k)) then ! end of region
                  ktop = k+1
                  in_region = .false.
                  call clean_region
               end if
            else ! not in region
               if (check(k)) then ! start of region
                  kbot = k
                  in_region = .true.
               end if
            end if
         end do

         if (in_region) then
            ktop = 1
            call clean_region
         end if


         contains


         subroutine clean_region
            integer :: kbot1, ktop1, kk
            include 'formats'
            if (dbg) write(*,3) 'clean_region semiconvective', kbot, ktop
            ! move top to below top convective region
            do while (s% mixing_type(ktop) == convective_mixing)
               ktop = ktop + 1
               if (ktop >= kbot) return
            end do
            if (dbg) write(*,2) 'new ktop 1', ktop
            ! move top to below top semiconvective region
            do while (s% mixing_type(ktop) == semiconvective_mixing)
               ktop = ktop + 1
               if (ktop >= kbot) return
            end do
            if (dbg) write(*,2) 'new ktop 2', ktop
            ! move bot to above bottom convective region
            do while (s% mixing_type(kbot) == convective_mixing)
               kbot = kbot - 1
               if (ktop >= kbot) return
            end do
            if (dbg) write(*,2) 'new kbot 1', kbot
            ! move bot to above bottom semiconvective region
            do while (s% mixing_type(kbot) == semiconvective_mixing)
               kbot = kbot - 1
               if (ktop >= kbot) return
            end do
            if (dbg) write(*,2) 'new kbot 2', kbot
            ! convert any semiconvective region between kbot and ktop
            kbot1 = kbot
            do while (kbot1 > ktop)
               ! move kbot1 to bottom of next semiconvective region
               do while (s% mixing_type(kbot1) == convective_mixing)
                  kbot1 = kbot1 - 1
                  if (kbot1 <= ktop) return
               end do
               ktop1 = kbot1
               ! move ktop1 to top of semiconvective region
               do while (s% mixing_type(ktop1) == semiconvective_mixing)
                  ktop1 = ktop1 - 1
                  if (ktop1 <= ktop) return
               end do
               s% D_mix(ktop1+1:kbot1) = s% D_mix(ktop1)
               s% mixing_type(ktop1+1:kbot1) = convective_mixing
               if (.not. s% conv_vel_flag) s% conv_vel(ktop1+1:kbot1) = s% conv_vel(ktop1)
               if (dbg) write(*,3) 'merge semiconvective island', kbot1, ktop1+1
               kbot1 = ktop1
            end do
         end subroutine clean_region


         logical function check(k)
            integer, intent(in) :: k
            check = &
               (s% mixing_type(k) == semiconvective_mixing) .or. &
               (s% mixing_type(k) == convective_mixing)
         end function check

      end subroutine remove_embedded_semiconvection


      subroutine do_mix_envelope(s)
         type (star_info), pointer :: s
         real(dp) :: T_mix_limit
         integer :: j, k, i, nz

         include 'formats'

         T_mix_limit = s% T_mix_limit
         !write(*,1) 'T_mix_limit', T_mix_limit
         if (T_mix_limit <= 0) return
         nz = s% nz
         j = 0
         do k = 1, nz ! search inward until find T >= T_mix_limit
            if (s% T(k) >= T_mix_limit) then
               j = k; exit
            end if
         end do
         if (j==0) j=nz ! all T < T_mix_limit
         ! find base of innermost convection that has T < T_mix_limit
         i = 0
         do k = j, 1, -1
            if (s% mixing_type(k) == convective_mixing) then
               i = k; exit
            end if
         end do
         if (i == 0) then
            return ! no convection in region with T < T_mix_limit
         end if
         ! extend convection region to surface
         j = maxloc(s% cdc(1:i), dim=1)
         if (.not. s% conv_vel_flag) s% conv_vel(1:i) = s% conv_vel(j)
         s% cdc(1:i) = s% cdc(j)
         do k = 1,i
            s% D_mix(k) = s% D_mix(j)
         end do
         s% mixing_type(1:i) = convective_mixing

      end subroutine do_mix_envelope


      subroutine get_convection_sigmas(s, dt, ierr)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer :: nz, k, j, species, ktop, kbot, bdy
         real(dp) :: sig_term_limit ! sig_term_limit is used to help convergence

         real(dp) :: siglim, xmstar, dq00, dqm1, cdcterm, dmavg, rho_face, &
            cdc, max_sig, D, xm1, x00, xp1, dm, dX, X, cushion, limit, &
            Tc, full_off, full_on, X_lim, dX_lim, qbot, qtop, &
            f1, f, df_dlnd00, df_dlndm1, df_dlnR, df_d_rho_face, alfa, beta
         logical :: in_convective_region
         real(dp), dimension(:), pointer :: sig, D_mix

         include 'formats'

         sig_term_limit = s% sig_term_limit

         ierr = 0
         nz = s% nz
         xmstar = s% xmstar
         sig => s% sig
         D_mix => s% D_mix

         sig(1) = 0d0
         do k = 2, nz
            D = D_mix(k)
            if (D == 0) then
               sig(k) = 0d0
               cycle
            end if
            rho_face = (s% dq(k-1)*s% rho(k) + s% dq(k)*s% rho(k-1))/ &
                        (s% dq(k-1) + s% dq(k))
            f1 = 4*pi*s% r(k)*s% r(k)*rho_face
            f = f1*f1
            cdc = D*f
            cdcterm = s% mix_factor*cdc
            dq00 = s% dq(k)
            dqm1 = s% dq(k-1)
            dmavg = xmstar*(dq00+dqm1)/2
            sig(k) = cdcterm/dmavg
            if (is_bad_num(sig(k))) then
               if (s% stop_for_bad_nums) then
                  write(*,2) 'sig(k)', k, sig(k)
                  stop 'get_convection_sigmas'
               end if
               sig(k) = 1d99
            end if
         end do

         ! can get numerical problems unless limit sig
         max_sig = maxval(sig(1:nz))
         if (max_sig < 1) return
         do k = 1, nz
            s% sig_raw(k) = sig(k)
            if (sig(k) == 0) cycle
            if (k > 1) then
               siglim = sig_term_limit*xmstar*min(s% dq(k),s% dq(k-1))/dt
            else
               siglim = sig_term_limit*xmstar*s% dq(k)/dt
            end if
            if (sig(k) > siglim) sig(k) = siglim
         end do

         species = s% species
         ! limit sigma to avoid negative mass fractions
         cushion = 10d0
         Tc = s% T(s% nz)
         full_off = s% Tcenter_max_for_sig_min_factor_full_off
         if (Tc <= full_off) return
         full_on = s% Tcenter_min_for_sig_min_factor_full_on
         limit = s% sig_min_factor_for_high_Tcenter
         if (Tc < full_on) then
            alfa = (Tc - full_off)/(full_on - full_off)
            ! limit is full on for alfa = 1 and limit is 1 for alfa = 0
            limit = limit*alfa + (1d0 - alfa)
         end if
         if (limit >= 1d0 .or. s% num_mix_regions == 0) return
         ! boundaries are in order from center to surface
         ! no bottom boundary at loc=nz included if center is mixed
         ! however, do include top boundary at loc=1 if surface is mixed
         if (s% top_mix_bdy(1)) then
            kbot = s% nz
            qbot = 0d0
            bdy = 1
         else
            kbot = s% mix_bdy_loc(1)
            qbot = s% q(kbot)
            bdy = 2
         end if
         if (.not. s% top_mix_bdy(bdy)) then
            stop 'bad mix bdy info 1'
         end if
         ktop = s% mix_bdy_loc(bdy)
         qtop = s% q(ktop)
         call do1_region
         do while (bdy < s% num_mix_boundaries)
            bdy = bdy+1
            if (s% top_mix_bdy(bdy)) then
               stop 'bad mix bdy info 2'
            end if
            kbot = s% mix_bdy_loc(bdy)
            qbot = s% q(kbot)
            bdy = bdy+1
            if (.not. s% top_mix_bdy(bdy)) then
               stop 'bad mix bdy info 3'
            end if
            ktop = s% mix_bdy_loc(bdy)
            qtop = s% q(ktop)
            call do1_region
         end do


         contains


         subroutine do1_region
            real(dp) :: delta_m, max_lim, alfa, beta, delta_m_upper, delta_m_lower
            delta_m = s% m(ktop) - s% m(kbot)
            delta_m_upper = s% delta_m_upper_for_sig_min_factor
            delta_m_lower = s% delta_m_lower_for_sig_min_factor
            if (delta_m >= delta_m_upper) then
               max_lim = 1d0
            else if (delta_m <= delta_m_lower) then
               max_lim = limit
            else
               alfa = (delta_m - delta_m_lower)/(delta_m_upper - delta_m_lower)
               beta = 1d0 - alfa
               max_lim = alfa + beta*limit
            end if
            do k=max(2,ktop),kbot
               call do1(k, max_lim, &
                  min(s% q(k) - qbot, qtop - s% q(k))*s% xmstar/Msun)
            end do
         end subroutine do1_region


         subroutine do1(k, max_lim, delta_m_to_bdy)
            integer, intent(in) :: k
            real(dp), intent(in) :: max_lim, delta_m_to_bdy
            real(dp) :: lim, max_delta_m_to_bdy
            include 'formats'
            siglim = sig(k)
            if (siglim == 0d0) return
            ! okay to increase limit up to max_lim
            max_delta_m_to_bdy = s% max_delta_m_to_bdy_for_sig_min_factor
            if (delta_m_to_bdy >= max_delta_m_to_bdy) return ! no change in sig
            lim = limit + (max_lim - limit)*delta_m_to_bdy/max_delta_m_to_bdy
            if (lim >= 1d0) return
            do j=1,species
               xm1 = s% xa(j,k-1)
               x00 = s% xa(j,k)
               if (xm1 > x00) then
                  X = xm1
                  dX = xm1 - x00
                  dm = s% dm(k-1)
               else
                  X = x00
                  dX = x00 - xm1
                  dm = s% dm(k)
                  dX = -dX
               end if
               if (X < 1d-5) cycle
               if (cushion*dt*dX*siglim > dm*X) then
                  siglim = dm*X/(dt*dX*cushion)
               end if
            end do
            if (siglim > lim*sig(k)) then
               sig(k) = siglim
            else
               sig(k) = lim*sig(k)
            end if
         end subroutine do1


      end subroutine get_convection_sigmas


      subroutine update_rotation_mixing_info(s, ierr)
         use rotation_mix_info, only: set_rotation_mixing_info
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz, k0
         logical :: set_min_am_nu_non_rot, okay
         real(dp) :: &
            am_nu_visc_factor, &
            am_nu_DSI_factor, &
            am_nu_SH_factor, &
            am_nu_SSI_factor, &
            am_nu_ES_factor, &
            am_nu_GSF_factor, &
            am_nu_ST_factor, &
            f, lgT, full_off, full_on
         real(dp), dimension(:), pointer :: & ! work vectors for tridiagonal solve
            sig, rhs, d, du, dl, bp, vp, xp, x

         include 'formats'

         ierr = 0
         nz = s% nz

         call set_rotation_mixing_info(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'update_rotation_mixing_info failed in call to set_rotation_mixing_info'
            return
         end if
         
         call check('after set_rotation_mixing_info')
         if (s% D_omega_flag) call check_D_omega('check_D_omega after set_rotation_mixing_info')

         ! include rotation part for mixing abundances
         full_on = s% D_mix_rotation_max_logT_full_on
         full_off = s% D_mix_rotation_min_logT_full_off
         do k = 2, nz
            lgT = s% lnT(k)/ln10
            if (lgT <= full_on) then
               f = 1d0
            else if (lgT >= full_off) then
               f = 0d0
            else ! lgT > full_on and < full_off
               f = (lgT - full_on) / (full_off - full_on)
            end if
            if (s% D_omega_flag) then
               s% D_mix_rotation(k) = f*s% am_D_mix_factor*s% D_omega(k)
            else
               s% D_mix_rotation(k) = &
                  f*s% am_D_mix_factor *(  &
                     ! note: have dropped viscosity from mixing
                     s% D_DSI_factor * s% D_DSI(k)  + &
                     s% D_SH_factor  * s% D_SH(k)   + &
                     s% D_SSI_factor * s% D_SSI(k)  + &
                     s% D_ES_factor  * s% D_ES(k)   + &
                     s% D_GSF_factor * s% D_GSF(k)  + &
                     s% D_ST_factor  * s% D_ST(k))
            end if
            s% D_mix(k) = s% D_mix_non_rotation(k) + s% D_mix_rotation(k)
         end do
         
         call check('after include rotation part for mixing abundances')

         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do k=2,nz
               write(*,2) 's% D_visc(k)', k, s% D_visc(k)
               write(*,2) 's% D_DSI(k)', k, s% D_DSI(k)
               write(*,2) 's% D_SH(k)', k, s% D_SH(k)
               write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
               write(*,2) 's% D_ES(k)', k, s% D_ES(k)
               write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
               write(*,2) 's% D_ST(k)', k, s% D_ST(k)
            end do
         end if

         if (s% model_number == -1) then
            k = 3
            write(*,2) 's% D_visc(k)', k, s% D_visc(k)
            write(*,2) 's% D_DSI(k)', k, s% D_DSI(k)
            write(*,2) 's% D_SH(k)', k, s% D_SH(k)
            write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
            write(*,2) 's% D_ES(k)', k, s% D_ES(k)
            write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
            write(*,2) 's% D_ST(k)', k, s% D_ST(k)
            write(*,2) 's% D_mix_non_rotation(k)', k, s% D_mix_non_rotation(k)
            write(*,2) 's% D_mix(k)', k, s% D_mix(k)
         end if

         am_nu_DSI_factor = s% am_nu_DSI_factor
         am_nu_SH_factor = s% am_nu_SH_factor
         am_nu_SSI_factor = s% am_nu_SSI_factor
         am_nu_ES_factor = s% am_nu_ES_factor
         am_nu_GSF_factor = s% am_nu_GSF_factor
         am_nu_ST_factor = s% am_nu_ST_factor
         am_nu_visc_factor = s% am_nu_visc_factor
         
         if ((.not. s% am_nu_rot_flag) .and. &
               (s% D_omega_flag .and. .not. s% job% use_D_omega_for_am_nu_rot)) then
            ! check for any am_nu factors > 0 and not same as for D_omega
            okay = .true.
            if (am_nu_DSI_factor >= 0 .and. am_nu_DSI_factor /= s% D_DSI_factor) okay = .false.
            if (am_nu_SH_factor >= 0 .and. am_nu_SH_factor /= s% D_SH_factor) okay = .false.
            if (am_nu_SSI_factor >= 0 .and. am_nu_SSI_factor /= s% D_SSI_factor) okay = .false.
            if (am_nu_DSI_factor >= 0 .and. am_nu_DSI_factor /= s% D_DSI_factor) okay = .false.
            if (am_nu_ES_factor >= 0 .and. am_nu_ES_factor /= s% D_ES_factor) okay = .false.
            if (am_nu_GSF_factor >= 0 .and. am_nu_GSF_factor /= s% D_GSF_factor) okay = .false.
            if (am_nu_ST_factor >= 0 .and. am_nu_ST_factor /= s% D_ST_factor) okay = .false.
            if (.not. okay) then
               write(*,*) 'have an am_nu factor >= 0 and not same as corresponding D_omega factor'
               write(*,*) 'so if want smoothing like D_omega, must set am_nu_rot_flag true'
               write(*,*) 'please fix this'
               ierr = -1
               return
            end if
         end if
         
         ! If am_nu_..._factor < -1, use the D_..._factor
         if (am_nu_DSI_factor < 0) am_nu_DSI_factor = s% D_DSI_factor
         if (am_nu_SH_factor < 0) am_nu_SH_factor = s% D_SH_factor
         if (am_nu_SSI_factor < 0) am_nu_SSI_factor = s% D_SSI_factor
         if (am_nu_ES_factor < 0) am_nu_ES_factor = s% D_ES_factor
         if (am_nu_GSF_factor < 0) am_nu_GSF_factor = s% D_GSF_factor
         if (am_nu_ST_factor < 0) am_nu_ST_factor = s% D_ST_factor
         if (am_nu_visc_factor < 0) am_nu_visc_factor = s% D_visc_factor

         ! set set_min_am_nu_non_rot to s% set_min_am_nu_non_rot if
         ! ye > min_center_ye_for_min_am_nu_non_rot
         ! and s% set_min_am_nu_non_rot
         set_min_am_nu_non_rot = &
            s% set_min_am_nu_non_rot .and. &
            s% ye(nz) >= s% min_center_Ye_for_min_am_nu_non_rot .and. &
            (.not. s% set_uniform_am_nu_non_rot)

         if (s% am_nu_rot_flag) then
            call set_am_nu_rot(ierr)
            if (ierr /= 0) return
         end if

         do k=1,nz
            if (s% set_uniform_am_nu_non_rot) then
               s% am_nu_non_rot(k) = s% uniform_am_nu_non_rot
            else
               s% am_nu_non_rot(k) = s% am_nu_factor* &
                  s% am_nu_non_rotation_factor*s% D_mix_non_rotation(k)
            end if
            if (set_min_am_nu_non_rot) &
               s% am_nu_non_rot(k) = max(s% min_am_nu_non_rot, s% am_nu_non_rot(k))
            if (s% am_nu_rot_flag) then
               ! already have am_nu_rot(k) from calling set_am_nu_rot above
            else if (s% D_omega_flag .and. s% job% use_D_omega_for_am_nu_rot) then
               s% am_nu_rot(k) = s% am_nu_factor*s% D_omega(k)
            else
               s% am_nu_rot(k) = s% am_nu_factor * ( &
                  am_nu_visc_factor*s% D_visc(k) + &
                  am_nu_DSI_factor*s% D_DSI(k) + &
                  am_nu_SH_factor*s% D_SH(k) + &
                  am_nu_SSI_factor*s% D_SSI(k) + &
                  am_nu_ES_factor*s% D_ES(k) + &
                  am_nu_GSF_factor*s% D_GSF(k) + &
                  am_nu_ST_factor*s% nu_ST(k))
            end if
            if (s% am_nu_rot(k) < 0d0) s% am_nu_rot(k) = 0d0
            s% am_nu_omega(k) = &
               s% am_nu_omega_non_rot_factor*s% am_nu_non_rot(k) + &
               s% am_nu_omega_rot_factor*s% am_nu_rot(k)
            if (s% am_nu_omega(k) < 0) then
               ierr = -1
               if (s% report_ierr) write(*,3) &
                  'update_rotation_mixing_info: am_nu_omega(k) < 0', k, nz, s% am_nu_omega(k), &
                  s% am_nu_non_rot(k), s% am_nu_rot(k), s% D_omega(k)
               return
            end if
            s% am_nu_j(k) = &
               s% am_nu_j_non_rot_factor*s% am_nu_non_rot(k) + &
               s% am_nu_j_rot_factor*s% am_nu_rot(k)
            if (s% am_nu_j(k) < 0) then
               ierr = -1
               if (s% report_ierr) write(*,3) &
                  'update_rotation_mixing_info: am_nu_j(k) < 0', k, nz, s% am_nu_j(k), &
                  s% am_nu_non_rot(k), s% am_nu_rot(k), s% D_omega(k)
               return
            end if
         end do

         if (s% use_other_am_mixing) then
            call s% other_am_mixing(s% id, ierr)
            if (ierr /= 0) return
         end if
         
         contains

         subroutine check(str)
            character(len=*) :: str
            integer :: k
            include 'formats'
            do k = 2, s% nz
               if (is_bad_num(s% D_mix(k))) then
                  ierr = -1
                  if (s% report_ierr) then
                     write(*,3) trim(str) // ' mixing_type, D_mix', k, s% mixing_type(k), s% D_mix(k)
                     if (s% rotation_flag) then
                        if (is_bad_num(s% D_mix_non_rotation(k))) &
                           write(*,2) 's% D_mix_non_rotation(k)', k, s% D_mix_non_rotation(k)
                        if (is_bad_num(s% D_visc(k))) write(*,2) 's% D_visc(k)', k, s% D_visc(k)
                        if (is_bad_num(s% D_DSI(k))) write(*,2) 's% D_DSI(k)', k, s% D_DSI(k)
                        if (is_bad_num(s% D_SH(k))) write(*,2) 's% D_SH(k)', k, s% D_SH(k)
                        if (is_bad_num(s% D_SSI(k))) write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
                        if (is_bad_num(s% D_ES(k))) write(*,2) 's% D_ES(k)', k, s% D_ES(k)
                        if (is_bad_num(s% D_GSF(k))) write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
                        if (is_bad_num(s% D_ST(k))) write(*,2) 's% D_ST(k)', k, s% D_ST(k)
                     end if
                  end if
                  if (s% stop_for_bad_nums) stop 'set mixing info'
               end if
            end do
         end subroutine check

         subroutine check_D_omega(str)
            character(len=*) :: str
            integer :: k
            include 'formats'
            do k = 2, s% nz
               if (is_bad_num(s% D_omega(k))) then
                  ierr = -1
                  if (s% report_ierr) then
                     write(*,3) trim(str) // ' mixing_type, D_omega', k, s% mixing_type(k), s% D_omega(k)
                     write(*,*) 's% doing_finish_load_model', s% doing_finish_load_model
                     if (s% rotation_flag) then
                        if (is_bad_num(s% D_mix_non_rotation(k))) &
                           write(*,2) 's% D_mix_non_rotation(k)', k, s% D_mix_non_rotation(k)
                        if (is_bad_num(s% D_visc(k))) write(*,2) 's% D_visc(k)', k, s% D_visc(k)
                        if (is_bad_num(s% D_DSI(k))) write(*,2) 's% D_DSI(k)', k, s% D_DSI(k)
                        if (is_bad_num(s% D_SH(k))) write(*,2) 's% D_SH(k)', k, s% D_SH(k)
                        if (is_bad_num(s% D_SSI(k))) write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
                        if (is_bad_num(s% D_ES(k))) write(*,2) 's% D_ES(k)', k, s% D_ES(k)
                        if (is_bad_num(s% D_GSF(k))) write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
                        if (is_bad_num(s% D_ST(k))) write(*,2) 's% D_ST(k)', k, s% D_ST(k)
                     end if
                  end if
                  if (s% stop_for_bad_nums) stop 'set mixing info'
               end if
            end do
         end subroutine check_D_omega
         
         subroutine set_am_nu_rot(ierr)
            use alloc
            use rotation_mix_info, only: smooth_for_rotation
            integer, intent(out) :: ierr
            integer :: i, k, nz
            real(dp) :: &
               dt, rate, d_ddt_dm1, d_ddt_d00, d_ddt_dp1, m, &
               d_dt, d_dt_in, d_dt_out, am_nu_rot_source
            include 'formats'
         
            ierr = 0
            nz = s% nz
            dt = s% dt
         
            if (s% am_nu_rot_flag .and. s% doing_finish_load_model) then
               do k=1,nz
                  s% am_nu_rot(k) = 0d0
               end do
            else if (s% am_nu_rot_flag) then
                     
               do k=1,nz
                  if (s% q(k) <= s% max_q_for_nu_omega_zero_in_convection_region .and. &
                      s% mixing_type(k) == convective_mixing) then
                     s% am_nu_rot(k) = 0d0
                     cycle
                  end if
                  am_nu_rot_source = s% am_nu_factor * ( &
                     am_nu_visc_factor*s% D_visc(k) + &
                     am_nu_DSI_factor*s% D_DSI(k) + &
                     am_nu_SH_factor*s% D_SH(k) + &
                     am_nu_SSI_factor*s% D_SSI(k) + &
                     am_nu_ES_factor*s% D_ES(k) + &
                     am_nu_GSF_factor*s% D_GSF(k) + &
                     am_nu_ST_factor*s% nu_ST(k))
                  if (is_bad(am_nu_rot_source)) then
                     write(*,2) 'am_nu_rot_source', k, am_nu_rot_source
                     stop 'set am_nu_rot'
                  end if
                  s% am_nu_rot(k) = am_nu_rot_source
                  if (is_bad(s% am_nu_rot(k))) then
                     write(*,2) 's% am_nu_rot(k)', k, s% am_nu_rot(k)
                     write(*,2) 'am_nu_rot_source', k, am_nu_rot_source
                     stop 'set am_nu_rot'
                  end if
               end do
               
               if (s% smooth_am_nu_rot > 0 .or. &
                    (s% nu_omega_mixing_rate > 0d0 .and. s% dt > 0)) then
                  
                  call do_alloc(ierr)
                  if (ierr /= 0) return

                  if (s% smooth_am_nu_rot > 0) then
                     call smooth_for_rotation(s, s% am_nu_rot, s% smooth_am_nu_rot, sig)
                  end if
            
                  if (s% nu_omega_mixing_rate > 0d0 .and. s% dt > 0) then ! mix am_nu_rot
            
                     rate = min(s% nu_omega_mixing_rate, 1d0/dt)
                     do k=2,nz-1
                        if (s% am_nu_rot(k) == 0 .or. s% am_nu_rot(k+1) == 0) then
                           sig(k) = 0
                        else if ((.not. s% nu_omega_mixing_across_convection_boundary) .and. &
                           s% mixing_type(k) /= convective_mixing .and. &
                               (s% mixing_type(k-1) == convective_mixing .or. &
                                s% mixing_type(k+1) == convective_mixing)) then
                            sig(k) = 0
                        else
                           sig(k) = rate*dt
                        end if       
                     end do
                     sig(1) = 0
                     sig(nz) = 0
            
                     do k=1,nz
                        if (k < nz) then
                           d_dt_in = sig(k)*(s% am_nu_rot(k+1) - s% am_nu_rot(k))
                        else
                           d_dt_in = -sig(k)*s% am_nu_rot(k)
                        end if
                        if (k > 1) then
                           d_dt_out = sig(k-1)*(s% am_nu_rot(k) - s% am_nu_rot(k-1))
                           d_ddt_dm1 = sig(k-1)
                           d_ddt_d00 = -(sig(k-1) + sig(k))
                        else
                           d_dt_out = 0
                           d_ddt_dm1 = 0
                           d_ddt_d00 = -sig(k)
                        end if
                        d_dt = d_dt_in - d_dt_out
                        d_ddt_dp1 = sig(k)
                        rhs(k) = d_dt
                        d(k) = 1d0 - d_ddt_d00
                        if (k < nz) then
                           du(k) = -d_ddt_dp1
                        else
                           du(k) = 0
                        end if
                        if (k > 1) dl(k-1) = -d_ddt_dm1               
                     end do
                     dl(nz) = 0
            
                     ! solve tridiagonal
                     bp(1) = d(1)
                     vp(1) = rhs(1)
                     do i = 2,nz
                        if (bp(i-1) == 0) then
                           write(*,*) 'failed in set_am_nu_rot', s% model_number
                           stop 'mix_am_nu_rot'
                           ierr = -1
                           return
                        end if
                        m = dl(i-1)/bp(i-1)
                        bp(i) = d(i) - m*du(i-1)
                        vp(i) = rhs(i) - m*vp(i-1)
                     end do
                     xp(nz) = vp(nz)/bp(nz)
                     x(nz) = xp(nz)
                     do i = nz-1, 1, -1
                        xp(i) = (vp(i) - du(i)*xp(i+1))/bp(i)
                        x(i) = xp(i)
                     end do
            
                     do k=2,nz
                        if (is_bad(x(k))) then
                           return
                           write(*,3) 's% am_nu_rot(k) prev, x', k, &
                              s% model_number, s% am_nu_rot(k), x(k), bp(i)
                           stop 'mix_am_nu_rot'
                        end if
                     end do
            
                     ! update am_nu_rot
                     do k=2,nz
                        s% am_nu_rot(k) = s% am_nu_rot(k) + x(k)
                        if (is_bad(s% am_nu_rot(k))) then
                           write(*,3) 's% am_nu_rot(k)', k, s% model_number, s% am_nu_rot(k)
                           stop 'mix_am_nu_rot'
                        end if
                        if (s% am_nu_rot(k) < 0d0) s% am_nu_rot(k) = 0d0
                     end do
                     s% am_nu_rot(1) = 0d0
                  
                  end if
                  
                  call dealloc

               end if
            
            end if
         
            if (s% am_nu_rot_flag) then ! check
               do k=1,nz
                  if (is_bad(s% am_nu_rot(k))) then
                     write(*,2) 'before return s% am_nu_rot(k)', k, s% am_nu_rot(k)
                     stop 'set_am_nu_rot'
                  end if
                  if (s% am_nu_rot(k) < 0d0) s% am_nu_rot(k) = 0d0
               end do
            end if         
         
         end subroutine set_am_nu_rot


         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
                sig, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                rhs, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                d, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                du, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                dl, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                bp, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                vp, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                xp, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                x, nz, nz_alloc_extra, 'mix_am_nu_rot', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays
         

      end subroutine update_rotation_mixing_info


      subroutine set_RTI_mixing_info(s, ierr)
         use chem_def, only: ih1, ihe4
         use star_utils, only: get_shock_info
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: &
            C, alpha_face, f, v, &
            am1, a00, ap1, min_dm, alfa0, alfa, cs, r, shock_mass_start, &
            log_max_boost, m_full_boost, m_no_boost, max_boost, &
            dm_for_center_eta_nondecreasing, min_eta
         integer :: k, nz, i_h1
         include 'formats'
         ierr = 0
         if (.not. s% RTI_flag) return
 
         nz = s% nz

         s% eta_RTI(1:nz) = 0d0
         s% etamid_RTI(1:nz) = 0d0
         s% boost_for_eta_RTI(1:nz) = 0d0
         s% sig_RTI(1:nz) = 0d0
         s% sigmid_RTI(1:nz) = 0d0

         if (s% RTI_C <= 0d0) return
             
         i_h1 = s% net_iso(ih1)
         
         shock_mass_start = 1d99
         do k = 1, nz
            if (s% u_flag) then
               v = s% u(k)
            else
               v = s% v(k)
            end if
            if (v > s% csound(k)) then
               if (k > 1) shock_mass_start = s% m(k) ! skip this after breakout
               exit
            end if
         end do
          
         min_dm = s% RTI_min_dm_behind_shock_for_full_on*Msun         
         log_max_boost = s% RTI_log_max_boost
         max_boost = exp10(log_max_boost)
         m_full_boost = s% RTI_m_full_boost*Msun
         m_no_boost = s% RTI_m_no_boost*Msun
         
         min_eta = -1d0
         dm_for_center_eta_nondecreasing = Msun*s% RTI_dm_for_center_eta_nondecreasing
         
         do k=1,nz
            
            f = max(0d0, s% X(k) - s% RTI_C_X0_frac*s% surface_h1)
            C = s% RTI_C*(1d0 + f*f*s% RTI_C_X_factor)
            
            if (s% m(k) < m_no_boost) then
               if (s% m(k) <= m_full_boost) then
                  C = C*max_boost
               else
                  alfa = (m_no_boost - s% m(k))/(m_no_boost - m_full_boost)
                  C = C*exp10(alfa*log_max_boost)
               end if
            end if

            if (shock_mass_start - s% m(k) < min_dm) &
               C = C*(shock_mass_start - s% m(k))/min_dm

            if (k > 1) then
               alpha_face = &
                  (s% dq(k-1)*s% alpha_RTI(k) + s% dq(k)*s% alpha_RTI(k-1))/ &
                     (s% dq(k-1) + s% dq(k))
               cs = s% csound_face(k)
               r = s% r(k)
               s% eta_RTI(k) = C*alpha_face*cs*r
               
               if (is_bad(s% eta_RTI(k))) then
                  ierr = -1
                  return
                  if (s% stop_for_bad_nums) then
                     write(*,2) 's% eta_RTI(k)', k, s% eta_RTI(k)
                     stop 'set_RTI_mixing_info'
                  end if
               end if
            
               if (s% m(k) - s% M_center <= dm_for_center_eta_nondecreasing) then
                  if (min_eta < 0d0) then
                     min_eta = s% eta_RTI(k)
                  else if (min_eta > s% eta_RTI(k)) then
                     s% eta_RTI(k) = min_eta
                  end if
               end if

            end if
            
            s% etamid_RTI(k) = max(min_eta, C*s% alpha_RTI(k)*s% csound(k)*s% rmid(k))
            s% boost_for_eta_RTI(k) = C/s% RTI_C
            
            if (is_bad(s% etamid_RTI(k))) then
               ierr = -1
               return
            end if

         end do
         
         call get_shock_info(s)

         ! sig_RTI(k) is mixing flow across face k in (gm sec^1)
         call get_RTI_sigmas(s, s% sig_RTI, s% eta_RTI, &
            s% rho_face, s% r, s% dm_bar, s% dt, ierr)
         if (ierr /= 0) return
         
         if (s% v_flag) then
            ! sigmid_RTI(k) is mixing flow at center k in (gm sec^1)
            call get_RTI_sigmas(s, s% sigmid_RTI, s% etamid_RTI, &
               s% rho, s% rmid, s% dm, s% dt, ierr)
            if (ierr /= 0) return
         end if

      end subroutine set_RTI_mixing_info


      subroutine get_RTI_sigmas(s, sig, eta, rho, r, dm_bar, dt, ierr)
         type (star_info), pointer :: s
         real(dp), pointer, dimension(:) :: sig, eta, rho, r, dm_bar
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer :: nz, k
         real(dp) :: D, f1, f, cdcterm, max_sig, &
            siglim, sig_term_limit, xmstar

         include 'formats'

         ierr = 0
         nz = s% nz
         xmstar = s% xmstar
         sig_term_limit = s% sig_term_limit

         sig(1) = 0d0
         max_sig = 0d0
         do k = 2, nz
            D = eta(k)
            if (D <= 0 .or. &
                  (s% q(k) > s% shock_q .and. s% shock_q > 0d0)) then
               sig(k) = 0d0
               cycle
            end if
            f1 = 4*pi*r(k)*r(k)*rho(k)
            f = f1*f1
            cdcterm = s% mix_factor*D*f
            sig(k) = cdcterm/dm_bar(k)
            if (is_bad(sig(k))) then
               if (s% stop_for_bad_nums) then
                  write(*,2) 'sig(k)', k, sig(k)
                  stop 'get_RTI_sigmas'
               end if
               sig(k) = 1d99
            end if
            if (sig(k) > max_sig) max_sig = sig(k)
         end do

         if (max_sig < 1) return
         do k = 1, nz
            if (sig(k) == 0) cycle
            siglim = sig_term_limit*dm_bar(k)/dt
            if (sig(k) > siglim) sig(k) = siglim
         end do

      end subroutine get_RTI_sigmas


      subroutine set_dPdr_dRhodr_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: rho, r00, alfa00, beta00, P_face00, rho_face00, &
            rp1, alfap1, betap1, dr_m1, dr_00, &
            c, d, am1, a00, ap1, v, rmid
         real(dp), pointer, dimension(:) :: dPdr, drhodr, P_face, rho_face
         integer :: k, nz, width
         logical, parameter :: do_slope_limiting = .false.
         include 'formats'
         ierr = 0
         if (.not. s% RTI_flag) return
         if (s% dPdr_dRhodr_info(1) >= 0d0) then
            if (is_bad(s% dPdr_dRhodr_info(1))) stop 'set_dPdr_dRhodr_info'
            return ! already set this step
         end if

         nz = s% nz
         
         call do_alloc(ierr)
         if (ierr /= 0) return

         do k=2,nz
            rho = s% rho(k)
            r00 = s% r(k)
            alfa00 = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta00 = 1d0 - alfa00
            P_face(k) = alfa00*s% P(k) + beta00*s% P(k-1)
            rho_face(k) = alfa00*s% rho(k) + beta00*s% rho(k-1)
         end do
         P_face(1) = s% P(1)
         rho_face(1) = s% rho(1)

         do k=1,nz
            if (s% u_flag) then
               v = s% u(k)
            else
               v = s% v(k)
            end if
            if (k == nz .or. k == 1 .or.  &
                  (s% alpha_RTI_start(k) == 0d0 .and. &
                     v/s% csound(k) < s% alpha_RTI_src_min_v_div_cs)) then
               drhodr(k) = 0d0
               dPdr(k) = 0d0
            else if (do_slope_limiting) then
               dr_m1 = s% r(k-1) - s% r(k)
               dr_00 = s% r(k) - s% r(k+1)
               dPdr(k) = slope_limit(P_face, k, dr_m1, dr_00)
               drhodr(k) = slope_limit(rho_face, k, dr_m1, dr_00)     
            else
               !dr_00 = s% r(k) - s% r(k+1)
               rmid = 0.5d0*(s% r(k) + s% r(k+1))
               dr_00 = s% dm(k)/(4d0*pi*rmid*rmid*s% rho(k)) ! don't subtract r's to get dr
               dPdr(k) = (P_face(k) - P_face(k+1))/dr_00
               drhodr(k) = (rho_face(k) - rho_face(k+1))/dr_00
            end if
         end do

         ! smooth dPdr and drhodr
         if (s% RTI_smooth_mass > 0d0 .and. s% RTI_smooth_iterations > 0) then
            call do_smoothing_by_mass( &
               s, s% RTI_smooth_mass, s% RTI_smooth_iterations, drhodr, ierr)
            if (ierr /= 0) return
            call do_smoothing_by_mass( &
               s, s% RTI_smooth_mass, s% RTI_smooth_iterations, dPdr, ierr)
            if (ierr /= 0) return
         else if (s% RTI_smooth_fraction < 1d0) then
            c = s% RTI_smooth_fraction
            d = 0.5d0*(1d0 - c)
            am1 = 0
            a00 = dPdr(1)
            ap1 = dPdr(2)
            do k=2,nz-1
               am1 = a00
               a00 = ap1
               ap1 = dPdr(k+1)
               dPdr(k) = c*a00 + d*(am1 + ap1)
            end do
            dPdr(nz) = c*a00 + 2*d*am1
            a00 = drhodr(1)
            ap1 = drhodr(2)
            do k=2,nz-1
               am1 = a00
               a00 = ap1
               ap1 = drhodr(k+1)
               if (am1 == 0d0 .and. a00 == 0d0) cycle
               drhodr(k) = c*a00 + d*(am1 + ap1)
            end do
            drhodr(nz) = c*a00 + 2*d*am1
         end if

         ! set
         do k=1,nz
            s% dPdr_info(k) = dPdr(k)
            s% dRhodr_info(k) = drhodr(k)
            s% dPdr_dRhodr_info(k) = min(0d0,dPdr(k)*drhodr(k))
         end do
         
         call dealloc

         contains

         real(dp) function slope_limit(v, k, dr_m1, dr_00)
            real(dp), intent(in) :: v(:), dr_m1, dr_00
            integer, intent(in) :: k
            real(dp) :: s_00, s_p1
            s_00 = (v(k-1) - v(k))/dr_m1
            s_p1 = (v(k) - v(k+1))/dr_00
            if (s_00*s_p1 < 0d0) then
               slope_limit = 0d0
            else if(abs(s_00) > abs(s_p1)) then
               slope_limit = s_p1
            else
               slope_limit = s_00
            end if
         end function slope_limit

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               P_face, nz, 0, 'set_dPdr_dRhodr_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               rho_face, nz, 0, 'set_dPdr_dRhodr_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dPdr, nz, 0, 'set_dPdr_dRhodr_info', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               drhodr, nz, 0, 'set_dPdr_dRhodr_info', ierr)
         if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine set_dPdr_dRhodr_info


      subroutine do_smoothing_by_mass( &
            s, smooth_mass, number_iterations, val, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: smooth_mass
         integer, intent(in) :: number_iterations
         real(dp), pointer :: val(:)
         integer, intent(out) :: ierr
         integer :: nz, iter, k_center, k_inner, k_outer, j, k
         real(dp) :: mlo, mhi, mmid, smooth_m, v, dm_half, mtotal, mass_from_cell
         real(dp), allocatable :: work(:)
         include 'formats'
         ierr = 0
         nz = s% nz
         allocate(work(nz))
         do k = 1, nz
            work(k) = val(k)
         end do
         smooth_m = smooth_mass*Msun
         dm_half = 0.5d0*smooth_m
         do iter = 1, number_iterations
            do k_center = nz-1, 2, -1
               mmid = s% m(k_center) - 0.5d0*s% dm(k_center)
               mlo = mmid - dm_half
               mhi = mmid + dm_half
               k_inner = k_center
               v = 0d0
               mtotal = 0d0
               do k=k_center, nz, -1
                  if (s% m(k) - 0.5d0*s% dm(k) <= mlo) then
                     k_inner = k
                     mass_from_cell = s% m(k) - mlo
                     v = v + val(k)*mass_from_cell
                     mtotal = mtotal + mass_from_cell
                     exit
                  end if
               end do
               k_outer = k_center
               do k=k_center,1,-1
                  if (s% m(k) >= mhi) then
                     k_outer = k
                     mass_from_cell = mhi - (s% m(k) - s% dm(k))
                     v = v + val(k)*mass_from_cell
                     mtotal = mtotal + mass_from_cell
                     exit
                  end if
               end do
               do k=k_outer+1,k_inner-1
                  v = v + val(k)*s% dm(k)
                  mtotal = mtotal + s% dm(k)
               end do
               if(mtotal>0d0) then
                  work(k_center) = v/mtotal
               else
                  work(k_center) = 0d0
               end if
            end do
            do k = 2, nz-1
               val(k) = work(k)
            end do
         end do
      end subroutine do_smoothing_by_mass


      subroutine set_dxdt_mix(s)

         type (star_info), pointer :: s

         real(dp) :: x00, xp1, xm1, dx00, dxp1, dm, sig00, sigp1, &
              flux00, dflux00_dxm1, dflux00_dx00, &
              fluxp1, dfluxp1_dx00, dfluxp1_dxp1

         integer :: j, k
         
         include 'formats'

         do k = 1, s% nz

            dm = s% dm(k)
            sig00 = s% sig(k)

            if (k < s% nz) then
               sigp1 = s% sig(k+1)
            else
               sigp1 = 0
            end if

            if (k > 1) then
               dflux00_dxm1 = -sig00
               dflux00_dx00 = sig00
            else
               dflux00_dxm1 = 0
               dflux00_dx00 = 0
            end if

            if (k < s% nz) then
               dfluxp1_dx00 = -sigp1
               dfluxp1_dxp1 = sigp1
            else
               dfluxp1_dx00 = 0
               dfluxp1_dxp1 = 0
            end if

            s% d_dxdt_mix_dx00(k) = (dfluxp1_dx00 - dflux00_dx00)/dm
            s% d_dxdt_mix_dxm1(k) = -dflux00_dxm1/dm
            s% d_dxdt_mix_dxp1(k) = dfluxp1_dxp1/dm


            do j=1, s% species
               x00 = s% xa(j,k)
               if (k > 1) then
                  xm1 = s% xa(j,k-1)
                  dx00 = xm1 - x00
                  flux00 = -sig00*dx00
               else
                  flux00 = 0
               end if
               if (k < s% nz) then
                  xp1 = s% xa(j,k+1)
                  dxp1 = x00 - xp1
                  fluxp1 = -sigp1*dxp1
               else
                  fluxp1 = 0
               end if
               s% dxdt_mix(j,k) = (fluxp1 - flux00)/dm
            end do

         end do

      end subroutine set_dxdt_mix



      subroutine add_RTI_turbulence(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: nz, k
         real(dp) :: coeff, D, P_face, r_face, q_face, cdc, Hp, vc, &
            alfa, beta, rho_face, lg_coeff, full_off, full_on, min_m
         logical :: dbg

         include 'formats'

         dbg = .false.

         ierr = 0
         nz = s% nz

         if (.not. s% RTI_flag) return

         coeff = s% composition_RTI_diffusion_factor
         if (coeff <= 0) return
         lg_coeff = log10(coeff)
         full_off = s% min_M_RTI_factors_full_off*Msun
         full_on = s% max_M_RTI_factors_full_on*Msun
         min_m = s% RTI_min_m_for_D_mix_floor*Msun
         do k = 2, nz
            D = s% eta_RTI(k)
            if (s% m(k) <= full_on) then
               D = D*coeff
            else if (s% m(k) < full_off) then
               D = D*exp10( &
                        lg_coeff*(full_off - s% m(k))/(full_off - full_on))
            ! else full off, i.e., coeff = 1
            end if
            if (is_bad(D)) then
               write(*,2) 'eta_RTI', k, D
               if (s% stop_for_bad_nums) stop 'add_RTI_turbulence'
               ierr = -1
               cycle
            end if
            if (D < s% RTI_D_mix_floor .and. s% m(k) >= min_m) &
               D = s% RTI_D_mix_floor
            if (D > s% D_mix(k)) &
               s% mixing_type(k) = rayleigh_taylor_mixing
            D = D + s% D_mix(k)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1 - alfa
            rho_face = alfa*s% rho(k) + beta*s% rho(k-1)
            P_face = alfa*s% P(k) + beta*s% P(k-1)
            r_face = s% r(k)
            q_face = s% q(k)
            cdc = pow2(pi4*s% r(k)*s% r(k)*rho_face)*D ! gm^2/sec
            if (s% cgrav(k) <= 0) then
               Hp = s% r(k)
            else
               Hp = P_face/(rho_face*s% cgrav(k)* &
                  (s% M_center + s% xmstar*q_face)/(r_face*r_face))
            end if
            vc = 3*D/(s% alpha_mlt(k)*Hp)
            s% cdc(k) = cdc
            s% D_mix(k) = D
            s% conv_vel(k) = vc
         end do

      end subroutine add_RTI_turbulence


      end module mix_info





