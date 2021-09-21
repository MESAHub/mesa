! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
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

      module astero_lib
      ! library for calculating neutrino losses from non-nuclear-burning sources
      ! neutrino losses that occur during nuclear reactions are included in the nuclear library
      ! the data interface for the library is defined in neu_def
      
      use const_def, only: dp
      use gyre_support, only: GYRE_IS_ENABLED

      use astero_support, only: &
            astero_get_cubic_all_freq_corr => get_cubic_all_freq_corr, &
            astero_get_combined_all_freq_corr => get_combined_all_freq_corr, &
            astero_get_power_law_all_freq_corr => get_power_law_all_freq_corr, &
            astero_get_sonoi_all_freq_corr => get_sonoi_all_freq_corr
      implicit none
      
      logical, parameter :: astero_gyre_is_enabled = GYRE_IS_ENABLED

      contains ! the procedure interface for the library
      ! client programs should only call these routines.


      subroutine run_star_astero( &
            extras_controls, inlist_astero_search_controls_fname)
         use astero_run_support, only: do_run_star_astero
         use astero_def, only: init_astero_def
         interface
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         character (len=256) :: inlist_astero_search_controls_fname
         optional inlist_astero_search_controls_fname
         call init_astero_def
         call do_run_star_astero( &
            extras_controls, inlist_astero_search_controls_fname)
      end subroutine run_star_astero
      
      
      ! this can be called from user run_star_extras check model routine
      subroutine adipls_get_one_el_info( &
            s, l, nu1, nu2, iscan, R, G, M, &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, store_for_adipls, &
            save_mode_info, order_to_save_in, save_mode_filename_in, &
            num, l_freq, l_inertia, l_order, l_em, ierr)
         use star_def, only: star_info
         use adipls_support, only: do_adipls_get_one_el_info
         type (star_info), pointer :: s
         integer, intent(in) :: l, iscan
         real(dp), intent(in) :: nu1, nu2, R, G, M
         logical, intent(in) :: &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, store_for_adipls, save_mode_info
         integer, intent(in) :: order_to_save_in
         character (len=*), intent(in) :: save_mode_filename_in
         integer, intent(out) :: num
         real(dp), pointer, dimension(:) :: l_freq, l_inertia
         integer, pointer, dimension(:) :: l_order, l_em
         integer, intent(out) :: ierr
         call do_adipls_get_one_el_info( &
            s, l, nu1, nu2, iscan, R, G, M, &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, store_for_adipls, &
            save_mode_info, order_to_save_in, save_mode_filename_in, &
            num, l_freq, l_inertia, l_order, l_em, ierr)
      end subroutine adipls_get_one_el_info
      
      
      subroutine astero_gyre_get_modes(id, el, store_model, ierr)
         use star_def, only: star_ptr, star_info
         use gyre_support, only: do_gyre_get_modes
         integer, intent(in)       :: id
         integer, intent(in)       :: el
         logical, intent(in)       :: store_model
         integer, intent(out)      :: ierr
         type (star_info), pointer :: s
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call do_gyre_get_modes(s, el, store_model, ierr)
      end subroutine astero_gyre_get_modes
      
      
      ! for surface_effects test case


      subroutine astero_get_one_el_info( &
            s, l, nu1, nu2, iscan, i1, i2, store_model, code, ierr)
         use astero_support, only: get_one_el_info
         use star_def, only: star_info
         type (star_info), pointer :: s
         integer, intent(in) :: l, iscan, i1, i2
         real(dp), intent(in) :: nu1, nu2
         logical, intent(in) :: store_model
         character (len=*), intent(in) :: code
         integer, intent(out) :: ierr
         call get_one_el_info( &
            s, l, nu1, nu2, iscan, i1, i2, store_model, code, ierr)
      end subroutine astero_get_one_el_info
      
      
      real(dp) function astero_interpolate_l0_inertia(freq)
         use astero_support, only: interpolate_l0_inertia
         real(dp), intent(in) :: freq
         astero_interpolate_l0_inertia = interpolate_l0_inertia(freq)
      end function astero_interpolate_l0_inertia
      
      
      subroutine astero_get_kjeldsen_radial_freq_corr( &
            a_div_r, b, nu_max, correction_factor, check_obs, &
            nl0, l0_obs, l0_freq, l0_freq_corr, l0_inertia)
         use astero_support, only: get_kjeldsen_radial_freq_corr
         real(dp), intent(in) :: a_div_r, b, nu_max, correction_factor
         logical, intent(in) :: check_obs ! if false, then l0_obs is not used
         integer, intent(in) :: nl0
         real(dp), intent(in), dimension(:) :: &
            l0_obs, l0_freq, l0_inertia
         real(dp), intent(inout) :: l0_freq_corr(:)
         call get_kjeldsen_radial_freq_corr( &
            a_div_r, b, nu_max, correction_factor, check_obs, &
            nl0, l0_obs, l0_freq, l0_freq_corr, l0_inertia)
      end subroutine astero_get_kjeldsen_radial_freq_corr
      
      
      end module astero_lib

