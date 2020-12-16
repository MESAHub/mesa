! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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

      module evolve_support

      use star_private_def
      use const_def

      implicit none

      private
      public :: set_current_to_old, new_generation, output, output_to_file


      contains


      subroutine new_generation(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp), pointer :: tmp(:,:), tmp1(:)
         integer, pointer :: itmp1(:)
         integer :: nz, i, k, i_u

         include 'formats'

         ierr = 0
         nz = s% nz
         
         s% nz_old = s% nz

         s% mstar_old = s% mstar

         s% xmstar_old = s% xmstar

         s% M_center_old = s% M_center

         s% v_center_old = s% v_center

         s% R_center_old = s% R_center

         s% L_center_old = s% L_center

         s% total_radiation_old = s% total_radiation

         s% total_angular_momentum_old = s% total_angular_momentum

         s% Lsurf_m_old = s% Lsurf_m
         
         if (.not. s% rsp_flag) then

            call flip(s% q, s% q_old, ierr)
            if (ierr /= 0) return

            call flip(s% dq, s% dq_old, ierr)
            if (ierr /= 0) return
         
            if (.not. s% conv_vel_flag) then
               call flip(s% conv_vel, s% conv_vel_old, ierr)
               if (ierr /= 0) return
            end if
            
            ! keep old values for rotation variables in order to do time smoothing in D_omega

            call flip(s% nu_ST, s% nu_ST_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_ST, s% D_ST_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_DSI, s% D_DSI_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_SH, s% D_SH_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_SSI, s% D_SSI_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_ES, s% D_ES_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_GSF, s% D_GSF_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_mix, s% D_mix_old, ierr)
            if (ierr /= 0) return

            call flip(s% dPdr_dRhodr_info, s% dPdr_dRhodr_info_old, ierr)
            if (ierr /= 0) return

            call flip(s% omega, s% omega_old, ierr)
            if (ierr /= 0) return

            call flip(s% j_rot, s% j_rot_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_omega, s% D_omega_old, ierr)
            if (ierr /= 0) return

            call flip(s% am_nu_rot, s% am_nu_rot_old, ierr)
            if (ierr /= 0) return

            call flip(s% D_smooth, s% D_smooth_old, ierr)
            if (ierr /= 0) return

            tmp => s% xh_old
            s% xh_old => s% xh
            call enlarge_if_needed_2(tmp,s% nvar_hydro,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(tmp)
            s% xh => tmp

            tmp => s% xa_old
            s% xa_old => s% xa
            call enlarge_if_needed_2(tmp,s% species,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(tmp)
            s% xa => tmp
         
         end if

         s% time_old = s% time

         s% L_nuc_burn_total_old = s% L_nuc_burn_total
         s% power_nuc_burn_old = s% power_nuc_burn
         s% power_h_burn_old = s% power_h_burn
         s% power_he_burn_old = s% power_he_burn
         s% power_c_burn_old = s% power_c_burn
         s% power_photo_old = s% power_photo
         s% power_z_burn_old = s% power_z_burn
         s% power_nuc_neutrinos_old = s% power_nuc_neutrinos
         s% power_nonnuc_neutrinos_old = s% power_nonnuc_neutrinos
         s% power_neutrinos_old = s% power_neutrinos

         do i=1,num_categories
            s% L_by_category_old(i) = s% L_by_category(i)
         end do

         s% L_phot_old = s% L_phot

         s% mstar_dot_old = s% mstar_dot

         s% L_surf_old = s% L_surf

         s% v_surf_old = s% v_surf

         s% gradT_excess_alpha_old = s% gradT_excess_alpha

         s% h1_czb_mass_old = s% h1_czb_mass_prev

         s% he_core_mass_old = s% he_core_mass

         s% c_core_mass_old = s% c_core_mass

         s% Teff_old = s% Teff

         s% center_eps_nuc_old = s% center_eps_nuc

         s% Lrad_div_Ledd_avg_surf_old = s% Lrad_div_Ledd_avg_surf

         s% w_div_w_crit_avg_surf_old = s% w_div_w_crit_avg_surf

         s% n_conv_regions_old = s% n_conv_regions
         s% cz_bot_mass_old(:) = s% cz_bot_mass(:)
         s% cz_top_mass_old(:) = s% cz_top_mass(:)         

         s% model_number_old = s% model_number

         s% dt_limit_ratio_old = s% dt_limit_ratio

         s% revised_net_name_old = s% revised_net_name

         s% revised_max_yr_dt_old = s% revised_max_yr_dt

         s% astero_revised_max_yr_dt_old = s% astero_revised_max_yr_dt

         s% total_internal_energy_old = s% total_internal_energy

         s% total_gravitational_energy_old = s% total_gravitational_energy

         s% total_radial_kinetic_energy_old = s% total_radial_kinetic_energy

         s% total_turbulent_energy_old = s% total_turbulent_energy

         s% total_rotational_kinetic_energy_old = s% total_rotational_kinetic_energy

         s% total_energy_old = s% total_energy

         s% cumulative_work_outward_at_surface_old = s% cumulative_work_outward_at_surface

         s% cumulative_work_inward_at_center_old = s% cumulative_work_inward_at_center

         s% cumulative_energy_error_old = s% cumulative_energy_error

         s% cumulative_eps_grav_old = s% cumulative_eps_grav

         s% cumulative_L_center_old = s% cumulative_L_center

         s% cumulative_L_surf_old = s% cumulative_L_surf

         s% cumulative_extra_heating_old = s% cumulative_extra_heating

         s% cumulative_irradiation_heating_old = s% cumulative_irradiation_heating

         s% cumulative_WD_sedimentation_heating_old = s% cumulative_WD_sedimentation_heating

         s% cumulative_nuclear_heating_old = s% cumulative_nuclear_heating

         s% cumulative_non_nuc_neu_cooling_old = s% cumulative_non_nuc_neu_cooling

         s% cumulative_sources_and_sinks_old = s% cumulative_sources_and_sinks

         s% log_P_center_old = s% log_P_center

         s% min_kap_floor_old = s% min_kap_floor

         s% was_in_implicit_wind_limit_old = s% was_in_implicit_wind_limit

         s% model_number_for_last_retry_old = s% model_number_for_last_retry

         s% dt_old = s% dt

         do i = 1, s% len_extra_work
            s% extra_work_old(i) = s% extra_work(i)
         end do

         do i = 1, s% len_extra_iwork
            s% extra_iwork_old(i) = s% extra_iwork(i)
         end do

         s% ixtra_old = s% ixtra
         s% xtra_old = s% xtra
         s% lxtra_old = s% lxtra

         call s% other_new_generation(s% id, ierr)
         
         if (s% fill_arrays_with_NaNs) s% need_to_setvars = .true.

         contains

         subroutine flip(ptr, ptr_old, ierr)
            real(dp), pointer, dimension(:) :: ptr, ptr_old
            integer, intent(out) :: ierr
            logical :: first_time
            ierr = 0
            tmp1 => ptr_old
            ptr_old => ptr
            first_time = (.not. associated(tmp1))
            call realloc_if_needed_1(tmp1,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) then
               call fill_with_NaNs(tmp1)
            else if (s% zero_when_allocate) then
               tmp1(:) = 0
            else if (first_time) then
               tmp1(1:nz) = -9d99
            end if
            ptr => tmp1
         end subroutine flip

      end subroutine new_generation


      subroutine set_current_to_old(s)
         type (star_info), pointer :: s
         real(dp), pointer :: p1(:)
         integer :: i, k

         include 'formats'

         s% nz = s% nz_old
         s% mstar = s% mstar_old
         s% xmstar = s% xmstar_old
         s% M_center = s% M_center_old
         s% v_center = s% v_center_old
         s% R_center = s% R_center_old
         s% L_center = s% L_center_old

         s% total_radiation = s% total_radiation_old

         s% total_angular_momentum = s% total_angular_momentum_old
         s% Lsurf_m = s% Lsurf_m_old
         s% L_nuc_burn_total = s% L_nuc_burn_total_old
         s% L_phot = s% L_phot_old
         
         s% power_nuc_burn = s% power_nuc_burn_old
         s% power_h_burn = s% power_h_burn_old
         s% power_he_burn = s% power_he_burn_old
         s% power_c_burn = s% power_c_burn_old
         s% power_photo = s% power_photo_old
         s% power_z_burn = s% power_z_burn_old
         s% power_nuc_neutrinos = s% power_nuc_neutrinos_old
         s% power_nonnuc_neutrinos = s% power_nonnuc_neutrinos_old
         s% power_neutrinos = s% power_neutrinos_old
         
         s% mstar_dot = s% mstar_dot_old
         s% gradT_excess_alpha = s% gradT_excess_alpha_old
         s% L_surf = s% L_surf_old
         s% v_surf = s% v_surf_old

         s% he_core_mass = s% he_core_mass_old
         s% c_core_mass = s% c_core_mass_old
         s% Teff = s% Teff_old
         s% center_eps_nuc = s% center_eps_nuc_old
         s% Lrad_div_Ledd_avg_surf = s% Lrad_div_Ledd_avg_surf_old
         s% w_div_w_crit_avg_surf = s% w_div_w_crit_avg_surf_old
         s% n_conv_regions = s% n_conv_regions_old
         s% cz_bot_mass(:) = s% cz_bot_mass_old(:)
         s% cz_top_mass(:) = s% cz_top_mass_old(:)         
         s% dt_limit_ratio = s% dt_limit_ratio_old
         s% revised_net_name = s% revised_net_name_old
         s% revised_max_yr_dt = s% revised_max_yr_dt_old
         s% astero_revised_max_yr_dt = s% astero_revised_max_yr_dt_old

         s% total_internal_energy = s% total_internal_energy_old
         s% total_gravitational_energy = s% total_gravitational_energy_old
         s% total_radial_kinetic_energy = s% total_radial_kinetic_energy_old
         s% total_turbulent_energy = s% total_turbulent_energy_old
         s% total_rotational_kinetic_energy = s% total_rotational_kinetic_energy_old
         s% total_energy = s% total_energy_old

         s% cumulative_work_outward_at_surface = s% cumulative_work_outward_at_surface_old
         s% cumulative_work_inward_at_center = s% cumulative_work_inward_at_center_old
         s% cumulative_energy_error = s% cumulative_energy_error_old
         s% cumulative_eps_grav = s% cumulative_eps_grav_old
         s% cumulative_L_center = s% cumulative_L_center_old
         s% cumulative_L_surf = s% cumulative_L_surf_old
         s% cumulative_extra_heating = s% cumulative_extra_heating_old
         s% cumulative_irradiation_heating = s% cumulative_irradiation_heating_old
         s% cumulative_WD_sedimentation_heating = s% cumulative_WD_sedimentation_heating_old
         s% cumulative_nuclear_heating = s% cumulative_nuclear_heating_old
         s% cumulative_non_nuc_neu_cooling = s% cumulative_non_nuc_neu_cooling_old
         s% cumulative_sources_and_sinks = s% cumulative_sources_and_sinks_old
         s% log_P_center = s% log_P_center_old
         s% min_kap_floor = s% min_kap_floor_old

         s% was_in_implicit_wind_limit = s% was_in_implicit_wind_limit_old
         
         s% model_number_for_last_retry = s% model_number_for_last_retry_old

         do i = 1, s% len_extra_work
            s% extra_work(i) = s% extra_work_old(i)
         end do

         do i = 1, s% len_extra_iwork
            s% extra_iwork(i) = s% extra_iwork_old(i)
         end do

         s% ixtra = s% ixtra_old
         s% xtra = s% xtra_old
         s% lxtra = s% lxtra_old

         call s% other_set_current_to_old(s% id)

      end subroutine set_current_to_old


      subroutine output(id, ierr)
         use star_utils, only: get_name_for_restart_file
         interface
            subroutine save_restart_info(iounit, id, ierr)
               integer, intent(in) :: iounit
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine save_restart_info
         end interface
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         character (len=strlen) :: filename, num_str, fstring
         type (star_info), pointer :: s
         integer :: num_digits

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call get_name_for_restart_file(s% model_number, s% photo_digits, num_str)
         filename = trim(s% photo_directory) // '/' // trim(num_str)
         call output_to_file(filename, id, ierr)
         if (ierr /= 0) return

         write(*, '(a)', advance='no') 'save ' // trim(filename)
         num_digits = 1 + log10(dble(max(1,s% model_number)))
         write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(*,fstring) ' for model ', s% model_number

      end subroutine output


      subroutine output_to_file(filename, id, ierr)
         use photo_out, only: output_star_photo
         character (len=*) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         integer :: iounit, k
         type (star_info), pointer :: s
         character(len=strlen) :: iomsg

         include 'formats'

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         open(newunit=iounit, file=trim(filename), action='write', &
            status='replace', iostat=ierr, iomsg=iomsg, form='unformatted')
         if (ierr == 0) then
            s% most_recent_photo_name = trim(filename)
            call output_star_photo(s, iounit, ierr)
            close(iounit)
         else
            write(*,*) trim(iomsg)
         endif

      end subroutine output_to_file


      end module evolve_support


