! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module do_one_utils
      
      use star_private_def
      use const_def
      use utils_lib, only: is_bad

      implicit none
      
      private
      public :: do_one_check_model, do_one_finish, &
         write_terminal_header, do_bare_bones_check_model, do_check_limits, &
         do_show_log_description, do_show_terminal_header, do_terminal_summary
      
      ! model log priorities
      integer, parameter :: delta_priority = 1
      integer, parameter :: phase_priority = 2

      
      
      contains
      
      
      logical function model_is_okay(s)
         type (star_info), pointer :: s
         ! for now, just check for valid number in the final dynamic timescale
         model_is_okay = ((s% dynamic_timescale - s% dynamic_timescale) .eq. 0d0) &
                        .and. ((s% dynamic_timescale + 1d0) > 1d0)
      end function model_is_okay

      
      subroutine do_one_finish(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine do_one_finish
      
      
      subroutine set_save_profiles_info(s, model_priority)
         type (star_info), pointer :: s
         integer, intent(in) :: model_priority
         s% need_to_save_profiles_now = .true.
         s% save_profiles_model_priority = model_priority
      end subroutine set_save_profiles_info
      
      
      subroutine write_terminal_header(s)
         type (star_info), pointer :: s
         if (s% model_number <= s% recent_log_header) return
         if (s% just_wrote_terminal_header) return
         s% recent_log_header = s% model_number
         call do_show_terminal_header(s)
         s% just_wrote_terminal_header = .true.
      end subroutine write_terminal_header
      
      
      subroutine do_show_log_description(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         write(*,*)
         write(*,'(a)') " The terminal output contains the following information"
         write(*,*)
         write(*,'(a)') "      'step' is the number of steps since the start of the run,"
         write(*,'(a)') "      'lg_dt' is log10 timestep in years,"
         write(*,'(a)') "      'age_yr' is the simulated years since the start run,"
         write(*,'(a)') "      'lg_Tcntr' is log10 center temperature (K),"
         write(*,'(a)') "      'lg_Dcntr' is log10 center density (g/cm^3),"
         write(*,'(a)') "      'lg_Tmax' is log10 max temperature (K),"
         write(*,'(a)') "      'Teff' is the surface temperature (K),"
         write(*,'(a)') "      'lg_R' is log10 surface radius (Rsun),"
         write(*,'(a)') "      'lg_L' is log10 surface luminosity (Lsun),"
         write(*,'(a)') "      'lg_LH' is log10 total PP and CNO hydrogen burning power (Lsun),"
         write(*,'(a)') "      'lg_L3a' is log10 total triple-alpha helium burning power (Lsun),"
         write(*,'(a)') "      'lg_gsurf' is log10 surface gravity,"
         write(*,'(a)') "      'lg_LNuc' is log10 nuclear power (Lsun),"
         write(*,'(a)') "      'lg_LNeu' is log10 total neutrino power (Lsun),"
         write(*,'(a)') "      'lg_Lphoto' is log10 total photodisintegration (Lsun),"
         write(*,'(a)') "      'Mass' is the total stellar baryonic mass (Msun),"
         write(*,'(a)') "      'lg_Mdot' is log10 magnitude of rate of change of mass (Msun/year),"
         write(*,'(a)') "      'lg_Dsurf' is log10 surface density (g/cm^3),"
         write(*,'(a)') "      'H_env' is the amount of mass where H is the most abundant iso,"
         write(*,'(a,e9.2)') "      'He_core' is the largest mass where He is most abundant iso."
         write(*,'(a,e9.2)') "      'C_core' is the largest mass where C is most abundant iso."
         write(*,'(a)') "      'H_cntr' is the center H1 mass fraction,"
         write(*,'(a)') "      'He_cntr' is the center He4 mass fraction,"
         write(*,'(a)') "      'C_cntr' is the center C12 mass fraction,"
         write(*,'(a)') "      'N_cntr' is the center N14 mass fraction,"
         write(*,'(a)') "      'O_cntr' is the center O16 mass fraction,"
         write(*,'(a)') "      'Ne_cntr' is the center Ne20 mass fraction,"
         write(*,'(a)') "      'gam_cntr' is the center plasma interaction parameter,"
         write(*,'(a)') "      'eta_cntr' is the center electron degeneracy parameter,"
         write(*,'(a)') "      'zones' is the number of zones in the current model,"
         write(*,'(a)') "      'iters' is the number of solver iterations for the current step,"
         write(*,'(a)') "      'retry' is the number of step retries required during the run,"
         write(*,'(a)') "      'dt_limit' is an indication of what limited the timestep."
         write(*,*)
         write(*,'(a)') " All this and more are saved in the LOGS directory during the run."
      end subroutine do_show_log_description
      
      
      subroutine do_show_terminal_header(s)
         type (star_info), pointer :: s
         integer :: id, ierr, io
         call output_terminal_header(s,terminal_iounit)
         if (len_trim(s% extra_terminal_output_file) > 0) then
            ierr = 0
            open(newunit=io, file=trim(s% extra_terminal_output_file), &
               action='write', position='append', iostat=ierr)
            if (ierr == 0) then
               call output_terminal_header(s,io)
               close(io)
            else
               write(*,*) 'failed to open extra_terminal_output_file ' // &
                  trim(s% extra_terminal_output_file)
            end if
         end if
      end subroutine do_show_terminal_header
      
      
      subroutine output_terminal_header(s,io)
         use chem_def, only: isi28
         type (star_info), pointer :: s
         integer, intent(in) :: io
         character (len=5) :: iters
         iters = 'iters'
         
         include 'formats'
         
         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'
         write(io,*)
         write(io,'(a)') &
            '       step    lg_Tmax     Teff     lg_LH      lg_Lnuc     Mass       ' // &
            'H_rich     H_cntr     N_cntr     Y_surf   eta_cntr   zones  retry'
            
         ! note that if the age is in days, then the timestep is automatically in seconds.
         if (trim(s% terminal_show_timestep_units) == 'seconds' .or. &
             trim(s% terminal_show_timestep_units) == 'secs') then
            if (s% terminal_show_log_dt) then
               write(io,'(a)',advance='no') '  lg_dt_sec'
            else
               write(io,'(a)',advance='no') '     dt_sec'
            end if
         else if (trim(s% terminal_show_timestep_units) == 'days') then
            if (s% terminal_show_log_dt) then         
               write(io,'(a)',advance='no') ' lg_dt_days'
            else
               write(io,'(a)',advance='no') '    dt_days'
            end if
         else if (trim(s% terminal_show_timestep_units) == 'years' .or. &
                  trim(s% terminal_show_timestep_units) == 'yrs') then
            if (s% terminal_show_log_dt) then         
               write(io,'(a)',advance='no') '  lg_dt_yrs'
            else
               write(io,'(a)',advance='no') '     dt_yrs'
            end if
         else
            write(*,*) 'unrecognized option for terminal_show_timestep_units ' // trim(s% terminal_show_timestep_units)
            return
         end if
         
         if (s% initial_z >= 1d-5) then
            write(io,'(a)') &
               '    lg_Tcntr    lg_R     lg_L3a     lg_Lneu     lg_Mdot    ' // &
               'He_core    He_cntr    O_cntr     Z_surf   gam_cntr   ' // iters // '  '
         else
            write(io,'(a)') &
               '    lg_Tcntr    lg_R     lg_L3a     lg_Lneu     lg_Mdot    ' // &
               'He_core    He_cntr    O_cntr     lg_Z_surf gam_cntr  ' // iters // '  '
         end if

         if (trim(s% terminal_show_age_units) == 'seconds' .or. &
             trim(s% terminal_show_age_units) == 'secs') then
            if (s% terminal_show_log_age) then
               write(io,'(a)',advance='no') 'lg_age_secs'
            else
               write(io,'(a)',advance='no') '   age_secs'
            end if
         else if (trim(s% terminal_show_age_units) == 'days') then
            if (s% terminal_show_log_age) then
               write(io,'(a)',advance='no') 'lg_age_days'
            else
               write(io,'(a)',advance='no') '   age_days'
            end if
         else if (trim(s% terminal_show_age_units) == 'years' .or. &
                  trim(s% terminal_show_age_units) == 'yrs') then
            if (s% terminal_show_log_age) then
               write(io,'(a)',advance='no') ' lg_age_yrs'
            else
               write(io,'(a)',advance='no') '    age_yrs'
            end if
         else
            write(*,*) 'unrecognized option for terminal_show_age_units ' // trim(s% terminal_show_age_units)
            return
         end if
         
         if (s% net_iso(isi28) == 0) then
            write(io,'(a)') &
               '    lg_Dcntr    lg_L     lg_LZ      lg_Lphoto   lg_Dsurf   ' // &
               'CO_core    C_cntr     Ne_cntr    Z_cntr   v_div_cs       dt_limit'
         else
            write(io,'(a)') &
               '    lg_Dcntr    lg_L     lg_LZ      lg_Lphoto   lg_Dsurf   ' // &
               'CO_core    C_cntr     Ne_cntr    Si_cntr  v_div_cs       dt_limit'
         end if
         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'
         write(io,*)
         
      end subroutine output_terminal_header
      
      
      subroutine do_terminal_summary(s)
         type (star_info), pointer :: s
         integer :: id, ierr, io
         call output_terminal_summary(s,terminal_iounit)
         if (len_trim(s% extra_terminal_output_file) > 0) then
            ierr = 0
            open(newunit=io, file=trim(s% extra_terminal_output_file), &
               action='write', position='append', iostat=ierr)
            if (ierr == 0) then
               call output_terminal_summary(s,io)
               close(io)
            else
               write(*,*) 'failed to open extra_terminal_output_file ' // &
                  trim(s% extra_terminal_output_file)
            end if
         end if
      end subroutine do_terminal_summary
      
      
      subroutine output_terminal_summary(s,io)
         use num_def, only:banded_matrix_type
         use const_def, only:secyer
         use chem_def
         use rates_def, only: i_rate
         use star_utils, only:eval_current_y, eval_current_z
         type (star_info), pointer :: s
         integer, intent(in) :: io
         
         real(dp) :: time_step, age, dt, Xmax, v, vsurf_div_csound, tmp, &
            sum_Lnuc, sum_LH, sum_LHe, sum_Lz, sum_Lphoto
         integer :: model, ierr, nz, iters
         character (len=3) :: id_str
         character (len=32) :: why
         character (len=90) :: fmt, fmt1, fmt2, fmt3, fmt4, fmt5
         
         include 'formats'
         
         age = s% star_age ! in years            
         if (trim(s% terminal_show_age_units) == 'seconds' .or. &
             trim(s% terminal_show_age_units) == 'secs') then
            age = age*secyer
         else if (trim(s% terminal_show_age_units) == 'days') then
            age = age*secyer/(24*60*60)
         end if
            
         time_step = s% time_step ! in years
         if (trim(s% terminal_show_timestep_units) == 'seconds' .or. &
             trim(s% terminal_show_timestep_units) == 'secs') then
            time_step = time_step*secyer
         else if (trim(s% terminal_show_timestep_units) == 'days') then
            time_step = time_step*secyer/(24*60*60)
         end if
        
         if (s% terminal_show_log_age) age = safe_log10(age)
         if (s% terminal_show_log_dt) time_step = safe_log10(time_step)

         model = s% model_number
         nz = s% nz

         ierr = 0         
         
         Xmax = dot_product(s% dq(1:nz), s% xa(s% species,1:nz))
         
         if (s% u_flag) then
            v = s% u(1)
         else if (s% v_flag) then
            v = s% v(1)
         else
            v = 0 ! s% r(1)*s% dlnR_dt(1)
         end if
         vsurf_div_csound = v / s% csound(1)

         dt = s% time_step*secyer
         
         sum_Lnuc = s% power_nuc_burn
         sum_LH = s% power_h_burn
         sum_LHe = s% power_he_burn
         sum_Lphoto = abs(s% power_photo)
         sum_Lz = s% power_z_burn 
         
         if (s% id == 1) then
            id_str = ''
         else
            write(id_str,'(i3)') s% id
         end if
         
         fmt1 = '(a3,i8,f11.6,'
         
         if (s% Teff < 1d4) then
            fmt2 = 'f11.3,'
         else
            fmt2 = '1pe11.3,0p,'
         end if
         
         if (s% star_mass >= 1d2) then
            fmt3 = '2f11.6,2(1pe11.3),0p,'
         else
            fmt3 = '4f11.6,'
         end if
         
         if (s% eta(s% nz) >= 1d3) then
            fmt4 = '3f11.6,e11.3,'
         else
            fmt4 = '3f11.6,f11.6,'
         end if
         fmt5 = '2i7)'
         
         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3) // trim(fmt4) // trim(fmt5)
         !write(*,*) 'fmt line1 ' // trim(fmt)
         write(io,fmt=fmt) &
            id_str, model, &
            s% log_max_temperature, &   ! fmt1
            s% Teff, &   ! fmt2
            safe_log10(sum_LH), & ! fmt3
            safe_log10(sum_Lnuc), &
            s% star_mass, &            
            s% star_mass - max(s% he_core_mass, s% co_core_mass), &
            s% center_h1, & ! fmt4
            s% center_n14, &
            s% surface_he3 + s% surface_he4, &
            s% eta(s% nz), &
            s% nz, & ! fmt5
            s% num_retries
         
         tmp = max(0d0, min(1d0, 1 - (s% surface_h1 + s% surface_he3 + s% surface_he4)))
         if (s% initial_z >= 1d-5) then
            fmt1 = '(1pe11.4, 0p, 9f11.6, '
         else
            tmp = safe_log10(tmp)
            fmt1 = '(1pe11.4, 0p, 8f11.6, e11.2, '
         end if
         if (s% gam(s% nz) >= 1d3) then
            fmt2 = 'e11.3, '
         else
            fmt2 = 'f11.6, '
         end if
         fmt3 = ' i7)'
         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3)
         !write(*,*) 'fmt line2 ' // trim(fmt)
         iters = s% num_solver_iterations
         write(io,fmt=fmt) &
            time_step,  &
            s% log_center_temperature, &
            s% log_surface_radius, &
            safe_log10(sum_LHe), &
            safe_log10(abs(s% power_neutrinos)), &
            safe_log10(abs(s% star_mdot)), &
            s% he_core_mass, &
            s% center_he3 + s% center_he4, &
            s% center_o16, &
            tmp, &
            s% gam(s% nz), &
            iters
         
         if (s% why_Tlim <= 0) then
            why = ''
         else
            why = dt_why_str(min(numTlim,s% why_Tlim))
            if (s% why_Tlim == Tlim_dX .and. s% Tlim_dX_species > 0 &
                     .and. s% Tlim_dX_species <= s% species) then
               why = trim(dt_why_str(s% why_Tlim)) // ' ' // &
                  trim(chem_isos% name(s% chem_id(s% Tlim_dX_species)))
            else if (s% why_Tlim == Tlim_dX_div_X .and. s% Tlim_dX_div_X_species > 0 &
                     .and. s% Tlim_dX_div_X_species <= s% species) then
               why = trim(dt_why_str(s% why_Tlim)) // ' ' // &
                  trim(chem_isos% name(s% chem_id(s% Tlim_dX_div_X_species)))
            else if (s% why_Tlim ==  Tlim_dX_nuc_drop .and. s% dX_nuc_drop_max_j > 0 &
                     .and. s% dX_nuc_drop_max_j <= s% species) then
               why = trim(dt_why_str(s% why_Tlim)) // ' ' // &
                  trim(chem_isos% name(s% chem_id(s% dX_nuc_drop_max_j)))
            else if (s% why_Tlim ==  Tlim_dlgL_nuc_cat) then 
               if (s% Tlim_dlgL_nuc_category > 0 &
                     .and. s% Tlim_dlgL_nuc_category <= num_categories ) then
                  why = trim(category_name(s% Tlim_dlgL_nuc_category))
               else
                  why = '???'
               end if
            end if         
         end if         
         
         if (s% net_iso(isi28) == 0) then
            tmp = 1 - (s% center_h1 + s% center_he3 + s% center_he4)
            fmt = '(1pe11.4, 0p, 5f11.6, 0p4f11.6, 0p, e11.3, a14)'
            !write(*,*) 'fmt line3 ' // trim(fmt)
            write(io,fmt) &
               age, &
               s% log_center_density, &
               s% log_surface_luminosity, &
               safe_log10(sum_Lz), &
               safe_log10(abs(s% power_photo)), &
               s% lnd(1)/ln10, &
               s% co_core_mass, &
               s% center_c12, &
               s% center_ne20, &
               tmp, &
               vsurf_div_csound, &
               trim(why)
         else
            tmp = s% center_si28
            fmt = '(1pe11.4, 0p, 5f11.6, 0p4f11.6, 0p, e11.3, a14)'
            !write(*,*) 'fmt line3 ' // trim(fmt)
            write(io,fmt) &
               age, &
               s% log_center_density, &
               s% log_surface_luminosity, &
               safe_log10(sum_Lz), &
               safe_log10(abs(s% power_photo)), &
               s% lnd(1)/ln10, &
               s% co_core_mass, &
               s% center_c12, &
               s% center_ne20, &
               tmp, &
               vsurf_div_csound, &
               trim(why)
         end if
         
         call show_trace_history_values(max(0, s% num_trace_history_values))
         write(io,*)
         
         s% just_wrote_terminal_header = .false.
         
         
         contains

         
         subroutine show_trace_history_values(num)
            use history, only: get_history_specs, get_history_values, get1_hist_value
            integer, intent(in) :: num
            real(dp) :: values(num)
            integer :: int_values(num), specs(num)
            logical :: is_int_value(num)
            logical :: failed_to_find_value(num)
            real(dp) :: val
            integer :: i
            include 'formats'
            call get_history_specs(s, num, s% trace_history_value_name, specs, .false.)
            call get_history_values( &
               s, num, specs, is_int_value, int_values, values, failed_to_find_value)
            do i = 1, num
               if (failed_to_find_value(i)) then
                  if (.not. get1_hist_value(s, s% trace_history_value_name(i), val)) then
                     cycle
                  end if
                  values(i) = val
                  if (is_bad(values(i))) then   
                     stop 'show_trace_history_values bad from get1_hist_value'
                  end if
               else if (is_int_value(i)) then
                  write(io,'(a40,i14)') &
                     trim(s% trace_history_value_name(i)), int_values(i)
                  cycle
               end if
               if ((values(i) == 0) .or. &
                        (abs(values(i)) > 1d-4 .and. abs(values(i)) < 1d4)) then
                  write(io,'(a40,99(f26.16))') &
                     trim(s% trace_history_value_name(i)), values(i)
               else
                  write(io,'(a40,99(1pd26.16))') &
                     trim(s% trace_history_value_name(i)), values(i)
                  if (is_bad(values(i))) then   
                     stop 'show_trace_history_values'
                  end if
               end if
            end do
         end subroutine show_trace_history_values


      end subroutine output_terminal_summary
      
      
      logical function get_history_info(s, do_write)
         type (star_info), pointer :: s
         logical, intent(in) :: do_write                  
         integer :: model
         logical :: write_history, write_terminal         
         include 'formats'
         model = s% model_number
         if (s% history_interval > 0) then
            write_history = (mod(model, s% history_interval) == 0) .or. do_write
         else
            write_history = .false.
         end if
         if (s% terminal_interval > 0) then
            write_terminal = (mod(model, s% terminal_interval) == 0) .or. do_write
         else
            write_terminal = .false.
         end if
         get_history_info = write_history .or. write_terminal         
         if (.not. get_history_info) return
         if (s% write_header_frequency*s% terminal_interval > 0) then
            if ( mod(model, s% write_header_frequency*s% terminal_interval) .eq. 0 &
                 .and. .not. s% doing_first_model_of_run) then
               write(*,*)
               call write_terminal_header(s)
            endif
         end if         
         if (write_terminal) call do_terminal_summary(s)  
         if (write_history) s% need_to_update_history_now = .true.   
               
      end function get_history_info
      
        
      integer function do_bare_bones_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         logical :: logged
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            do_bare_bones_check_model = terminate
            return
         end if
         logged = get_history_info( s, .false. )
         do_bare_bones_check_model = keep_going
      end function do_bare_bones_check_model
      
      
      subroutine save_profile(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call set_save_profiles_info(s, phase_priority)
      end subroutine save_profile

        
      integer function do_check_limits(id)
         use rates_def
         use chem_def
         use chem_lib, only: chem_get_iso_id
         use star_utils
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr, i, j, k, cid, k_burn, k_omega, nz, max_abs_vel_loc, &
            period_number, max_period_number
         real(dp) :: log_surface_gravity, v_div_csound_max, bound_mass, &
            power_nuc_burn, power_h_burn, power_he_burn, power_z_burn, logQ, max_logQ, min_logQ, &
            envelope_fraction_left, avg_x, v_surf, csound_surf, delta_nu, v_surf_div_v_esc, &
            ratio, dt_C, peak_burn_vconv_div_cs, min_pgas_div_p, v_surf_div_v_kh, GREKM_avg_abs, &
            max_omega_div_omega_crit, omega_div_omega_crit, log_Teff, Lnuc_div_L, max_abs_vel, &
            species_mass_for_min_limit, species_mass_for_max_limit
            
         include 'formats'
         
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            do_check_limits = terminate
            return
         end if
         
         if (s% RSP_flag) then
            max_period_number = s% RSP_max_num_periods
            period_number = s% RSP_num_periods
            GREKM_avg_abs = s% rsp_GREKM_avg_abs
         else
            max_period_number = 0
            period_number = -1
            GREKM_avg_abs = 1d99
         end if

         nz = s% nz
         do_check_limits = keep_going
         
         species_mass_for_min_limit = get_species_mass(s% star_species_mass_min_limit_iso)
         species_mass_for_max_limit = get_species_mass(s% star_species_mass_max_limit_iso)

         csound_surf = s% csound_face(1)
         if (s% u_flag) then
            v_surf =  s% u(1)
            v_div_csound_max = maxval(abs(s% u(1:nz)/s% csound(1:nz)))
         else if (s% v_flag) then
            v_surf = s% v(1)
            v_div_csound_max = maxval(abs(s% v(1:nz)/s% csound_face(1:nz)))
         else
            v_surf = abs(s% r(1) * s% dlnR_dt(1))
            v_div_csound_max = 0d0
         end if
         
         bound_mass = get_bound_mass(s)/Msun
         
         if(s%u_flag) then
            max_abs_vel_loc = maxloc(abs(s%u(1:nz)),dim=1)
            max_abs_vel = s%u(max_abs_vel_loc)
         else if (s%v_flag) then
            max_abs_vel_loc = maxloc(abs(s%v(1:nz)),dim=1)
            max_abs_vel = s%v(max_abs_vel_loc)
         else
            max_abs_vel_loc = -1
            max_abs_vel = 0d0
         end if
         
         
         if (s% photosphere_r > 0d0) then
            v_surf_div_v_kh = abs(v_surf)*s% kh_timescale/s% photosphere_r
            if (s% cgrav(1) <= 0d0) then
               v_surf_div_v_esc = 0d0
            else
               v_surf_div_v_esc = v_surf/sqrt(2*s% cgrav(1)*s% m(1)/(s% photosphere_r*Rsun))
            end if
         else
            v_surf_div_v_kh = 0d0
            v_surf_div_v_esc = 0d0
         end if
         
         log_surface_gravity = s% log_surface_gravity
         power_nuc_burn = s% power_nuc_burn
         power_h_burn = s% power_h_burn
         power_he_burn = s% power_he_burn
         power_z_burn = s% power_z_burn
         log_Teff = safe_log10(s% Teff)
         if (s% L_phot > 0d0) then
            Lnuc_div_L = s% L_nuc_burn_total / s% L_phot
            if (.not. s% get_delta_nu_from_scaled_solar) then
               delta_nu = 1d6/(2*s% photosphere_acoustic_r) ! microHz
            else
               delta_nu = &
                  s% delta_nu_sun*sqrt(s% star_mass)*pow3(s% Teff/s% Teff_sun) / &
                     pow(s% L_phot,0.75d0)
            end if
         else
            Lnuc_div_L = 0d0
            delta_nu = 0d0
         end if
         
         if (s% dxdt_nuc_factor > 0d0) then
            k = maxloc(s% eps_nuc(1:nz), dim=1)
            peak_burn_vconv_div_cs = s% conv_vel(k)/s% csound(k)
         else
            peak_burn_vconv_div_cs = 0d0
         end if
         
         if (s% initial_mass > s% he_core_mass) then
            envelope_fraction_left = &
               (s% star_mass - s% he_core_mass)/(s% initial_mass - s% he_core_mass)
         else
            envelope_fraction_left = 1
         end if
         
         max_logQ = -99
         min_logQ = 100
         do k = 1, s% nz
            logQ = s% lnd(k)/ln10 - 2*s% lnT(k)/ln10 + 12
            if (s% lnT(k)/ln10 < 5.5d0) then ! only worry about lower T cases
               if (logQ > max_logQ) max_logQ = logQ
            end if
            if (logQ < min_logQ) min_logQ = logQ
         end do
         
         min_pgas_div_p = 1d99
         do k = s% nz, 1, -1
            if (s% q(k) > s% Pgas_div_P_limit_max_q) exit
            if (s% pgas(k)/s% p(k) < min_pgas_div_p) min_pgas_div_p = s% pgas(k)/s% p(k)
         end do
         
         max_omega_div_omega_crit = 0; k_omega = 0
         if (s% rotation_flag .and. s% omega_div_omega_crit_limit > 0) then
            do k = 1, s% nz
               omega_div_omega_crit = abs(s% omega(k))/omega_crit(s,k) 
               if (omega_div_omega_crit > max_omega_div_omega_crit) then
                  k_omega = k
                  max_omega_div_omega_crit = omega_div_omega_crit
               end if
            end do
         end if
         
         if(s% max_abs_rel_run_E_err > 0d0) then
            if (abs(s% cumulative_energy_error/s% total_energy) > s% max_abs_rel_run_E_err &
                  .and. .not. s% doing_relax) then
               write(*, '(/,a,/, 2e20.10)') &
                  'stop because abs rel_run_E_err exceeds limit (max_abs_rel_run_E_err) ', &
                  abs(s% cumulative_energy_error/s% total_energy), s% max_abs_rel_run_E_err
               do_check_limits = terminate
               s% termination_code = t_max_abs_rel_run_E_err
               s% result_reason = abs_rel_run_E_err
               return
            end if
         end if
         
         if (max_abs_vel > clight) then
            write(*, '(/,a,/, I5,1X,2e20.10)') &
               'retry because maximum velocity exceeds speed of light ',max_abs_vel_loc,max_abs_vel,max_abs_vel/clight
            do_check_limits = retry
            return
         end if

         if (peak_burn_vconv_div_cs > 0.75d0*s% peak_burn_vconv_div_cs_limit) then
            write(*,1) 'peak_burn_vconv_div_cs: ', &
               peak_burn_vconv_div_cs / s% peak_burn_vconv_div_cs_limit, &
               peak_burn_vconv_div_cs, s% peak_burn_vconv_div_cs_limit
            k = maxloc(s% eps_nuc(1:nz), dim=1)
            write(*,2) 'maxloc eps_nuc', k, s% conv_vel(k), s% csound(k), s% eps_nuc(k)
            stop 'test do_one_utils'
         end if
         
         if (s% fe_core_infall < s% fe_core_infall_limit .and. &
             s% fe_core_infall > 0.99d0*s% fe_core_infall_limit) &
            write(*,1) 'nearing fe_core_infall limit', &
               s% fe_core_infall, s% fe_core_infall_limit
         
         if (s% non_fe_core_infall < s% non_fe_core_infall_limit .and. &
             s% non_fe_core_infall > 0.99d0*s% non_fe_core_infall_limit) &
            write(*,1) 'nearing non_fe_core_infall limit', &
               s% non_fe_core_infall, s% non_fe_core_infall_limit
         
         if (s% non_fe_core_rebound > 0.99d0*s% non_fe_core_rebound_limit) &
            write(*,1) 'nearing non_fe_core_rebound limit', &
               s% non_fe_core_rebound, s% non_fe_core_rebound_limit
         
         if (max_omega_div_omega_crit > 0.75d0*s% omega_div_omega_crit_limit .and. &
               s% omega_div_omega_crit_limit > 0 .and. k_omega > 0) &
            write(*,2) 'omega_div_omega_crit', k_omega, &
               max_omega_div_omega_crit, s% omega_div_omega_crit_limit, &
               s% m(k_omega)/Msun, s% r_equatorial(k_omega)/Rsun, &
               s% omega(k_omega), &
               sqrt(s% cgrav(k_omega)*s% m(k_omega)/ pow3(s% r_equatorial(k_omega)))
         
         if (s% star_age >= s% max_age .and. s% max_age > 0) then 
            call compare_to_target('star_age >= max_age', s% star_age, s% max_age, &
                  t_max_age)
                  
         else if (s% time >= s% max_age_in_seconds .and. s% max_age_in_seconds > 0) then 
            call compare_to_target('time >= max_age_in_seconds', &
               s% time, s% max_age_in_seconds, t_max_age)
            
         else if (max_omega_div_omega_crit >= s% omega_div_omega_crit_limit .and. &
               s% omega_div_omega_crit_limit > 0) then 
            write(*, '(/,a,/, 2e20.10)') &
               'stop max_omega_div_omega_crit >= omega_div_omega_crit_limit', &
               max_omega_div_omega_crit, s% omega_div_omega_crit_limit
            do_check_limits = terminate
            s% termination_code = t_max_omega_div_omega_crit
            s% result_reason = result_reason_normal
                  
         else if (peak_burn_vconv_div_cs >= s% peak_burn_vconv_div_cs_limit) then 
            write(*, '(/,a,/, 2e20.10)') &
               'stop peak_burn_vconv_div_cs >= peak_burn_vconv_div_cs_limit', &
               peak_burn_vconv_div_cs, s% peak_burn_vconv_div_cs_limit
            do_check_limits = terminate
            s% termination_code = t_peak_burn_vconv_div_cs_limit
            s% result_reason = result_reason_normal
            
         else if (s% model_number >= s% max_model_number .and. s% max_model_number >= 0) then 
            write(*, '(/,a,/, 2i9)') 'stop because model_number >= max_model_number', &
               s% model_number, s% max_model_number
            do_check_limits = terminate
            s% termination_code = t_max_model_number
            s% result_reason = result_reason_normal
            
         else if (period_number >= max_period_number .and. max_period_number >= 0) then 
            write(*, '(/,a,/, 2i9)') 'stop because period_number >= max_period_number', &
               period_number, max_period_number
            do_check_limits = terminate
            s% termination_code = t_max_period_number
            s% result_reason = result_reason_normal
            
         else if (GREKM_avg_abs < s% RSP_GREKM_avg_abs_limit &
                  .and. s% RSP_GREKM_avg_abs_limit >= 0 &
                  .and. period_number >= 10) then 
            write(*, '(/,a,/, 2e20.10)') &
               'stop because GREKM_avg_abs < RSP_GREKM_avg_abs_limit', &
               GREKM_avg_abs, s% RSP_GREKM_avg_abs_limit
            do_check_limits = terminate
            s% termination_code = 0
            s% result_reason = result_reason_normal
           
         else if (s% center_degeneracy >= s% eta_center_limit) then 
            call compare_to_target('center_degeneracy >= eta_center_limit', &
               s% center_degeneracy, s% eta_center_limit, t_eta_center_limit)
            
         else if (s% log_center_temperature >= s% log_center_temp_limit) then 
            call compare_to_target('log_center_temperature >= log_center_temp_limit', &
               s% log_center_temperature, s% log_center_temp_limit, t_log_center_temp_limit)
            
         else if (s% log_center_temperature <= s% log_center_temp_lower_limit) then 
            call compare_to_target('log_center_temperature <= log_center_temp_lower_limit', &
               s% log_center_temperature, s% log_center_temp_lower_limit, &
               t_log_center_temp_lower_limit)
            
         else if (s% max_entropy >= s% max_entropy_limit) then 
            call compare_to_target('max_entropy >= max_entropy_limit', &
               s% max_entropy, s% max_entropy_limit, t_max_entropy_limit)
            
         else if (s% max_entropy <= s% max_entropy_lower_limit) then 
            call compare_to_target('max_entropy <= max_entropy_lower_limit', &
               s% max_entropy, s% max_entropy_lower_limit, &
               t_max_entropy_lower_limit)
            
         else if (s% center_entropy >= s% center_entropy_limit) then 
            call compare_to_target('center_entropy >= center_entropy_limit', &
               s% center_entropy, s% center_entropy_limit, t_center_entropy_limit)
            
         else if (s% center_entropy <= s% center_entropy_lower_limit) then 
            call compare_to_target('center_entropy <= center_entropy_lower_limit', &
               s% center_entropy, s% center_entropy_lower_limit, &
               t_center_entropy_lower_limit)
            
         else if (s% log_center_density <= s% log_center_density_lower_limit) then 
            call compare_to_target('log_center_density <= log_center_density_lower_limit', &
               s% log_center_density, s% log_center_density_lower_limit, &
               t_log_center_density_lower_limit)
            
         else if (s% log_center_density >= s% log_center_density_limit) then 
            call compare_to_target('log_center_density >= log_center_density_limit', &
               s% log_center_density, s% log_center_density_limit, t_log_center_density_limit)
            
         else if (s% center_gamma > s% gamma_center_limit) then 
            call compare_to_target('center_gamma > gamma_center_limit', &
               s% center_gamma, s% gamma_center_limit, t_gamma_center_limit)
            
         else if (s% log_max_temperature >= s% log_max_temp_upper_limit) then 
            call compare_to_target('log_max_temperature >= log_max_temp_upper_limit', &
               s% log_max_temperature, s% log_max_temp_upper_limit, t_log_max_temp_upper_limit)
            
         else if (s% log_max_temperature <= s% log_max_temp_lower_limit) then 
            call compare_to_target('log_max_temperature <= log_max_temp_lower_limit', &
               s% log_max_temperature, s% log_max_temp_lower_limit, t_log_max_temp_lower_limit)
            
         else if (s% center_he4 < s% HB_limit .and. s% center_h1 < 1d-4) then 
            call compare_to_target('center he4 < HB_limit', s% center_he4, s% HB_limit, t_HB_limit)
            
         else if (s% star_mass_min_limit > 0 .and. s% star_mass <= s% star_mass_min_limit) then 
            call compare_to_target('star_mass <= star_mass_min_limit', &
               s% star_mass, s% star_mass_min_limit, t_star_mass_min_limit)
            
         else if (s% star_mass_max_limit > 0 .and. s% star_mass >= s% star_mass_max_limit) then 
            call compare_to_target('star_mass >= star_mass_max_limit', &
               s% star_mass, s% star_mass_max_limit, t_star_mass_max_limit)
            
         else if (s% bound_mass_min_limit > 0 .and. bound_mass <= s% bound_mass_min_limit) then 
            call compare_to_target('bound_mass <= bound_mass_min_limit', &
               bound_mass, s% bound_mass_min_limit, t_bound_mass_min_limit)
            
         else if (s% bound_mass_max_limit > 0 .and. bound_mass >= s% bound_mass_max_limit) then 
            call compare_to_target('bound_mass >= bound_mass_max_limit', &
               bound_mass, s% bound_mass_max_limit, t_bound_mass_max_limit)
            
         else if (species_mass_for_min_limit >= 0 .and. &
               species_mass_for_min_limit <= s% star_species_mass_min_limit) then 
            call compare_to_target( &
               trim(s% star_species_mass_min_limit_iso) // ' total mass <= star_species_mass_min_limit', &
               species_mass_for_min_limit, s% star_species_mass_min_limit, t_star_species_mass_min_limit)
            
         else if (species_mass_for_max_limit >= s% star_species_mass_max_limit) then 
            call compare_to_target( &
               trim(s% star_species_mass_max_limit_iso) // ' total mass >= star_species_mass_max_limit', &
               species_mass_for_max_limit, s% star_species_mass_max_limit, t_star_species_mass_max_limit)
                        
         else if (s% xmstar_min_limit > 0 .and. s% xmstar <= s% xmstar_min_limit) then 
            call compare_to_target('xmstar <= xmstar_min_limit', &
               s% xmstar, s% xmstar_min_limit, t_xmstar_min_limit)
            
         else if (s% xmstar_max_limit > 0 .and. s% xmstar >= s% xmstar_max_limit) then 
            call compare_to_target('xmstar >= xmstar_max_limit', &
               s% xmstar, s% xmstar_max_limit, t_xmstar_max_limit)
            
         else if (s% star_mass - s% he_core_mass < s% envelope_mass_limit) then 
            call compare_to_target('envelope mass < envelope_mass_limit', &
               s% star_mass - s% he_core_mass, s% envelope_mass_limit, &
               t_envelope_mass_limit)
            
         else if (envelope_fraction_left < s% envelope_fraction_left_limit) then 
            call compare_to_target('envelope_fraction_left < limit', &
               envelope_fraction_left, s% envelope_fraction_left_limit, &
               t_envelope_fraction_left_limit)
            
         else if (s% he_core_mass >= s% he_core_mass_limit) then 
            call compare_to_target('he_core_mass >= he_core_mass_limit', &
               s% he_core_mass, s% he_core_mass_limit, t_he_core_mass_limit)
            
         else if (s% co_core_mass >= s% co_core_mass_limit) then 
            call compare_to_target('co_core_mass >= co_core_mass_limit', &
               s% co_core_mass, s% co_core_mass_limit, t_co_core_mass_limit)
            
         else if (s% one_core_mass >= s% one_core_mass_limit) then 
            call compare_to_target('one_core_mass >= one_core_mass_limit', &
               s% one_core_mass, s% one_core_mass_limit, t_one_core_mass_limit)
            
         else if (s% fe_core_mass >= s% fe_core_mass_limit) then 
            call compare_to_target('fe_core_mass >= fe_core_mass_limit', &
               s% fe_core_mass, s% fe_core_mass_limit, t_fe_core_mass_limit)
            
         else if (s% neutron_rich_core_mass >= s% neutron_rich_core_mass_limit) then 
            call compare_to_target('neutron_rich_core_mass >= neutron_rich_core_mass_limit', &
               s% neutron_rich_core_mass, s% neutron_rich_core_mass_limit, t_neutron_rich_core_mass_limit)
            
         else if ( &
               s% he_core_mass >= s% co_core_mass .and. &
               s% co_core_mass > 0 .and. &
               s% center_he4 < 1d-4 .and. &
               s% he_core_mass - s% co_core_mass < s% he_layer_mass_lower_limit) then 
            call compare_to_target('he layer mass < he_layer_mass_lower_limit', &
               s% he_core_mass - s% co_core_mass, s% he_layer_mass_lower_limit, &
               t_he_layer_mass_lower_limit)
            
         else if (abs(safe_log10(power_h_burn) - s% log_surface_luminosity) <= &
                  s% abs_diff_lg_LH_lg_Ls_limit &
                  .and. s% abs_diff_lg_LH_lg_Ls_limit > 0) then 
            call compare_to_target('abs(lg_LH - lg_Ls) <= limit', &
               abs(safe_log10(power_h_burn) - s% log_surface_luminosity), &
                  s% abs_diff_lg_LH_lg_Ls_limit, t_abs_diff_lg_LH_lg_Ls_limit)

         else if (s% Teff <= s% Teff_lower_limit) then 
            call compare_to_target('Teff <= Teff_lower_limit', &
               s% Teff, s% Teff_lower_limit, t_Teff_lower_limit)
               
         else if (s% Teff >= s% Teff_upper_limit) then 
            call compare_to_target('Teff >= Teff_upper_limit', &
               s% Teff, s% Teff_upper_limit, t_Teff_upper_limit)

         else if (delta_nu <= s% delta_nu_lower_limit .and. s% delta_nu_lower_limit > 0) then 
            call compare_to_target('delta_nu <= delta_nu_lower_limit', &
               delta_nu, s% delta_nu_lower_limit, t_delta_nu_lower_limit)
               
         else if (delta_nu >= s% delta_nu_upper_limit .and. s% delta_nu_upper_limit > 0) then 
            call compare_to_target('delta_nu >= delta_nu_upper_limit', &
               delta_nu, s% delta_nu_upper_limit, t_delta_nu_upper_limit)

         else if (s% delta_Pg <= s% delta_Pg_lower_limit .and. s% delta_Pg_lower_limit > 0) then 
            call compare_to_target('delta_Pg <= delta_Pg_lower_limit', &
               s% delta_Pg, s% delta_Pg_lower_limit, t_delta_Pg_lower_limit)
               
         else if (s% delta_Pg >= s% delta_Pg_upper_limit .and. s% delta_Pg_upper_limit > 0) then 
            call compare_to_target('delta_Pg >= delta_Pg_upper_limit', &
               s% delta_Pg, s% delta_Pg_upper_limit, t_delta_Pg_upper_limit)
               
         else if (s% shock_mass >= s% shock_mass_upper_limit .and. &
               s% shock_mass_upper_limit > 0) then 
            call compare_to_target('shock_mass >= shock_mass_upper_limit', &
               s% shock_mass, s% shock_mass_upper_limit, t_shock_mass_upper_limit)
               
         else if (s% outer_mach1_mass >= s% mach1_mass_upper_limit .and. &
               s% mach1_mass_upper_limit > 0) then 
            call compare_to_target('mach1_mass >= mach1_mass_upper_limit', &
               s% outer_mach1_mass, s% mach1_mass_upper_limit, t_mach1_mass_upper_limit)

         else if (s% photosphere_m - s% M_center/Msun <= s% photosphere_m_sub_M_center_limit) then 
            call compare_to_target( &
               'photosphere_m - M_center/Msun <= photosphere_m_sub_M_center_limit', &
               s% photosphere_m - s% M_center/Msun, &
               s% photosphere_m_sub_M_center_limit, &
               t_photosphere_m_sub_M_center_limit)

         else if (s% photosphere_m <= s% photosphere_m_lower_limit) then 
            call compare_to_target('photosphere_m <= photosphere_m_lower_limit', &
               s% photosphere_m, s% photosphere_m_lower_limit, t_photosphere_m_lower_limit)
               
         else if (s% photosphere_m >= s% photosphere_m_upper_limit) then 
            call compare_to_target('photosphere_m >= photosphere_m_upper_limit', &
               s% photosphere_m, s% photosphere_m_upper_limit, t_photosphere_m_upper_limit)

         else if (s% photosphere_r <= s% photosphere_r_lower_limit) then 
            call compare_to_target('photosphere_r <= photosphere_r_lower_limit', &
               s% photosphere_r, s% photosphere_r_lower_limit, t_photosphere_r_lower_limit)
               
         else if (s% photosphere_r >= s% photosphere_r_upper_limit) then 
            call compare_to_target('photosphere_r >= photosphere_r_upper_limit', &
               s% photosphere_r, s% photosphere_r_upper_limit, t_photosphere_r_upper_limit)

         else if (log_Teff <= s% log_Teff_lower_limit) then 
            call compare_to_target('log_Teff <= log_Teff_lower_limit', &
               log_Teff, s% log_Teff_lower_limit, t_log_Teff_lower_limit)
               
         else if (log_Teff >= s% log_Teff_upper_limit) then 
            call compare_to_target('log_Teff >= log_Teff_upper_limit', &
               log_Teff, s% log_Teff_upper_limit, t_log_Teff_upper_limit)

         else if (s% log_surface_temperature <= s% log_Tsurf_lower_limit) then 
            call compare_to_target('log_surface_temperature <= log_Tsurf_lower_limit', &
               s% log_surface_temperature, s% log_Tsurf_lower_limit, t_log_Tsurf_lower_limit)
               
         else if (s% log_surface_temperature >= s% log_Tsurf_upper_limit) then 
            call compare_to_target('log_surface_temperature >= log_Tsurf_upper_limit', &
               s% log_surface_temperature, s% log_Tsurf_upper_limit, t_log_Tsurf_upper_limit)

         else if (s% log_surface_radius <= s% log_Rsurf_lower_limit) then 
            call compare_to_target('log_surface_radius <= log_Rsurf_lower_limit', &
               s% log_surface_radius, s% log_Rsurf_lower_limit, t_log_Rsurf_lower_limit)
               
         else if (s% log_surface_radius >= s% log_Rsurf_upper_limit) then 
            call compare_to_target('log_surface_radius >= log_Rsurf_upper_limit', &
               s% log_surface_radius, s% log_Rsurf_upper_limit, t_log_Rsurf_upper_limit)

         else if (s% log_surface_pressure <= s% log_Psurf_lower_limit) then 
            call compare_to_target('log_surface_pressure <= log_Psurf_lower_limit', &
               s% log_surface_pressure, s% log_Psurf_lower_limit, t_log_Psurf_lower_limit)
               
         else if (s% log_surface_pressure >= s% log_Psurf_upper_limit) then 
            call compare_to_target('log_surface_pressure >= log_Psurf_upper_limit', &
               s% log_surface_pressure, s% log_Psurf_upper_limit, t_log_Psurf_upper_limit)

         else if (s% log_surface_density <= s% log_Dsurf_lower_limit) then 
            call compare_to_target('log_surface_density <= log_Dsurf_lower_limit', &
               s% log_surface_density, s% log_Dsurf_lower_limit, t_log_Dsurf_lower_limit)
               
         else if (s% log_surface_density >= s% log_Dsurf_upper_limit) then 
            call compare_to_target('log_surface_density >= log_Dsurf_upper_limit', &
               s% log_surface_density, s% log_Dsurf_upper_limit, t_log_Dsurf_upper_limit)

         else if (s% log_surface_luminosity <= s% log_L_lower_limit) then 
            call compare_to_target('log_surface_luminosity <= log_L_lower_limit', &
               s% log_surface_luminosity, s% log_L_lower_limit, t_log_L_lower_limit)
               
         else if (s% log_surface_luminosity >= s% log_L_upper_limit) then 
            call compare_to_target('log_surface_luminosity >= log_L_upper_limit', &
               s% log_surface_luminosity, s% log_L_upper_limit, t_log_L_upper_limit)

         else if (log_surface_gravity <= s% log_g_lower_limit) then 
            call compare_to_target('log_surface_gravity <= log_g_lower_limit', &
               log_surface_gravity, s% log_g_lower_limit, t_log_g_lower_limit)
            
         else if (log_surface_gravity >= s% log_g_upper_limit) then 
            call compare_to_target('log_surface_gravity >= log_g_upper_limit', &
               log_surface_gravity, s% log_g_upper_limit, t_log_g_upper_limit)

         else if (power_nuc_burn >= s% power_nuc_burn_upper_limit) then 
            call compare_to_target('power_nuc_burn >= power_nuc_burn_upper_limit', &
               power_nuc_burn, s% power_nuc_burn_upper_limit, t_power_nuc_burn_upper_limit)

         else if (power_h_burn >= s% power_h_burn_upper_limit) then 
            call compare_to_target('power_h_burn >= power_h_burn_upper_limit', &
               power_h_burn, s% power_h_burn_upper_limit, t_power_h_burn_upper_limit)

         else if (power_he_burn >= s% power_he_burn_upper_limit) then 
            call compare_to_target('power_he_burn >= power_he_burn_upper_limit', &
               power_he_burn, s% power_he_burn_upper_limit, t_power_he_burn_upper_limit)

         else if (power_z_burn >= s% power_z_burn_upper_limit) then 
            call compare_to_target('power_z_burn >= power_z_burn_upper_limit', &
               power_z_burn, s% power_z_burn_upper_limit, t_power_z_burn_upper_limit)

         else if (power_nuc_burn < s% power_nuc_burn_lower_limit) then 
            call compare_to_target('power_nuc_burn < power_nuc_burn_lower_limit', &
               power_nuc_burn, s% power_nuc_burn_lower_limit, t_power_nuc_burn_lower_limit)

         else if (power_h_burn < s% power_h_burn_lower_limit) then 
            call compare_to_target('power_h_burn < power_h_burn_lower_limit', &
               power_h_burn, s% power_h_burn_lower_limit, t_power_h_burn_lower_limit)

         else if (power_he_burn < s% power_he_burn_lower_limit) then 
            call compare_to_target('power_he_burn < power_he_burn_lower_limit', &
               power_he_burn, s% power_he_burn_lower_limit, t_power_he_burn_lower_limit)

         else if (power_z_burn < s% power_z_burn_lower_limit) then 
            call compare_to_target('power_z_burn < power_z_burn_lower_limit', &
               power_z_burn, s% power_z_burn_lower_limit, t_power_z_burn_lower_limit)

         else if (s% center_Ye < s% center_Ye_lower_limit) then 
            call compare_to_target('center_Ye < center_Ye_lower_limit', &
               s% center_Ye, s% center_Ye_lower_limit, t_center_Ye_lower_limit)

         else if (s% R_center < s% center_R_lower_limit) then 
            call compare_to_target('R_center < center_R_lower_limit', &
               s% R_center, s% center_R_lower_limit, t_center_R_lower_limit)

         else if (s% fe_core_infall > s% fe_core_infall_limit) then 
            if (abs(s% error_in_energy_conservation/s% total_energy_end) < &
                  s% hard_limit_for_rel_error_in_energy_conservation) then
               do_check_limits = terminate
               s% result_reason = result_reason_normal
               s% termination_code = t_fe_core_infall_limit
               write(*, '(/,a,/, 99e20.10)') &
                  'stop because fe_core_infall > fe_core_infall_limit', &
                  s% fe_core_infall, s% fe_core_infall_limit
            else
               write(*,2) 'rel_E_err too large for fe_core_infall termination', &
                  s% model_number, s% error_in_energy_conservation/abs(s% total_energy_end)
            end if

         else if (s% non_fe_core_infall > s% non_fe_core_infall_limit) then 
            call compare_to_target('non_fe_core_infall > non_fe_core_infall_limit', &
               s% non_fe_core_infall, s% non_fe_core_infall_limit, t_non_fe_core_infall_limit)

         else if (s% non_fe_core_rebound > s% non_fe_core_rebound_limit) then 
            call compare_to_target('non_fe_core_rebound > non_fe_core_rebound_limit', &
               s% non_fe_core_rebound, s% non_fe_core_rebound_limit, t_non_fe_core_rebound_limit)

         else if (v_surf/csound_surf > s% v_div_csound_surf_limit) then 
            call compare_to_target('v_surf/csound_surf > v_div_csound_surf_limit', &
               v_surf/csound_surf, s% v_div_csound_surf_limit, t_v_div_csound_surf_limit)

         else if (v_div_csound_max > s% v_div_csound_max_limit) then 
            call compare_to_target('v_div_csound_max > v_div_csound_max_limit', &
               v_div_csound_max, s% v_div_csound_max_limit, t_v_div_csound_max_limit)

         else if (s% min_gamma1 < s% gamma1_limit) then 
            call compare_to_target('min_gamma1 < gamma1_limit', &
               s% min_gamma1, s% gamma1_limit, t_gamma1_limit)            

         else if (min_pgas_div_p < s% Pgas_div_P_limit) then 
            call compare_to_target('min_pgas_div_p < Pgas_div_P_limit', &
               min_pgas_div_p, s% Pgas_div_P_limit, t_Pgas_div_P_limit)            

         else if (Lnuc_div_L <= s% Lnuc_div_L_lower_limit) then 
            call compare_to_target('Lnuc_div_L <= Lnuc_div_L_lower_limit', &
               Lnuc_div_L, s% Lnuc_div_L_lower_limit, t_Lnuc_div_L_lower_limit)
               
         else if (Lnuc_div_L >= s% Lnuc_div_L_upper_limit) then 
            call compare_to_target('Lnuc_div_L >= Lnuc_div_L_upper_limit', &
               Lnuc_div_L, s% Lnuc_div_L_upper_limit, t_Lnuc_div_L_upper_limit)

         else if (v_surf_div_v_kh <= s% v_surf_div_v_kh_lower_limit) then 
            call compare_to_target('v_surf_div_v_kh <= v_surf_div_v_kh_lower_limit', &
               v_surf_div_v_kh, s% v_surf_div_v_kh_lower_limit, t_v_surf_div_v_kh_lower_limit)
               
         else if (v_surf_div_v_kh >= s% v_surf_div_v_kh_upper_limit) then 
            call compare_to_target('v_surf_div_v_kh >= v_surf_div_v_kh_upper_limit', &
               v_surf_div_v_kh, s% v_surf_div_v_kh_upper_limit, t_v_surf_div_v_kh_upper_limit)
               
         else if (v_surf_div_v_esc >= s% v_surf_div_v_esc_limit) then 
            call compare_to_target('v_surf_div_v_esc >= v_surf_div_v_esc_limit', &
               v_surf_div_v_esc, s% v_surf_div_v_esc_limit, t_v_surf_div_v_esc_limit)
               
         else if (v_surf*1d-5 >= s% v_surf_kms_limit) then 
            call compare_to_target('v_surf_kms >= v_surf_kms_limit', &
               v_surf*1d-5, s% v_surf_kms_limit, t_v_surf_kms_limit)
               
         else if (s% cumulative_extra_heating >= &
                  s% stop_when_reach_this_cumulative_extra_heating .and.&
                  s% stop_when_reach_this_cumulative_extra_heating > 0) then 
            call compare_to_target( &
               'cumulative_extra_heating >= limit', &
               s% cumulative_extra_heating, &
               s% stop_when_reach_this_cumulative_extra_heating, &
               t_cumulative_extra_heating_limit)
               
         else if (s% stop_near_zams .and. &
                  Lnuc_div_L >= s% Lnuc_div_L_zams_limit) then
            do_check_limits = terminate
            s% termination_code = t_Lnuc_div_L_zams_limit
            s% result_reason = result_reason_normal
            write(*, '(/,a,/, 99e20.10)') &
               'stop because Lnuc_div_L >= Lnuc_div_L_zams_limit', Lnuc_div_L, s% Lnuc_div_L_zams_limit
               
         else if (s% stop_at_phase_PreMS .and. s% phase_of_evolution == phase_PreMS) then
            do_check_limits = terminate
            s% termination_code = t_phase_PreMS
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_PreMS'
               
         else if (s% stop_at_phase_ZAMS .and. s% phase_of_evolution == phase_ZAMS) then
            do_check_limits = terminate
            s% termination_code = t_phase_ZAMS
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_ZAMS'
               
         else if (s% stop_at_phase_IAMS .and. s% phase_of_evolution == phase_IAMS) then
            do_check_limits = terminate
            s% termination_code = t_phase_IAMS
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_IAMS'
               
         else if (s% stop_at_phase_TAMS .and. s% phase_of_evolution == phase_TAMS) then
            do_check_limits = terminate
            s% termination_code = t_phase_TAMS
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_TAMS'
               
         else if (s% stop_at_phase_He_Burn .and. s% phase_of_evolution == phase_He_Burn) then
            do_check_limits = terminate
            s% termination_code = t_phase_He_Burn
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_He_Burn'
               
         else if (s% stop_at_phase_ZACHeB .and. s% phase_of_evolution == phase_ZACHeB) then
            do_check_limits = terminate
            s% termination_code = t_phase_ZACHeB
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_ZACHeB'
               
         else if (s% stop_at_phase_TACHeB .and. s% phase_of_evolution == phase_TACHeB) then
            do_check_limits = terminate
            s% termination_code = t_phase_TACHeB
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_TACHeB'
               
         else if (s% stop_at_phase_TP_AGB .and. s% phase_of_evolution == phase_TP_AGB) then
            do_check_limits = terminate
            s% termination_code = t_phase_TP_AGB
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_TP_AGB'
               
         else if (s% stop_at_phase_C_Burn .and. s% phase_of_evolution == phase_C_Burn) then
            do_check_limits = terminate
            s% termination_code = t_phase_C_Burn
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_C_Burn'
               
         else if (s% stop_at_phase_Ne_Burn .and. s% phase_of_evolution == phase_Ne_Burn) then
            do_check_limits = terminate
            s% termination_code = t_phase_Ne_Burn
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_Ne_Burn'
               
         else if (s% stop_at_phase_O_Burn .and. s% phase_of_evolution == phase_O_Burn) then
            do_check_limits = terminate
            s% termination_code = t_phase_O_Burn
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_O_Burn'
               
         else if (s% stop_at_phase_Si_Burn .and. s% phase_of_evolution == phase_Si_Burn) then
            do_check_limits = terminate
            s% termination_code = t_phase_Si_Burn
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_Si_Burn'
               
         else if (s% stop_at_phase_WDCS .and. s% phase_of_evolution == phase_WDCS) then
            do_check_limits = terminate
            s% termination_code = t_phase_WDCS
            s% result_reason = result_reason_normal
            write(*, '(/,a,/)') 'stop because phase_of_evolution == phase_WDCS'

         end if
                  
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_central_limits
            if (s% xa_central_lower_limit(j) <= 0) cycle
            if (len_trim(s% xa_central_lower_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_central_lower_limit_species(j))
            if (cid <= 0) then
               write(*,*)
               write(*,2) '<' // trim(s% xa_central_lower_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_central_lower_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = center_avg_x(s,i)
            if (avg_x < s% xa_central_lower_limit(j)) then
               call compare_to_target('have dropped below central lower limit for ' // &
                     trim(s% xa_central_lower_limit_species(j)), &
                     avg_x, s% xa_central_lower_limit(j), t_xa_central_lower_limit)
               exit
            end if
         end do
         
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_central_limits
            if (s% xa_central_upper_limit(j) <= 0) cycle
            if (s% xa_central_upper_limit(j) >= 1) cycle
            if (len_trim(s% xa_central_upper_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_central_upper_limit_species(j))
            if (cid <= 0) then
               !cycle
               write(*,*)
               write(*,2) '<' // trim(s% xa_central_upper_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_central_upper_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = center_avg_x(s,i)
            if (avg_x > s% xa_central_upper_limit(j)) then
               call compare_to_target('have risen above central upper limit for ' // &
                     trim(s% xa_central_upper_limit_species(j)), &
                     avg_x, s% xa_central_upper_limit(j), t_xa_central_upper_limit)
               exit
            end if
         end do
         
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_surface_limits
            if (s% xa_surface_lower_limit(j) <= 0) cycle
            if (len_trim(s% xa_surface_lower_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_surface_lower_limit_species(j))
            if (cid <= 0) then
               write(*,*)
               write(*,2) '<' // trim(s% xa_surface_lower_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_surface_lower_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = surface_avg_x(s,i)
            if (avg_x < s% xa_surface_lower_limit(j)) then
               call compare_to_target('have dropped below surface lower limit for ' // &
                     trim(s% xa_surface_lower_limit_species(j)), &
                     avg_x, s% xa_surface_lower_limit(j), t_xa_surface_lower_limit)
               exit
            end if
         end do
         
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_surface_limits
            if (s% xa_surface_upper_limit(j) <= 0) cycle
            if (s% xa_surface_upper_limit(j) >= 1) cycle
            if (len_trim(s% xa_surface_upper_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_surface_upper_limit_species(j))
            if (cid <= 0) then
               !cycle
               write(*,*)
               write(*,2) '<' // trim(s% xa_surface_upper_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_surface_upper_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = surface_avg_x(s,i)
            if (avg_x > s% xa_surface_upper_limit(j)) then
               call compare_to_target('have risen above surface upper limit for ' // &
                     trim(s% xa_surface_upper_limit_species(j)), &
                     avg_x, s% xa_surface_upper_limit(j), t_xa_surface_upper_limit)
               exit
            end if
         end do
         
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_average_limits
            if (s% xa_average_lower_limit(j) <= 0) cycle
            if (len_trim(s% xa_average_lower_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_average_lower_limit_species(j))
            if (cid <= 0) then
               !cycle
               write(*,*)
               write(*,2) '<' // trim(s% xa_average_lower_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_average_lower_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = dot_product(s% dq(1:nz), s% xa(i,1:nz))
            if (avg_x < s% xa_average_lower_limit(j)) then
               call compare_to_target('have dropped below average lower limit for ' // &
                     trim(s% xa_average_lower_limit_species(j)), &
                     avg_x, s% xa_average_lower_limit(j), t_xa_average_lower_limit)
               exit
            end if
         end do
         
         if (do_check_limits /= keep_going) return
         
         do j=1,num_xa_average_limits
            if (s% xa_average_upper_limit(j) <= 0) cycle
            if (s% xa_average_upper_limit(j) >= 1) cycle
            if (len_trim(s% xa_average_upper_limit_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_average_upper_limit_species(j))
            if (cid <= 0) then
               write(*,*)
               write(*,2) '<' // trim(s% xa_average_upper_limit_species(j)) // '>'
               write(*,2) 'is invalid for xa_average_upper_limit_species', j
               write(*,*)
               do_check_limits = terminate
               return
            end if
            i = s% net_iso(cid)
            if (i == 0) cycle
            avg_x = dot_product(s% dq(1:nz), s% xa(i,1:nz))
            if (avg_x > s% xa_average_upper_limit(j)) then
               call compare_to_target('have risen above average upper limit for ' // &
                     trim(s% xa_average_upper_limit_species(j)), &
                     avg_x, s% xa_average_upper_limit(j), t_xa_average_upper_limit)
               exit
            end if
         end do
         
         
         contains
         
         
         subroutine compare_to_target(str, value, target_value, termination_code)
            character (len=*), intent(in) :: str
            real(dp), intent(in) :: value, target_value
            integer, intent(in) :: termination_code
            real(dp) :: err
            include 'formats'
            err = abs(value - target_value)/ &
               (s% when_to_stop_atol + s% when_to_stop_rtol*max(abs(value),abs(target_value)))
            if (err > 1) then
               do_check_limits = redo
               s% dt = 0.5d0*s% dt
               write(*,'(/,a,5e20.10)') &
                  'redo with smaller timestep to get closer to stopping target '  // trim(str), &
                  value, target_value
            else
               do_check_limits = terminate
               s% result_reason = result_reason_normal
               s% termination_code = termination_code
               write(*, '(/,a,/, 99e20.10)') 'stop because ' // trim(str), value, target_value
            end if
         end subroutine compare_to_target
         
         
         real(dp) function get_species_mass(str) ! Msun
            use chem_lib, only: chem_get_iso_id
            character(len=*), intent(in) :: str
            integer :: id, j
            get_species_mass = -1d0
            id = chem_get_iso_id(str)
            if (id > 0) then
               j = s% net_iso(id)
               if (j > 0) then
                  get_species_mass = dot_product(s% dm(1:s% nz),s% xa(j,1:s% nz))/Msun
               end if
            end if
         end function get_species_mass
         
         
      end function do_check_limits

        
      integer function do_one_check_model(id)
         use rates_def, only: i_rate
         use chem_def, only: i_burn_c
         use star_utils, only: update_time, total_times
         integer, intent(in) :: id
         
         logical :: must_do_profile
         real(dp), parameter :: log_he_temp = 7.8d0
         real(dp), parameter :: d_tau_min = 1d-2, d_tau_max = 1d0
         real(dp), parameter :: little_step_factor = 10d0, little_step_size = 10d0
         real(dp) :: v, surf_grav, power_he_burn, power_z_burn, &
            power_neutrinos
         integer :: model, profile_priority, ierr
         integer, parameter :: tau_ramp = 50
         type (star_info), pointer :: s
         logical :: logged
         integer :: nz
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            do_one_check_model = terminate
            return
         end if

         nz = s% nz
         must_do_profile = time_to_profile(s) ! also updates phase_of_evolution
         profile_priority = delta_priority
         model = s% model_number
         do_one_check_model = keep_going
         
         do_one_check_model = do_check_limits(id)
         if (do_one_check_model /= keep_going) then
            if (dbg) write(*,*) 'do_check_limits /= keep_going'
            write(*,*)
            must_do_profile = .true.
         end if
         
         if (s% u_flag) then
            v = s% u(1)
         else if (s% v_flag) then
            v = s% v(1)
         else
            v = s% r(1) * s% dlnR_dt(1)
         end if
         
         power_he_burn = s% power_he_burn
         power_z_burn = s% power_z_burn
         power_neutrinos = s% power_neutrinos
         
         if (must_do_profile) profile_priority = phase_priority
         
         logged = get_history_info(s, must_do_profile)

         if (logged .and. s% write_profiles_flag) then
            if (s% model_number .eq. s% profile_model &
               .or. (s% profile_interval > 0 .and. &
                     (s% doing_first_model_of_run .or. &
                     mod(s% model_number,s% profile_interval) == 0))) then
               if (s% write_profiles_flag) must_do_profile = .true.
               if (s% model_number == s% profile_model .or.&
                   s% doing_first_model_of_run .or. &
                   (mod(s% model_number, s% priority_profile_interval) == 0)) then
                  profile_priority = phase_priority
               end if
            end if
            if ( must_do_profile ) then
               if (dbg) write(*,*) 'do_one_check_model: call set_save_profiles_info'
               call set_save_profiles_info(s, profile_priority)
            end if
         end if
         
      end function do_one_check_model
      
      
      logical function time_to_profile(s)
         use chem_def, only: ih1
         use star_utils
         type (star_info), pointer :: s
         ! end-of-run and helium break-even are always done. 
         ! this function decides on other models to be profiled.
         real(dp), parameter :: center_he_drop = 1d-2, surface_t_drop = 4d-2
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         time_to_profile = .false.
         
!         select case ( s% phase_of_evolution )
!         case ( phase_starting )
!            if ( arrived_main_seq(s) ) then
!               if (dbg) write(*,*) 'arrived_main_seq'
!               time_to_profile = .true.
!               s% prev_tsurf = s% log_surface_temperature
!               if (s% center_h1 > center_h_going) then
!                  !s% phase_of_evolution = phase_early_main_seq
!               else if ( s% center_h1 > center_h_gone ) then
!                  !s% phase_of_evolution = phase_mid_main_seq
!               else if ( s% center_he4 > center_he_going ) then
!                  !s% phase_of_evolution = phase_he_ignition_over
!               else
!                  !s% phase_of_evolution = phase_helium_burning
!               end if
!            end if
!         case ( phase_early_main_seq )
!            if ( s% center_h1 < center_h_going ) then
!               time_to_profile = .true.
!               s% prev_tsurf = s% log_surface_temperature
!               !s% phase_of_evolution = phase_mid_main_seq
!            end if
!         case ( phase_mid_main_seq )
!            if ( s% center_h1 < center_h_gone &
!                  .and. s% log_surface_temperature < s% prev_tsurf-surface_t_drop ) then
!               time_to_profile = .true.
!               !s% phase_of_evolution = phase_wait_for_he
!            end if
!         case ( phase_wait_for_he )
!         case ( phase_he_igniting ) ! for non-flash ignition of helium core
!            if ( s% center_he4 <= s% ignition_center_xhe-center_he_drop &
!                  .and. s% log_surface_luminosity > s% prev_luminosity ) then
!               time_to_profile = .true.
!               !s% phase_of_evolution = phase_he_ignition_over
!               s% prev_tcntr2 = s% prev_tcntr1; s% prev_age2 = s% prev_age1
!               s% prev_tcntr1 = s% log_center_temperature; s% prev_age1 = s% star_age
!               if ( s% log_surface_luminosity > s% he_luminosity_limit ) &
!                  s% he_luminosity_limit = s% log_surface_luminosity
!            end if
!            s% prev_luminosity = s% log_surface_luminosity
!         case ( phase_he_ignition_over )
!            if ( s% center_he4 < center_he_going ) then
!               time_to_profile = .true.
!               !s% phase_of_evolution = phase_helium_burning
!            end if      
!            s% prev_tcntr2 = s% prev_tcntr1; s% prev_age2 = s% prev_age1
!            s% prev_tcntr1 = s% log_center_temperature; s% prev_age1 = s% star_age
!         case ( phase_carbon_burning )
!         case ( phase_helium_burning )
!         end select
!         
      end function time_to_profile
      
      
      subroutine dummy_before_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = id ! so that we use that arg
         ierr = 0
      end subroutine dummy_before_evolve


      
      end module do_one_utils
      
