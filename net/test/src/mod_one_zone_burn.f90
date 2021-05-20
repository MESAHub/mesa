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

      module mod_one_zone_support
      use chem_def
      use chem_lib
      use math_lib
      use net_def
      use net_lib
      use const_def, only: Qconv, secyer, kerg, avo, ln10
      use rates_def
      use utils_lib, only: mesa_error
      
      implicit none

      character(len=256) :: net_name

      character(len=256):: final_abundances_filename
      logical :: save_final_abundances

      character(len=256):: initial_abundances_filename
      logical :: read_initial_abundances

      character(len=256):: T_Rho_history_filename
      logical :: read_T_Rho_history

      
      character(len=256):: burn_filename
      real(dp) :: burn_tend, burn_rho, burn_temp, &
         burn_rtol, burn_atol, burn_P, burn_xmin, burn_xmax, &
         eps, odescal, stptry
      logical :: trace, burn_dbg, use_pivoting
      real(dp) :: min_for_show_peak_abundances
      integer :: max_num_for_show_peak_abundances
      
      integer, parameter :: max_num_burn_isos_to_show = 1000
      character(len=iso_name_length) :: names_of_isos_to_show(max_num_burn_isos_to_show)
      integer :: num_names_of_isos_to_show

      integer, parameter :: max_num_isos_for_Xinit = 1000
      character(len=iso_name_length) :: names_of_isos_for_Xinit(max_num_isos_for_Xinit)
      real(dp) :: values_for_Xinit(max_num_isos_for_Xinit)
      integer :: num_isos_for_Xinit
      logical :: uniform_Xinit

      integer :: screening_mode
      integer, pointer :: which_rates(:)

      integer, parameter :: io_out = 35
      real(dp) :: data_output_min_t


      character (len=32) :: which_solver
      integer :: decsol_switch
      character (len=32) :: small_mtx_decsol, large_mtx_decsol

      logical :: show_net_reactions_info
      integer :: which_rates_choice
      
      real(dp) :: rattab_logT_lower_bound, rattab_logT_upper_bound

      character(len=256):: data_filename, data_heading_line
      character (len=64) :: net_file, cache_suffix
      
      integer :: handle, eos_handle, net_handle
      type (Net_General_Info), pointer  :: g
      integer :: species, num_reactions
      integer, pointer :: reaction_id(:)
      integer, dimension(:), pointer :: net_iso, chem_id

      real(dp) :: z, abar, zbar, z2bar, z53bar, ye, eps_neu_total
      real(dp) :: eta, d_eta_dlnT, d_eta_dlnRho
      real(dp), dimension(:), pointer :: &
            xin, xin_copy, d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT
      real(dp), pointer :: d_dxdt_dx(:, :)  
      
      real(dp) :: weak_rate_factor

      integer :: max_steps ! maximal number of allowed steps.

      real(dp), pointer :: x_initial(:) ! (species)

      real(dp) :: burn_lnE, burn_lnS
      real(dp) :: burn_logT, burn_logRho, &
         burn_eta, burn_deta_dlnT, burn_Cv, burn_d_Cv_dlnT
      
      real(dp) :: T_prev, time_prev, eps_nuc_prev, eps_neu_prev, cp_prev
      real(dp), pointer :: x_previous(:) ! (species)

      real(dp), dimension(:), pointer :: peak_abundance, peak_time

      integer :: num_times_for_burn
      integer, parameter :: max_num_times = 10000
      real(dp), dimension(max_num_times) :: &
         times_for_burn, log10Ts_for_burn, log10Rhos_for_burn, &
         etas_for_burn, log10Ps_for_burn

      logical :: burn_at_constant_P, clip, show_peak_x_and_time, &
         burn_at_constant_density
      real(dp) :: starting_temp, pressure

      real(dp) :: max_step_size ! maximal step size.

      real(dp), pointer :: rate_factors(:) ! (num_reactions)
      integer, pointer :: net_reaction_ptr(:) 
      
      integer, parameter :: max_num_reactions_to_track = 100
      integer :: num_reactions_to_track
      character(len=maxlen_reaction_Name) :: &
         reaction_to_track(max_num_reactions_to_track)
      integer :: index_for_reaction_to_track(max_num_reactions_to_track)
      
      integer, parameter :: max_num_special_rate_factors = 100
      integer :: num_special_rate_factors
      real(dp) :: special_rate_factor(max_num_special_rate_factors)
      character(len=maxlen_reaction_Name) :: &
         reaction_for_special_factor(max_num_special_rate_factors)

      character (len=16) :: set_rate_c12ag, set_rate_n14pg, set_rate_3a, &
         set_rate_1212
      
      logical :: show_Qs, quiet, complete_silence_please, &
         show_ye_stuff
      
      real(dp) :: starting_logT
      
      logical :: dbg = .false.
      
      
      contains
      
      
      integer function burn_isos_for_Xinit(i)
         integer, intent(in) :: i 
         burn_isos_for_Xinit = chem_get_iso_id(names_of_isos_for_Xinit(i))
      end function burn_isos_for_Xinit
      
      
      integer function burn_isos_to_show(i)
         integer, intent(in) :: i 
         burn_isos_to_show = chem_get_iso_id(names_of_isos_to_show(i))
      end function burn_isos_to_show

      
      subroutine Do_One_Zone_Burn(net_file_in)
         use num_lib, only: solver_option
         use mtx_lib, only: decsol_option
         use eos_lib
         use mtx_def
         use interp_1d_lib, only: interp_pm
         use interp_1d_def, only: pm_work_size
         use test_net_support, only: Setup_eos
         use net_lib, only: get_net_reaction_table_ptr
         use rates_lib, only: rates_reaction_id
         use utils_lib, only: set_nan
         
         character (len=*), intent(in) :: net_file_in
         
         character (len=256) :: net_file
         real(dp) :: logRho, logT, Rho, T, xsum, &
           eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT
         real(dp), target :: eps_nuc_categories(num_categories)
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         integer :: i, j, ierr, iounit, decsol_choice, solver_choice, num_times

         real(dp), dimension(:), pointer :: times, dxdt_source_term
         real(dp), dimension(:), pointer :: &
            log10Ts_f1, log10Rhos_f1, etas_f1, log10Ps_f1
         real(dp), dimension(:,:), pointer :: &
            log10Ts_f, log10Rhos_f, etas_f, log10Ps_f
         
         ! args to control the solver -- see num/public/num_isolve.dek
         real(dp) :: h 
         ! absolute and relative error tolerances
         real(dp), pointer :: rtol(:) ! relative error tolerance (species)
         real(dp), pointer :: atol(:) ! absolute error tolerance (species)
         integer :: itol ! switch for rtol and atol
         
         real(dp), pointer :: burn_work_array(:), net_work_array(:)
         
         real(dp), pointer :: ending_x(:) ! (species)
         integer :: nfcn    ! number of function evaluations
         integer :: njac    ! number of jacobian evaluations
         integer :: nstep   ! number of computed steps
         integer :: naccpt  ! number of accepted steps
         integer :: nrejct  ! number of rejected steps
         integer :: max_order_used
         
         integer :: iout, caller_id, cid, ir, burn_lwork, net_lwork
         integer(8) :: time0, time1, clock_rate

         real(dp) :: ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS, dt
         real(dp) :: ending_log10T, starting_log10T, avg_eps_nuc, ending_eps_neu_total
         real(dp) :: time_doing_net, time_doing_eos, told
         real(dp) :: xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction

         integer, parameter :: nwork = pm_work_size
         real(dp), pointer :: pm_work(:)
         
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo
         type (Net_General_Info), pointer :: g
         
         include 'formats'
         
         ierr = 0
         told = 0
         
         net_file = net_file_in
         netinfo => netinfo_target

         call test_net_setup(net_file)

         net_handle = handle
         call get_net_ptr(net_handle, g, ierr)
         if (ierr /= 0) return

         call Setup_eos(eos_handle)
!         g% max_rate_times_dt = max_rate_times_dt
         
         logT = burn_logT
         T = burn_temp     
         logRho = burn_logRho
         Rho = burn_rho
         
         if (read_T_Rho_history) then
            num_times = max_num_times
         else if (num_times_for_burn <= 0) then
            num_times = 1
         else
            num_times = num_times_for_burn
         end if

         if (num_names_of_isos_to_show < 0) num_names_of_isos_to_show = species
         
         allocate( &
            rate_factors(num_reactions), rtol(species), atol(species), &
            x_initial(species), x_previous(species), ending_x(species), times(num_times), &
            log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
            etas_f1(4*num_times), log10Ps_f1(4*num_times), &
            peak_abundance(species), peak_time(species), &
            stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed for Do_One_Zone_Burn'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get_net_reaction_table_ptr(net_handle, net_reaction_ptr, ierr)
         if (ierr /= 0) then
            write(*,*) 'bad net?  get_net_reaction_table_ptr failed'
            return
         end if

         rate_factors(:) = 1

         if (num_special_rate_factors > 0) then
            do i=1,num_special_rate_factors
               if (len_trim(reaction_for_special_factor(i)) == 0) cycle
               ir = rates_reaction_id(reaction_for_special_factor(i))
               j = 0
               if (ir > 0) j = net_reaction_ptr(ir)
               if (j <= 0) cycle
               rate_factors(j) = special_rate_factor(i)
               write(*,1) 'set special rate factor for ' // &
                     trim(reaction_for_special_factor(i)), special_rate_factor(i)
            end do
         end if
         
         if (num_reactions_to_track > 0) then
            do i=1,num_reactions_to_track
               index_for_reaction_to_track(i) = 0
               if (len_trim(reaction_to_track(i)) == 0) cycle
               ir = rates_reaction_id(reaction_to_track(i))
               j = 0
               if (ir > 0) j = net_reaction_ptr(ir)
               if (j <= 0) cycle
               index_for_reaction_to_track(i) = j
               write(*,1) 'track rate ' // trim(reaction_to_track(i))
            end do
         end if
         
         log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
         log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
         etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)
         log10Ps_f(1:4,1:num_times) => log10Ps_f1(1:4*num_times)
         
         peak_abundance(:) = 0

         xin = 0
         eta = 0
         
         iout = 1
         itol = 0         
         
         rtol(:) = burn_rtol
         atol(:) = burn_atol
         
         xin = 0
         if (read_initial_abundances) then
            call read_X(ierr)
            if (ierr /= 0) return
         else if (uniform_Xinit) then
            xin(:) = 0.5d0/(species-1)
            j = net_iso(ih1)
            if (j <= 0) stop 'where is the h?'
            xin(j) = 0.5d0
         else
            do i = 1, num_isos_for_Xinit
               cid = burn_isos_for_Xinit(i)
               j = net_iso(cid)
               if (j <= 0 .or. j > species) then
                  write(*,*) 'bad names_of_isos_for_Xinit ' // &
                        trim(names_of_isos_for_Xinit(i))
                  call mesa_error(__FILE__,__LINE__)
               end if
               xin(j) = values_for_Xinit(i)
            end do
         end if
         
         !xin(:) = xin(:)/sum(xin(:))
         
         if (read_T_Rho_history) then
            call do_read_T_Rho_history(ierr)
            if (ierr /= 0) return
         end if
         
         if (num_times_for_burn <= 0) then
            times(1) = burn_tend
            log10Ts_f(1,1) = logT
            log10Ts_f(2:4,1) = 0
            log10Rhos_f(1,1) = logRho
            log10Rhos_f(2:4,1) = 0
            etas_f(1:1,1) = burn_eta
            etas_f(2:4,1) = 0
         else
            times(1:num_times) = times_for_burn(1:num_times)
            log10Ts_f(1,1:num_times) = log10Ts_for_burn(1:num_times)
            log10Ts_f(2:4,1:num_times) = 0
            log10Rhos_f(1,1:num_times) = log10Rhos_for_burn(1:num_times)
            log10Rhos_f(2:4,1:num_times) = 0
            etas_f(1,1:num_times) = etas_for_burn(1:num_times)
            etas_f(2:4,1:num_times) = 0
            logRho = log10Rhos_for_burn(1)
            burn_rho = exp10(logRho)
            logT = log10Ts_for_burn(1)
            burn_temp = exp10(logT)
            burn_tend = times_for_burn(num_times)
            if (.false.) then
               write(*,*)
               do i=1,num_times
                  write(*,2) 'history', i, times_for_burn(i), log10Ts_for_burn(i), &
                     log10Rhos_for_burn(i), etas_for_burn(i)
               end do
               write(*,*)
            end if
         end if
         starting_logT = logT
         
         h = 1d-2*burn_tend ! 1d-14
         !write(*,1) 'h', h
         !stop
         
         x_initial(1:species) = xin(1:species)
         x_previous(1:species) = xin(1:species)
         caller_id = 0
         dxdt_source_term => null()
         
         if (.not. quiet) then
            write(*,*)
            write(*,*)
            write(*,1) 'one zone burn ' // trim(net_file)
            write(*,*)
            write(*,2) 'number of species', species
            write(*,*)
            write(*,1) 'temp', burn_temp
            write(*,1) 'logT', burn_logT
            write(*,1) 'rho', burn_rho
            write(*,1) 'logRho', burn_logRho
            write(*,1) 'eta', burn_eta
            write(*,*)
            write(*,1) 'tend', burn_tend
            write(*,1) 'tend/secyer', burn_tend/secyer
            write(*,*)
            write(*,1) 'initial abundances'
            call show_X(xin,.false.,.false.)
         end if
            
! data_heading_line was not set and writing out nulls. change it - fxt 
!         write(io_out,'(a)') trim(data_heading_line)
         write(data_heading_line,'(99(a,1pe14.6))') 'temp =',burn_temp,' rho =',burn_rho
         write(io_out,'(a)') trim(data_heading_line) 

         write(io_out,'(a7,99(a26,1x))',advance='no') &
            'i', &
            'signed_log_avg_eps_nuc', &
            'avg_eps_nuc', &
            'signed_log_eps_nuc', &
            'eps_nuc', &
            'eps_neu', &
            'delta_logT', &
            'logT', &
            'logRho', &
            'logPgas', &
            'logP', &
            'abar', &
            'zbar', &
            'entropy', &
            'logE', &
            'time', &
            'lg_time', &
            'lg_yrs', &
            'dt', &
            'lg_dt', &
            'ye', &
            'xsum_sub_1'
      
         do i=1,num_names_of_isos_to_show
            if (num_names_of_isos_to_show < species) then
               cid = burn_isos_to_show(i)
            else
               if (i > species) exit
               cid = chem_id(i)
            endif
            j = net_iso(cid)
            if (j == 0) cycle
            write(io_out,'(a26,1x)',advance='no') 'lg_' // trim(chem_isos% name(cid))
         end do
         do i=1,num_names_of_isos_to_show
            if (num_names_of_isos_to_show < species) then
               cid = burn_isos_to_show(i)
            else
               if (i > species) exit
               cid = chem_id(i)
            endif
            j = net_iso(cid)
            if (j == 0) cycle
            write(io_out,'(a26,1x)',advance='no') trim(chem_isos% name(cid))
         end do
         do i=1,num_reactions_to_track
            j = index_for_reaction_to_track(i)
            if (j == 0) cycle
            write(io_out,'(a26,1x)',advance='no')  'eps_' // trim(reaction_to_track(i))
            write(io_out,'(a26,1x)',advance='no')  'raw_' // trim(reaction_to_track(i))
            write(io_out,'(a26,1x)',advance='no')  'scrn_' // trim(reaction_to_track(i))
         end do
         write(io_out,*) 
         
         if (show_net_reactions_info) then
            write(*,'(a)') ' species'
            do j=1,species
               cid = chem_id(j)
               write(*,'(a)') chem_isos% name(cid)
            end do
            do j=1,species
               cid = chem_id(j)
               write(*,'(3i5,3x,a)') j, &
                  chem_isos% Z(cid), &
                  chem_isos% N(cid), &
                  chem_isos% name(cid)
            end do
            write(*,*)

            call show_net_reactions_and_info(handle, 6, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in show_net_reactions'
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*,*)
         end if
         
         if (.not. quiet) then
            write(*,1) 'h', h
            write(*,1) 'max_step_size', max_step_size
            write(*,2) 'max_steps', max_steps
            write(*,2) 'screening_mode', screening_mode
            write(*,*)
         end if
         
         if (species >= decsol_switch) then
            decsol_choice = decsol_option(large_mtx_decsol, ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR: unknown large_mtx_decsol ' // trim(large_mtx_decsol)
               return
            end if
         else
            decsol_choice = decsol_option(small_mtx_decsol, ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR: unknown small_mtx_decsol ' // trim(small_mtx_decsol)
               return
            end if
         end if
         
         solver_choice = solver_option(which_solver, ierr)
         if (ierr /= 0) then
            write(*,*) 'ERROR: unknown value for which_solver ' // trim(which_solver)
            return
         end if
         
         nullify(pm_work)
         
         call system_clock(time0,clock_rate)
         time_doing_net = -1
         time_doing_eos = -1
            
         if (burn_at_constant_density) then
            
            burn_lwork = net_1_zone_burn_const_density_work_size(handle,ierr)
            if (ierr /= 0) return
            net_lwork = net_work_size(handle,ierr)
            if (ierr /= 0) return
            allocate(net_work_array(net_lwork), burn_work_array(burn_lwork))
            
            call set_nan(burn_work_array) ! TESTING
            
            starting_log10T = burn_logT
            logT = burn_logT
            logRho = burn_logRho
            
            call net_1_zone_burn_const_density( &
               net_handle, eos_handle, species, num_reactions, &
               0d0, burn_tend, xin, logT, logRho, &
               get_eos_info_for_burn_at_const_density, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode,  &
               stptry, max_steps, eps, odescal, &
               use_pivoting, trace, burn_dbg, burn_finish_substep, &
               burn_lwork, burn_work_array, net_lwork, net_work_array, &
               ending_x, eps_nuc_categories, ending_log10T, &
               avg_eps_nuc, ending_eps_neu_total, &
               nfcn, njac, nstep, naccpt, nrejct, ierr)
            if (ierr == 0 .and. .not. quiet) then
               write(*,*)
               write(*,1) 'constant log10Rho', burn_logRho
               write(*,1) 'starting log10T', starting_log10T
               write(*,1) 'ending log10T', ending_log10T
               write(*,1) 'change in log10T', ending_log10T - starting_log10T
               write(*,1) 'change in T', exp10(ending_log10T/starting_log10T)
               write(*,1) 'avg_eps_nuc', avg_eps_nuc
               write(*,1) 'ending_eps_neu_total', ending_eps_neu_total
            end if
            deallocate(burn_work_array)

         else if (burn_at_constant_P) then
            if (num_times_for_burn <= 0) then
               log10Ps_f(1,1) = log10(pressure)
               log10Ps_f(2:4,1:num_times) = 0
            else
               log10Ps_f(1,1:num_times) = log10Ps_for_burn(1:num_times)
               log10Ps_f(2:4,1:num_times) = 0
               pressure = log10Ps_f(1,1)
               starting_temp = log10Ts_f(1,1)
            end if
            if (.not. quiet) then
               write(*,1) 'pressure', pressure
               write(*,1) 'starting_temp', starting_temp
            end if
            if (num_times > 1) then ! create interpolant
               allocate(pm_work(num_times*nwork))
               call interp_pm(times, num_times, log10Ps_f1, nwork, pm_work, 'net_1_zone_burn', ierr)
               if (ierr /= 0) stop 'failed in interp for logTs'
            end if
            
            call net_1_zone_burn_const_P( &
               net_handle, eos_handle, species, num_reactions, &
               solver_choice, starting_temp, xin, clip, &
               num_times, times, log10Ps_f1, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, screening_mode, &
               h, max_step_size, max_steps, rtol, atol, itol, burn_xmin, burn_xmax, &
               decsol_choice, caller_id, burn_solout, iout, &
               ending_x, ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS, &
               nfcn, njac, nstep, naccpt, nrejct, time_doing_net, time_doing_eos, ierr)
            if (ierr == 0 .and. .not. quiet) then
               write(*,*)
               write(*,1) 'pressure', pressure
               write(*,1) 'starting_temp', starting_temp
               write(*,1) 'ending_temp', ending_temp
               write(*,1) 'initial_rho', initial_rho
               write(*,1) 'ending_rho', ending_rho
               write(*,1) 'initial_entropy', exp(initial_lnS)/(avo*kerg)
               write(*,1) 'ending_entropy', exp(ending_lnS)/(avo*kerg)
            end if
         else
            if (num_times > 1) then ! create interpolants
               allocate(pm_work(num_times*nwork))
               call interp_pm(times, num_times, log10Ts_f1, nwork, pm_work, 'net_1_zone_burn', ierr)
               if (ierr /= 0) stop 'failed in interp for logTs'
               call interp_pm(times, num_times, log10Rhos_f1, nwork, pm_work, 'net_1_zone_burn', ierr)
               if (ierr /= 0) stop 'failed in interp for logRhos'
               call interp_pm(times, num_times, etas_f1, nwork, pm_work, 'net_1_zone_burn', ierr)
               if (ierr /= 0) stop 'failed in interp for etas'
            end if
            
            burn_lwork = net_1_zone_burn_work_size(handle,ierr)
            if (ierr /= 0) return
            net_lwork = net_work_size(handle,ierr)
            if (ierr /= 0) return
            allocate(net_work_array(net_lwork), burn_work_array(burn_lwork))
            
            call set_nan(burn_work_array) ! TESTING
            
            call net_1_zone_burn( &
               net_handle, eos_handle, species, num_reactions, 0d0, burn_tend, xin, &
               num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode,  & 
               stptry, max_steps, eps, odescal, &
               use_pivoting, trace, burn_dbg, burn_finish_substep, &
               burn_lwork, burn_work_array, & 
               net_lwork, net_work_array, & 
               ending_x, eps_nuc_categories, &
               avg_eps_nuc, eps_neu_total, &
               nfcn, njac, nstep, naccpt, nrejct, ierr)
               
            deallocate(net_work_array, burn_work_array)
            
         end if
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate

         if (ierr /= 0) then
            write(*,*) 'net_1_zone_burn returned ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (.not. quiet) then
            write(*,*)
            write(*,*)
         end if
         
         if (.not. complete_silence_please) then
            if (nstep >= max_steps) then
               write(*,2) 'hit max number of steps', nstep
               stop 'burn'
            end if
            write(*,2) 'number of species', species
            write(*,1) 'large final abundances', min_for_show_peak_abundances
            call show_X(ending_x(:),show_peak_x_and_time,.true.)
            if (show_ye_stuff) then
               call basic_composition_info( &
                     species, chem_id, ending_x, xh, xhe, z, abar, zbar, z2bar, z53bar, ye, &
                     mass_correction, xsum)
               write(*,1) 'changes in abundances'
               do i=1,species
                  x_previous(i) = ending_x(i) - x_initial(i)
               end do
               min_for_show_peak_abundances = -1d99
               call show_X(x_previous(:),.false.,.true.)
               write(*,*)
               write(*,1) 'ye', ye
               write(*,*)
               if (burn_at_constant_density) then
                  write(*,*)
                  write(*,1) 'constant log10Rho', burn_logRho
                  write(*,1) 'starting log10T', starting_log10T
                  write(*,1) 'ending log10T', ending_log10T
                  write(*,1) 'change in log10T', ending_log10T - starting_log10T
                  write(*,1) 'change in T', exp10(ending_log10T) - exp10(starting_log10T)
                  write(*,1) 'ending_eps_neu_total', ending_eps_neu_total
                  write(*,1) 'avg_eps_nuc', avg_eps_nuc
                  write(*,*)
                  write(*,1) 'tolerance eps', eps
                  write(*,1) 'odescal', odescal
               else if (.not. burn_at_constant_P) then
                  write(*,1) 'ending eps_neu_total', eps_neu_total
                  write(*,1) 'avg_eps_nuc', avg_eps_nuc
                  write(*,*)
                  write(*,1) 'logT', logT
                  write(*,1) 'logRho', logRho
                  write(*,*)
                  write(*,1) 'burn_temp', burn_temp
                  write(*,1) 'burn_rho', burn_rho
                  write(*,*)
                  write(*,1) 'burn_rtol', burn_rtol
                  write(*,1) 'burn_atol', burn_atol
               else
                  write(*,1) 'burn_temp', burn_temp
                  write(*,1) 'burn_rho', burn_rho
                  write(*,*)
                  write(*,1) 'burn_rtol', burn_rtol
                  write(*,1) 'burn_atol', burn_atol
               end if
               write(*,*)
               write(*,1) 'burn_tend (seconds)', burn_tend
               write(*,*)
               write(*,1) trim(net_name)
            end if
            write(*,*)
         end if

         if (.not. quiet) then
            write(*,2) 'nfcn', nfcn
            write(*,2) 'njac', njac
            write(*,2) 'naccpt', naccpt
            write(*,2) 'nrejct', nrejct
            write(*,*)
            write(*,2) 'number of steps', nstep
            write(*,*)
            write(*,*)
            write(*,1) 'output file ' // trim(data_filename)
            write(*,*)
            write(*,'(/,a30,99f18.3,/)') 'runtime (seconds)', dt
            write(*,*)
         end if
         
         if (save_final_abundances) then
            if (.not. quiet) write(*,*) 'save final abundances to ' // trim(final_abundances_filename)
            open(newunit=iounit, file=trim(final_abundances_filename), iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open final_abundances_filename ' // trim(final_abundances_filename)
            else
               write(iounit,2) 'species', species
               do j = 1, species
                  cid = chem_id(j)
                  write(iounit,3) trim(chem_isos% name(cid)), &
                     chem_isos% Z(cid), &
                     chem_isos% N(cid), &
                     max(0d0,ending_x(j))
               end do
               close(iounit)
            end if
         end if
         
         if (associated(pm_work)) deallocate(pm_work)
         deallocate( &
            rate_factors, rtol, atol, &
            ending_x, x_initial, x_previous, times, &
            log10Ts_f1, log10Rhos_f1, &
            etas_f1, log10Ps_f1, &
            peak_abundance, peak_time)
         
         
         contains
         
         
         subroutine get_eos_info_for_burn_at_const_density( &
               eos_handle, species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
            use eos_lib, only: eosDT_get
            use eos_def
            integer, intent(in) :: eos_handle, species
            integer, pointer :: chem_id(:) ! maps species to chem id
            integer, pointer :: net_iso(:) ! maps chem id to species number
            real(dp), intent(in) :: &
               xa(:), rho, logRho, T, logT
            real(dp), intent(out) :: &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT
            integer, intent(out) :: ierr

            real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
            real(dp) :: d_dxa(num_eos_d_dxa_results,species)
            
            include 'formats'
            ierr = 0
            
            Cv = burn_Cv
            d_Cv_dlnT = burn_d_Cv_dlnT
            eta = burn_eta
            d_eta_dlnT = burn_deta_dlnT
            
            if (burn_Cv <= 0d0 .or. burn_eta <= 0d0) then
               call eosDT_get( &
                  eos_handle, species, chem_id, net_iso, xa, &
                  Rho, logRho, T, logT, &
                  res, d_dlnd, d_dlnT, d_dxa, ierr)

               if (ierr /= 0) then
                  write(*,*) 'failed in eosDT_get'
                  return
               end if
               if (burn_Cv <= 0d0) then
                  Cv = res(i_cv)
                  d_Cv_dlnT = d_dlnT(i_cv)
               end if
               if (burn_eta <= 0d0) then
                  eta = res(i_eta)
                  d_eta_dlnT = d_dlnT(i_eta)
               end if
            end if
         
         end subroutine get_eos_info_for_burn_at_const_density


         subroutine burn_finish_substep(step, time, y, ierr)
            integer,intent(in) :: step
            real(dp), intent(in) :: time, y(:)
            real(dp), dimension(size(y)) :: x
            integer, intent(out) :: ierr 
            include 'formats'
            ierr = 0            
            !if (burn_at_constant_density) then
            !write(*,2) 'finish_substep time xh logT T', &
            !   step, time, y(1), y(species+1)/ln10, exp(y(species+1))
            !return
            !end if
            
            do i=1,species
               cid = chem_id(i)
               x(i) = y(i)*dble(chem_isos% Z_plus_N(cid))
            end do
            
            if (burn_at_constant_density) then
                x(species+1) = y(species+1)
                logT = x(species+1)/ln10
            end if
            call burn_solout1( &
               step, told, time, logT, logRho, species, x, ierr)
            told = time
         end subroutine burn_finish_substep
         

         real(dp) function interp_y(i, s, rwork_y, iwork_y, ierr)
            use const_def, only: dp
            integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
            real(dp), intent(in) :: s ! interpolation x value (between xold and x).
            real(dp), intent(inout), target :: rwork_y(*)
            integer, intent(inout), target :: iwork_y(*)
            integer, intent(out) :: ierr
            ierr = 0
            interp_y = 0
         end function interp_y


         subroutine do_read_T_Rho_history(ierr)
            integer, intent(out) :: ierr
            character (len=256) :: buffer, string
            integer :: i, n, iounit, t, num_isos, id, k
            
            include 'formats'
            
            ierr = 0
            write(*,*) 'read T Rho history from ' // trim(T_Rho_history_filename)
            open(newunit=iounit, file=trim(T_Rho_history_filename), &
               action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open'
               return
            end if
            
            read(iounit,*,iostat=ierr) num_times
            if (ierr /= 0) then
               write(*,*) 'first line should have num_times'
               close(iounit)
               return
            end if
            
            if (num_times > max_num_times) then
               write(*,3) 'num_times > max_num_times', num_times, max_num_times
               close(iounit)
               return
            end if
            
            if (.false.) write(*,2) 'num_times', num_times
            
            do i=1,num_times
               read(iounit,*,iostat=ierr) times_for_burn(i), log10Ts_for_burn(i), &
                  log10Rhos_for_burn(i), etas_for_burn(i)
               if (ierr /= 0) then
                  write(*,2) 'failed reading line', i+1
                  write(*,'(a)') 'line should have values for time, logT, logRho, and eta'
                  exit
               end if
               if (.false.) write(*,2) 'history', i, times_for_burn(i), log10Ts_for_burn(i), &
                  log10Rhos_for_burn(i), etas_for_burn(i)
            end do
         
            close(iounit)
            
            num_times_for_burn = num_times

         end subroutine do_read_T_Rho_history
         
         
         subroutine read_X(ierr)
            use utils_def
            use utils_lib
            integer, intent(out) :: ierr
            character (len=256) :: buffer, string
            integer :: i, n, iounit, t, num_isos, id, k
            
            include 'formats'
            
            write(*,*) 'read initial abundances from ' // trim(initial_abundances_filename)
            open(newunit=iounit, file=trim(initial_abundances_filename), &
               action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open'
               return
            end if
            
            n = 0
            i = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token .or. string /= 'species') then
               write(*,*) 'expect to find specification of number of species at start of file'
               ierr = -1; return
            end if
            t = token(iounit, n, i, buffer, string)
            read(string,fmt=*,iostat=ierr) num_isos
            if (t /= name_token .or. ierr /= 0) then
               write(*,*) 'expect to find specification of number of species at start of file'
               ierr = -1; return
            end if
            if (num_isos /= species) then
               write(*,2) 'expect to find number of species equal to those in current net', species
               ierr = -1; return
            end if
            do k = 1, species
               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  write(*,*) 'failed to find iso name at start of line: ' // trim(string)
                  ierr = -1; return
               end if
               id = get_nuclide_index(string)
               if (id <= 0) then
                  write(*,*) 'failed to recognize name of iso ' // trim(string)
                  ierr = -1
                  return
               end if
               j = net_iso(id)
               if (j <= 0 .or. j > species) then
                  write(*,*) 'iso not in current net ' // trim(string)
                  ierr = 1
                  return
               end if
               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  write(*,*) 'failed to read iso abundance: ' // &
                     trim(chem_isos% name(id)) // ' ' // trim(string)
                  ierr = -1; return
               end if
               read(string,fmt=*,iostat=ierr) xin(j)
               if (ierr /= 0) then
                  write(*,*) 'failed to read iso abundance: ' &
                     // trim(chem_isos% name(id)) // ' ' // trim(string)
               end if
            end do
            close(iounit)
         end subroutine read_X
         
   
         subroutine show_X(X,show_peak,do_sort)
            use num_lib, only: qsort
            real(dp) :: X(:)
            logical, intent(in) :: show_peak, do_sort 
            real(dp), target :: v_t(species) 
            integer, target :: index_t(species) 
            real(dp), pointer :: v(:) 
            integer, pointer :: index(:) 
            integer :: j
            real(dp) :: xsum
            include 'formats'
            v => v_t
            index => index_t

            if (do_sort) then
            
               v(1:species) = abs(x(1:species))
               call qsort(index, species, v)
               do i=1,species
                  if (i > max_num_for_show_peak_abundances) exit
                  j = index(species+1-i)
                  if (x(j) > min_for_show_peak_abundances) &
                     write(*,2) trim(chem_isos% name(chem_id(j))), i, x(j)
                  if (x(j) > 1.1d0 .or. x(j) < -0.1d0) then
                     write(*,1) 'bad x for ' // trim(chem_isos% name(chem_id(j))), x(j)
                     call mesa_error(__FILE__,__LINE__)
                  end if
               end do
               
            else

               do j=1, species
                  if (x(j) > min_for_show_peak_abundances) &
                     write(*,2) trim(chem_isos% name(chem_id(j))), j, x(j)
                  if (x(j) > 1.1d0 .or. x(j) < -0.1d0) then
                     write(*,1) 'bad x for ' // trim(chem_isos% name(chem_id(j))), x(j)
                     call mesa_error(__FILE__,__LINE__)
                  end if

                  !write(*,1) 'xin(net_iso(i' // trim(chem_isos% name(chem_id(j))) // ')=', x(j)

               end do

            end if
            !stop

            write(*,*)
            xsum = sum(x(1:species))
            write(*,1) 'xsum', xsum
            write(*,*)
            if (.not. show_peak) return
            write(*,*)
            write(*,1) 'peak x and time'
            do j=1, species
               if (peak_abundance(j) >= min_for_show_peak_abundances) &
                  write(*,1) trim(chem_isos% name(chem_id(j))), &
                     peak_abundance(j), peak_time(j)
            end do
            write(*,*)
         end subroutine show_X


         subroutine burn_solout( &
               step, told, time, n, x, rwork_y, iwork_y, interp_y, &
               lrpar, rpar, lipar, ipar, irtrn)
            use const_def
            use eos_def
            use eos_lib
            integer, intent(in) :: step, n, lrpar, lipar
            real(dp), intent(in) :: told, time
            real(dp), intent(inout) :: x(:)
            real(dp), intent(inout), target :: rwork_y(*)
            integer, intent(inout), target :: iwork_y(*)
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer :: i, cid
            interface
               include 'num_interp_y.dek'
            end interface
            integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
            
            call burn_solout1( &
               step, told, time, rpar(r_burn_prev_lgT), rpar(r_burn_prev_lgRho), n, x, irtrn)
         end subroutine burn_solout


         subroutine burn_solout1( &
               step, told, time, logT_in, logRho_in, n, x, irtrn)
            use const_def
            use eos_def
            use eos_lib
            use chem_lib, only: get_Q
            integer, intent(in) :: step, n
            real(dp), intent(in) :: told, time, logT_in, logRho_in
            real(dp), intent(in) :: x(:)
            integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.

            real(dp) :: logT, logRho, lgPgas, Pgas, Prad, lgP, avg_eps_nuc
            real(dp) :: eps_neu_total, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, &
              eps_nuc_categories(num_categories)
            real(dp), dimension(:), pointer :: &
               rate_raw, rate_raw_dT, rate_raw_dRho, &
               rate_screened, rate_screened_dT, rate_screened_dRho

            real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
            real(dp) :: dlnRho_dlnT_const_P, d_epsnuc_dlnT_const_P, d_Cv_dlnT
            real(dp) :: res(num_eos_basic_results)
            real(dp) :: d_dlnRho_const_T(num_eos_basic_results) 
            real(dp) :: d_dlnT_const_Rho(num_eos_basic_results) 
            real(dp) :: d_dabar_const_TRho(num_eos_basic_results) 
            real(dp) :: d_dzbar_const_TRho(num_eos_basic_results) 

            real(dp) :: Rho, T, xsum, d_eps_nuc_dx(species), dx, enuc, &
                  dt, energy, entropy, burn_ergs, &
                  xh, xhe, Z, mass_correction
         
            integer :: i, j, lwork, adjustment_iso, cid, ierr, max_j
            real(dp), dimension(species) :: dabar_dx, dzbar_dx, eps_nuc_dx, dmc_dx
            real(dp), pointer :: work(:), actual_Qs(:), actual_neuQs(:)
            logical, pointer :: from_weaklib(:)
            type (Net_General_Info), pointer  :: g
            logical :: skip_jacobian

            include 'formats'
         
            irtrn = 0
            if (time == 0) return

            logT = logT_in
            logRho = logRho_in
            
            if ((.not. quiet) .and. step > 1 .and. mod(step,50) == 0) then
               max_j = maxloc(x(1:species),dim=1)
               write(*,2) 'step, time, logT, logRho, ' // trim(chem_isos% name(chem_id(max_j))), &
                  step, time, logT, logRho, x(max_j)
               if (.false.) then
                  do j=1,species
                     write(*,2) trim(chem_isos% name(chem_id(j))), j, x(j)
                  end do
                  write(*,1) 'sum(x)', sum(x(1:species))
                  write(*,*)
               end if
            end if

            ierr = 0
            lwork = net_work_size(handle, ierr) 
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
            allocate(work(lwork), &
                  actual_Qs(num_reactions), actual_neuQs(num_reactions), &
                  from_weaklib(num_reactions), &
                  stat=ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
            xin(1:species) = x(1:species)

            call composition_info( &
                  species, chem_id, xin(1:species), xh, xhe, z, abar, zbar, z2bar, z53bar, ye, &
                  mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
            Z = max(0d0, min(1d0, 1d0 - (xh + xhe)))
         
            if (burn_at_constant_P) then

               logT = x(n)/ln10
               T = exp10(logT)
               Prad = Radiation_Pressure(T)
               Pgas = pressure - Prad
               lgPgas = log10(Pgas)
         
               call eosPT_get( &
                  eos_handle, Z, xh, abar, zbar, &
                  species, chem_id, net_iso, x, &
                  Pgas, lgPgas, T, logT, &
                  Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
                  res, d_dlnRho_const_T, d_dlnT_const_Rho,  &
                  d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
                  !ierr)
               if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

            else ! this is okay for burn_at_constant_density as well as constant T and Rho
         
              ! logT = rpar(r_burn_prev_lgT)
              ! logRho = rpar(r_burn_prev_lgRho)
               T = exp10(logT)
               Rho = exp10(logRho)
         
               call eosDT_get_legacy( &
                  eos_handle, Z, xh, abar, zbar, &
                  species, chem_id, net_iso, x, &
                  Rho, logRho, T, logT, &
                  res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
                  !Pgas, Prad, energy, entropy, ierr)
               if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            
            end if
         
            lgPgas = res(i_lnPgas)/ln10
            Pgas = exp10(lgPgas)
            Prad = Radiation_Pressure(T)
            lgP = log10(Pgas + Prad)
            skip_jacobian = .false.
            eta = res(i_eta)
            d_eta_dlnT = 0d0
            d_eta_dlnRho = 0d0
         
            call net_get_with_Qs( &
                  handle, skip_jacobian, netinfo, species, num_reactions, &
                  xin(1:species), T, logT, Rho, logRho, &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
                  screening_mode, &
                  eps_nuc_categories, eps_neu_total, &
                  lwork, work, actual_Qs, actual_neuQs, from_weaklib, ierr)
            if (ierr /= 0) then
               write(*,*)
               write(*,1) 'logT', logT
               write(*,1) 'logRho', logRho
               write(*, *) 'bad return from net_get'
               call mesa_error(__FILE__,__LINE__)
            end if

            call get_net_ptr(handle, g, ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
            dt = time - told

            ! set burn_ergs according to change from initial abundances
            eps_nuc = 0
            xsum = 0
            burn_ergs = 0
            do i=1,species
               cid = chem_id(i)
               
               dx = x(i) - x_initial(i)
               xsum = xsum + x(i)
               burn_ergs = burn_ergs + &
                  (get_Q(chem_isos,cid))*dx/chem_isos% Z_plus_N(cid)
                     
               dx = x(i) - x_previous(i)
               eps_nuc = eps_nuc + &
                  (get_Q(chem_isos,cid))*dx/chem_isos% Z_plus_N(cid)
                     
            end do
            avg_eps_nuc = burn_ergs*Qconv/time - eps_neu_total
            eps_nuc = eps_nuc*Qconv/dt - eps_neu_total
            burn_logT = logT
            burn_logRho = logRho
            burn_lnS = res(i_lnS)
            burn_lnE = res(i_lnE)
         
            x_previous(1:species) = x(1:species)
         
            if (time >= data_output_min_t) then
         
               call get_net_rate_ptrs(g% handle, &
                  rate_screened, rate_screened_dT, rate_screened_dRho, &
                  rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
                  ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in get_net_rate_ptrs'
                  call mesa_error(__FILE__,__LINE__)
               end if

               write(io_out,'(i7,99(1pe26.16,1x))',advance='no') &
                  step, &
                  sign(1d0,avg_eps_nuc)*log10(max(1d0,abs(avg_eps_nuc))), &
                  avg_eps_nuc, &
                  sign(1d0,eps_nuc)*log10(max(1d0,abs(eps_nuc))), &
                  eps_nuc, &
                  eps_neu_total, &
                  burn_logT - starting_logT, &
                  burn_logT, &
                  burn_logRho, &
                  lgPgas, &
                  lgP, &
                  abar, &
                  zbar, &
                  exp(burn_lnS)/(avo*kerg), &
                  burn_lnE/ln10, &
                  time, safe_log10(time), safe_log10(time/secyer), &
                  time - told, safe_log10(time - told), ye, xsum-1
               do i=1,num_names_of_isos_to_show
                  if (num_names_of_isos_to_show < species) then
                     cid = burn_isos_to_show(i)
                  else
                     if (i > species) exit
                     cid = chem_id(i)
                  endif
                  j = net_iso(cid)
                  if (j == 0) cycle

! output mass fractions, not abundances - fxt
!                  write(io_out,'(1pe26.16,1x)',advance='no') safe_log10(x(j))
                  write(io_out,'(1pe26.16,1x)',advance='no') safe_log10(x(j)*chem_isos% Z_plus_N(cid))
               end do
               do i=1,num_names_of_isos_to_show
                  if (num_names_of_isos_to_show < species) then
                     cid = burn_isos_to_show(i)
                  else
                     if (i > species) exit
                     cid = chem_id(i)
                  endif
                  j = net_iso(cid)
                  if (j == 0) cycle

! output mass fractions, not abundances - fxt
!                  write(io_out,'(1pe26.16,1x)',advance='no') x(j)
                  write(io_out,'(1pe26.16,1x)',advance='no') x(j)*chem_isos% Z_plus_N(cid)
               end do
               do i=1,num_reactions_to_track
                  j = index_for_reaction_to_track(i)
                  if (j == 0) cycle
                  write(io_out,'(1pe26.16,1x)',advance='no') rate_raw(j)
                  write(io_out,'(1pe26.16,1x)',advance='no') rate_screened(j)
               end do
               write(io_out,*) 
               do j=1, species
                  if (x(j) > peak_abundance(j)) then
                     peak_abundance(j) = x(j)
                     peak_time(j) = time
                  end if
               end do
            end if

            if (show_Qs) then
               write(*,*)
               write(*,1) 'logT', logT
               write(*,1) 'logRho', logRho
               write(*,*)
               write(*,'(30x,4a20)') 'Q total', 'Q neutrino', 'Q total-neutrino'
               do i = 1, num_reactions
                  if (from_weaklib(i)) then
                     write(*,'(a30,99f20.10)') 'weaklib ' // trim(reaction_Name(reaction_id(i))), &
                        actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
                  else
                     write(*,'(a30,99f20.10)') trim(reaction_Name(reaction_id(i))), &
                        actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
                  end if
               end do
               write(*,*)
               stop 'show_Qs'
            end if
         
         
            deallocate(work, actual_Qs, actual_neuQs, from_weaklib)

         end subroutine burn_solout1


      end subroutine Do_One_Zone_Burn
      
      
      subroutine test_net_setup(net_file_in)
         character (len=*), intent(in) :: net_file_in
         integer :: ierr, i
         
         include 'formats'
         
         net_file = net_file_in

         call net_init(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_handle failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_start_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_start_def failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call read_net_file(net_file, handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_net_file failed ', trim(net_file)
            call mesa_error(__FILE__,__LINE__)
         end if

         call net_finish_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_finish_def failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_ptr failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         species = g% num_isos
         num_reactions = g% num_reactions

         allocate(which_rates(rates_reaction_id_max), reaction_id(num_reactions))
         which_rates(:) = which_rates_choice

         call set_which_rates(ierr)
         if (ierr /= 0) then
            write(*,*) 'set_which_rates failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_set_which_rates(handle, which_rates, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_set_which_rates failed'
            call mesa_error(__FILE__,__LINE__)
         end if

         call net_setup_tables(handle, cache_suffix, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_setup_tables failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get_chem_id_table(handle, species, chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_chem_id_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get_net_iso_table(handle, net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_net_iso_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if

         call get_reaction_id_table(handle, num_reactions, reaction_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_reaction_id_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         allocate( &
               xin(species), xin_copy(species), d_eps_nuc_dx(species), &
               dxdt(species), d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species, species))
     
      end subroutine test_net_setup
      
      
      subroutine set_which_rates(ierr)
         use rates_def
         use rates_lib
         integer, intent(out) :: ierr
         integer :: which_rate
         
         ierr = 0
         
         if (len_trim(set_rate_c12ag) > 0) then
            if (set_rate_c12ag == 'NACRE') then
               which_rate = use_rate_c12ag_NACRE
            else if (set_rate_c12ag == 'Buchmann') then
               which_rate = use_rate_c12ag_JR
            else if (set_rate_c12ag == 'Kunz') then
               which_rate = use_rate_c12ag_Kunz
            else
               write(*,*) 'invalid string for set_rate_c12ag ' // trim(set_rate_c12ag)
               ierr = -1
               return
            end if
            call set_which_rate_c12ag(which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_3a) > 0) then
            if (set_rate_3a == 'NACRE') then
               which_rate = use_rate_3a_NACRE
            else if (set_rate_3a == 'Fynbo') then
               which_rate = use_rate_3a_JR
            else if (set_rate_3a == 'CF88') then
               which_rate = use_rate_3a_CF88
            else if (set_rate_3a == 'FL87') then
               which_rate = use_rate_3a_FL87
            else
               write(*,*) 'invalid string for set_rate_3a ' // trim(set_rate_3a)
               ierr = -1
               return
            end if
            call set_which_rate_3a(which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_1212) > 0) then
            if (set_rate_1212 == 'CF88_1212') then
               stop 'fix use_rate_1212_CF88'
               !which_rate = use_rate_1212_CF88
            else
               write(*,*) 'invalid string for set_rate_1212 ' // trim(set_rate_1212)
               ierr = -1
               return
            end if
            call set_which_rate_1212(which_rates, which_rate)
         end if
         
      end subroutine set_which_rates


      end module mod_one_zone_support



      module mod_one_zone_burn
      use chem_lib
      use net_def
      use net_lib
      use rates_lib, only: rates_init
      use rates_def
      use const_lib
      use utils_lib
      use mtx_def

      use mod_one_zone_support
      
      implicit none
      
      integer :: ierr, unit
         
      namelist /one_zone/ &
         mesa_dir, net_name, quiet, show_ye_stuff, num_names_of_isos_to_show, names_of_isos_to_show, &
         final_abundances_filename, save_final_abundances, show_peak_x_and_time, &
         initial_abundances_filename, read_initial_abundances, &
         read_T_Rho_history, T_Rho_history_filename, &
         num_isos_for_Xinit, names_of_isos_for_Xinit, values_for_Xinit, uniform_Xinit, &
         burn_tend, burn_rho, burn_temp, burn_logRho, burn_logT, burn_eta, burn_deta_dlnT, &
         burn_Cv, burn_d_Cv_dlnT, &
         eps, odescal, stptry, trace, burn_dbg, use_pivoting, &
         burn_rtol, burn_atol, burn_xmin, burn_xmax, weak_rate_factor, &
         min_for_show_peak_abundances, max_num_for_show_peak_abundances, &
         data_output_min_t, data_filename, &
         which_solver, screening_mode, which_rates_choice, &
         data_heading_line, show_net_reactions_info, &
         rattab_logT_lower_bound, rattab_logT_upper_bound, max_steps, max_step_size, &
         decsol_switch, small_mtx_decsol, large_mtx_decsol, &
         burn_at_constant_P, burn_at_constant_density, &
         starting_temp, pressure, cache_suffix, &
         num_times_for_burn, times_for_burn, log10Ts_for_burn, &
         log10Rhos_for_burn, etas_for_burn, log10Ps_for_burn, &
         set_rate_c12ag, set_rate_n14pg, set_rate_3a, set_rate_1212, &
         show_Qs, num_reactions_to_track, reaction_to_track,  &
         num_special_rate_factors, reaction_for_special_factor, special_rate_factor

      contains
      
      
      subroutine do_one_burn(filename, qt)
      character(len=*) :: filename
      logical, intent(in) :: qt
      
      include 'formats'
      
      ! set defaults
      
      mesa_dir = '../..'         
      net_name = 'test.net'
      quiet = .false.
      show_ye_stuff = .false.
      cache_suffix = '0'
      final_abundances_filename = ''
      save_final_abundances = .false.
      initial_abundances_filename = ''
      read_initial_abundances = .false.
      read_T_Rho_history = .false.
      T_Rho_history_filename = ''
      burn_tend = 10 ! seconds
      burn_rho = -1d99
      burn_temp = -1d99
      burn_logT = -1d99
      burn_logRho = -1d99
      burn_eta = -1d99
      burn_deta_dlnT = -1d99
      burn_Cv = -1d99
      burn_d_Cv_dlnT = -1d99
      burn_rtol = 1d-8
      burn_atol = 1d-9
      burn_xmin = -1d-10
      burn_xmax = 1d0+1d-10
      eps = 1d-8
      odescal = 1d-12
      stptry = 0d0
      trace = .false.
      burn_dbg = .false.
      use_pivoting = .true.
      show_net_reactions_info = .false.
      show_Qs = .false.
      decsol_switch = 50
         ! if current number of species <= switch,
            ! then use small_mtx_decsol,
            ! else use large_mtx_decsol.
      small_mtx_decsol = 'lapack'
      large_mtx_decsol = ''
      which_solver = 'rodas4_solver'
      rattab_logT_lower_bound = -1
      rattab_logT_upper_bound = -1
      max_steps = 50000
      uniform_Xinit = .false.
      burn_at_constant_P = .false.
      burn_at_constant_density = .false.
      starting_temp = -1
      pressure = -1
      num_times_for_burn = 0 ! <= 0 means don't use the arrays
      times_for_burn = 0
      log10Ts_for_burn = 0
      log10Rhos_for_burn = 0
      etas_for_burn = 0
      log10Ps_for_burn = 0
      max_step_size = 0
      
      min_for_show_peak_abundances = 1d-3 ! show if peak is > this
      max_num_for_show_peak_abundances = 21
      show_peak_x_and_time = .true.
      
      data_filename = 'one_zone_burn.data'
      data_output_min_t = -99
      
      num_names_of_isos_to_show = -1

      num_isos_for_Xinit = 4
      names_of_isos_for_Xinit(1:num_isos_for_Xinit) = (/ &
         'he4', 'c12', 'n14', 'o16' /)
      values_for_Xinit(1:num_isos_for_Xinit) = (/ &
         0.95d0, 0.005d0, 0.035d0, 0.010d0 /)
      
      screening_mode = extended_screening
      which_rates_choice = rates_NACRE_if_available

      num_special_rate_factors = 0 ! must be <= max_num_special_rate_factors
      reaction_for_special_factor(:) = ''
      special_rate_factor(:) = 1
      
      num_reactions_to_track = 0
      reaction_to_track(:) = ''

      set_rate_c12ag = '' ! empty string means ignore this control
         ! one of 'NACRE', 'Buchmann', or 'Kunz'
      set_rate_n14pg = '' ! empty string means ignore this control
         ! one of 'NACRE', 'Imbriani', or 'CF88'
      set_rate_3a = '' ! empty string means ignore this control
         ! one of 'NACRE', 'Fynbo', 'CF88', or 'FL87'
            ! FL87 is Fushiki and Lamb, Apj, 317, 368-388, 1987
            ! and includes both strong screening and pyconuclear
      set_rate_1212 = '' ! empty string means ignore this control
         ! one of 'CF88_basic_1212', 'CF88_multi_1212'
         ! CF88_basic_1212 is the single rate approximation from CF88.
         ! CF88_multi_1212 combines the rates for the n, p, and a channels.
            ! c12(c12,n)mg23, c12(c12,p)na23, and c12(c12,a)ne20
            ! uses neutron branching from dayras, switkowski, and woosley, 1976.
      
      weak_rate_factor = 1
      
      ! read inlist
      
      open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
      if (ierr /= 0) then
         write(*, *) 'Failed to open control namelist file ', trim(filename)
         call mesa_error(__FILE__,__LINE__)
      else
         read(unit, nml=one_zone, iostat=ierr)  
         close(unit)
         if (ierr /= 0) then
            write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
            write(*, '(a)') &
               'The following runtime error message might help you find the problem'
            write(*, *) 
            open(unit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            read(unit, nml=one_zone)
            call mesa_error(__FILE__,__LINE__)
         end if  
      end if

      ! do initialization
      
      if(burn_temp<0.d0 .and. burn_logT<0.d0) then
         call mesa_error(__FILE__,__LINE__,"Must set either burn_temp or burn_logT")
         stop
      end if
      if (burn_temp < 0) burn_temp = exp10(burn_logT)
      if (burn_logT < 0) burn_logT = log10(burn_temp)
      if(burn_rho<0.d0 .and. burn_logRho<0.d0) then
         call mesa_error(__FILE__,__LINE__,"Must set either burn_rho or burn_logRho")
         stop
      end if
      if (burn_rho < 0) burn_rho = exp10(burn_logRho)
      if (burn_logRho < 0) burn_logRho = log10(burn_rho)

      starting_temp = burn_temp

      allocate(net_iso(num_chem_isos), chem_id(num_chem_isos))
      
      !reaclib_filename = 'jina_reaclib_results_20130213default2'
      !write(*,*) 'changing reaclib_filename'
      
      open(unit=io_out, file=trim(data_filename), action='write', iostat=ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      complete_silence_please = qt
      call Do_One_Zone_Burn(net_name)
         
      open(unit=io_out, file=trim(data_filename), action='write', iostat=ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      end subroutine do_one_burn
      
      
      end module mod_one_zone_burn
