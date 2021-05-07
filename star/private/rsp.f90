! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton, Radek Smolec & The MESA Team
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

      module rsp
      use star_def, only: star_ptr, star_info
      use rsp_eval_eos_and_kap, only: init_for_rsp_eos_and_kap
      use rsp_def
      use rsp_step, only: calculate_energies, init_HYD, HYD

      implicit none
      
      private
      public :: rsp_setup_part1, rsp_setup_part2, rsp_one_step, &
         build_rsp_model, rsp_total_energy_integrals, do1_rsp_build
      
      contains
      
      
      subroutine do1_rsp_build(s,ierr) 
         ! call from other_rsp_build_model after changing params.
         ! can change rsp_* params; but cannot change nz or net.
         ! multiple calls are ok to search.
         use rsp_build, only: do_rsp_build
         use hydro_vars, only: set_vars, set_Teff
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         include 'formats'
         call init_def(s)
         call init_for_rsp_eos_and_kap(s)  
         s% rsp_period = 0d0
         call do_rsp_build(s,ierr)
         if (ierr /= 0) return
         do k=1,NZN
            s% v(k) = 0d0
            s% xh(s% i_v,k) = s% v(k)
         end do
         ierr = 0
         call finish_build_rsp_model(s,ierr)
         if (ierr /= 0) return
         s% doing_finish_load_model = .true.
         call set_vars(s, 0d0, ierr)
         if (ierr /= 0) return
         s% doing_finish_load_model = .false.
         call set_Teff(s, ierr)
         if (ierr /= 0) return
      end subroutine do1_rsp_build


      subroutine build_rsp_model(s,ierr)
         ! called by model_builder in place of loading a model
         use alloc, only: allocate_star_info_arrays, set_RSP_flag
         use rsp_build, only: do_rsp_build
         use net, only: do_micro_change_net
         use const_def, only: Lsun, Rsun, Msun
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: i, j, k
         include 'formats'
         NSTART = 1
         s% nz = s% RSP_nz         
         if (s% job% change_initial_net) then
            call do_micro_change_net(s, s% job% new_net_name, ierr)
         else
            call do_micro_change_net(s, 'o18_and_ne22.net', ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed in do_micro_change_net'
            return
         end if         
         s% tau_factor = s% RSP_surface_tau/s% tau_base
         call init_def(s)
         call init_allocate(s,s% nz)         
         call allocate_star_info_arrays(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in allocate_star_info_arrays'
            return
         end if
         call set_RSP_flag(s% id, .true., ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_RSP_flag'
            return
         end if
         call init_for_rsp_eos_and_kap(s)  
         s% rsp_period = 0d0
         if (s% RSP_use_atm_grey_with_kap_for_Psurf) then
            s% tau_factor = s% RSP_tau_surf_for_atm_grey_with_kap/s% tau_base
            write(*,1) 'set s% tau_factor using RSP_tau_surf', s% tau_factor
         end if
         if (s% use_other_rsp_build_model) then
            call s% other_rsp_build_model(s% id,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in other_rsp_build_model'
               return
            end if
         else
            call do_rsp_build(s,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_rsp_build'
               return
            end if            
         end if
         if (.not. s% use_RSP_new_start_scheme) &
            call set_random_velocities(s)
         rsp_tau_factor = s% tau_factor
         PERIODLIN = s% rsp_period
         if (s% rsp_period > 0d0) then
            s% rsp_dt = s% rsp_period/dble(s% rsp_target_steps_per_cycle)
         else
            s% rsp_dt = 1d0
         end if
         s% dt_next = s% rsp_dt
         s% dt = s% rsp_dt
         if (is_bad(s% dt)) then
            write(*,1) 'dt', s% dt
            return
         end if
         call finish_build_rsp_model(s, ierr)      
         write(*,2) 'nz', s%nz
         write(*,1) 'T(nz)', s% T(s%nz)             
         write(*,1) 'L_center/Lsun', s% L_center/Lsun           
         write(*,1) 'R_center/Rsun', s% R_center/Rsun           
         write(*,1) 'M_center/Msun', s% M_center/Msun           
         write(*,1) 'L(1)/Lsun', s% L(1)/Lsun           
         write(*,1) 'R(1)/Rsun', s% r(1)/Rsun           
         write(*,1) 'M(1)/Msun', s% m(1)/Msun           
         write(*,1) 'v(1)/1d5', s% v(1)/1d5       
         write(*,1) 'tau_factor', s% tau_factor   
         write(*,1) 'tau_base', s% tau_base   
         write(*,*) 
      end subroutine build_rsp_model
      

      subroutine finish_build_rsp_model(s,ierr)
         use star_utils, only: &
            normalize_dqs, set_qs, set_m_and_dm, set_dm_bar, &
            store_rho_in_xh, store_T_in_xh, store_r_in_xh
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: i, k, j
         include 'formats'
         do i=1,NZN
            k = NZN+1 - i
            s% Prad(k) = f_Edd_isotropic*s% erad(k)/s% Vol(k)
            if (is_bad(s% Prad(k))) then
               write(*,2) 's% Prad(k)', k, s% Prad(k)
               write(*,2) 's% erad(k)', k, s% erad(k)
               write(*,2) 's% Vol(k)', k, s% Vol(k)
               stop 'build_rsp_model'
            end if
            s% dq(k) = s% dm(k)/s% xmstar
            if (is_bad(s% Vol(k)) .or. s% Vol(k) <= 0d0) then
               write(*,2) 's% Vol(k)', I, s% Vol(k)
               stop 'build_rsp_model'
            end if
            call store_rho_in_xh(s, k, 1d0/s% Vol(k))
            call store_T_in_xh(s, k, s% T(k))
            call store_r_in_xh(s, k, s% r(k))
            s% xh(s% i_Et_RSP,k) = s% RSP_w(k)*s% RSP_w(k)
            do j=1,s% species
               s% xa(j,k) = xa(j)
            end do
            s% xh(s% i_v,k) = s% v(k)
         end do
         s% dq(s% nz) = (s% m(NZN) - s% M_center)/s% xmstar
         if (.not. s% do_normalize_dqs_as_part_of_set_qs) then
            call normalize_dqs(s, NZN, s% dq, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'normalize_dqs failed in finish_build_rsp_model'
               return
            end if
         end if
         call set_qs(s, NZN, s% q, s% dq, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_qs'
            stop 'build_rsp_model'
         end if
         call set_m_and_dm(s)
         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)
      end subroutine finish_build_rsp_model

      
      subroutine set_random_velocities(s)
         use star_utils, only: rand
         type (star_info), pointer :: s
         integer :: i, k
         if (s% RSP_Avel /= 0d0 .or. s% RSP_Arnd /= 0d0) then
            do i=1,NZN
               k = NZN+1-i
               s% v(k) = s% v(k) + &
                  1d5*dble(i)/dble(NZN)*(s% RSP_Avel + s% RSP_Arnd*(2d0*rand(s) - 1d0))
            end do
            write(*,*) 'set random velocities'
            s% RSP_have_set_velocities = .true.
         end if      
      end subroutine set_random_velocities
      

      subroutine rsp_setup_part1(s,restart,ierr)
         ! called by finish_load_model before set_vars
         use const_def, only: crad
         use rsp_eval_eos_and_kap, only: &
            restart_rsp_eos_and_kap, get_surf_P_T_kap
         use alloc, only: resize_star_info_arrays
         use star_utils, only: get_XYZ
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), target :: copy_info
         type (star_info), pointer :: c, prv
         real(dp) :: tau_surf, kap_guess, T_surf, Psurf, kap_surf, Teff_atm, Y
         integer :: i, k
         logical :: okay
         include 'formats'
         ierr = 0
         if (s% RSP_use_atm_grey_with_kap_for_Psurf .and. &
              (s% atm_option /= 'T_tau' .or. &
               s% atm_T_tau_relation /= 'Eddington' .or. &
               s% atm_T_tau_opacity /= 'iterated')) then
            write(*,*)
            write(*,*) 'when RSP_use_atm_grey_with_kap_for_Psurf, must set:'
            write(*,*) '   atm_option = ''T_tau'''
            write(*,*) '   atm_T_tau_relation = ''Eddington'''
            write(*,*) '   atm_T_tau_opacity = ''iterated'''
            ierr = -1
            stop 'rsp_setup_part1'
            return
         end if
         if (restart) then
            NSTART = 2
            call init_def(s)
            call restart_rsp_eos_and_kap(s)
         else
            prev_cycle_run_num_steps = 0
            run_num_iters_prev_period = 0
            run_num_retries_prev_period = 0
            NSTART = 1
            if (s% job% load_saved_model_for_RSP .and. &
                  .not. s% job% create_RSP_model) then
               call init_def(s)
               call init_allocate(s,s% nz)
               call get_XYZ(s, s% xa(:,1), s% RSP_X, Y, s% RSP_Z)               
               call init_for_rsp_eos_and_kap(s)
               IWORK=0
               NZN = s% nz
               ELSTA = s% L(1)
               RSTA = s% r(1)  
               s% rsp_dt = s% dt_next
               if (s% max_timestep > 0d0 .and. s% rsp_dt > s% max_timestep) &
                  s% rsp_dt = s% max_timestep
               rsp_tau_factor = s% tau_factor
               s% rsp_period = s% rsp_dt*dble(s% RSP_target_steps_per_cycle)               
               s% RSP_have_set_velocities = .true.
               call copy_from_xh_to_rsp(s,-1)
               do k=1,NZN
                  s% L_start(k) = 0d0
               end do
               if (s% RSP_use_atm_grey_with_kap_for_Psurf) then
                  tau_surf = s% RSP_tau_surf_for_atm_grey_with_kap
                  kap_guess = 1d-2
                  call get_surf_P_T_kap(s, &
                     s% m(1), s% r(1), s% L(1), tau_surf, kap_guess, &
                     T_surf, Psurf, kap_surf, Teff_atm, ierr)
                  if (ierr /= 0) stop 'failed in get_surf_P_T_kap'
               else if (s% RSP_use_Prad_for_Psurf) then
                  Psurf = crad*s% T(1)**4/3d0
               else
                  Psurf = 0d0
               end if
               Psurf_from_atm = Psurf
            else
               s% dt_next = s% rsp_dt
               s% dt = s% rsp_dt
            end if
            rsp_min_dr_div_cs = 1d99
            rsp_min_rad_diff_time = 1d99
            call begin_calculation(s,restart,ierr)  
         end if  
         s% tau_factor = rsp_tau_factor
      end subroutine rsp_setup_part1
         

      subroutine rsp_setup_part2(s, restart, want_rsp_model, is_rsp_model, ierr)
         use hydro_vars, only: set_Teff
         use rsp_step
         ! called by finish_load_model after set_vars
         type (star_info), pointer :: s
         logical, intent(in) :: restart, want_rsp_model, is_rsp_model
         integer, intent(out) :: ierr
         integer :: i, j, k, species, nz, op_err
         real(dp), allocatable :: w_avg(:)
         real(dp) :: &
            dFr_dr_out, dFr_dr_in, dFr_dr_00, &
            dFr_dT_out, dFr_dT_00, dFr_dVol_00, &
            Lr, Lc, POM
         logical :: okay
         include 'formats'
         if (restart) then
            call finish_rsp_photo_in(s)
            return
         end if
         ierr = 0
         nz = s% nz
         if (want_rsp_model .and. .not. is_rsp_model) then
            ! have loaded a file that is not an rsp model and now must convert model
            do k=1,nz
               s% erad(k) = crad*s% T(k)**4/s% rho(k)
               s% xh(s% i_erad_RSP,k) = s% erad(k)
            end do
            !$OMP PARALLEL DO PRIVATE(I,op_err) SCHEDULE(dynamic,2)
            do i = 1,nz
               call do1_specific_volume(s,i)
               call do1_eos_and_kap(s,i,op_err)
               if (op_err /= 0) ierr = op_err
               call calc_Prad(s,i)
            end do
            !$OMP END PARALLEL DO
            if (ierr /= 0) return
            !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
            do k=1,nz ! set Fr using updated eos and kap vars
               call T_form_of_calc_Fr(s, nz+1-k, s% Fr(k), &
                  dFr_dr_out, dFr_dr_in, dFr_dr_00, &
                  dFr_dT_out, dFr_dT_00, dFr_dVol_00)
               s% xh(s% i_Fr_RSP,k) = s% Fr(k)
            end do
            !$OMP END PARALLEL DO
            if (ALFAC == 0d0 .or. ALFAS == 0d0) then
               s% RSP_w(1:nz) = 0d0
               s% RSP_Et(1:nz) = 0d0
               s% xh(s% i_Et_RSP,1:nz) = 0d0
            else
               !$OMP PARALLEL DO PRIVATE(I) SCHEDULE(dynamic,2)
               do i = 1,nz
                  call calc_Hp_face(s,i)
                  call calc_Y_face(s,i)
                  call calc_PII_face(s,i)
               end do
               !$OMP END PARALLEL DO
               allocate(w_avg(nz))
               do k=1,nz
                  Lr = 4d0*pi*s% r(k)**2*s% Fr(k)
                  Lc = s% L(k) - Lr
                  POM = P4*(s% r(k)**2)*(ALFAC/ALFAS)*s% T(k)/s% Vol(k)
                  w_avg(k) = Lc/(POM*s% PII(k))
               end do
               s% RSP_w(1) = 0d0
               s% RSP_w(nz) = 0d0
               do k=2,nz-1
                  s% RSP_w(k) = 0.5d0*(w_avg(k) + w_avg(k+1))
                  if (s% RSP_w(k) < 0d0) s% RSP_w(k) = 0d0
                  s% RSP_Et(k) = s% RSP_w(k)**2
                  s% xh(s% i_Et_RSP,k) = s% RSP_Et(k)               
                  s% RSP_w(k) = sqrt(s% xh(s% i_Et_RSP,k))
               end do   
            end if         
         end if
         call finish_after_build_model(s)
         call copy_results(s)  
         call set_Teff(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'rsp_setup_part2: set_Teff returned ierr', ierr
            return
         end if
         TEFF = s% Teff
         if (s% rsp_period > 0d0) then
            s% rsp_dt = s% rsp_period/dble(s% rsp_target_steps_per_cycle)
            s% dt_next = s% rsp_dt*s% RSP_initial_dt_factor
            s% dt = s% dt_next
         end if
         if (s% use_other_rsp_build_model .and. &
               s% set_RSP_Psurf_to_multiple_of_initial_P1 > 0d0) then
            s% RSP_Psurf = s% Peos(1)*s% set_RSP_Psurf_to_multiple_of_initial_P1
            write(*,1) 'rsp_setup_part2 set RSP_Psurf', s% RSP_Psurf
         end if
      end subroutine rsp_setup_part2


      subroutine get_LINA_info(s,ierr)
         use rsp_lina, only: do_LINA
         type (star_info), pointer :: s      
         integer, intent(out) :: ierr  
         
         real(dp), allocatable :: VEL(:,:)
         real(dp), allocatable, dimension(:) :: &
            M, DM, DM_BAR, R, Vol, T, RSP_Et, Lr
         integer :: NMODES, I, k, sz
         real(dp) :: amix1, amix2, velkm
         
         include 'formats'
         ierr = 0
         
         if (s% RSP_kick_vsurf_km_per_sec == 0d0) then
            write(*,*) 'skip calling LINA since RSP_kick_vsurf_km_per_sec = 0'
            return
         end if
         
         sz = NZN+1
         
         allocate(VEL(sz,15), &
            M(sz), DM(sz), DM_BAR(sz), R(sz), Vol(sz), T(sz), RSP_Et(sz), Lr(sz))
            
         do i=1,NZN
            k = NZN+1-i 
            M(i) = s% m(k)
            DM(i) = s% dm(k)
            DM_BAR(i) = s% dm_bar(k)
            R(i) = s% r(k)
            Vol(i) = s% Vol(k)
            T(i) = s% T(k)
            RSP_Et(i) = s% RSP_Et(k)
            Lr(i) = 4d0*pi*s% r(k)**2*s% Fr(k)
         end do                    
        
         NMODES = s% RSP_nmodes

         call do_LINA(s, s% RSP_L*SUNL, NZN, NMODES, VEL, &
            s% rsp_LINA_periods, s% rsp_LINA_growth_rates, &
            M, DM, DM_BAR, R, Vol, T, RSP_Et, Lr, ierr)
         if (ierr /= 0) return

         write(*,'(a)') '            P(days)         growth'
         do I=1,NMODES
            write(*,'(I3,2X,99e16.5)') I-1, &
               s% rsp_LINA_periods(I)/86400.d0, &
               s% rsp_LINA_growth_rates(I)
         enddo
         
         s% rsp_period = &
            s% rsp_LINA_periods(s% RSP_mode_for_setting_PERIODLIN + 1)
         
         amix1 = s% RSP_fraction_1st_overtone
         amix2 = s% RSP_fraction_2nd_overtone
         velkm = s% RSP_kick_vsurf_km_per_sec
         s% v_center = 0d0
         do i=1,NZN
            k = NZN+1-i    
            s% v(k)=1.0d5*VELKM* &
               ((1.0d0-AMIX1-AMIX2)*vel(i,1)+AMIX1*vel(i,2)+AMIX2*vel(i,3))
         end do
            
         s% RSP_have_set_velocities = .true.
         
      end subroutine get_LINA_info  


      subroutine rsp_one_step(s,ierr)
         use brunt, only: do_brunt_N2
         use rsp_step, only: rsp_set_Teff, &
            turn_off_time_weighting, turn_on_time_weighting
         type (star_info), pointer :: s      
         integer, intent(out) :: ierr    
         integer :: k, j, k_max_abs_rel_hse_err
         real(dp) :: rand, hse_err, max_abs_rel_hse_err
         logical :: restart
         
         include 'formats'

         ierr = 0
         s% RSP_just_set_velocities = .false.
         if (.not. s% RSP_have_set_velocities) then
         
            max_abs_rel_hse_err = 0d0
            k_max_abs_rel_hse_err = 0
            do k=2,s% nz
               hse_err = &
                  (s% Peos(k-1) - s% Peos(k))/(-s% cgrav(k)*s% m(k)*s% dm_bar(k)/(4d0*pi*s% r(k)**4)) - 1d0
               if (abs(hse_err) >= max_abs_rel_hse_err) then
                  max_abs_rel_hse_err = abs(hse_err)
                  k_max_abs_rel_hse_err = k
               end if
            end do

            s% need_to_save_profiles_now = .true.
            s% RSP_just_set_velocities = .true.        
         
            write(*,3) 'relaxation max_abs_rel_hse_err, days', s% model_number, k_max_abs_rel_hse_err, &
               max_abs_rel_hse_err, s% time/(24*3600)
            
            if (.not. s% use_other_RSP_linear_analysis) then    
               call get_LINA_info(s,ierr)               
            else               
               ! must set gradT before calling since gyre needs it.
               ! Y_face is superadiabatic gradient
               do k=1,NZN
                  s% gradT_sub_grada(k) = s% Y_face(k)
                  if (k > 1) s% gradT(k) = &
                     s% Y_face(k) + 0.5d0*(s% grada(k-1) + s% grada(k))
                  s% abar(k) = abar
                  s% zbar(k) = zbar
                  s% X(k) = X
                  s% Z(k) = Z
                  do j=1,s% species
                     s% xa(j,k) = xa(j)
                  end do
               end do
               s% gradT(1) = s% gradT(2)
               s% calculate_Brunt_N2 = .true.
               call do_brunt_N2(s, 1, NZN, ierr)
               if (ierr /= 0) return
               restart = .false.
               call s% other_rsp_linear_analysis(s% id, restart, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) &
                     write(*,*) 'other_rsp_linear_analysis ierr', ierr
                  return
               end if  
               s% RSP_have_set_velocities = .true.                
            end if                              
            
            PERIODLIN = s% rsp_period
            s% rsp_dt = s% rsp_period/dble(s% rsp_target_steps_per_cycle)
            s% dt = s% rsp_dt  
            
            s% cumulative_energy_error_old = 0d0
            s% time = 0d0
            s% time_old = 0d0
            write(*,*) 'automatically resets age and cumulative energy error info when sets velocities'
            s% need_to_save_profiles_now = .true.
            
            call set_random_velocities(s)
            
         end if 
         
         s% do_history_file = s% RSP_have_set_velocities ! don't write history entries until set velocities
         !call turn_on_time_weighting(s)
         
         if (s% dt > s% RSP_max_dt .and. s% RSP_max_dt > 0d0) then
            s% dt = s% RSP_max_dt
         end if
         
         call do1_step(s,ierr) 
         if (ierr /= 0) return
         
         call copy_results(s)  
         call rsp_set_Teff(s)
         if (s% RSP_write_map) call add_to_map

         
         contains
                           
         subroutine add_to_map
            use profile_getval, only: get_profile_val
            integer :: i, k, NPCH1, NPCH2, IP, n, io
            real(dp) :: ph_x, FASE
            character (len=256) :: fname
            include 'formats'
            NPCH1 = s% RSP_map_first_period
            NPCH2 = s% RSP_map_last_period
            IP = s% RSP_num_periods
            io = 77
            if (IP+1.ge.NPCH1.and.IP+1.le.NPCH2) then
               if(.not. writing_map) then
                  call read_map_specs(s,ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in read_map_specs'
                     return
                  end if
                  if (len_trim(s% RSP_map_filename) == 0) &
                     s% RSP_map_filename = 'map.data'
                  fname = trim(s% log_directory) // '/' // trim(s% RSP_map_filename)
                  open(io,file=trim(fname),status='unknown')
                  write(*,*) 'writing ' // trim(fname)
                  write(io,'(i18,1x,i4)',advance='no') 1, 2
                  do n=1,num_map_cols
                     write(io,'(1x,i18)',advance='no') n+2
                  end do
                  write(io,*)
                  write(io,'(a18,1x,a4)',advance='no') 'phase', 'zone'
                  do n=1,num_map_cols
                     write(io,'(1x,a18)',advance='no') trim(map_col_names(n))
                  end do
                  write(io,*)
                  writing_map = .true.
                  done_writing_map = .false.
                  s% need_to_set_history_names_etc = .true.
                  s% star_history_name = s% RSP_map_history_filename
                  FASE0 = s% time
               endif     
               FASE=(s% time-FASE0)/s% rsp_period
               !write(*,4) 'add to map', s% model_number, IP, NPCH2, FASE
               do k=1,NZN,s% RSP_map_zone_interval ! gnuplot pm3d map
                  I = NZN+1 - k
                  if(I.gt.IBOTOM.and.I.lt.NZN) then
                     write(io,'(d18.10,1x,i4)',advance='no') FASE, k
                     do n=1,num_map_cols
                        write(io,'(1x,d18.10)',advance='no') &
                           get_profile_val(s,map_ids(n),k)
                     end do
                     write(io,*)
                     !write(io,778) FASE,I,s% T(k), &
                     !    s% Hp_face(k),s% Y_face(k),s% PII(k),s% Lc(k)/s% L(k), &
                     !    4d0*pi*s% r(k)**2*s% Fr(k)/s% L(k),s% Lt(k)/s% L(k), &
                     !    s% RSP_w(k)**2,s% egas(k)+s% erad(k),s% csound(k), &
                     !    s% Pt(k),s% Pgas(k)+s% Prad(k),s% Eq(k), &
                     !    s% COUPL(k)
                  endif
               enddo
               !write(io,*)
            end if
            if(IP.eq.NPCH2 .and. .not. done_writing_map) then
               close(io)
               fname = trim(s% log_directory) // '/' // trim(s% RSP_map_filename)
               write(*,*) '  close ' // trim(fname)
               done_writing_map = .true.
            end if
 778  format(d16.10,1x,i3,14(1x,d16.10))                     
         end subroutine add_to_map
         
      end subroutine rsp_one_step
      
      
      subroutine read_map_specs(s,ierr)
         use utils_lib
         use utils_def
         use profile_getval, only: get_profile_id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: iounit, n, i, t, id, col
         character (len=strlen) :: buffer, string, filename

         include 'formats'

         ierr = 0
            
         filename = s% RSP_map_columns_filename
         if (len_trim(filename) == 0) filename = 'map_columns.list'
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            return
         end if

         n = 0
         i = 0
         col = 0
         
         num_map_cols = 0

         do
            t = token(iounit, n, i, buffer, string)
            if (t == eof_token) exit
            if (t /= name_token) then
               write(*,*) 'error reading map columns list ' // trim(filename)
               stop 'read_map_specs'
               call error; return
            end if
            id = get_profile_id(s, string)
            if (id < 0) then
               write(*,*) 'error: <' // trim(string) // '> in map columns is not in your profile list'
               write(*,*) 'please add it to your profile columns list and try again'
               write(*,*) 'also, replace any TAB characters by spaces in ' // trim(filename)
               stop 'read_map_specs'
               call error; return
            end if
            col = col+1
            if (col > max_map_cols) then
               write(*,*) 'error: ' // trim(filename) // ' has too many map columns'
               write(*,*) 'the max is currently fixed at ', max_map_cols
               write(*,*) 'please delete some and try again'
               call error; return
            end if
            map_col_names(col) = trim(string)
            map_ids(col) = id
         end do
         
         num_map_cols = col

         close(iounit)
         
         contains

         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error
      
      end subroutine read_map_specs
      

      subroutine do1_step(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: i, k
         real(dp) :: dr_div_cs, r_in, r_00, max_dt, target_dt, total_radiation

         include 'formats'
         
         ID=ID+1
      
         !call rsp_dump_for_debug(s)
         !stop 'rsp_dump_for_debug'
         target_dt = min( &
            s% rsp_period/dble(s% RSP_target_steps_per_cycle), &
            s% dt*s% max_timestep_factor)
         if (s% rsp_max_dt > 0) target_dt = s% rsp_max_dt ! force the timestep
         s% dt = target_dt

         if (is_bad(s% dt)) then
            write(*,1) 'dt', s% dt
            write(*,1) 'rsp_period', s% rsp_period
            write(*,2) 'RSP_target_steps_per_cycle', s% RSP_target_steps_per_cycle
            write(*,1) 'max_timestep_factor', s% max_timestep_factor
            stop 'do1_step 1'
         end if

         max_dt = rsp_min_dr_div_cs*s% RSP_max_dt_times_min_dr_div_cs
         if (s% RSP_max_dt_times_min_rad_diff_time > 0d0 .and. rsp_min_rad_diff_time > 0d0) then
            if (rsp_min_rad_diff_time*s% RSP_max_dt_times_min_rad_diff_time < max_dt) then
               max_dt = rsp_min_rad_diff_time*s% RSP_max_dt_times_min_rad_diff_time
               if (s% dt > max_dt) then
                  write(*,3) 'dt limited by rad diff time', NZN+1-i_min_rad_diff_time, s% model_number, &
                     s% dt, rsp_min_rad_diff_time, s% RSP_max_dt_times_min_rad_diff_time
                  !stop 'rsp'
               end if
            end if
         end if
         if (s% dt > max_dt) then
            if (s% RSP_report_limit_dt) &
               write(*,4) 'limit dt to max_dt', s% model_number
            s% dt = max_dt
         end if
         
         if (is_bad(s% dt) .or. s% dt <= 0d0) then
            write(*,1) 'dt', s% dt
            write(*,1) 'RSP_max_dt_times_min_dr_div_cs', s% RSP_max_dt_times_min_dr_div_cs
            write(*,1) 'rsp_min_dr_div_cs', rsp_min_dr_div_cs
            write(*,1) 'rsp_min_rad_diff_time', rsp_min_rad_diff_time
            stop 'do1_step 2'
         end if
         
         ierr = 0
         call HYD(s,ierr) 
         if (ierr /= 0) return
         ! s% dt might have been reduced by retries in HYD
         s% time = s% time_old + s% dt
         s% rsp_dt = s% dt ! will be used to set dt for next step
         
         ! set this here for use in next step. to avoid restart problems.
         rsp_min_dr_div_cs = 1d99
         i_min_dr_div_cs = -1
         r_00 = s% R_center
         do i = 1,nzn
            k = NZN+1-i
            r_in = r_00
            r_00 = s% r(k)
            if (abs(s% v(k))/s% csound(k) < s% RSP_v_div_cs_threshold_for_dt_limit) cycle
            k = nzn+1-i
            dr_div_cs = (r_00 - r_in)/s% csound(k)
            if (dr_div_cs < rsp_min_dr_div_cs) then
               rsp_min_dr_div_cs = dr_div_cs
               i_min_dr_div_cs = i
            end if
         end do
         
         rsp_min_rad_diff_time = 1d99
         i_min_rad_diff_time = -1
         if (s% RSP_max_dt_times_min_rad_diff_time > 0d0) then
            rsp_min_rad_diff_time = dt_for_radiative_diffusion(i_min_rad_diff_time)
         end if
         
         call calculate_work_integrals(s)      
         call calculate_energies(s,total_radiation)
         call gather_pulse_statistics(s)
         if (s% RSP_max_num_periods < 0 .or. &
             s% rsp_num_periods < s% RSP_max_num_periods) return
         call get_GRPDV(s)
                  
         contains

         real(dp) function dt_for_radiative_diffusion(i_min_rad_diff_time)
            integer, intent(out) :: i_min_rad_diff_time
            real(dp) :: min_dt, dt_cell, l, D, dr
            integer :: k, k_min_dt, nz
            include 'formats'
            min_dt = 1d99
            nz = s% nz
            k_min_dt = -1
            i_min_rad_diff_time = -1
            do k=1,nz
               l = s% V(k)/s% opacity(k) ! photon mean free path
               D = clight*l/3d0 ! diffusion coefficient, clight/(3*opacity*rho)
               if (k < nz) then
                  dr = s% r(k) - s% r(k+1)
               else
                  dr = s% r(k) - s% r_center
               end if
               ! if curious, ask Lars about the Pgas/Prad term
               dt_cell = dr**2/D*(1d0 + s% Pgas(k)/s% Prad(k))
               if (dt_cell < min_dt) then
                  min_dt = dt_cell
                  k_min_dt = k
               end if
            end do
            i_min_rad_diff_time = NZN-k_min_dt+1
            dt_for_radiative_diffusion = min_dt
         end function dt_for_radiative_diffusion         
         
      end subroutine do1_step
      
      
      subroutine gather_pulse_statistics(s) ! assumes have set EKMAX and EKMIN
         ! updates LMAX, LMIN, RMAX, RMIN, 
         !     s% rsp_GREKM, s% rsp_GREKM_avg_abs, s% rsp_DeltaR, s% rsp_DeltaMAG
         type (star_info), pointer :: s
         logical :: cycle_complete
         integer :: i, k
         include 'formats'
         if(s% L(1)/SUNL.gt.LMAX) LMAX=s% L(1)/SUNL
         if(s% L(1)/SUNL.lt.LMIN) LMIN=s% L(1)/SUNL      
         INSIDE=0      
         call check_cycle_completed(s,cycle_complete)
         ULL=UN
         TE_start=s% time
         if (cycle_complete) then
            if (s% rsp_num_periods > 1) then
               s% rsp_GREKM = (EKMAX-EKMAXL)/(EKMAX+EKMAXL)*2.d0
               s% rsp_GREKM_avg_abs = s% rsp_GREKM_avg_abs_frac_new*abs(s% rsp_GREKM) + &
                  (1d0 - s% rsp_GREKM_avg_abs_frac_new)*s% rsp_GREKM_avg_abs
            end if
            EKDEL = EKMAX-EKMIN
            EKMAXL = EKMAX
            EKMAX =-10.d50
            EKMIN = -EKMAX
            s% rsp_DeltaR = RMAX-RMIN
            RMAX = -SUNR
            RMIN = -RMAX
            s% rsp_DeltaMAG = 2.5d0*log10(LMAX/LMIN)
            LMIN = 10.d10
            LMAX = -LMIN
            VMAX = 0
         end if
      end subroutine gather_pulse_statistics


      subroutine get_GRPDV(s)
         type (star_info), pointer :: s
         integer :: I, k
         PDVWORK=0.d0
         do I=1,NZN
            k = NZN+1-i
            WORK(I)=  WORK(I)+(VV0(I)-s% Vol(k))*s% dm(k)* &
                 (THETA*(PPP0(I)+PPQ0(I))+ &
                 THETA1*((s% Pgas(k)+s% Prad(k))+s% Pvsc(k))) &
                 -s% dt*s% dm(k)*s% Eq(k)
            WORKQ(I)=  WORKQ(I)+(VV0(I)-s% Vol(k))*s% dm(k)* &
                 (THETA*PPQ0(I)+THETA1*s% Pvsc(k))
            PDVWORK=PDVWORK+WORK(i)
         enddo
         s% rsp_GRPDV=PDVWORK/EKDEL
         if (is_bad(s% rsp_GRPDV)) s% rsp_GRPDV=0d0
      end subroutine get_GRPDV
      
      
      subroutine begin_calculation(s,restart,ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         real(dp) :: total_radiation
         include 'formats'
         ierr = 0
         FIRST  = 0 
         TT1    = 0.d0
         EKMAX  = -10.d50
         EKMIN  = -EKMAX
         EKMAXL = EKMAX
         s% rsp_num_periods =-1
         ID = 0   !number of timesteps done in one period
         INSIDE = 1 ! for initial call
         !s% mstar = M(1)
         call set_star_vars(s,ierr)
         if(s% rsp_num_periods.eq.1)s% rsp_GREKM=0.d0
         EKDEL  = EKMAX-EKMIN
         EKMAXL = EKMAX
         EKMAX  =-10.d50
         EKMIN  =-EKMAX
         LMIN   = 10.d10
         LMAX   =-LMIN
         RMAX   = -SUNR
         RMIN   = -RMAX
         VMAX   = 0
         s% rsp_GREKM = 0
         s% rsp_GREKM_avg_abs = 0
         s% rsp_GRPDV = 0
         s% rsp_DeltaR = 0
         s% rsp_DeltaMAG = 0
         ID   = 0
         E0   = 0.d0
         call init_HYD(s,ierr)         
         if (ierr /= 0) return
         s% rsp_num_periods = 0
         call calculate_energies(s,total_radiation)
         E0 = EDE_start 
         call calculate_work_integrals(s)   
      end subroutine begin_calculation
      
      
      subroutine calculate_work_integrals(s)
         type (star_info), pointer :: s
         integer :: i, k
         real(dp) :: dt, dm, dVol, P_tw, Pvsc_tw, Ptrb_tw
         character (len=256) :: fname
         dt = s% dt
         ! LAST STEP OF PdV
         if(INSIDE.eq.1.and.IWORK.eq.1) then  
            IWORK=0
            do I=1,NZN
               k = NZN+1-i
               dm = s% dm(k)
               dVol = VV0(I) - s% Vol_start(k)
               P_tw = THETA*PPP0(I) + &
                  THETA1*(s% Pgas_start(k) + s% Prad_start(k))
               Pvsc_tw = THETA*PPQ0(I) + THETA1*s% Pvsc_start(k)
               Ptrb_tw = THETAT*PPT0(I) + THETAT1*s% Ptrb_start(k)
               WORK(I) = WORK(I) + &
                    dVol*s% dm(k)*(P_tw + Pvsc_tw + Ptrb_tw) &
                  - dt*dm*s% Eq(k)
               WORKQ(I) = WORKQ(I) + dVol*dm*Pvsc_tw
               WORKT(I) = WORKT(I) + dVol*dm*Ptrb_tw
               WORKC(I) = WORKC(I) - dt*dm*s% Eq(k)
            enddo
            if (s% rsp_num_periods == s% RSP_work_period) then
               fname = trim(s% log_directory) // '/' // trim(s% RSP_work_filename)
               write(*,*) 'write work integral data to ' // trim(fname)
               open(78,file=trim(fname),status='unknown')
               do I=1,NZN
                  k = NZN+1-i
                  write(78,'(I3,tr1,F7.4,4(tr1,d16.10))') &
                     I, log10(sum(s% dm(1:k))), &
                     WORK(I)/EKDEL, WORKQ(I)/EKDEL, &
                     WORKT(I)/EKDEL, WORKC(I)/EKDEL
               enddo
               close(78)
            end if
            PDVWORK=0.d0
            do I=1,NZN
               k = NZN+1-i
               PDVWORK=PDVWORK+WORK(I)
               WORK(I)=0.d0    
               WORKQ(I)=0.d0
               WORKT(I)=0.d0
               WORKC(I)=0.d0
            enddo
            s% rsp_GRPDV=PDVWORK/EKDEL
         endif

         ! INITIAL STEP OF PdV:
         if((INSIDE.eq.1.and.IWORK.eq.0).or. &
            (s% rsp_num_periods.eq.0.and.IWORK.eq.0))then
            IWORK=1 
            do I=1,NZN
               k = NZN+1-i
               VV0(I) = s% Vol_start(k)
               PPP0(I) = s% Pgas_start(k) + s% Prad_start(k)
               PPQ0(I) = s% Pvsc_start(k)
               PPT0(I) = s% Ptrb_start(k)
               PPC0(I) = s% Chi_start(k)
            enddo      
         endif

         ! FIRST AND NEXT STEPS of PdV:
         if(IWORK.eq.1)then
            do I=1,NZN
               k = NZN+1-i
               dm = s% dm(k)
               dVol = s% Vol(k) - s% Vol_start(k)
               P_tw = THETA*(s% Pgas(k) + s% Prad(k)) &
                  + THETA1*(s% Pgas_start(k) + s% Prad_start(k))
               Pvsc_tw = THETA*s% Pvsc(k) + THETA1*s% Pvsc_start(k)
               Ptrb_tw = THETAT*s% Ptrb(k) + THETAT1*s% Ptrb_start(k)
               WORK(I) = WORK(I) + &
                  dm*(dVol*(P_tw + Pvsc_tw + Ptrb_tw) - dt*s% Eq(k))
               WORKQ(I)=  WORKQ(I) + dm*dVol*Pvsc_tw
               WORKT(I)=  WORKT(I) + dm*dVol*Ptrb_tw
               WORKC(I)=  WORKC(I) - dt*dm*s% Eq(k)
            enddo       
         endif
      end subroutine calculate_work_integrals


      subroutine rsp_total_energy_integrals(s, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total, total_radiation)
         type (star_info), pointer :: s
         real(dp), intent(out) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total, total_radiation
         include 'formats'
         call calculate_energies(s,total_radiation)
         total_internal_energy = ETHE
         total_gravitational_energy = EGRV
         total_radial_kinetic_energy = EKIN
         total_rotational_kinetic_energy = 0d0
         total_turbulent_energy = ECON
         sum_total = ETOT
      end subroutine rsp_total_energy_integrals

      
      subroutine check_cycle_completed(s,cycle_complete)
         ! uses ULL, FIRST, s% RSP_min_max_R_for_periods, PERIODLIN, s% RSP_min_PERIOD_div_PERIODLIN
         ! depends on s% RSP_have_set_velocities = .true.
         ! sets s% rsp_num_periods, s% rsp_period
         
         type (star_info), pointer :: s
         logical, intent(out) :: cycle_complete
         real(dp) :: TET, min_PERIOD
         include 'formats'
         TET = s% time
         cycle_complete = .false.
         UN=s% v(1)
         if(UN.gt.0.d0.and.ULL.le.0.d0) then
            RMIN=s% r(1)/SUNR
         end if
         if (s% model_number.eq.1) return
         if (.not. s% RSP_have_set_velocities) return
         if (s% RSP_min_max_R_for_periods > 0d0 .and. &
               s% r(1)/SUNR < s% RSP_min_max_R_for_periods) return
         if (UN/s% csound(1) > VMAX) then
            VMAX = UN/s% csound(1)
         end if
         if(UN*ULL.gt.0.0d0.or.UN.gt.0.d0) return
         T0=TET
         min_PERIOD = PERIODLIN*s% RSP_min_PERIOD_div_PERIODLIN
         if (abs(UN-ULL).gt.1.0d-10) T0=TE_start-(TE_start-TET)*ULL/(ULL-UN)
         if (min_PERIOD > 0d0 .and. T0-TT1 < min_PERIOD) return
         if(FIRST.eq.1)then   
            cycle_complete = .true.
            s% rsp_num_periods=s% rsp_num_periods+1
            s% rsp_period=T0-TT1
            if (is_bad(s% rsp_period)) then
               write(*,2) 'PERIOD', s% model_number, s% rsp_period
               write(*,2) 'TET', s% model_number, TET
               write(*,2) 'T0', s% model_number, T0
               write(*,2) 'TT1', s% model_number, TT1
               stop 'check_cycle_completed'
            end if
            RMAX=s% r(1)/SUNR !(NOT INTERPOLATED)
            write(*,'(a7,i7,f11.5,a9,f11.5,a14,f9.5,a9,i3,a7,i6,a16,f9.5,a6,i10,a6,f10.3)')  &
               'period', s% rsp_num_periods, s% rsp_period/(24*3600), &
               'delta R', RMAX - RMIN, &
               'max vsurf/cs', VMAX, &
               'retries',  &
                  s% num_retries - run_num_retries_prev_period,     &
               'steps',  &
                  s% model_number - prev_cycle_run_num_steps, &
               'avg iters/step',  &
                  dble(s% total_num_solver_iterations - run_num_iters_prev_period)/ &
                  dble(s% model_number - prev_cycle_run_num_steps), &
               'step', s% model_number, 'days', s% time/(24*3600)
            prev_cycle_run_num_steps = s% model_number
            run_num_iters_prev_period = s% total_num_solver_iterations
            run_num_retries_prev_period = s% num_retries
            TT1=T0
            INSIDE=1 
            VMAX = 0d0
         else
             write(*,*) 'first maximum radius, period calculations start at model, day', &
                s% model_number, s% time/(24*3600)
             TT1=T0
             FIRST=1
             ID=0
         endif
      end subroutine check_cycle_completed
      
 
      end module rsp
      
