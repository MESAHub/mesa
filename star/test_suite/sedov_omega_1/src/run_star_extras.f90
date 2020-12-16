! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      include "test_suite_extras_def.inc"
      
      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns 
         
         if (.not. s% x_logical_ctrl(1)) return
         
         call create_Sedov_test_model(id, s, ierr)
         if (ierr /= 0) stop 'failed in create_Sedov_test_model'
         
      end subroutine extras_controls
      
      
      subroutine create_Sedov_test_model(id, s, ierr)
         use eos_lib
         use chem_lib, only: basic_composition_info
         use chem_def, only: ih1, ihe4
         use utils_lib, only: is_bad
         
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         logical :: log_mesh
         integer :: i, k, nz, h1, he4, k_inject
         real(dp) :: &
            xh, xhe, z, abar, zbar, z53bar, z2bar, ye, mass_correction, sumx, &
            gamma, R_inject, R_min, R_max, dr, rho, T, lnd, rho_0, omega, P_floor, &
            P, P_base, energy, lnE, r00, rmid, rp1, lnR, lnR00, r_for_rho, &
            sum_dq, cell_vol, cell_mass, sum_cell_mass, M_star, lnR_min, lnR_max, dlnR, &
            E_inject, dm_inject, E_left_to_do, E_done, E_cell, energy_inject
         
         include 'formats'
         
         ierr = 0
         
         ! set basic controls
         !s% use_ODE_var_eqn_pairing = .true.
         s% u_flag = .true.
         s% RTI_flag = .false.
         s% constant_L = .true.

         gamma = 1.4d0
         s% gamma_law_hydro = gamma
         
         log_mesh = s% split_merge_amr_log_zoning
         nz = s% split_merge_amr_nz_baseline
         s% nz = nz      
         E_inject = s% x_ctrl(1)
         R_inject = s% x_ctrl(2)   
         R_min = s% x_ctrl(3)
         R_max = s% x_ctrl(4)
         rho_0 = s% x_ctrl(5)
         P_floor = s% x_ctrl(6)
         omega = s% x_ctrl(7)   
         P_base = s% x_ctrl(8)
         r_for_rho = s% x_ctrl(9)
         
         write(*,1) 'R_inject', R_inject
         write(*,1) 'omega', omega
         
         lnR_min = log(R_min)
         lnR_max = log(R_max)
         dr = (R_max - R_min)/(nz-1)
         dlnR = (lnR_max - lnR_min)/(nz-1)

         call star_set_net(id, 'basic_plus_fe56.net', ierr)
         if (ierr /= 0) stop 'failed in star_set_net'    
         
         call star_set_var_info(id, ierr)
         if (ierr /= 0) stop 'failed in star_set_var_info'     
         
         call star_set_chem_names(id, ierr)
         if (ierr /= 0) stop 'failed in star_set_chem_names'  
         
         call star_allocate_arrays(id, ierr)
         if (ierr /= 0) stop 'failed in star_allocate_arrays' 
         
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)

! for density gradient rho = rho_0 * radius**(-omega)
         
         s% R_center = 0d0
         E_done = 0
         k_inject = 0
         
         do k = nz, 1, -1
            if (k==nz) then
               r00 = R_min
               rp1 = s% R_center
               rmid = 0.5d0*r00
               lnR00 = lnR_min
            else 
               if (log_mesh) then
                  lnR00 = lnR00 + dlnR
                  rp1 = r00
                  r00 = exp(lnR00)
                  dr = r00 - rp1
               else
                  rp1 = r00
                  r00 = rp1 + dr
                  lnR00 = log(r00)
               end if
               rmid = r00-0.5d0*dr
            end if
            if (rmid > R_inject) then
               P = P_floor
            else
               k_inject = k
               P = P_floor + P_base*0.5d0*(1d0 + cos(pi*rmid/R_inject))
            end if
            s% r(k) = r00
            s% lnR(k) = lnR00
            s% xa(1:s% species,k) = 0d0
            s% xa(h1,k) = 1d0
            s% u(k) = 0
            s% alpha_RTI(k) = 0d0
            rho = rho_0/(pow(r_for_rho, omega) + pow(rmid, omega))
            s% rho(k) = rho
            s% lnd(k) = log(rho)
            s% dm(k) = rho*(4*pi/3)*(r00*r00*r00 - rp1*rp1*rp1)
            if (k == nz) then
               s% m(k) = s% M_center + s% dm(k)
            else
               s% m(k) = s% m(k+1) + s% dm(k)
            end if
            s% L(k) = 0d0
            s% P(k) = P
            s% lnP(k) = log(P)
            s% xa(1:s% species,k) = 0d0
            if (k == nz) then
               s% xa(he4,k) = 1d0 ! serves as marker of how much the innermost cell has expanded
            else
               s% xa(h1,k) = 1d0
            end if
            call basic_composition_info( &
               s% species, s% chem_id, s% xa(:,k), xh, xhe, z, &
               abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
            call eos_gamma_DP_get_ET( &
               abar, rho, P, gamma, s% energy(k), s% T(k), ierr)
            if (ierr /= 0) then
               stop 'failed in eos_gamma_DP_get_ET'
               return
            end if
            s% abar(k) = abar
            s% lnE(k) = log(s% energy(k))
            s% lnT(k) = log(s% T(k))
            if (k >= k_inject) then
               E_done = E_done + s% dm(k)*s% energy(k)
               !write(*,2) 'energy E_done', k, s% energy(k), E_done/E_inject, E_done, E_inject
            end if
         end do
         
         write(*,2) 'nz_inject r', nz-k_inject+1, s% r(k_inject)
                  
         sum_cell_mass = 0
         do k=1,nz ! sum small numbers 1st
            sum_cell_mass = sum_cell_mass + s% dm(k)
         end do
         
         write(*,1) 'sum_cell_mass/Msun', sum_cell_mass/Msun
         write(*,1) 'E_done', E_done
         write(*,1) 'E want', E_inject
         write(*,1) '(E_done - E)/E', (E_done - E_inject)/E_inject
         
         s% L_center = 0d0
         s% star_mass = s% m(1)/Msun
         s% mstar = s% m(1)
         s% xmstar = s% m(1) - s% M_center
         s% M_center = 0d0
         s% q(1) = 1d0
         sum_dq = 0
         do k=1,nz-1
            s% dq(k) = s% dm(k)/s% xmstar
            sum_dq = sum_dq + s% dq(k)
            s% q(k+1) = s% q(k) - s% dq(k)
         end do
         s% dq(nz) = s% q(nz)
         
         write(*,1) 'mass', s% m(1)
         write(*,2) 's% dq(nz)', nz, s% dq(k)
            
         s% model_number = 0
         s% star_age = 0
         s% initial_z = 0.02d0
         
         call star_write_model(id, 'sedov_start.mod', ierr)
         write(*,2) 'sedov_start.mod', nz
         if (ierr /= 0) stop 'failed in star_write_model' 

         write(*,*)
         write(*,*) 'finished create_Sedov_test_model'
         stop
         
      end subroutine create_Sedov_test_model
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         integer :: k
         logical :: AOK
         character (len=strlen) :: test
         real(dp) :: total_energy, r_shock, rho_shock, u_shock, energy_shock, p_shock
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% time < 1d0) return
         k = maxloc(s% u(1:s% nz), dim=1)         
         total_energy = 1.46687d0
         r_shock = 1.000000D+00
         rho_shock = 6.000000D+00
         u_shock = 4.166667D-01
         energy_shock = 8.680555D-02
         p_shock = 2.083333D-01
         
         AOK = &
            abs(s% time - 1d0) < 1d-14 .and. &
            abs(s% total_energy - total_energy) < 0.001d0 .and. &
            abs(s% r(k) - r_shock) < 0.005d0 .and. &
            abs(s% rho(k) - rho_shock) < 0.5d0 .and. &
            abs(s% u(k) - u_shock) < 0.01d0 .and. &
            abs(s% energy(k) - energy_shock) < 0.002d0 .and. &
            abs(s% P(k) - p_shock) < 0.01d0
            
         if (.not. AOK) then
            write(*,*) 'abs(s% time - 1d0) < 1d-14', abs(s% time - 1d0) < 1d-14
            write(*,*) 'abs(s% total_energy - total_energy) < 0.001', &
               abs(s% total_energy - total_energy) < 0.001d0, s% total_energy, total_energy
            write(*,*) 'abs(s% r(k) - r_shock) < 0.005', &
               abs(s% r(k) - r_shock) < 0.005d0, s% r(k), r_shock
            write(*,*) 'abs(s% rho(k) - rho_shock) < 0.5', &
               abs(s% rho(k) - rho_shock) < 0.5d0, s% rho(k), rho_shock
            write(*,*) 'abs(s% u(k) - u_shock) < 0.01', &
               abs(s% u(k) - u_shock) < 0.01d0, s% u(k), u_shock
            write(*,*) 'abs(s% energy(k) - energy_shock) < 0.002', &
               abs(s% energy(k) - energy_shock) < 0.002d0, s% energy(k), energy_shock
            write(*,*) 'abs(s% P(k) - p_shock) < 0.01', &
               abs(s% P(k) - p_shock) < 0.01d0, s% P(k), p_shock
         end if

         write(*,2) 'omega', k, s% x_ctrl(7)
         write(*,2) 's% time - 1d0', k, s% time - 1d0, s% time
         write(*,2) '(s% total_energy - total_energy)/total_energy', k, &
               (s% total_energy - total_energy)/total_energy, s% total_energy, total_energy
         write(*,2) '(s% r(k) - r_shock)/r_shock', k, &
               (s% r(k) - r_shock)/r_shock, s% r(k), r_shock
         write(*,2) '(s% rho(k) - rho_shock)/rho_shock', k, &
               (s% rho(k) - rho_shock)/rho_shock, s% rho(k), rho_shock
         write(*,2) '(s% u(k) - u_shock)/u_shock', k, &
               (s% u(k) - u_shock)/u_shock, s% u(k), u_shock
         write(*,2) '(s% energy(k) - energy_shock)/energy_shock', k, &
               (s% energy(k) - energy_shock)/energy_shock, s% energy(k), energy_shock
         write(*,2) '(s% P(k) - p_shock)/p_shock', k, &
               (s% P(k) - p_shock)/p_shock, s% P(k), p_shock

         if (AOK) then
            write(*,*) 'All tests are okay'
         else
            write(*,*) 'failed to satisfy some tests'
         end if

         call test_suite_after_evolve(s, ierr)

      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: t
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (n == 0) return
         
         names(1) = 'Pmax_P'
         names(2) = 'Pmax_r_1m13'
         names(3) = 'Pmax_v'
         names(4) = 'Pmax_rho'
         names(5) = 'Pmax_T'
         names(6) = 'log_Pmax_P'
         names(7) = 'log_Pmax_r'
         names(8) = 'log_Pmax_v'
         names(9) = 'log_Pmax_rho'
         names(10) = 'log_Pmax_T'
         names(11) = 'Pmax_r_div_t'
         names(12) = 'Pmax_m_div_Msun'
         
         k = maxloc(s% P(1:s% nz), dim=1)
         if (k == s% nz) k = s% nz-1
         t = s% time
         
         !write(*,2) 'Pmax k r*1d-13', k, 0.5d0*(s% r(k)+s% r(k+1))*1d-13
                  
         vals(1) = s% P(k) ! Pmax_P
         vals(2) = 0.5d0*(s% r(k)+s% r(k+1))*1d-13 ! Pmax_r_1m13
         vals(3) = 0.5d0*(s% u(k)+s% v(k+1)) ! Pmax_v
         vals(4) = s% rho(k) ! Pmax_rho
         vals(5) = s% T(k) ! Pmax_T
         vals(6) = log10(s% P(k)) ! Pmax_P
         vals(7) = log10(0.5d0*(s% r(k)+s% r(k+1))) ! Pmax_r
         vals(8) = log10(0.5d0*(s% u(k)+s% v(k+1))) ! Pmax_v
         vals(9) = log10(s% rho(k)) ! Pmax_rho
         vals(10) = log10(s% T(k)) ! Pmax_T
         vals(11) = 0.5d0*(s% r(k)+s% r(k+1))/s% time ! Pmax_r_div_t
         vals(12) = s% m(k)/Msun ! Pmax_m_div_Msun
         
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp, avo, kerg
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: gamma, energy
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (n == 0) return
         names(1) = 'Pressure_x_600'
         names(2) = 'Rho_x_10'
         names(3) = 'velocity_x_100'
         names(4) = 'energy_div_100'
         gamma = 1.4d0
         do k=1,s% nz
            vals(k,1) = s% P(k)*600d0
            vals(k,2) = s% rho(k)*10d0
            vals(k,3) = s% u(k)*100d0
            energy = s% P(k)/s% rho(k)/(gamma - 1d0)
            vals(k,4) = energy*1d-2
         end do
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      


      end module run_star_extras
      
