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
         
         call create_Noh_test_model(id, s, ierr)
         if (ierr /= 0) stop 'failed in create_Noh_test_model'
         
      end subroutine extras_controls
      
      
      subroutine create_Noh_test_model(id, s, ierr)
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
            xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx, &
            gamma, R_min, R_max, dr, rho, T, lnd, rho_0, v_0, &
            P, lnP, P_base, energy, lnE, r00, rmid, rp1, lnR, lnR00, &
            sum_dq, cell_vol, cell_mass, sum_cell_mass, M_star, lnR_min, lnR_max, dlnR, &
            E_inject, dm_inject, E_left_to_do, E_done, E_cell, energy_inject
         
         include 'formats'
         
         ierr = 0
         
         ! set basic controls
         s% u_flag = .true.
         s% RTI_flag = .false.
         s% constant_L = .true.

         gamma = 1.6667d0
         s% gamma_law_hydro = gamma
         
         log_mesh = s% split_merge_amr_log_zoning
         nz = s% split_merge_amr_nz_baseline
         s% nz = nz      

         v_0 = s% x_ctrl(1)
         R_min = s% x_ctrl(3)
         R_max = s% x_ctrl(4)
         rho_0 = s% x_ctrl(5)
         P = s% x_ctrl(6)
         lnP = log(P)
         
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
         
         M_star = 4*pi*rho_0*R_max*R_max*R_max/3
         s% R_center = 0d0
         
         do k = nz, 1, -1
            if (k==nz) then
               r00 = R_min
               rp1 = s% R_center
               lnR00 = lnR_min
            else 
               if (log_mesh) then
                  lnR00 = lnR00 + dlnR
                  rp1 = r00
                  r00 = exp(lnR00)
               else
                  rp1 = r00
                  r00 = rp1 + dr
                  lnR00 = log(r00)
               end if
            end if
            s% r(k) = r00
            s% lnR(k) = lnR00
            s% xa(1:s% species,k) = 0d0
            s% xa(h1,k) = 1d0
            s% u(k) = v_0
            s% alpha_RTI(k) = 0d0
            rho = rho_0
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
            s% lnP(k) = lnP
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
            s% zbar(k) = zbar
            s% z53bar(k) = z53bar
            s% lnE(k) = log(s% energy(k))
            s% lnT(k) = log(s% T(k))
         end do
         
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
         
         call star_write_model(id, 'noh_start.mod', ierr)
         write(*,2) 'noh_start.mod', nz
         if (ierr /= 0) stop 'failed in star_write_model' 

         write(*,*)
         write(*,*) 'finished create_Noh_test_model'
         stop
         
      end subroutine create_Noh_test_model
      
      
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
         real(dp) :: dt, min_dr
         character (len=strlen) :: test
         integer :: k, min_k
         logical :: okay
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         okay = .true.
         min_dr = 1
         min_k = 0
         do k=1, s% nz
            if (abs(s% r(k) - 0.1d0) < min_dr) then
               min_dr = abs(s% r(k) - 0.1d0)
               min_k = k
            end if
         end do
         do k = min(s% nz,min_k + 10), s% nz
            if (s% r(k) < 0.09d0) then
               if (abs(s% u(k)) > 0.001d0) then
                  okay = .false.
                  write(*,2) 'u bad', k, s% u(k), s% r(k)
               end if
               exit
            end if
         end do
         do k = max(1,min_k - 10), 1, -1
            if (s% r(k) > 0.11d0) then
               if (abs(s% u(k) + 1d0) > 0.001d0) then
                  okay = .false.
                  write(*,2) 'u bad', k, s% u(k), s% r(k)
               end if
               exit
            end if
         end do
         if (okay) write(*,*) 'velocities look reasonable near r=0.1'

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
      
