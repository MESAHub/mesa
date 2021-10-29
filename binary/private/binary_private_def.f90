! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module binary_private_def
      
      use binary_def

      implicit none

   ! history column options

      integer, parameter :: bh_model_number = 1
      integer, parameter :: bh_age = bh_model_number + 1
      integer, parameter :: bh_donor_index = bh_age + 1
      integer, parameter :: bh_period_days = bh_donor_index + 1
      integer, parameter :: bh_period_hr = bh_period_days + 1
      integer, parameter :: bh_period_minutes = bh_period_hr + 1
      integer, parameter :: bh_lg_separation = bh_period_minutes + 1
      integer, parameter :: bh_binary_separation = bh_lg_separation + 1
      integer, parameter :: bh_eccentricity =  bh_binary_separation + 1
      integer, parameter :: bh_star_1_radius = bh_eccentricity + 1
      integer, parameter :: bh_star_2_radius = bh_star_1_radius + 1
      integer, parameter :: bh_rl_1 = bh_star_2_radius + 1
      integer, parameter :: bh_rl_2 = bh_rl_1 + 1
      integer, parameter :: bh_rl_overflow_1 = bh_rl_2 + 1
      integer, parameter :: bh_rl_overflow_2 = bh_rl_overflow_1 + 1
      integer, parameter :: bh_rl_relative_overflow_1 = bh_rl_overflow_2 + 1
      integer, parameter :: bh_rl_relative_overflow_2 = bh_rl_relative_overflow_1 + 1
      integer, parameter :: bh_P_rot_div_P_orb_1 = bh_rl_relative_overflow_2 + 1
      integer, parameter :: bh_P_rot_div_P_orb_2 = bh_P_rot_div_P_orb_1 + 1
      integer, parameter :: bh_lg_t_sync_1 = bh_P_rot_div_P_orb_2 + 1
      integer, parameter :: bh_lg_t_sync_2 = bh_lg_t_sync_1 + 1
      integer, parameter :: bh_star_1_mass = bh_lg_t_sync_2 + 1
      integer, parameter :: bh_lg_star_1_mass = bh_star_1_mass + 1
      integer, parameter :: bh_star_2_mass = bh_lg_star_1_mass + 1
      integer, parameter :: bh_lg_star_2_mass = bh_star_2_mass + 1
      integer, parameter :: bh_sum_of_masses = bh_lg_star_2_mass + 1
      integer, parameter :: bh_lg_mtransfer_rate = bh_sum_of_masses + 1
      integer, parameter :: bh_lg_mstar_dot_1 = bh_lg_mtransfer_rate + 1
      integer, parameter :: bh_lg_mstar_dot_2 = bh_lg_mstar_dot_1 + 1
      integer, parameter :: bh_lg_system_mdot_1 = bh_lg_mstar_dot_2 + 1
      integer, parameter :: bh_lg_system_mdot_2 = bh_lg_system_mdot_1 + 1
      integer, parameter :: bh_lg_wind_mdot_1 = bh_lg_system_mdot_2 + 1
      integer, parameter :: bh_lg_wind_mdot_2 = bh_lg_wind_mdot_1 + 1
      integer, parameter :: bh_star_1_div_star_2_mass = bh_lg_wind_mdot_2 + 1
      integer, parameter :: bh_delta_star_1_mass = bh_star_1_div_star_2_mass + 1
      integer, parameter :: bh_delta_star_2_mass = bh_delta_star_1_mass + 1
      integer, parameter :: bh_fixed_xfer_fraction = bh_delta_star_2_mass + 1
      integer, parameter :: bh_eff_xfer_fraction = bh_fixed_xfer_fraction + 1
      integer, parameter :: bh_lg_mdot_edd = bh_eff_xfer_fraction + 1
      integer, parameter :: bh_mdot_edd_eta = bh_lg_mdot_edd + 1
      integer, parameter :: bh_lg_accretion_luminosity = bh_mdot_edd_eta + 1
      integer, parameter :: bh_bh_spin = bh_lg_accretion_luminosity + 1
      integer, parameter :: bh_v_orb_1 = bh_bh_spin + 1
      integer, parameter :: bh_v_orb_2 = bh_v_orb_1 + 1
      integer, parameter :: bh_lg_F_irr = bh_v_orb_2 + 1
      integer, parameter :: bh_J_orb = bh_lg_F_irr + 1
      integer, parameter :: bh_J_spin_1 = bh_J_orb + 1
      integer, parameter :: bh_J_spin_2 = bh_J_spin_1 + 1
      integer, parameter :: bh_J_total = bh_J_spin_2 + 1
      integer, parameter :: bh_Jdot = bh_J_total + 1
      integer, parameter :: bh_jdot_mb = bh_Jdot + 1
      integer, parameter :: bh_jdot_gr = bh_jdot_mb + 1
      integer, parameter :: bh_jdot_ml = bh_jdot_gr + 1
      integer, parameter :: bh_jdot_ls = bh_jdot_ml + 1
      integer, parameter :: bh_jdot_missing_wind = bh_jdot_ls + 1
      integer, parameter :: bh_extra_jdot = bh_jdot_missing_wind + 1
      integer, parameter :: bh_accretion_mode = bh_extra_jdot + 1
      integer, parameter :: bh_acc_am_div_kep_am = bh_accretion_mode + 1
      integer, parameter :: bh_edot =  bh_acc_am_div_kep_am + 1
      integer, parameter :: bh_edot_tidal =  bh_edot + 1
      integer, parameter :: bh_edot_enhance =  bh_edot_tidal + 1
      integer, parameter :: bh_extra_edot =  bh_edot_enhance + 1
      integer, parameter :: bh_point_mass_index =  bh_extra_edot + 1
      integer, parameter :: bh_ignore_rlof_flag =  bh_point_mass_index + 1
      integer, parameter :: bh_model_twins_flag =  bh_ignore_rlof_flag + 1
      integer, parameter :: bh_CE_flag =  bh_model_twins_flag + 1
      integer, parameter :: bh_CE_lambda1 =  bh_CE_flag + 1
      integer, parameter :: bh_CE_lambda2 =  bh_CE_lambda1 + 1
      integer, parameter :: bh_CE_Ebind1 =  bh_CE_lambda2 + 1
      integer, parameter :: bh_CE_Ebind2 =  bh_CE_Ebind1 + 1
      integer, parameter :: bh_CE_num1 =  bh_CE_Ebind2 + 1
      integer, parameter :: bh_CE_num2 =  bh_CE_num1 + 1
      
      integer, parameter :: bh_col_id_max = bh_CE_num2
      
      character (len=maxlen_binary_history_column_name) :: binary_history_column_name(bh_col_id_max)
      
      contains
      
      
      subroutine binary_history_column_names_init(ierr)
         integer, intent(out) :: ierr
         
         integer :: i, cnt
         ierr = 0
         cnt = 0
         binary_history_column_name(:) = ''

         binary_history_column_name(bh_model_number) = 'model_number'
         binary_history_column_name(bh_age) = 'age'
         binary_history_column_name(bh_donor_index) = 'donor_index'
         binary_history_column_name(bh_period_days) = 'period_days'
         binary_history_column_name(bh_period_hr) = 'period_hr'
         binary_history_column_name(bh_period_minutes) = 'period_minutes'
         binary_history_column_name(bh_lg_separation) = 'lg_separation'
         binary_history_column_name(bh_binary_separation) = 'binary_separation'
         binary_history_column_name(bh_eccentricity) = 'eccentricity'
         binary_history_column_name(bh_star_1_radius) = 'star_1_radius'
         binary_history_column_name(bh_star_2_radius) = 'star_2_radius'
         binary_history_column_name(bh_rl_1) = 'rl_1'
         binary_history_column_name(bh_rl_2) = 'rl_2'
         binary_history_column_name(bh_rl_overflow_1) = 'rl_overflow_1'
         binary_history_column_name(bh_rl_overflow_2) = 'rl_overflow_2'
         binary_history_column_name(bh_rl_relative_overflow_1) = 'rl_relative_overflow_1'
         binary_history_column_name(bh_rl_relative_overflow_2) = 'rl_relative_overflow_2'
         binary_history_column_name(bh_P_rot_div_P_orb_1) = 'P_rot_div_P_orb_1'
         binary_history_column_name(bh_P_rot_div_P_orb_2) = 'P_rot_div_P_orb_2'
         binary_history_column_name(bh_lg_t_sync_1) = 'lg_t_sync_1'
         binary_history_column_name(bh_lg_t_sync_2) = 'lg_t_sync_2'
         binary_history_column_name(bh_star_1_mass) = 'star_1_mass'
         binary_history_column_name(bh_lg_star_1_mass) = 'lg_star_1_mass'
         binary_history_column_name(bh_star_2_mass) = 'star_2_mass'
         binary_history_column_name(bh_lg_star_2_mass) = 'lg_star_2_mass'
         binary_history_column_name(bh_sum_of_masses) = 'sum_of_masses'
         binary_history_column_name(bh_lg_mtransfer_rate) = 'lg_mtransfer_rate'
         binary_history_column_name(bh_lg_mstar_dot_1) = 'lg_mstar_dot_1'
         binary_history_column_name(bh_lg_mstar_dot_2) = 'lg_mstar_dot_2'
         binary_history_column_name(bh_lg_system_mdot_1) = 'lg_system_mdot_1'
         binary_history_column_name(bh_lg_system_mdot_2) = 'lg_system_mdot_2'
         binary_history_column_name(bh_lg_wind_mdot_1) = 'lg_wind_mdot_1'
         binary_history_column_name(bh_lg_wind_mdot_2) = 'lg_wind_mdot_2'
         binary_history_column_name(bh_star_1_div_star_2_mass) = 'star_1_div_star_2_mass'
         binary_history_column_name(bh_delta_star_1_mass) = 'delta_star_1_mass'
         binary_history_column_name(bh_delta_star_2_mass) = 'delta_star_2_mass'
         binary_history_column_name(bh_fixed_xfer_fraction) = 'fixed_xfer_fraction'
         binary_history_column_name(bh_eff_xfer_fraction) = 'eff_xfer_fraction'
         binary_history_column_name(bh_lg_mdot_edd) = 'lg_mdot_edd'
         binary_history_column_name(bh_mdot_edd_eta) = 'mdot_edd_eta'
         binary_history_column_name(bh_lg_accretion_luminosity) = 'lg_accretion_luminosity'
         binary_history_column_name(bh_bh_spin) = 'bh_spin'
         binary_history_column_name(bh_v_orb_1) = 'v_orb_1'
         binary_history_column_name(bh_v_orb_2) = 'v_orb_2'
         binary_history_column_name(bh_lg_F_irr) = 'lg_F_irr'
         binary_history_column_name(bh_J_orb) = 'J_orb'
         binary_history_column_name(bh_J_spin_1) = 'J_spin_1'
         binary_history_column_name(bh_J_spin_2) = 'J_spin_2'
         binary_history_column_name(bh_J_total) = 'J_total'
         binary_history_column_name(bh_Jdot) = 'Jdot'
         binary_history_column_name(bh_jdot_mb) = 'jdot_mb'
         binary_history_column_name(bh_jdot_gr) = 'jdot_gr'
         binary_history_column_name(bh_jdot_ml) = 'jdot_ml'
         binary_history_column_name(bh_jdot_ls) = 'jdot_ls'
         binary_history_column_name(bh_jdot_missing_wind) = 'jdot_missing_wind'
         binary_history_column_name(bh_extra_jdot) = 'extra_jdot'
         binary_history_column_name(bh_accretion_mode) = 'accretion_mode'
         binary_history_column_name(bh_acc_am_div_kep_am) = 'acc_am_div_kep_am'
         binary_history_column_name(bh_edot) = 'edot'
         binary_history_column_name(bh_edot_tidal) = 'edot_tidal'
         binary_history_column_name(bh_edot_enhance) = 'edot_enhance'
         binary_history_column_name(bh_extra_edot) = 'extra_edot'
         binary_history_column_name(bh_point_mass_index) = 'point_mass_index'
         binary_history_column_name(bh_ignore_rlof_flag) = 'ignore_rlof_flag'
         binary_history_column_name(bh_model_twins_flag) = 'model_twins_flag'
         binary_history_column_name(bh_CE_flag) = 'CE_flag'
         binary_history_column_name(bh_CE_lambda1) = 'CE_lambda1'
         binary_history_column_name(bh_CE_lambda2) = 'CE_lambda2'
         binary_history_column_name(bh_CE_Ebind1) = 'CE_Ebind1'
         binary_history_column_name(bh_CE_Ebind2) = 'CE_Ebind2'
         binary_history_column_name(bh_CE_num1) = 'CE_num1'
         binary_history_column_name(bh_CE_num2) = 'CE_num2'
                  
         cnt = 0
         do i=1,bh_col_id_max
            if (len_trim(binary_history_column_name(i)) == 0) then
               write(*,*) 'missing name for log column id', i
               if (i > 1) write(*,*) 'following ' // trim(binary_history_column_name(i-1))
               write(*,*) 
               cnt = cnt+1
            end if
         end do

         if (cnt > 0) then
            ierr = -1
            return
         end if

      end subroutine binary_history_column_names_init         

      subroutine binary_private_def_init
         use num_def
         use utils_lib, only: get_compiler_version, get_mesasdk_version

         integer :: i      
         logical :: okay
         integer :: ierr
         
         include 'formats'
         
         okay = .true.
         ierr = 0

         binary_dt_why_str(1:b_numTlim) = ''
         
         binary_dt_why_str(b_Tlim_comp) = 'b_companion'
         binary_dt_why_str(b_Tlim_roche) = 'b_RL'
         binary_dt_why_str(b_Tlim_jorb) = 'b_jorb'
         binary_dt_why_str(b_Tlim_env) = 'b_envelope'
         binary_dt_why_str(b_Tlim_sep) = 'b_separation'
         binary_dt_why_str(b_Tlim_ecc) = 'b_eccentricity'
         binary_dt_why_str(b_Tlim_dm) = 'b_deltam'
         
         do i=1,b_numTlim
            if (len_trim(binary_dt_why_str(i)) == 0) then
               if (i > 1) then
                  write(*,2) 'missing binary_dt_why_str following ' // trim(binary_dt_why_str(i-1)), i
               else
                  write(*,2) 'missing binary_dt_why_str 1'
               end if
               okay = .false.
            end if
         end do
         
         if (.not. okay) call mesa_error(__FILE__,__LINE__,'binary_private_def_init')

         
         !here we store useful information about the compiler and SDK
         call get_compiler_version(compiler_name,compiler_version_name)
         call get_mesasdk_version(mesasdk_version_name,ierr)
         call date_and_time(date=date)
         
      end subroutine binary_private_def_init         
      
      integer function alloc_binary(ierr)
         integer, intent(out) :: ierr
         integer :: i
         type (binary_info), pointer :: b
         
         ierr = 0
         alloc_binary = -1
!$omp critical (binary_handle)
         if (.not. have_initialized_binary_handles) then
            do i = 1, max_binary_handles
               binary_handles(i)% binary_id = i
               binary_handles(i)% in_use = .false.
            end do
            have_initialized_binary_handles = .true.
         end if
         do i = 1, max_binary_handles
            if (.not. binary_handles(i)% in_use) then
               binary_handles(i)% in_use = .true.
               alloc_binary = i
               exit
            end if
         end do
!$omp end critical (binary_handle)
         if (alloc_binary == -1) then
            ierr = -1
            return
         end if
         if (binary_handles(alloc_binary)% binary_id /= alloc_binary) then
            ierr = -1
            return
         end if
         b => binary_handles(alloc_binary)
         
      end function alloc_binary
      
      
      subroutine free_binary(b)
         type (binary_info), pointer :: b
         binary_handles(b% binary_id)% in_use = .false.
      end subroutine free_binary
      

      end module binary_private_def

