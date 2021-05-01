! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton & The MESA Team
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


      module hydro_temperature

      use star_private_def
      use const_def
      use utils_lib, only: mesa_error, is_bad
      use auto_diff
      use auto_diff_support

      implicit none

      private
      public :: do1_alt_dlnT_dm_eqn, do1_gradT_eqn, do1_dlnT_dm_eqn

      contains
      


      ! just relate L_rad to T gradient.
      ! d_P_rad/dm = -<opacity_face>*L_rad/(clight*area^2) -- see, e.g., K&W (5.12)
      ! P_rad = (1/3)*crad*T^4
      ! d_P_rad/dm = (crad/3)*(T(k-1)^4 - T(k)^4)/dm_bar
      ! L_rad = L - L_non_rad, L_non_rad = L_start - L_rad_start
      ! L_rad_start = (-d_P_rad/dm_bar*clight*area^2/<opacity_face>)_start
      subroutine do1_alt_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         real(dp) :: alfa, beta, scale, dm_bar
         type(auto_diff_real_star_order1) :: L_ad, r_00, area, area2, Lrad_ad, &
            kap_00, kap_m1, kap_face, d_P_rad_expected_ad, T_m1, T4_m1, T_00, T4_00, &
            P_rad_m1, P_rad_00, d_P_rad_actual_ad, resid
         
         integer :: i_equL, i
         logical :: dbg
         logical :: test_partials

         include 'formats'
         ierr = 0
         i_equL = s% i_equL
         if (i_equL == 0) return

         if (.not. s% use_dPrad_dm_form_of_T_gradient_eqn) then
            ierr = -1
            return
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         dbg = .false.

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa

         scale = s% energy_start(k)*s% rho_start(k)
         dm_bar = s% dm_bar(k)
         L_ad = wrap_L_00(s,k)
         r_00 = wrap_r_00(s,k)
         area = pi4*pow2(r_00); area2 = pow2(area)

         if ((check_flag_and_val(s% conv_vel_flag, s% conv_vel, k)) .or. &
               (.not. s% conv_vel_flag .and. s% lnT(k)/ln10 <= s% max_logT_for_mlt &
               .and. s% mixing_type(k) == convective_mixing .and. s% gradr(k) > 0d0 &
               .and. abs(s% gradr(k) - s% gradT(k)) > abs(s% gradr(k))*1d-5)) then
            Lrad_ad = L_ad*s% gradT_ad(k)/s% gradr_ad(k) ! C&G 14.109
         else
            Lrad_ad = L_ad
         end if

         kap_00 = wrap_kap_00(s,k)
         kap_m1 = wrap_kap_m1(s,k)
         kap_face = alfa*kap_00 + beta*kap_m1
         if (kap_face%val < s% min_kap_for_dPrad_dm_eqn) &
            kap_face = s% min_kap_for_dPrad_dm_eqn
                  
         ! calculate expected d_P_rad from current L_rad
         d_P_rad_expected_ad = -dm_bar*kap_face*Lrad_ad/(clight*area2)
         
         ! calculate actual d_P_rad in current model
         T_m1 = wrap_T_m1(s,k); T4_m1 = pow4(T_m1)
         T_00 = wrap_T_00(s,k); T4_00 = pow4(T_00)

         !d_P_rad_expected = d_P_rad_expected*s% gradr_factor(k) !TODO(Pablo): check this

         P_rad_m1 = (crad/3d0)*T4_m1
         P_rad_00 = (crad/3d0)*T4_00
         d_P_rad_actual_ad = P_rad_m1 - P_rad_00
         
         ! residual
         resid = (d_P_rad_expected_ad - d_P_rad_actual_ad)/scale 
         s% equ(i_equL, k) = resid%val
         
         if (is_bad(resid%val)) then
!$OMP critical (star_alt_dlntdm_bad_num)
            write(*,2) 'resid%val', k, resid%val
            if (s% stop_for_bad_nums) stop 'do1_alt_dlnT_dm_eqn'
!$OMP end critical (star_alt_dlntdm_bad_num)
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% gradT(k)
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_alt_dlnT_dm_eqn', ierr)

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_alt_dlnT_dm_eqn', s% solver_test_partials_var
         end if

         contains 
         
         logical function check_flag_and_val(flag, array, index)
            logical,intent(in) :: flag
            real(dp), dimension(:),intent(in) :: array
            integer, intent(in) :: index

            check_flag_and_val = .false.
            if(flag) then
               if(array(index)>0d0) check_flag_and_val = .true.
            end if
         end function check_flag_and_val

      end subroutine do1_alt_dlnT_dm_eqn


      subroutine do1_gradT_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            resid, gradT, dlnT, dlnP
         integer :: i_equL
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_equL = s% i_equL
         if (i_equL == 0) return

         gradT = s% gradT_ad(k)
         dlnT = wrap_lnT_m1(s,k) - wrap_lnT_00(s,k)
         dlnP = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)

         resid = gradT*dlnP - dlnT
         s% equ(i_equL, k) = resid%val

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) stop 'do1_gradT_eqn'
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'gradT', k, gradT
            write(*,2) 'dlnT', k, dlnT
            write(*,2) 'dlnP', k, dlnP
            stop 'do1_gradT_eqn'
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_equL,k)
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_gradT_eqn', ierr)
         
         !call set_xtras
            
         contains
         
         subroutine set_xtras
            use auto_diff_support
            use star_utils, only: get_Lrad
            type(auto_diff_real_star_order1) :: &
               T4m1, T400, kap_m1, kap_00, alfa, beta, kap_face, &
               diff_T4_div_kap
            T4m1 = pow4(wrap_T_m1(s,k))
            T400 = pow4(wrap_T_00(s,k))
            kap_m1 = wrap_kap_m1(s,k)
            kap_00 = wrap_kap_00(s,k)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
            kap_face = alfa*kap_00 + beta*kap_m1
            diff_T4_div_kap = (T4m1 - T400)/kap_face
            s% xtra1_array(k) = s% T_start(k)
            s% xtra2_array(k) = T4m1%val - T400%val
            s% xtra3_array(k) = kap_face%val
            s% xtra4_array(k) = diff_T4_div_kap%val
            s% xtra5_array(k) = get_Lrad(s,k) 
            s% xtra6_array(k) = 1
         end subroutine set_xtras
         
      end subroutine do1_gradT_eqn


      subroutine do1_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: resid, &
            dlnPdm, Ppoint, gradT, dlnTdm, T00, Tm1, dT, Tpoint, lnTdiff
         real(dp) :: delm, alfa
         integer :: i_equL
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_equL = s% i_equL
         if (i_equL == 0) return
         
         if (s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn .or. &
               (s% X(k) <= s% max_X_for_gradT_eqn .and. s% max_X_for_gradT_eqn > 0d0)) then
            call do1_gradT_eqn(s, k, skip_partials, nvar, ierr)            
            return
         end if

         if (s% use_dPrad_dm_form_of_T_gradient_eqn .or. s% conv_vel_flag) then
            call do1_alt_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)            
            return
         end if
         
         ! dT/dm = dP/dm * T/P * grad_T, grad_T = dlnT/dlnP from MLT.
         ! but use hydrostatic value for dP/dm in this.
         ! this is because of limitations of MLT for calculating grad_T.
         ! (MLT assumes hydrostatic equilibrium)
         ! see comment in K&W chpt 9.1.
         
         call eval_dlnPdm_qhse(s, k, dlnPdm, Ppoint, ierr)
         if (ierr /= 0) return

         gradT = s% gradT_ad(k)
         dlnTdm = dlnPdm*gradT

         Tm1 = wrap_T_m1(s,k)
         T00 = wrap_T_00(s,k)
         dT = Tm1 - T00
         alfa = s% dm(k-1)/(s% dm(k-1) + s% dm(k))
         Tpoint = alfa*T00 + (1d0 - alfa)*Tm1
         lnTdiff = dT/Tpoint ! use this in place of lnT(k-1)-lnT(k)
         delm = (s% dm(k) + s% dm(k-1))/2
         
         resid = delm*dlnTdm - lnTdiff
         s% equ(i_equL, k) = resid%val

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) stop 'hydro eqns'
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'lnTdiff', k, lnTdiff
            write(*,2) 'delm', k, delm
            write(*,2) 'dlnPdm', k, dlnPdm
            write(*,2) 'gradT', k, gradT
            stop 'i_equL'
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_equL,k)
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_dlnT_dm_eqn', ierr)

      end subroutine do1_dlnT_dm_eqn



      ! only used for dlnT_dm equation
      subroutine eval_dlnPdm_qhse(s, k, & ! calculate the expected dlnPdm for HSE
            dlnPdm_qhse, Ppoint, ierr)
         use hydro_momentum, only: expected_HSE_grav_term
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dlnPdm_qhse, Ppoint
         integer, intent(out) :: ierr

         real(dp) :: alfa
         type(auto_diff_real_star_order1) :: grav, area, P00, Pm1
         include 'formats'

         ierr = 0

         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by Ppoint to make it unitless

         ! for rotation, multiply gravity by factor fp.  MESA 2, eqn 22.
         
         call expected_HSE_grav_term(s, k, grav, area, ierr)
         if (ierr /= 0) return
         
         P00 = wrap_Peos_00(s,k)
         if (s% using_velocity_time_centering) P00 = 0.5d0*(P00 + s% Peos_start(k))
         
         if (k == 1) then
            Pm1 = 0d0
            Ppoint = P00
         else
            Pm1 = wrap_Peos_m1(s,k)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*P00 + (1d0-alfa)*Pm1
         end if
         
         dlnPdm_qhse = grav/(area*Ppoint) ! note that expected_HSE_grav_term is negative
                  
         if (is_bad(dlnPdm_qhse%val)) then
            ierr = -1
            s% retry_message = 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
            if (s% report_ierr) then
!$OMP critical (hydro_vars_crit1)
               write(*,*) 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
               stop
!$OMP end critical (hydro_vars_crit1)
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'dlnPdm_qhse', k, dlnPdm_qhse
               stop 'eval_dlnPdm_qhse'
            end if
            return
         end if

      end subroutine eval_dlnPdm_qhse      
      
      end module hydro_temperature

