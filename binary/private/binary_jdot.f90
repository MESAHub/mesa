! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton, Pablo Marchant & The MESA Team
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


      module binary_jdot

      use const_def
      use star_lib
      use star_def
      use math_lib
      use binary_def

      implicit none

      contains

      real(dp) function get_jdot(b)
         type (binary_info), pointer :: b

         integer :: ierr
         
         ! calculate jdot from gravitational wave radiation
         if (.not. b% do_jdot_gr) then
             b% jdot_gr = 0d0
         else if (.not. b% use_other_jdot_gr) then
             call default_jdot_gr(b% binary_id, ierr)
         else
             call b% other_jdot_gr(b% binary_id, ierr)
         end if
            
         ! calculate jdot for mass ejected from system
         if (.not. b% do_jdot_ml) then
             b% jdot_ml = 0d0
         else if (.not. b% use_other_jdot_ml) then
             call default_jdot_ml(b% binary_id, ierr)
         else
             call b% other_jdot_ml(b% binary_id, ierr)
         end if

         ! solve jdot due to L-S coupling
         if (.not. b% do_jdot_ls) then
             b% jdot_ls = 0d0
         else if (.not. b% use_other_jdot_ls) then
             call default_jdot_ls(b% binary_id, ierr)
         else
             call b% other_jdot_ls(b% binary_id, ierr)
         end if

         ! solve jdot due to "missing wind" (see binary_controls.defaults)
         if (.not. b% do_jdot_missing_wind) then
             b% jdot_missing_wind = 0d0
         else if (.not. b% use_other_jdot_missing_wind) then
             call default_jdot_missing_wind(b% binary_id, ierr)
         else
             call b% other_jdot_missing_wind(b% binary_id, ierr)
         end if

         ! calculate jdot from magnetic braking
         if (.not. b% do_jdot_mb) then
             b% jdot_mb = 0d0
         else if (.not. b% use_other_jdot_mb) then
             call default_jdot_mb(b% binary_id, ierr)
         else
             call b% other_jdot_mb(b% binary_id, ierr)
         end if
         
         ! calculate extra jdot
         if (.not. b% use_other_extra_jdot) then
             b% extra_jdot = 0
         else 
             call b% other_extra_jdot(b% binary_id, ierr)
         end if
         
         get_jdot = (b% jdot_mb + b% jdot_gr + b% jdot_ml + b% jdot_missing_wind + &
            b% extra_jdot) * b% jdot_multiplier + b% jdot_ls
         
      end function get_jdot

      subroutine default_jdot_gr(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: bs4, clight5, cgrav3
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         bs4 = pow4(b% separation)
         clight5 = pow5(clight)
         cgrav3 = standard_cgrav*standard_cgrav*standard_cgrav
         b% jdot_gr = -32d0 * cgrav3 * b% m(b% a_i) * b% m(b% d_i) * (b% m(b% a_i) + b% m(b% d_i)) / &
             (5d0 * clight5 * bs4) * b% angular_momentum_j
      end subroutine default_jdot_gr

      subroutine default_jdot_ml(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: alfa
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         !mass lost from vicinity of donor
         b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             pow2(b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)*2*pi/b% period *&
             sqrt(1 - pow2(b% eccentricity))
         !mass lost from vicinity of accretor
         b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             pow2(b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)*2*pi/b% period *&
             sqrt(1 - pow2(b% eccentricity))
         !mass lost from circumbinary coplanar toroid
         b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
             sqrt(standard_cgrav * (b% m(1) + b% m(2)) * b% separation)
      end subroutine default_jdot_ml

      subroutine default_jdot_ls(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: delta_J
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_ls = 0
         ! ignore in first step, or if not doing rotation
         if (b% doing_first_model_of_run) &
            return
         ! bulk change in spin angular momentum takes tides into account
         delta_J = b% s_donor% total_angular_momentum_old - &
             b% s_donor% total_angular_momentum
         ! ignore angular momentum lost through winds
         if (b% s_donor% mstar_dot < 0) &
            delta_J = delta_J - b% s_donor% angular_momentum_removed * &
               abs(b% mdot_system_wind(b% d_i) / b% s_donor% mstar_dot)
         b% jdot_ls = b% jdot_ls + delta_J

         ! Repeat for accretor
         if (b% point_mass_i == 0) then
            delta_J = b% s_accretor% total_angular_momentum_old - &
               b% s_accretor% total_angular_momentum
            if (b% s_accretor% mstar_dot < 0) then
               ! all AM lost from the accretor is lost from the system
               delta_J = delta_J - b% s_accretor% angular_momentum_removed
            end if
            b% jdot_ls = b% jdot_ls + delta_J
         else if (b% model_twins_flag) then
            b% jdot_ls = b% jdot_ls + b% jdot_ls
         end if

         b% jdot_ls = b% jdot_ls / b% s_donor% dt
      end subroutine default_jdot_ls

      subroutine default_jdot_missing_wind(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_missing_wind = 0
         if (b% point_mass_i /= 0) return

         s => b% s_accretor

         if (s% mstar_dot < 0) then
            b% jdot_missing_wind = b% mtransfer_rate * b% fixed_xfer_fraction
         else
            b% jdot_missing_wind = b% mdot_system_wind(b% a_i)
         end if
         b% jdot_missing_wind = b% jdot_missing_wind * s% j_rot(1)

      end subroutine default_jdot_missing_wind

      subroutine check_jdot_mb_conditions(b, s, apply_jdot_mb, qconv_env)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         logical, intent(out) :: apply_jdot_mb
         real(dp), intent(out) :: qconv_env
         
         real(dp) :: qrad_core
         integer :: i, k, id

         include 'formats.inc'

         ! calculate how much of inner region is convective
         qrad_core = 0d0
         do k = s% nz, 1, -1
            if (s% q(k) > b% jdot_mb_qlim_for_check_rad_core .and. & 
               (qrad_core == 0d0 .or. s% mixing_type(k) /= convective_mixing)) exit
            if (s% mixing_type(k) == convective_mixing) &
               qrad_core = qrad_core + s% dq(k)
         end do

         ! calculate how much of the envelope
         qconv_env = 0d0
         do k = 1, s% nz
            if (s% q(k) < b% jdot_mb_qlim_for_check_conv_env .and. &
               (qconv_env == 0d0 .or. s% mixing_type(k) /= convective_mixing)) exit
            if (s% mixing_type(k) == convective_mixing) &
               qconv_env = qconv_env + s% dq(k)
         end do

         apply_jdot_mb = .true.
         if (qconv_env < b% jdot_mb_min_qconv_env) then
            apply_jdot_mb = .false.
            return
         end if

         if (qconv_env > b% jdot_mb_max_qconv_env) then
            apply_jdot_mb = .false.
            return
         end if

         if (qrad_core > b% jdot_mb_max_qrad_core) then
            apply_jdot_mb = .false.
            return
         end if
         
     end subroutine check_jdot_mb_conditions

      subroutine default_jdot_mb(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         logical :: apply_mdot_mb
         real(dp) :: rsun4,two_pi_div_p3, qconv_env, jdot_scale
         logical :: apply_jdot_mb
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_mb = 0
         rsun4 = pow4(rsun)
         two_pi_div_p3 = (2.0d0*pi/b% period)*(2.0d0*pi/b% period)*(2.0d0*pi/b% period)

         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
         call check_jdot_mb_conditions(b, b% s_donor, apply_jdot_mb, qconv_env)
         if (apply_jdot_mb .or. b% keep_mb_on) then
            jdot_scale = 1d0
            if (b% jdot_mb_scale_for_low_qconv_env) then
               !scale jdot for tiny convective envelope, from Podsiadlowski et al. (2002)
               !The Astrophysical Journal, Volume 565, Issue 2, pp. 1107-1133
               if (qconv_env > b% jdot_mb_mass_frac_for_scale) then
                  jdot_scale = 1d0
               else
                  jdot_scale = exp(-b% jdot_mb_mass_frac_for_scale/max(1d-99,qconv_env)+1)
               end if
            end if
            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
                           pow(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3*jdot_scale
            write(*,*) "check jdot_scale", 1, jdot_scale, b% jdot_mb
            b% using_jdot_mb(b% d_i) = .true.
            if ((apply_jdot_mb .or. b% keep_mb_on) .and. .not. b% using_jdot_mb_old(b% d_i)) then
               write(*,*) 'turn on magnetic braking for star ', b% d_i
            end if
         else if (.not. (apply_jdot_mb .or. b% keep_mb_on) .and. b% using_jdot_mb_old(b% d_i)) then
            ! required mdot for the implicit scheme may drop drastically,
            ! so its neccesary to increase change factor to avoid implicit 
            ! scheme from getting stuck
            b% change_factor = b% max_change_factor
            b% using_jdot_mb(b% d_i) = .false.
            write(*,*) 'turn off magnetic braking for star ', b% d_i
         end if

         if (b% point_mass_i == 0 .and. b% include_accretor_mb) then
            call check_jdot_mb_conditions(b, b% s_accretor, apply_jdot_mb, qconv_env)
            if (apply_jdot_mb .or. b% keep_mb_on) then
               jdot_scale = 1d0
               if (b% jdot_mb_scale_for_low_qconv_env) then
                  !scale jdot for tiny convective envelope, from Podsiadlowski et al. (2002)
                  !The Astrophysical Journal, Volume 565, Issue 2, pp. 1107-1133
                  if (qconv_env > b% jdot_mb_mass_frac_for_scale) then
                     jdot_scale = 1d0
                  else
                     jdot_scale = exp(-b% jdot_mb_mass_frac_for_scale/max(1d-99,qconv_env)+1)
                  end if
               end if
               b% jdot_mb = b% jdot_mb - &
                           3.8d-30*b% m(b% a_i)*rsun4* &
                           pow(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3*jdot_scale
               b% using_jdot_mb(b% a_i) = .true.
               if ((apply_jdot_mb .or. b% keep_mb_on) .and. .not. b% using_jdot_mb_old(b% a_i)) then
                  write(*,*) 'turn on magnetic braking for star ', b% a_i
               end if
            else if (.not. (apply_jdot_mb .or. b% keep_mb_on) .and. b% using_jdot_mb_old(b% a_i)) then
               ! required mdot for the implicit scheme may drop drastically,
               ! so its neccesary to increase change factor to avoid implicit 
               ! scheme from getting stuck
               b% change_factor = b% max_change_factor
               b% using_jdot_mb(b% a_i) = .false.
               write(*,*) 'turn off magnetic braking for star ', b% a_i
            end if
         else if (b% model_twins_flag .and. b% include_accretor_mb) then
            b% jdot_mb = 2*b% jdot_mb
            b% using_jdot_mb(b% a_i) = .true.
         end if

      end subroutine default_jdot_mb

      end module binary_jdot
