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

      module eps_grav

      use star_private_def
      use const_def
      use chem_def, only: chem_isos
      use utils_lib, only: mesa_error, is_bad
      use auto_diff

      implicit none

      private
      public :: eval_eps_grav_and_partials, zero_eps_grav_and_partials

      contains


      subroutine eval_eps_grav_and_partials(s, k, ierr)
         use auto_diff_support, only: wrap
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: f, d_eps_grav_dlnR00, d_eps_grav_dlnRp1
         include 'formats'
         ierr = 0

         if (s% dt <= 0) then
            call zero_eps_grav_and_partials(s, k)
            return
         end if

         call eval1_eps_grav_and_partials(s, k, ierr)
         if (ierr /= 0) return

         if (s% use_other_eps_grav) then
            ! note: call this after 1st doing the standard calculation
            
            stop 'need to change other_eps_grav to auto_diff'
            
            call s% other_eps_grav(s% id, k, s% dt, ierr)
            if (ierr /= 0) return
         end if

         f = s% eps_grav_factor
         if (abs(f - 1d0) > 1d-12) then
            s% eps_grav(k) = f*s% eps_grav(k)
            s% d_eps_grav_dlndm1(k) = f*s% d_eps_grav_dlndm1(k)
            s% d_eps_grav_dlnd00(k) = f*s% d_eps_grav_dlnd00(k)
            s% d_eps_grav_dlndp1(k) = f*s% d_eps_grav_dlndp1(k)
            s% d_eps_grav_dlnTm1(k) = f*s% d_eps_grav_dlnTm1(k)
            s% d_eps_grav_dlnT00(k) = f*s% d_eps_grav_dlnT00(k)
            s% d_eps_grav_dlnTp1(k) = f*s% d_eps_grav_dlnTp1(k)
            s% d_eps_grav_dlnR00(k) = f*s% d_eps_grav_dlnR00(k)
            s% d_eps_grav_dL00(k) = f*s% d_eps_grav_dL00(k)
            s% d_eps_grav_dLp1(k) = f*s% d_eps_grav_dLp1(k)
            s% d_eps_grav_dlnPgas00_const_T(k) = f*s% d_eps_grav_dlnPgas00_const_T(k)
            s% d_eps_grav_dlnPgasm1_const_T(k) = f*s% d_eps_grav_dlnPgasm1_const_T(k)
            s% d_eps_grav_dlnPgasp1_const_T(k) = f*s% d_eps_grav_dlnPgasp1_const_T(k)
            s% d_eps_grav_dlnTm1_const_Pgas(k) = f*s% d_eps_grav_dlnTm1_const_Pgas(k)
            s% d_eps_grav_dlnT00_const_Pgas(k) = f*s% d_eps_grav_dlnT00_const_Pgas(k)
            s% d_eps_grav_dlnTp1_const_Pgas(k) = f*s% d_eps_grav_dlnTp1_const_Pgas(k)
            s% d_eps_grav_dlnRp1(k) = f*s% d_eps_grav_dlnRp1(k)
            s% d_eps_grav_dv00(k) = f*s% d_eps_grav_dv00(k)
            s% d_eps_grav_dvp1(k) = f*s% d_eps_grav_dvp1(k)
         end if
         
         ! this is an interim solution so that users of eps_grav can convert to the auto_diff form 
         ! when change eps_grav to use auto_diff, can delete the s% d_eps_grav_d... 
         
         d_eps_grav_dlnR00 = s% d_eps_grav_dlnR00(k)
         d_eps_grav_dlnRp1 = s% d_eps_grav_dlnRp1(k)
         if (.not. s% solver_use_lnR) then
            d_eps_grav_dlnR00 = d_eps_grav_dlnR00/s% r(k)
            if (k < s% nz) then
               d_eps_grav_dlnRp1 = d_eps_grav_dlnRp1/s% r(k+1)
            else
               d_eps_grav_dlnRp1 = 0d0
            end if
         end if
         
         call wrap(s% eps_grav_ad(k), s% eps_grav(k), &
            s% d_eps_grav_dlndm1(k), s% d_eps_grav_dlnd00(k), s% d_eps_grav_dlndp1(k), &
            s% d_eps_grav_dlnTm1(k), s% d_eps_grav_dlnT00(k), s% d_eps_grav_dlnTp1(k), &
            0d0, 0d0, 0d0, &
            0d0, d_eps_grav_dlnR00, d_eps_grav_dlnRp1, &
            0d0, s% d_eps_grav_dv00(k), s% d_eps_grav_dvp1(k), &
            0d0, s% d_eps_grav_dL00(k), s% d_eps_grav_dLp1(k), &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)

      end subroutine eval_eps_grav_and_partials


      subroutine eval1_eps_grav_and_partials(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         logical :: using_PC

         include 'formats'
         ierr = 0

         using_PC = (s% eos_frac_PC(k) .gt. 0)

         if (using_PC .and. s% gam_start(k) >= s% Gamma_lnS_eps_grav_full_on) then
            call do_eps_grav_with_lnS(s, k, ierr)
         else if (using_PC .and. s% gam_start(k) > s% Gamma_lnS_eps_grav_full_off) then
            call blend_with_lnS_form(s, k, ierr)
         else
            call do_eps_grav_with_lnd(s, k, ierr)
         end if

         if (ierr /= 0 .or. is_bad(s% eps_grav(k))) then
            ierr = -1
            s% retry_message = 'failed in eval_eps_grav_and_partials'
            if (s% report_ierr) then
               write(*,2) &
                  'failed in eval_eps_grav_and_partials', k, s% eps_grav(k)
            end if
            if (s% stop_for_bad_nums) then
               stop 'eval1_eps_grav_and_partials'
            end if
            return
         end if

      end subroutine eval1_eps_grav_and_partials


      subroutine do_eps_grav_with_lnd(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         include 'formats'

         call do_eps_grav_with_lnd_Lagrangian(s, k, ierr)
         if (ierr /= 0) return

         if (s% include_composition_in_eps_grav) call include_composition_in_eps_grav(s, k)
         
      end subroutine do_eps_grav_with_lnd
      
      
      subroutine do_eps_grav_with_lnd_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: dlnd_dt, dlnT_dt
         dlnd_dt = s% dxh_lnd(k) * s% dVARDOT_dVAR
         dlnT_dt = s% dxh_lnT(k) * s% dVARDOT_dVAR
         call do_lnd_eps_grav(s, k, dlnd_dt, dlnT_dt, ierr)
      end subroutine do_eps_grav_with_lnd_Lagrangian


      ! this uses the given args to calculate -T*ds/dt
      subroutine do_lnd_eps_grav(s, k, dlnd_dt, dlnT_dt, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT, i_latent_ddlnT, i_latent_ddlnRho
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dlnd_dt, dlnT_dt
         integer, intent(out) :: ierr

         real(dp) :: dT1_dlnTdot, a1, da1_dlnd, da1_dlnT, &
            T1, dT1_dlnd, dT1_dlnT00, dT1_dlnd_dt, dT1_d_dlnTdt, &
            a2, da2_dlnd, da2_dlnT, &
            T2, dT2_dlnT, dT2_dlnd00, &
            T3, dT3_dlnd, dT3_dlnT, &
            latent_heat, d_latent_heat_dlnT, d_latent_heat_dlnRho
         logical :: test_partials

         real(dp) :: theta

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (s% use_time_centered_eps_grav) then
            theta = 0.5_dp
         else
            theta = 1.0_dp
         end if

         call zero_eps_grav_and_partials(s, k)
         
         !s% eps_grav(k) = -s% T(k)*s% cp(k)* &
         !      ((1-s% grada(k)*s% chiT(k))*dlnT_dt &
         !        - s% grada(k)*s% chiRho(k)*dlnd_dt)

         a1 = 1 - s% grada(k)*s% chiT(k)
         da1_dlnd = -(s% d_eos_dlnd(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnd(i_chiT,k))
         da1_dlnT = -(s% d_eos_dlnT(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnT(i_chiT,k))

         T1 = dlnT_dt*a1
         dT1_dlnd = dlnT_dt*da1_dlnd

         dT1_d_dlnTdt = a1
         dT1_dlnT00 = s% dVARDOT_dVAR*a1 + dlnT_dt*da1_dlnT

         a2 = s% grada(k)*s% chiRho(k)
         da2_dlnd = s% d_eos_dlnd(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnd(i_chiRho,k)
         da2_dlnT = s% d_eos_dlnT(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnT(i_chiRho,k)

         T2 = dlnd_dt*a2
         dT2_dlnT = dlnd_dt*da2_dlnT

         dT2_dlnd00 = s% dVARDOT_dVAR*a2 + dlnd_dt*da2_dlnd

         T3 = -s% T(k)*s% cp(k)
         dT3_dlnd = -s% T(k)*s% d_eos_dlnd(i_Cp,k)
         dT3_dlnT = -s% T(k)*(s% cp(k) + s% d_eos_dlnT(i_Cp,k))

         ! eps_grav = T3*(T1-T2)
         s% eps_grav(k) = T3*(T1-T2)

         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnd00(k) = (T3*(dT1_dlnd - dT2_dlnd00) + dT3_dlnd*(T1-T2))

         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnT00(k) = (T3*(dT1_dlnT00 - dT2_dlnT) + dT3_dlnT*(T1-T2))


         ! for time centered version
         if (s% use_time_centered_eps_grav) then

            a1 = 1 - s% grada_start(k)*s% chiT_start(k)
            da1_dlnd = 0d0
            da1_dlnT = 0d0

            T1 = dlnT_dt*a1
            dT1_dlnd = dlnT_dt*da1_dlnd

            dT1_d_dlnTdt = a1
            dT1_dlnT00 = s% dVARDOT_dVAR*a1 + dlnT_dt*da1_dlnT

            a2 = s% grada_start(k)*s% chiRho_start(k)
            da2_dlnd = 0d0
            da2_dlnT = 0d0

            T2 = dlnd_dt*a2
            dT2_dlnT = dlnd_dt*da2_dlnT

            dT2_dlnd00 = s% dVARDOT_dVAR*a2 + dlnd_dt*da2_dlnd

            T3 = -s% T_start(k)*s% cp_start(k)
            dT3_dlnd = 0d0
            dT3_dlnT = 0d0

            ! eps_grav = T3*(T1-T2)
            s% eps_grav(k) = theta * s% eps_grav(k) + (1d0-theta) * T3*(T1-T2)

            s% d_eps_grav_dlndm1(k) = 0
            s% d_eps_grav_dlndp1(k) = 0
            s% d_eps_grav_dlnd00(k) = theta * s% d_eps_grav_dlnd00(k) + (1d0-theta) * (T3*(dT1_dlnd - dT2_dlnd00) + dT3_dlnd*(T1-T2))

            s% d_eps_grav_dlnTm1(k) = 0
            s% d_eps_grav_dlnTp1(k) = 0
            s% d_eps_grav_dlnT00(k) = theta * s% d_eps_grav_dlnT00(k) + (1d0-theta) * (T3*(dT1_dlnT00 - dT2_dlnT) + dT3_dlnT*(T1-T2))

         end if

         ! phase transition latent heat
         latent_heat = dlnd_dt * s%latent_ddlnRho(k) + dlnT_dt * s%latent_ddlnT(k)
         d_latent_heat_dlnT = s% dVARDOT_dVAR * s%latent_ddlnT(k) + &
                        dlnT_dt * s% d_eos_dlnT(i_latent_ddlnT, k) + dlnd_dt * s% d_eos_dlnT(i_latent_ddlnRho, k)
         d_latent_heat_dlnRho = s% dVARDOT_dVAR * s%latent_ddlnRho(k) + &
                        dlnT_dt * s% d_eos_dlnd(i_latent_ddlnT, k) + dlnd_dt * s% d_eos_dlnd(i_latent_ddlnRho, k)

         s% eps_grav(k) = s%eps_grav(k) - latent_heat
         s% d_eps_grav_dlnT00(k) = s% d_eps_grav_dlnT00(k) - d_latent_heat_dlnT
         s% d_eps_grav_dlnd00(k) = s% d_eps_grav_dlnd00(k) - d_latent_heat_dlnRho

         if (is_bad(s% eps_grav(k))) then
            ierr = -1
            s% retry_message = 'do_lnd_eps_grav -- bad value for eps_grav'
            if (s% report_ierr) &
               write(*,2) 'do_lnd_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            if (s% stop_for_bad_nums) stop 'do_lnd_eps_grav'
            return
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% chiT(k) ! s% grada(k) ! a1 ! s% eps_grav(k)
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = s% d_eos_dlnd(i_chiT,k) ! s% d_eos_dlnd(i_grad_ad,k) ! da1_dlnd ! s% d_eps_grav_dlnd00(k)
            write(*,*) 'do_lnd_eps_grav chiT', s% solver_test_partials_var
         end if

         if (k == s% trace_k) then
            write(*,5) 'do_lnd_eps_grav', &
               k, s% solver_iter, s% solver_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,2) 'dlnd_dt', k, dlnd_dt
            write(*,2) 'dlnT_dt', k, dlnT_dt
            write(*,2) 's% cp(k)', k, s% cp(k)
            write(*,2) 's% grada(k)', k, s% grada(k)
            write(*,2) 's% chiT(k)', k, s% chiT(k)
            write(*,2) 'a1', k, a1
            write(*,2) 'da1_dlnd', k, da1_dlnd
            write(*,2) 'da1_dlnT', k, da1_dlnT
            write(*,2) 's% chiRho(k)', k, s% chiRho(k)
            write(*,2) 'a2', k, a2
            write(*,2) 'da2_dlnd', k, da2_dlnd
            write(*,2) 'da2_dlnT', k, da2_dlnT
            write(*,2) 'T1', k, T1
            write(*,2) 'T2', k, T2
            write(*,2) 'T3', k, T3
            write(*,2) 'dT1_dlnd', k, dT1_dlnd
            write(*,2) 'dT2_dlnd00', k, dT2_dlnd00
            write(*,2) 'dT3_dlnd', k, dT3_dlnd
            write(*,2) 'dT1_dlnT00', k, dT1_dlnT00
            write(*,2) 'dT2_dlnT', k, dT2_dlnT
            write(*,2) 'dT3_dlnT', k, dT3_dlnT
         end if

      end subroutine do_lnd_eps_grav


      subroutine do_eps_grav_with_lnS(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         include 'formats'

         call do_eps_grav_with_lnS_Lagrangian(s, k, ierr)
         if (ierr /= 0) return

         if (s% include_composition_in_eps_grav) call include_composition_in_eps_grav(s, k)

      end subroutine do_eps_grav_with_lnS


      subroutine do_eps_grav_with_lnS_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         include 'formats'
         call do_lnS_eps_grav(s, k, s% lnS_start(k), ierr)
      end subroutine do_eps_grav_with_lnS_Lagrangian


      subroutine do_lnS_eps_grav(s, k, lnS_start, ierr)
         use eos_def, only: i_lnS
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: lnS_start
         integer, intent(out) :: ierr

         real(dp) :: entropy, entropy_start, T, dS_dlnT, dS_dlnd

         include 'formats'
         ierr = 0
         call zero_eps_grav_and_partials(s, k)

         entropy = exp(s% lnS(k))
         T = s% T(k)
         entropy_start = exp(lnS_start)

         s% eps_grav(k) = -T*(entropy - entropy_start)*s% dVARDOT_dVAR

         if (is_bad(s% eps_grav(k))) then
            ierr = -1
            s% retry_message = 'do_lnS_eps_grav -- bad value for eps_grav'
            if (s% report_ierr) &
               write(*,2) 'do_lnS_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            if (s% stop_for_bad_nums) stop 'do_lnS_eps_grav'
            return
         end if

         dS_dlnT = s% dS_dT_for_partials(k)*s% T(k) ! instead of   entropy*s% d_eos_dlnT(i_lnS,k)
         dS_dlnd = s% dS_drho_for_partials(k)*s% rho(k) ! instead of   entropy*s% d_eos_dlnd(i_lnS,k)
         s% d_eps_grav_dlnT00(k) = -T*dS_dlnT*s% dVARDOT_dVAR + s% eps_grav(k)
         s% d_eps_grav_dlnd00(k) = -T*dS_dlnd*s% dVARDOT_dVAR

         if (k == s% trace_k) then
            write(*,5) 'do_lnS_eps_grav', &
               k, s% solver_iter, s% solver_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 'entropy', k, entropy
            write(*,2) 'T', k, T
            write(*,2) 'entropy_start', k, entropy_start
         end if

      end subroutine do_lnS_eps_grav


      subroutine blend_with_lnS_form(s, k, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa, Gamma
         include 'formats'
         ierr = 0
         Gamma = s% gam_start(k)
         alfa = (Gamma - s% Gamma_lnS_eps_grav_full_off) / &
            (s% Gamma_lnS_eps_grav_full_on - s% Gamma_lnS_eps_grav_full_off)
         ! alfa is fraction of lnS form in result
         if (alfa >= 1 .or. alfa <= 0) then
            ierr = -1
            s% retry_message = 'eval_eps_grav_and_partials -- error in blend_with_lnS_form'
          if (s% report_ierr) write(*, *) s% retry_message
            return
         end if
         call combine_two_eps_gravs( &
            s, k, alfa, 1d0 - alfa, do_eps_grav_with_lnS, do_eps_grav_with_lnd, ierr)
      end subroutine blend_with_lnS_form


      recursive subroutine combine_two_eps_gravs( &
            s, k, alfa, beta, eps_grav_proc1, eps_grav_proc2, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: alfa, beta
         interface
            subroutine eps_grav_proc1(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc1
            subroutine eps_grav_proc2(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc2
         end interface
         integer, intent(out) :: ierr

         real(dp) :: &
            eps_grav, d_eps_grav_dlndm1, d_eps_grav_dlnd00, d_eps_grav_dlndp1, &
            d_eps_grav_dlnTm1, d_eps_grav_dlnT00, d_eps_grav_dlnTp1, &
            d_eps_grav_dlnR00, d_eps_grav_dlnRp1, d_eps_grav_dL00, d_eps_grav_dLp1, &
            d_eps_grav_dlnPgas00_const_T, &
            d_eps_grav_dlnPgasm1_const_T, d_eps_grav_dlnPgasp1_const_T, &
            d_eps_grav_dlnTm1_const_Pgas, d_eps_grav_dlnT00_const_Pgas, &
            d_eps_grav_dlnTp1_const_Pgas, d_eps_grav_dv00, d_eps_grav_dvp1

         include 'formats'
         ierr = 0

         ! alfa is multiplier of result from calling eps_grav_proc1
         ! beta is multiplier of result from calling eps_grav_proc2
         ! i.e., eps_grav = alfa*eps_grav1 + beta*eps_grav2

         if (alfa > 1d0 .or. alfa < 0d0 .or. beta > 1d0 .or. beta < 0d0) then
            s% retry_message = 'failed in combine_two_eps_gravs'
            if (s% report_ierr) write(*, *) s% retry_message
            ierr = -1
            return
         end if

         ! result is alfa*eps_grav_proc1 + beta*eps_grav_proc2

         if (alfa > 0d0) then
            call eps_grav_proc1(s, k, ierr)
            if (ierr /= 0) return
            if (beta == 0d0) return
            ! save results
            eps_grav = s% eps_grav(k)
            d_eps_grav_dlndm1 = s% d_eps_grav_dlndm1(k)
            d_eps_grav_dlnd00 = s% d_eps_grav_dlnd00(k)
            d_eps_grav_dlndp1 = s% d_eps_grav_dlndp1(k)
            d_eps_grav_dlnTm1 = s% d_eps_grav_dlnTm1(k)
            d_eps_grav_dlnT00 = s% d_eps_grav_dlnT00(k)
            d_eps_grav_dlnTp1 = s% d_eps_grav_dlnTp1(k)
            d_eps_grav_dlnR00 = s% d_eps_grav_dlnR00(k)
            d_eps_grav_dlnRp1 = s% d_eps_grav_dlnRp1(k)
            d_eps_grav_dL00 = s% d_eps_grav_dL00(k)
            d_eps_grav_dLp1 = s% d_eps_grav_dLp1(k)
            d_eps_grav_dlnPgas00_const_T = s% d_eps_grav_dlnPgas00_const_T(k)
            d_eps_grav_dlnPgasm1_const_T = s% d_eps_grav_dlnPgasm1_const_T(k)
            d_eps_grav_dlnPgasp1_const_T = s% d_eps_grav_dlnPgasp1_const_T(k)
            d_eps_grav_dlnTm1_const_Pgas = s% d_eps_grav_dlnTm1_const_Pgas(k)
            d_eps_grav_dlnT00_const_Pgas = s% d_eps_grav_dlnT00_const_Pgas(k)
            d_eps_grav_dlnTp1_const_Pgas = s% d_eps_grav_dlnTp1_const_Pgas(k)
            d_eps_grav_dv00 = s% d_eps_grav_dv00(k)
            d_eps_grav_dvp1 = s% d_eps_grav_dvp1(k)
         else ! not needed, but to keep the compiler happy we set these to 0
            eps_grav = 0
            d_eps_grav_dlndm1 = 0
            d_eps_grav_dlnd00 = 0
            d_eps_grav_dlndp1 = 0
            d_eps_grav_dlnTm1 = 0
            d_eps_grav_dlnT00 = 0
            d_eps_grav_dlnTp1 = 0
            d_eps_grav_dlnR00 = 0
            d_eps_grav_dlnRp1 = 0
            d_eps_grav_dL00 = 0
            d_eps_grav_dLp1 = 0
            d_eps_grav_dlnPgas00_const_T = 0
            d_eps_grav_dlnPgasm1_const_T = 0
            d_eps_grav_dlnPgasp1_const_T = 0
            d_eps_grav_dlnTm1_const_Pgas = 0
            d_eps_grav_dlnT00_const_Pgas = 0
            d_eps_grav_dlnTp1_const_Pgas = 0
            d_eps_grav_dv00 = 0
            d_eps_grav_dvp1 = 0
         end if

         call eps_grav_proc2(s, k, ierr)
         if (ierr /= 0) return
         if (alfa == 0d0) return
         
         ! combine results
         s% eps_grav(k) = alfa*eps_grav + beta*s% eps_grav(k)

         s% d_eps_grav_dlndm1(k) = alfa*d_eps_grav_dlndm1 + beta*s% d_eps_grav_dlndm1(k)
         s% d_eps_grav_dlnd00(k) = alfa*d_eps_grav_dlnd00 + beta*s% d_eps_grav_dlnd00(k)
         s% d_eps_grav_dlndp1(k) = alfa*d_eps_grav_dlndp1 + beta*s% d_eps_grav_dlndp1(k)

         s% d_eps_grav_dlnTm1(k) = alfa*d_eps_grav_dlnTm1 + beta*s% d_eps_grav_dlnTm1(k)
         s% d_eps_grav_dlnT00(k) = alfa*d_eps_grav_dlnT00 + beta*s% d_eps_grav_dlnT00(k)
         s% d_eps_grav_dlnTp1(k) = alfa*d_eps_grav_dlnTp1 + beta*s% d_eps_grav_dlnTp1(k)

         s% d_eps_grav_dlnPgas00_const_T(k) = &
            alfa*d_eps_grav_dlnPgas00_const_T + beta*s% d_eps_grav_dlnPgas00_const_T(k)
         s% d_eps_grav_dlnPgasm1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasm1_const_T + beta*s% d_eps_grav_dlnPgasm1_const_T(k)
         s% d_eps_grav_dlnPgasp1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasp1_const_T + beta*s% d_eps_grav_dlnPgasp1_const_T(k)

         s% d_eps_grav_dlnTm1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTm1_const_Pgas + beta*s% d_eps_grav_dlnTm1_const_Pgas(k)
         s% d_eps_grav_dlnT00_const_Pgas(k) = &
            alfa*d_eps_grav_dlnT00_const_Pgas + beta*s% d_eps_grav_dlnT00_const_Pgas(k)
         s% d_eps_grav_dlnTp1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTp1_const_Pgas + beta*s% d_eps_grav_dlnTp1_const_Pgas(k)

         s% d_eps_grav_dlnR00(k) = alfa*d_eps_grav_dlnR00 + beta*s% d_eps_grav_dlnR00(k)
         s% d_eps_grav_dlnRp1(k) = alfa*d_eps_grav_dlnRp1 + beta*s% d_eps_grav_dlnRp1(k)

         s% d_eps_grav_dL00(k) = alfa*d_eps_grav_dL00 + beta*s% d_eps_grav_dL00(k)
         s% d_eps_grav_dLp1(k) = alfa*d_eps_grav_dLp1 + beta*s% d_eps_grav_dLp1(k)

         s% d_eps_grav_dv00(k) = alfa*d_eps_grav_dv00 + beta*s% d_eps_grav_dv00(k)
         s% d_eps_grav_dvp1(k) = alfa*d_eps_grav_dvp1 + beta*s% d_eps_grav_dvp1(k)

      end subroutine combine_two_eps_gravs


      subroutine include_composition_in_eps_grav(s, k)

         use eos_support, only: get_eos !, get_peos
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas

         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
            e, e_start, de, d_de_dlnd, d_de_dlnT, &
            e_with_xa_start, d_e_with_xa_start_dlnd, d_e_with_xa_start_dlnT, &
            e_with_DT_start, d_e_with_DT_start_dlnd, d_e_with_DT_start_dlnT
         real(dp), dimension(num_eos_basic_results) :: &
            res, dres_dlnd, dres_dlnT
         real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
         integer :: j, ierr
         logical :: test_partials

         real(dp) :: Pgas_with_DT_start

         real(dp) :: theta

         include 'formats'
         ierr = 0

         ! for now, bail if in new material
         if (k < s% k_below_just_added) return

         if (s% use_time_centered_eps_grav) then
            theta = 0.5_dp
         else
            theta = 1.0_dp
         end if

         ! some EOSes have composition partials and some do not
         ! those currently without dx partials are PC & Skye
         ! however, composition partials of lnE & lnP were revised by the call to
         ! fix_d_eos_dxa_partials at the beginning of eval_equ_for_solver

         ! directly fill array with composition derivatives
         ! we only have lnE and lnP derivs, so use
         ! eps_grav = -de/dt + (P/rho) * dlnd/dt
         do j=1, s% species
            s% d_eps_grav_dx(j,k) = -s% energy(k) * s% dlnE_dxa_for_partials(j,k) * s% dVARDOT_dVAR + &
               (s% P(k) / s% Rho(k)) * s% dlnP_dxa_for_partials(j,k) * s% dxh_lnd(k) * s% dVARDOT_dVAR
         end do

         e = s% energy(k)
         call get_eos( &
            s, k, s% xa_start(:,k), &
            s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
            res, dres_dlnd, dres_dlnT, &
            dres_dxa, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start', k
            return
         end if

         e_with_xa_start = exp(res(i_lnE))
         d_e_with_xa_start_dlnd = dres_dlnd(i_lnE)*e_with_xa_start
         d_e_with_xa_start_dlnT = dres_dlnT(i_lnE)*e_with_xa_start

         de = (e - e_with_xa_start)
         d_de_dlnd = (s% dE_dRho_for_partials(k)*s% Rho(k) - d_e_with_xa_start_dlnd)
         d_de_dlnT = (s% Cv_for_partials(k)*s% T(k) - d_e_with_xa_start_dlnT)

         if (s% use_time_centered_eps_grav) then

            e_start = s% energy_start(k)

            call get_eos( &
               s, k, s% xa(:,k), &
               s% rho_start(k), s% lnd_start(k)/ln10, s% T_start(k), s% lnT_start(k)/ln10, &
               res, dres_dlnd, dres_dlnT, &
               dres_dxa, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start', k
               return
            end if

            e_with_DT_start = exp(res(i_lnE))
            d_e_with_DT_start_dlnd = dres_dlnd(i_lnE)*e_with_DT_start
            d_e_with_DT_start_dlnT = dres_dlnT(i_lnE)*e_with_DT_start

            de = theta * de + (1d0 - theta) * (e_with_DT_start - e_start)
            d_de_dlnd = theta * d_de_dlnd
            d_de_dlnT = theta * d_de_dlnT

            Pgas_with_DT_start = exp(res(i_lnPgas))

            ! combine with end-of-step composition derivatives
            ! until EOSes always provide composition derivatives, ignore this
            ! the end-of-step term should generally be similar enough
            ! do j=1, s% species
            !    s% d_eps_grav_dx(j,k) = theta * s% d_eps_grav_dx(j,k) + (1d0 - theta) *  &
            !       (-e_with_DT_start * dres_dxa(i_lnE, j) * s% dVARDOT_dVAR + &
            !       (Pgas_with_DT_start/s% rho_start(k)) * dres_dxa(i_lnPgas, j) * s% dxh_lnd(k) * s% dVARDOT_dVAR)
            ! end do

         end if

         s% eps_grav_composition_term(k) = -de * s% dVARDOT_dVAR
         s% eps_grav(k) = s% eps_grav(k) + s% eps_grav_composition_term(k)

         s% d_eps_grav_dlnd00(k) = s% d_eps_grav_dlnd00(k) - d_de_dlnd * s% dVARDOT_dVAR
         s% d_eps_grav_dlnT00(k) = s% d_eps_grav_dlnT00(k) - d_de_dlnT * s% dVARDOT_dVAR


         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = de
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_de_dlnT
            write(*,*) 'get_dedt', s% solver_test_partials_var
         end if

         if (k == s% trace_k) then
            write(*,5) 'include_composition_in_eps_grav', &
               k, s% solver_iter, s% solver_adjust_iter, &
               s% model_number, s% eps_grav(k)
         end if

         if (is_bad(s% eps_grav(k))) then
          if (s% report_ierr) write(*, *) s% retry_message
            if (s% report_ierr) then
               write(*,2) 'include_composition_in_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
               write(*,2) 'eps_grav_composition_term', k, s% eps_grav_composition_term(k)
               !stop 'include_composition_in_eps_grav'
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'include_composition_in_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
               stop 'include_composition_in_eps_grav'
            end if
            return
         end if

      end subroutine include_composition_in_eps_grav


      subroutine zero_eps_grav_and_partials(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% eps_grav_ad(k) = 0
         s% eps_grav(k) = 0
         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlnd00(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnT00(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnPgasm1_const_T(k) = 0
         s% d_eps_grav_dlnPgas00_const_T(k) = 0
         s% d_eps_grav_dlnPgasp1_const_T(k) = 0
         s% d_eps_grav_dlnTm1_const_Pgas(k) = 0
         s% d_eps_grav_dlnT00_const_Pgas(k) = 0
         s% d_eps_grav_dlnTp1_const_Pgas(k) = 0
         s% d_eps_grav_dlnR00(k) = 0
         s% d_eps_grav_dlnRp1(k) = 0
         s% d_eps_grav_dL00(k) = 0
         s% d_eps_grav_dLp1(k) = 0
         s% d_eps_grav_dv00(k) = 0
         s% d_eps_grav_dvp1(k) = 0
         s% d_eps_grav_dx(:,k) = 0
         s% eps_grav_composition_term(k) = 0
      end subroutine zero_eps_grav_and_partials


      end module eps_grav

