! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
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

      module hydro_rsp2

      use star_private_def
      use const_def
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use accurate_sum_auto_diff_star_order1
      use star_utils

      implicit none

      private
      public :: &
         do1_rsp2_L_eqn, do1_turbulent_energy_eqn, do1_rsp2_Hp_eqn, &
         compute_Eq_cell, compute_Uq_face, set_RSP2_vars, &
         Hp_face_for_rsp2_val, Hp_face_for_rsp2_eqn, set_etrb_start_vars, &
         RSP2_adjust_vars_before_call_solver, get_RSP2_alfa_beta_face_weights
      
      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, & ! Ptrb
         x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3, & ! PII_face and Lc
         x_ALFAC = (1.d0/2.d0)*sqrt_2_div_3, & ! Lc
         x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3, & ! DAMP
         x_GAMMAR = 2.d0*sqrt(3.d0) ! DAMPR

      contains
      
      
      subroutine set_RSP2_vars(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr    
         type(auto_diff_real_star_order1) :: x
         integer :: k, op_err
         include 'formats'      
         ierr = 0
         op_err = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=1,s%nz
            ! Hp_face(k) <= 0 means it needs to be set.  e.g., after read file
            if (s% Hp_face(k) <= 0) then
               s% Hp_face(k) = get_scale_height_face_val(s,k)
               s% xh(s% i_Hp,k) = s% Hp_face(k)
            end if
            x = compute_Y_face(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            x = compute_PII_face(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            !Pvsc           skip?
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 1', s% model_number
            return
         end if
         !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=1,s% nz
            x = compute_Chi_cell(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            x = compute_Eq_cell(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            x = compute_C(s, k, op_err) ! COUPL
            if (op_err /= 0) ierr = op_err
            x = compute_L_face(s, k, op_err) ! Lr, Lt, Lc
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 2', s% model_number
            return
         end if
         do k = 1, s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0; 
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
         do k = s% nz + 1 - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM) , s% nz
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0; 
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
      end subroutine set_RSP2_vars


      subroutine do1_rsp2_L_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr         
         type(auto_diff_real_star_order1) ::  &
            L_expected, L_actual,resid
         real(dp) :: scale, residual, L_start_max
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% ctrl% solver_test_partials_k)
         test_partials = .false.         
         if (.not. s% RSP2_flag) then
            ierr = -1
            return
         end if

         ierr = 0
         !L_expected = compute_L_face(s, k, ierr)
         !if (ierr /= 0) return        
         L_expected = s% Lr_ad(k) + s% Lc_ad(k) + s% Lt_ad(k)
         L_actual = wrap_L_00(s, k)  
         L_start_max = maxval(s% L_start(1:s% nz))
         scale = 1d0/L_start_max
         if (is_bad(scale)) then
            write(*,2) 'do1_rsp2_L_eqn scale', k, scale
            call mesa_error(__FILE__,__LINE__,'do1_rsp2_L_eqn')
         end if
         resid = (L_expected - L_actual)*scale         
      
         residual = resid%val
         s% equ(s% i_equL, k) = residual         
         if (test_partials) then
            s% solver_test_partials_val = residual 
         end if
         
         call save_eqn_residual_info(s, k, nvar, s% i_equL, resid, 'do1_rsp2_L_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_rsp2_L_eqn', s% solver_test_partials_var
         end if      
      end subroutine do1_rsp2_L_eqn
      

      subroutine do1_rsp2_Hp_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr         
         type(auto_diff_real_star_order1) ::  &
            Hp_expected, Hp_actual,resid
         real(dp) :: residual, Hp_start
         logical :: test_partials
         include 'formats'
         !test_partials = (k == s% ctrl% solver_test_partials_k)
         test_partials = .false.
         
         if (.not. s% RSP2_flag) then
            ierr = -1
            return
         end if

         ierr = 0
         Hp_expected = Hp_face_for_rsp2_eqn(s, k, ierr)
         if (ierr /= 0) return        
         Hp_actual = wrap_Hp_00(s, k)  
         Hp_start = s% Hp_face_start(k)
         resid = (Hp_expected - Hp_actual)/max(Hp_expected,Hp_actual)    
      
         residual = resid%val
         s% equ(s% i_equ_Hp, k) = residual         
         if (test_partials) then
            s% solver_test_partials_val = residual 
         end if
         
         if (residual > 1d3) then
         !$omp critical (hydro_rsp2_1)
            write(*,2) 'residual', k, residual
            write(*,2) 'Hp_expected', k, Hp_expected%val
            write(*,2) 'Hp_actual', k, Hp_actual%val
            call mesa_error(__FILE__,__LINE__,'do1_rsp2_Hp_eqn')
         !$omp end critical (hydro_rsp2_1)
         end if
         
         call save_eqn_residual_info(s, k, nvar, s% i_equ_Hp, resid, 'do1_rsp2_Hp_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_rsp2_Hp_eqn', s% solver_test_partials_var
         end if      
         
      end subroutine do1_rsp2_Hp_eqn
   
   
      real(dp) function Hp_face_for_rsp2_val(s, k, ierr) result(Hp_face) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face_ad
         ierr = 0
         Hp_face_ad = Hp_face_for_rsp2_eqn(s, k, ierr)
         if (ierr /= 0) return
         Hp_face = Hp_face_ad%val
      end function Hp_face_for_rsp2_val
      
   
      function Hp_face_for_rsp2_eqn(s, k, ierr) result(Hp_face) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face
         type(auto_diff_real_star_order1) :: &
            rho_face, area, dlnPeos, &
            r_00, Peos_00, d_00, Peos_m1, d_m1, Peos_div_rho, &
            d_face, Peos_face, alt_Hp_face, A
         real(dp) :: alfa, beta
         integer :: j
         include 'formats'         
         ierr = 0
         if (k > s% nz) then
            Hp_face = 1d0 ! not used
            return
         end if
         if (k > 1 .and. .not. s% ctrl% RSP2_assume_HSE) then
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
            area = 4d0*pi*pow2(wrap_r_00(s,k))
            dlnPeos = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
            Hp_face = -s% dm_bar(k)/(area*rho_face*dlnPeos)
         else
            r_00 = wrap_r_00(s, k) ! not time-centered in RSP
            d_00 = wrap_d_00(s, k)
            Peos_00 = wrap_Peos_00(s, k)
            if (k == 1) then
               Peos_div_rho = Peos_00/d_00
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
            else
               d_m1 = wrap_d_m1(s, k)
               Peos_m1 = wrap_Peos_m1(s, k)
               call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
               Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
               if (k==-104) then
                  write(*,3) 'RSP2 Hp P_div_rho Pdrho_00 Pdrho_m1', k, s% solver_iter, &
                     Hp_face%val, Peos_div_rho%val, Peos_00%val/d_00%val, Peos_m1%val/d_m1%val
                  !write(*,3) 'RSP2 Hp r2_div_Gm r_start r', k, s% solver_iter, &
                  !   Hp_face%val, pow2(r_00%val)/(s% cgrav(k)*s% m(k)), &
                  !   s% r_start(k), r_00%val
               end if
               if (s% ctrl% alt_scale_height_flag) then
                  call mesa_error(__FILE__,__LINE__,'Hp_face_for_rsp2_eqn: cannot use alt_scale_height_flag')
                  ! consider sound speed*hydro time scale as an alternative scale height
                  d_face = alfa*d_00 + beta*d_m1
                  Peos_face = alfa*Peos_00 + beta*Peos_m1
                  alt_Hp_face = sqrt(Peos_face/s% cgrav(k))/d_face
                  if (alt_Hp_face%val < Hp_face%val) then ! blend
                     A = pow2(alt_Hp_face/Hp_face) ! 0 <= A%val < 1
                     Hp_face = A*Hp_face + (1d0 - A)*alt_Hp_face
                  end if
               end if
            end if
         end if
      end function Hp_face_for_rsp2_eqn


      subroutine do1_turbulent_energy_eqn(s, k, nvar, ierr)
         use star_utils, only: set_energy_eqn_scal, save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr         
         integer :: j
         ! for OLD WAY
         type(auto_diff_real_star_order1) :: &
            d_turbulent_energy_ad, Ptrb_dV_ad, dt_C_ad, dt_Eq_ad
         type(auto_diff_real_star_order1) :: xi0, xi1, xi2, A0, Af, w_00
         type(auto_diff_real_star_order1) :: tst, resid_ad, dt_dLt_dm_ad
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         logical :: non_turbulent_cell, test_partials
         real(dp) :: residual, atol, rtol, scal
         include 'formats'
         !test_partials = (k == s% ctrl% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         w_00 = wrap_w_00(s,k)
         
         non_turbulent_cell = &
            s% ctrl% mixing_length_alpha == 0d0 .or. &
            k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)         
         if (.not. s% RSP2_flag) then           
            resid_ad = w_00 - s% w_start(k) ! just hold w constant when not using RSP2
         else if (non_turbulent_cell) then
            resid_ad = w_00/s% csound(k) ! make w = 0         
         else 
            call setup_d_turbulent_energy(ierr); if (ierr /= 0) return ! erg g^-1 = cm^2 s^-2
            call setup_Ptrb_dV_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_C_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call setup_dt_Eq_ad(ierr); if (ierr /= 0) return ! erg g^-1
            call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)         
            ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
            esum_ad = d_turbulent_energy_ad + Ptrb_dV_ad + dt_dLt_dm_ad - dt_C_ad - dt_Eq_ad ! erg g^-1
            resid_ad = esum_ad
            
            if (k==-35 .and. s% solver_iter == 1) then
                  write(*,3) 'RSP2 w dEt PdV dtC dtEq', k, s% solver_iter, &
                     w_00%val, d_turbulent_energy_ad%val, Ptrb_dV_ad%val, dt_C_ad%val, dt_Eq_ad%val
            end if

            resid_ad = resid_ad*scal/s%dt ! to make residual unitless, must cancel out the dt in scal
            
         end if

         residual = resid_ad%val
         s% equ(s% i_detrb_dt, k) = residual

         if (test_partials) then
            tst = residual
            s% solver_test_partials_val = tst%val
            if (s% solver_iter == 12) &
               write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k), tst%val
         end if
         
         call save_eqn_residual_info(s, k, nvar, s% i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)     ! xi0 good , xi1 partial 0, xi2 good.  Af horrible.'
            write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k)/ln10, tst%val
         end if      

         contains
         
         subroutine setup_d_turbulent_energy(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            ierr = 0
            d_turbulent_energy_ad = wrap_etrb_00(s,k) - get_etrb_start(s,k)
         end subroutine setup_d_turbulent_energy
         
         ! Ptrb_dV_ad = Ptrb_ad*dV_ad
         subroutine setup_Ptrb_dV_ad(ierr) ! erg g^-1
            use star_utils, only: calc_Ptrb_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Ptrb_ad, PT0, dV_ad, d_00
            call calc_Ptrb_ad_tw(s, k, Ptrb_ad, PT0, ierr)
            if (ierr /= 0) return
            d_00 = wrap_d_00(s,k)
            dV_ad = 1d0/d_00 - 1d0/s% rho_start(k)
            Ptrb_dV_ad = Ptrb_ad*dV_ad ! erg cm^-3 cm^-3 g^-1 = erg g^-1
         end subroutine setup_Ptrb_dV_ad

         subroutine setup_dt_dLt_dm_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr            
            type(auto_diff_real_star_order1) :: Lt_00, Lt_p1, dLt_ad
            real(dp) :: Lt_00_start, Lt_p1_start, L_theta
            include 'formats'
            ierr = 0
            if (s% using_velocity_time_centering .and. &
                     s% ctrl% include_L_in_velocity_time_centering) then
               L_theta = s% ctrl% L_theta_for_velocity_time_centering
            else
               L_theta = 1d0
            end if
            Lt_00 = L_theta*s% Lt_ad(k) + (1d0 - L_theta)*s% Lt_start(k)
            if (k == s% nz) then
               Lt_p1 = 0d0
            else
               Lt_p1 = L_theta*shift_p1(s% Lt_ad(k+1)) + (1d0 - L_theta)*s% Lt_start(k+1)
               if (ierr /= 0) return
            end if
            dt_dLt_dm_ad = (Lt_00 - Lt_p1)*s%dt/s%dm(k)
         end subroutine setup_dt_dLt_dm_ad
         
         subroutine setup_dt_C_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: C
            C = s% COUPL_ad(k) ! compute_C(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_C_ad = s%dt*C
         end subroutine setup_dt_C_ad
                  
         subroutine setup_dt_Eq_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Eq_cell
            Eq_cell = s% Eq_ad(k) ! compute_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_Eq_ad = s%dt*Eq_cell
         end subroutine setup_dt_Eq_ad
      
      end subroutine do1_turbulent_energy_eqn
      
      
      subroutine get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: alfa, beta
         ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
         if (k == 1) call mesa_error(__FILE__,__LINE__,'bad k==1 for get_RSP2_alfa_beta_face_weights')
         if (s% ctrl% RSP2_use_mass_interp_face_values) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
         else
            alfa = 0.5d0
            beta = 0.5d0
         end if
      end subroutine get_RSP2_alfa_beta_face_weights

      
      function compute_Y_face(s, k, ierr) result(Y_face) ! superadiabatic gradient [unitless]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face
         type(auto_diff_real_star_order1) :: Hp_face, Y1, Y2, QQ_div_Cp_face, &
            r_00, d_00, Peos_00, Cp_00, T_00, chiT_00, chiRho_00, QQ_00, lnT_00, &
            r_m1, d_m1, Peos_m1, Cp_m1, T_m1, chiT_m1, chiRho_m1, QQ_m1, lnT_m1, &
            dlnT_dlnP, grad_ad_00, grad_ad_m1, grad_ad_face, dlnT, dlnP, alt_Y_face
         real(dp) :: dm_bar, alfa, beta
         include 'formats'
         ierr = 0
         
         if (k > s% nz) then
            Y_face = 0d0
            return
         end if
         
         if (k == 1 .or. s% ctrl% mixing_length_alpha == 0d0) then
            Y_face = 0d0
            s% Y_face(k) = 0d0
            s% Y_face_ad(k) = 0d0
            return
         end if
         
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         
         if (s% ctrl% RSP2_use_RSP_eqn_for_Y_face) then
      
            dm_bar = s% dm_bar(k)
            Hp_face = wrap_Hp_00(s,k)      
            r_00 = wrap_r_00(s, k)
            d_00 = wrap_d_00(s, k)
            Peos_00 = wrap_Peos_00(s, k)
            Cp_00 = wrap_Cp_00(s, k)
            T_00 = wrap_T_00(s, k)
            chiT_00 = wrap_chiT_00(s, k)
            chiRho_00 = wrap_chiRho_00(s, k)
            QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            lnT_00 = wrap_lnT_00(s,k)
      
            r_m1 = wrap_r_m1(s, k)
            d_m1 = wrap_d_m1(s, k)
            Peos_m1 = wrap_Peos_m1(s, k)
            Cp_m1 = wrap_Cp_m1(s, k)
            T_m1 = wrap_T_m1(s, k)
            chiT_m1 = wrap_chiT_m1(s, k)
            chiRho_m1 = wrap_chiRho_m1(s, k)
            QQ_m1 = chiT_m1/(d_m1*T_m1*chiRho_m1)
            lnT_m1 = wrap_lnT_m1(s,k)
            QQ_div_Cp_face = alfa*QQ_00/Cp_00 + beta*QQ_m1/Cp_m1
            ! QQ units (g cm^-3 K)^-1 = g^-1 cm^3 K^-1
            ! Cp units erg g^-1 K^-1 = g cm^2 s^-2 g^-1 K^-1 = cm^2 s^-2 K^-1
            ! QQ/Cp units = (g^-1 cm^3 K^-1)/(cm^2 s^-2 K^-1)
            !  = g^-1 cm^3 K^-1 cm^-2 s^2 K
            !  = g^-1 cm s^2
            ! P units = erg cm^-3 = g cm^2 s^-2 cm^-3 = g cm^-1 s^-2
            ! QQ/Cp*P is unitless.
         
            Y1 = QQ_div_Cp_face*(Peos_m1 - Peos_00) - (lnT_m1 - lnT_00)
            ! Y1 unitless
         
            Y2 = 4d0*pi*pow2(r_00)*Hp_face*2d0/(1d0/d_00 + 1d0/d_m1)/dm_bar
            ! units = cm^2 cm / (cm^3 g^-1) / g
            !       = cm^2 cm cm^-3 g g^-1 = unitless
         
            Y_face = Y1*Y2 ! unitless
            
            if (k==-35) then
               write(*,3) 'RSP2 Y_face Y1 Y2', k, s% solver_iter, s% Y_face(k), Y1%val, Y2%val
               write(*,3) 'Peos', k, s% solver_iter, Peos_00%val
               write(*,3) 'Peos', k-1, s% solver_iter, Peos_m1%val
               write(*,3) 'QQ', k, s% solver_iter, QQ_00%val
               write(*,3) 'QQ', k-1, s% solver_iter, QQ_m1%val
               write(*,3) 'Cp', k, s% solver_iter, Cp_00%val
               write(*,3) 'Cp', k-1, s% solver_iter, Cp_m1%val
               write(*,3) 'lgT', k, s% solver_iter, lnT_00%val/ln10
               write(*,3) 'lgT', k-1, s% solver_iter, lnT_m1%val/ln10
               write(*,3) 'lgd', k, s% solver_iter, s% lnd(k)/ln10
               write(*,3) 'lgd', k-1, s% solver_iter, s% lnd(k-1)/ln10
               !call mesa_error(__FILE__,__LINE__,'compute_Y_face')
            end if

         else
         
            grad_ad_00 = wrap_grad_ad_00(s,k)
            grad_ad_m1 = wrap_grad_ad_m1(s,k)
            grad_ad_face = alfa*grad_ad_00 + beta*grad_ad_m1
            dlnT = wrap_lnT_m1(s,k) - wrap_lnT_00(s,k)
            dlnP = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
            dlnT_dlnP = dlnT/dlnP
            if (is_bad(dlnT_dlnP%val)) then
               alt_Y_face = 0d0
            else if (s% ctrl% use_Ledoux_criterion .and. s% ctrl% calculate_Brunt_B) then
               ! gradL = grada + gradL_composition_term
               alt_Y_face = dlnT_dlnP - (grad_ad_face + s% gradL_composition_term(k))
            else
               alt_Y_face = dlnT_dlnP - grad_ad_face
            end if
            if (is_bad(alt_Y_face%val)) alt_Y_face = 0
            Y_face = alt_Y_face
            
         end if

         s% Y_face_ad(k) = Y_face
         s% Y_face(k) = Y_face%val

      end function compute_Y_face
      
      
      function compute_PII_face(s, k, ierr) result(PII_face) ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Cp_00, Cp_m1, Cp_face, Y_face
         real(dp) :: ALFAS_ALFA, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            PII_face = 0d0
            return
         end if
         if (k == 1 .or. s% ctrl% mixing_length_alpha == 0d0 .or. &
               k == s% nz) then ! just skip k == nz to be like RSP
            PII_face = 0d0
            s% PII(k) = 0d0
            s% PII_ad(k) = 0d0
            return
         end if
         Y_face = s% Y_face_ad(k) ! compute_Y_face(s, k, ierr)
         if (ierr /= 0) return
         Cp_00 = wrap_Cp_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         Cp_face = alfa*Cp_00 + beta*Cp_m1 ! ergs g^-1 K^-1
         ALFAS_ALFA = x_ALFAS*s% ctrl% mixing_length_alpha
         PII_face = ALFAS_ALFA*Cp_face*Y_face
         s% PII(k) = PII_face%val
         s% PII_ad(k) = PII_face
         if (k == -2 .and. s% PII(k) < 0d0) then
            write(*,2) 's% PII(k)', k, s% PII(k)
            write(*,2) 'Cp_face', k, Cp_face%val
            write(*,2) 'Y_face', k, Y_face%val
            !write(*,2) 'PII_face%val', k, PII_face%val
            !write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            call mesa_error(__FILE__,__LINE__,'compute_PII_face')
         end if
      end function compute_PII_face
      
      
      function compute_d_v_div_r(s, k, ierr) result(d_v_div_r) ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         v_00 = wrap_v_00(s,k)
         v_p1 = wrap_v_p1(s,k)
         r_00 = wrap_r_00(s,k)
         r_p1 = wrap_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1
      end function compute_d_v_div_r
      
      
      function compute_d_v_div_r_opt_time_center(s, k, ierr) result(d_v_div_r) ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         v_00 = wrap_opt_time_center_v_00(s,k)
         v_p1 = wrap_opt_time_center_v_p1(s,k)
         r_00 = wrap_opt_time_center_r_00(s,k)
         r_p1 = wrap_opt_time_center_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1
      end function compute_d_v_div_r_opt_time_center


      function wrap_Hp_cell(s, k) result(Hp_cell) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Hp_cell
         Hp_cell = 0.5d0*(wrap_Hp_00(s,k) + wrap_Hp_p1(s,k))
      end function wrap_Hp_cell
      
   
      function Hp_cell_for_Chi(s, k, ierr) result(Hp_cell) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell
         type(auto_diff_real_star_order1) :: &
            d_00, Peos_00, rmid
         real(dp) :: mmid, cgrav_mid
         integer :: j
         include 'formats'         
         ierr = 0
         
         Hp_cell = wrap_Hp_cell(s, k)
         return
         
         d_00 = wrap_d_00(s, k)
         Peos_00 = wrap_Peos_00(s, k)
         if (k < s% nz) then
            rmid = 0.5d0*(wrap_r_00(s,k) + wrap_r_p1(s,k))
            mmid = 0.5d0*(s% m(k) + s% m(k+1))
            cgrav_mid = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
         else
            rmid = 0.5d0*(wrap_r_00(s,k) + s% r_center)
            mmid = 0.5d0*(s% m(k) + s% m_center)
            cgrav_mid = s% cgrav(k)
         end if
         Hp_cell = pow2(rmid)*Peos_00/(d_00*cgrav_mid*mmid)
         if (s% ctrl% alt_scale_height_flag) then
            call mesa_error(__FILE__,__LINE__,'Hp_cell_for_Chi: cannot use alt_scale_height_flag')
         end if
      end function Hp_cell_for_Chi
      
      
      function compute_Chi_cell(s, k, ierr) result(Chi_cell) 
         ! eddy viscosity energy (Kuhfuss 1986) [erg]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            rho2, r6_cell, d_v_div_r, Hp_cell, w_00, d_00, r_00, r_p1
         real(dp) :: f, ALFAM_ALFA
         include 'formats'
         ierr = 0
         ALFAM_ALFA = s% ctrl% RSP2_alfam*s% ctrl% mixing_length_alpha
         if (ALFAM_ALFA == 0d0 .or. &
               k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
               k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            Chi_cell = 0d0
            if (k >= 1 .and. k <= s% nz) then
               s% Chi(k) = 0d0
               s% Chi_ad(k) = 0d0
            end if
         else
            Hp_cell = Hp_cell_for_Chi(s, k, ierr)
            if (ierr /= 0) return
            d_v_div_r = compute_d_v_div_r(s, k, ierr)
            if (ierr /= 0) return
            w_00 = wrap_w_00(s,k)
            d_00 = wrap_d_00(s,k)
            f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
            rho2 = pow2(d_00)
            r_00 = wrap_r_00(s,k)
            r_p1 = wrap_r_p1(s,k)
            r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
            Chi_cell = f*rho2*r6_cell*d_v_div_r*Hp_cell*w_00
            ! units = g^-1 cm s^-1 g^2 cm^-6 cm^6 s^-1 cm
            !       = g cm^2 s^-2
            !       = erg            
         end if
         s% Chi(k) = Chi_cell%val
         s% Chi_ad(k) = Chi_cell

      end function compute_Chi_cell

      
      function compute_Eq_cell(s, k, ierr) result(Eq_cell) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: d_v_div_r, Chi_cell
         include 'formats'
         ierr = 0
         if (s% ctrl% mixing_length_alpha == 0d0 .or. &
             k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            Eq_cell = 0d0
            if (k >= 1 .and. k <= s% nz) s% Eq_ad(k) = 0d0
         else
            Chi_cell = s% Chi_ad(k) ! compute_Chi_cell(s,k,ierr)
            if (ierr /= 0) return
            d_v_div_r = compute_d_v_div_r_opt_time_center(s, k, ierr)
            if (ierr /= 0) return
            Eq_cell = 4d0*pi*Chi_cell*d_v_div_r/s% dm(k) ! erg s^-1 g^-1
         end if
         s% Eq(k) = Eq_cell%val
         s% Eq_ad(k) = Eq_cell
      end function compute_Eq_cell


      function compute_Uq_face(s, k, ierr) result(Uq_face) ! cm s^-2, acceleration
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Uq_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_00
         include 'formats'
         ierr = 0         
         if (s% ctrl% mixing_length_alpha == 0d0 .or. &
             k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            Uq_face = 0d0
         else
            r_00 = wrap_opt_time_center_r_00(s,k)
            Chi_00 = s% Chi_ad(k) ! compute_Chi_cell(s,k,ierr)
            if (k > 1) then
               !Chi_m1 = shift_m1(compute_Chi_cell(s,k-1,ierr))
               Chi_m1 = shift_m1(s% Chi_ad(k-1))
               if (ierr /= 0) return
            else
               Chi_m1 = 0d0
            end if
            Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*s% dm_bar(k))
            
            if (k==-56) then
               write(*,3) 'RSP2 Uq chi_m1 chi_00 r', k, s% solver_iter, &
                  Uq_face%val, Chi_m1%val, Chi_00%val, r_00%val
            end if
            
         end if
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
         s% Uq(k) = Uq_face%val
      end function compute_Uq_face


      function compute_Source(s, k, ierr) result(Source) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source
         ! source_div_w assumes RSP2_source_seed == 0
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Peos_00, Cp_00, chiT_00, chiRho_00, QQ_00, &
            Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell, &
            grad_ad_00, P_QQ_div_Cp
         include 'formats'
         ierr = 0
         w_00 = wrap_w_00(s, k)
         T_00 = wrap_T_00(s, k)                  
         d_00 = wrap_d_00(s, k)         
         Peos_00 = wrap_Peos_00(s, k)         
         Cp_00 = wrap_Cp_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            
         Hp_face_00 = wrap_Hp_00(s,k)
         PII_face_00 = s% PII_ad(k) ! compute_PII_face(s, k, ierr)
         if (ierr /= 0) return
         
         if (k == s% nz) then
            PII_div_Hp_cell = PII_face_00/Hp_face_00
         else
            Hp_face_p1 = wrap_Hp_p1(s,k)
            if (ierr /= 0) return
            !PII_face_p1 = shift_p1(compute_PII_face(s, k+1, ierr))
            PII_face_p1 = shift_p1(s% PII_ad(k+1))
            if (ierr /= 0) return
            PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         end if
         
         ! Peos_00*QQ_00/Cp_00 = grad_ad if all perfect.
         !grad_ad_00 = wrap_grad_ad_00(s, k)
         P_QQ_div_Cp = Peos_00*QQ_00/Cp_00 ! use this to be same as RSP
         Source = (w_00 + s% ctrl% RSP2_source_seed)*PII_div_Hp_cell*T_00*P_QQ_div_Cp
         
         ! PII units same as Cp = erg g^-1 K^-1
         ! P*QQ/Cp is unitless (see Y_face)
         ! Source units = (erg g^-1 K^-1) cm^-1 cm s^-1 K
         !     = erg g^-1 s^-1
         
         if (k==-109) then
            write(*,3) 'RSP2 Source w PII_div_Hp T_P_QQ_div_Cp', k, s% solver_iter, &
               Source%val, w_00%val, PII_div_Hp_cell%val, T_00%val*P_QQ_div_Cp% val
            !write(*,3) 'RSP2 PII_00 PII_p1 Hp_00 Hp_p1', k, s% solver_iter, &
            !   PII_face_00%val, PII_face_p1%val, Hp_face_00%val, Hp_face_p1%val
         end if
         s% SOURCE(k) = Source%val

      end function compute_Source


      function compute_D(s, k, ierr) result(D) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D
         type(auto_diff_real_star_order1) :: dw3, w_00
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell
         include 'formats'
         ierr = 0
         if (s% ctrl% mixing_length_alpha == 0d0) then
            D = 0d0
         else
            Hp_cell = wrap_Hp_cell(s,k)
            w_00 = wrap_w_00(s,k)
            dw3 = pow3(w_00) - pow3(s% ctrl% RSP2_w_min_for_damping)
            D = (s% ctrl% RSP2_alfad*x_CEDE/s% ctrl% mixing_length_alpha)*dw3/Hp_cell
            ! units cm^3 s^-3 cm^-1 = cm^2 s^-3 = erg g^-1 s^-1
         end if
         if (k==-50) then
            write(*,3) 'RSP2 DAMP w Hp_cell dw3', k, s% solver_iter, &
               D%val, w_00%val, Hp_cell%val, dw3% val
         end if
         s% DAMP(k) = D%val
      end function compute_D


      function compute_Dr(s, k, ierr) result(Dr) ! erg g^-1 s^-1 = cm^2 s^-3
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Cp_00, kap_00, Hp_cell, POM2
         real(dp) :: gammar, alpha, POM
         include 'formats'
         ierr = 0
         alpha = s% ctrl% mixing_length_alpha
         gammar = s% ctrl% RSP2_alfar*x_GAMMAR
         if (gammar == 0d0) then
            Dr = 0d0
            s% DAMPR(k) = 0d0
            return
         end if
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Hp_cell = wrap_Hp_cell(s,k)
         POM = 4d0*boltz_sigma*pow2(gammar/alpha) ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00) 
            ! K^3 / ((g cm^-3)^2 (erg g^-1 K^-1) (cm^2 g^-1))
            ! K^3 / (cm^-4 erg K^-1) = K^4 cm^4 erg^-1
         Dr = get_etrb(s,k)*POM*POM2/pow2(Hp_cell)
         ! (erg cm^-2 K^-4 s^-1) (K^4 cm^4 erg^-1) cm^2 s^-2 cm^-2
         ! cm^2 s^-3 = erg g^-1 s^-1
         s% DAMPR(k) = Dr%val
      end function compute_Dr


      function compute_C(s, k, ierr) result(C) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: C
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source, D, Dr
         if (s% ctrl% mixing_length_alpha == 0d0 .or. &
             k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            if (k >= 1 .and. k <= s% nz) then
               s% SOURCE(k) = 0d0
               s% DAMP(k) = 0d0
               s% DAMPR(k) = 0d0
               s% COUPL(k) = 0d0
               s% COUPL_ad(k) = 0d0
            end if
            C = 0d0
            return
         end if
         Source = compute_Source(s, k, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         s% COUPL(k) = C%val
         s% COUPL_ad(k) = C
      end function compute_C


      function compute_L_face(s, k, ierr) result(L_face) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr, Lc, Lt
         call compute_L_terms(s, k, L_face, Lr, Lc, Lt, ierr)
      end function compute_L_face


      subroutine compute_L_terms(s, k, L, Lr, Lc, Lt, ierr)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L, Lr, Lc, Lt
         real(dp) :: L_val
         integer, intent(out) :: ierr         
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            L = 0d0
            L%val = s% L_center
            Lr = 0d0
            Lc = 0d0
            Lt = 0d0
            return
         end if
         Lr = compute_Lr(s, k, ierr)
         if (ierr /= 0) return
         if (k == 1) then            
            Lc = 0d0
            Lt = 0d0
         else
            Lc = compute_Lc(s, k, ierr)
            if (ierr /= 0) return
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
         end if
         L = Lr + Lc + Lt
         s% Lr_ad(k) = Lr
         s% Lc_ad(k) = Lc
         s% Lt_ad(k) = Lt
      end subroutine compute_L_terms


      function compute_Lr(s, k, ierr) result(Lr) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_00, T400, Erad, T_m1, T4m1, &
            kap_00, kap_m1, kap_face, diff_T4_div_kap, BW, BK
         real(dp) :: alfa
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lr = s% L_center
         else
            r_00 = wrap_r_00(s,k) ! not time centered
            area = 4d0*pi*pow2(r_00)
            T_00 = wrap_T_00(s,k)
            T400 = pow4(T_00)
            if (k == 1) then ! Lr(1) proportional to Erad in cell(1)
               Erad = crad * T400
               Lr = s% ctrl% RSP2_Lsurf_factor * area * clight * Erad
               s% Lr(k) = Lr%val
               return
            end if
            T_m1 = wrap_T_m1(s,k)
            T4m1 = pow4(T_m1)            
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            kap_00 = wrap_kap_00(s,k)
            kap_m1 = wrap_kap_m1(s,k)
            kap_face = alfa*kap_00 + (1d0 - alfa)*kap_m1
            diff_T4_div_kap = (T4m1 - T400)/kap_face

            if (s% ctrl% RSP2_use_Stellingwerf_Lr) then ! RSP style
               BW = log(T4m1/T400)
               if (abs(BW%val) > 1d-20) then
                  BK = log(kap_m1/kap_00)
                  if (abs(1d0 - BK%val/BW%val) > 1d-15 .and. abs(BW%val - BK%val) > 1d-15) then
                     diff_T4_div_kap = (T4m1/kap_m1 - T400/kap_00)/(1d0 - BK/BW)
                  end if
               end if
            end if
            Lr = -crad*clight/3d0*diff_T4_div_kap*pow2(area)/s% dm_bar(k)       
            ! units (erg cm^-3 K^-4) (cm s^-1) (K^4 cm^-2 g cm^4) g^-1 = erg s^-1  
         
            !s% xtra1_array(k) = s% T_start(k)
            !s% xtra2_array(k) = T4m1%val - T400%val
            !s% xtra3_array(k) = kap_face%val
            !s% xtra4_array(k) = diff_T4_div_kap%val
            !s% xtra5_array(k) = Lr%val/Lsun   
            !s% xtra6_array(k) = 1

         end if
         s% Lr(k) = Lr%val
      end function compute_Lr


      function compute_Lc(s, k, ierr) result(Lc) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lc_div_w_face
         Lc = compute_Lc_terms(s, k, Lc_div_w_face, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc


      function compute_Lc_terms(s, k, Lc_div_w_face, ierr) result(Lc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc, Lc_div_w_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area, &
            T_m1, T_00, d_m1, d_00, w_m1, w_00, T_rho_face, PII_face, w_face
         real(dp) :: ALFAC, ALFAS, alfa, beta
         include 'formats'
         ierr = 0
         if (s% ctrl% mixing_length_alpha == 0d0 .or. &
             k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            Lc = 0d0
            Lc_div_w_face = 1
            return
         end if
         r_00 = wrap_r_00(s, k)
         area = 4d0*pi*pow2(r_00)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)         
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         w_m1 = wrap_w_m1(s, k)
         w_00 = wrap_w_00(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         PII_face = s% PII_ad(k) ! compute_PII_face(s, k, ierr)
         w_face = alfa*w_00 + beta*w_m1
         ALFAC = x_ALFAC
         ALFAS = x_ALFAS
         Lc_div_w_face = area*(ALFAC/ALFAS)*T_rho_face*PII_face
         ! units = cm^2 K g cm^-3 ergs g^-1 K^-1 = ergs cm^-1
         Lc = w_face*Lc_div_w_face
         ! units = cm s^-1 ergs cm^-1 = ergs s^-1
         if (k == -458) then
            write(*,2) 'Lc%val', k, Lc%val
            write(*,2) 'w_face%val', k, w_face%val
            write(*,2) 'Lc_div_w_face', k, Lc_div_w_face%val
            write(*,2) 'PII_face%val', k, PII_face%val
            write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            call mesa_error(__FILE__,__LINE__,'compute_Lc_terms')
         end if
      end function compute_Lc_terms


      function compute_Lt(s, k, ierr) result(Lt) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lt
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area2, d_m1, d_00, &
            rho2_face, Hp_face, w_m1, w_00, w_face, etrb_m1, etrb_00
         real(dp) :: alpha_alpha_t, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lt = 0d0
            return
         end if 
         alpha_alpha_t = s% ctrl% mixing_length_alpha*s% ctrl% RSP2_alfat
         if (alpha_alpha_t == 0d0 .or. &
             k <= s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% ctrl% RSP2_nz_div_IBOTOM)) then
            Lt = 0d0
            s% Lt(k) = 0d0
            return
         end if
         r_00 = wrap_r_00(s,k)   
         area2 = pow2(4d0*pi*pow2(r_00))
         d_m1 = wrap_d_m1(s,k)
         d_00 = wrap_d_00(s,k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         rho2_face = alfa*pow2(d_00) + beta*pow2(d_m1)
         w_m1 = wrap_w_m1(s,k)
         w_00 = wrap_w_00(s,k)
         w_face = alfa*w_00 + beta*w_m1
         etrb_m1 = wrap_etrb_m1(s,k)
         etrb_00 = wrap_etrb_00(s,k)
         Hp_face = wrap_Hp_00(s,k)      
         ! Ft = - alpha_t * rho_face * alpha * Hp_face * w_face * detrb/dr (thesis eqn 2.44)
         ! replace dr by dm_bar/(area*rho_face)
         ! Ft = - alpha_alpha_t * rho_face * Hp_face * w_face * (area*rho_face) * detrb/dm_bar
         ! Lt = area * Ft
         ! Lt = -alpha_alpha_t * (area*rho_face)**2 * Hp_face * w_face * (etrb(k-1) - etrb(k))/dm_bar
         Lt = - alpha_alpha_t * area2 * rho2_face * Hp_face * w_face * (etrb_m1 - etrb_00) / s% dm_bar(k)  
         ! units = (cm^4) (g^2 cm^-6) (cm) (cm s^-1) (ergs g^-1) g^-1 = erg s^-1
         s% Lt(k) = Lt%val      
      end function compute_Lt


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, op_err
         type(auto_diff_real_star_order1) :: Y_face, Lt
         include 'formats'
         ierr = 0
         do k=1,s%nz
            Y_face = compute_Y_face(s, k, ierr)
            if (ierr /= 0) return
            s% Y_face_start(k) = Y_face%val  
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_start(k) = Lt%val  
            s% w_start(k) = s% w(k)
            s% Hp_face_start(k) = s% Hp_face(k)
         end do    
      end subroutine set_etrb_start_vars

      
      subroutine RSP2_adjust_vars_before_call_solver(s, ierr) ! replaces check_omega in RSP
         ! JAK OKRESLIC OMEGA DLA PIERWSZEJ ITERACJI
         use micro, only: do_eos_for_cell
         type (star_info), pointer :: s
         integer, intent(out) :: ierr    
         real(dp) :: PII_div_Hp, QQ, SOURCE, Hp_cell, DAMP, POM, POM2, DAMPR, del, soln
         type(auto_diff_real_star_order1) :: x
         integer :: k, op_err
         include 'formats'         
         ierr = 0
         if (s% ctrl% mixing_length_alpha == 0d0) return
         
         !$OMP PARALLEL DO PRIVATE(k,PII_div_Hp,QQ,SOURCE,Hp_cell,DAMP,POM,POM2,DAMPR,del,soln) SCHEDULE(dynamic,2)
         do k=s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent+1, &
               s% nz - max(1,int(s% nz/s% ctrl% RSP_nz_div_IBOTOM))
               
            if (s% w(k) > s% ctrl% RSP2_w_min_for_damping) cycle
            
            PII_div_Hp = 0.5d0*(s% PII(k)/s% Hp_face(k) + s% PII(k+1)/s% Hp_face(k+1))
            QQ = s% chiT(k)/(s% rho(k)*s% T(k)*s% chiRho(k))         
            SOURCE = PII_div_Hp*s% T(k)*s% Peos(k)*QQ/s% Cp(k)
            
            Hp_cell = 0.5d0*(s% Hp_face(k) + s% Hp_face(k+1))
            DAMP = (s% ctrl% RSP2_alfad*x_CEDE/s% ctrl% mixing_length_alpha)/Hp_cell            
            
            POM = 4d0*boltz_sigma*pow2(s% ctrl% RSP2_alfar*x_GAMMAR/s% ctrl% mixing_length_alpha)
            POM2 = pow3(s% T(k))/(pow2(s% rho(k))*s% Cp(k)*s% opacity(k)) 
            DAMPR = POM*POM2/pow2(Hp_cell)
            
            del = pow2(DAMPR) + 4d0*DAMP*SOURCE
            
            if (k==-35) then
               write(*,2) 'del', k, del
               write(*,2) 'DAMPR', k, DAMPR
               write(*,2) 'DAMP', k, DAMP
               write(*,2) 'SOURCE', k, SOURCE
               write(*,2) 'POM', k, PII_div_Hp
               write(*,2) 'POM2', k, s% T(k)*s% Peos(k)*QQ/s% Cp(k)
               write(*,2) 's% Hp_face(k)', k, s% Hp_face(k)
               write(*,2) 's% Hp_face(k+1)', k+1, s% Hp_face(k+1)
               write(*,2) 's% PII(k)', k, s% PII(k)
               write(*,2) 's% PII(k+1)', k+1, s% PII(k+1)
               write(*,2) 's% Y_face(k)', k, s% Y_face(k)
               write(*,2) 's% Y_face(k+1)', k+1, s% Y_face(k+1)
            end if
            
            if (del < 0d0) cycle
            soln = (-DAMPR + sqrt(del))/(2d0*DAMP)
            if (k==-35) write(*,2) 'soln', k, soln
            if (soln > 0d0) then
               ! i tried soln = sqrt(soln) here. helps solver convergence, but hurts the model results.
               if (s% ctrl% RSP2_report_adjust_w) &
                  write(*,3) 'RSP2_adjust_vars_before_call_solver w', k, s% model_number, s% w(k), soln
               s% w(k) = soln
            end if

         end do
         !$OMP END PARALLEL DO
      end subroutine RSP2_adjust_vars_before_call_solver




      end module hydro_rsp2

