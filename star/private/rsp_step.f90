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

      module rsp_step
      use rsp_def
      use star_def, only: star_ptr, star_info
      use const_def, only: dp, qp, crad
      use utils_lib, only: is_bad

      implicit none
      
      private
      public :: calculate_energies, init_HYD, HYD, &
         turn_off_time_weighting, turn_on_time_weighting, &
         eval_vars, eval_eqns, calc_equations, save_start_vars, &
         do1_eos_and_kap, calc_Fr, do1_specific_volume, get_Psurf, &
         set_f_Edd, calc_Prad, calc_Hp_face, calc_Y_face, &
         calc_PII_face, calc_Pvsc, calc_Pturb, calc_Chi, calc_Eq, &
         calc_source_sink, acceleration_eqn, calc_cell_equations, &
         T_form_of_calc_Fr, calc_Lc, calc_Lt, check_omega, rsp_set_Teff
      
      logical, parameter :: call_is_bad = .false.
      
      real(dp) :: x_dbg
      
      integer, parameter :: i_var_Vol = 99 ! for remeshing tests with dfridr
      
      integer, parameter :: &
         i_var_T = 2, i_var_w = 3, i_var_er = 4, i_var_Fr = 5, i_var_R = 6, & ! R must be last
      
         i_r_dr_in2 = 1, &
         i_r_dT_in = 2, &
         i_r_dw_in = 3, &
         i_r_der_in = 4, &
         i_r_dFr_in = 5, &

         i_r_dr_in = i_r_dr_in2 + NV, &
         i_r_dr_00 = i_r_dr_in + NV, &
         i_r_dr_out = i_r_dr_00 + NV, &
         i_r_dr_out2 = i_r_dr_out + NV, &
         
         i_r_dT_00 = i_r_dT_in + NV, &
         i_r_dT_out = i_r_dT_00 + NV, &
         i_r_dT_out2 = i_r_dT_out + NV, &
         
         i_r_dw_00 = i_r_dw_in + NV, &
         i_r_dw_out = i_r_dw_00 + NV, &
         i_r_dw_out2 = i_r_dw_out + NV, &
         
         i_r_der_00 = i_r_der_in + NV, &
         i_r_der_out = i_r_der_00 + NV, &
         i_r_der_out2 = i_r_der_out + NV, &
         
         i_r_dFr_00 = i_r_dFr_in + NV, &
         i_r_dFr_out = i_r_dFr_00 + NV, &
         i_r_dFr_out2 = i_r_dFr_out + NV, &
         
         i_T_dT_in2 = 1, &
         i_T_dw_in2 = 2, &
         i_T_der_in2 = 3, &
         i_T_dFr_in2 = 4, &
         i_T_dr_in2 = 5, &
         
         i_T_dT_in = i_T_dT_in2 + NV, &
         i_T_dT_00 = i_T_dT_in + NV, &
         i_T_dT_out = i_T_dT_00 + NV, &
         i_T_dT_out2 = i_T_dT_out + NV, &
         
         i_T_dw_in = i_T_dw_in2 + NV, &
         i_T_dw_00 = i_T_dw_in + NV, &
         i_T_dw_out = i_T_dw_00 + NV, &
         
         i_T_der_in = i_T_der_in2 + NV, &
         i_T_der_00 = i_T_der_in + NV, &
         i_T_der_out = i_T_der_00 + NV, &
         
         i_T_dFr_in = i_T_dFr_in2 + NV, &
         i_T_dFr_00 = i_T_dFr_in + NV, &
         i_T_dFr_out = i_T_dFr_00 + NV, &
         
         i_T_dr_in = i_T_dr_in2 + NV, &
         i_T_dr_00 = i_T_dr_in + NV, &
         i_T_dr_out = i_T_dr_00 + NV, &
         
         i_w_dw_in2 = 1, &
         i_w_der_in2 = 2, &
         i_w_dFr_in2 = 3, &
         i_w_dr_in2 = 4, &
         i_w_dT_in = 5, &
         
         i_w_dw_in = i_w_dw_in2 + NV, &
         i_w_dw_00 = i_w_dw_in + NV, &
         i_w_dw_out = i_w_dw_00 + NV, &
         i_w_dw_out2 = i_w_dw_out + NV, &
         
         i_w_der_in = i_w_der_in2 + NV, &
         i_w_der_00 = i_w_der_in + NV, &
         i_w_der_out = i_w_der_00 + NV, &
         
         i_w_dFr_in = i_w_dFr_in2 + NV, &
         i_w_dFr_00 = i_w_dFr_in + NV, &
         i_w_dFr_out = i_w_dFr_00 + NV, &
         
         i_w_dr_in = i_w_dr_in2 + NV, &
         i_w_dr_00 = i_w_dr_in + NV, &
         i_w_dr_out = i_w_dr_00 + NV, &
         
         i_w_dT_00 = i_w_dT_in + NV, &
         i_w_dT_out = i_w_dT_00 + NV, &
         i_w_dT_out2 = i_w_dT_out + NV, &
         
         i_er_der_in2 = 1, &
         i_er_dFr_in2 = 2, &
         i_er_dr_in2 = 3, &
         i_er_dT_in = 4, &
         i_er_dw_in = 5, &

         i_er_der_in = i_er_der_in2 + NV, &
         i_er_der_00 = i_er_der_in + NV, &
         i_er_der_out = i_er_der_00 + NV, &
         i_er_der_out2 = i_er_der_out + NV, &
         
         i_er_dFr_in = i_er_dFr_in2 + NV, &
         i_er_dFr_00 = i_er_dFr_in + NV, &
         i_er_dFr_out = i_er_dFr_00 + NV, &

         i_er_dr_in = i_er_dr_in2 + NV, &
         i_er_dr_00 = i_er_dr_in + NV, &
         i_er_dr_out = i_er_dr_00 + NV, &
         
         i_er_dT_00 = i_er_dT_in + NV, &
         i_er_dT_out = i_er_dT_00 + NV, &
         i_er_dT_out2 = i_er_dT_out + NV, &
         
         i_er_dw_00 = i_er_dw_in + NV, &
         i_er_dw_out = i_er_dw_00 + NV, &
         i_er_dw_out2 = i_er_dw_out + NV, &
         
         i_Fr_dFr_in2 = 1, &
         i_Fr_dr_in2 = 2, &
         i_Fr_dT_in = 3, &
         i_Fr_dw_in = 4, &
         i_Fr_der_in = 5, &

         i_Fr_dFr_in = i_Fr_dFr_in2 + NV, &
         i_Fr_dFr_00 = i_Fr_dFr_in + NV, &
         i_Fr_dFr_out = i_Fr_dFr_00 + NV, &
         i_Fr_dFr_out2 = i_Fr_dFr_out + NV, &

         i_Fr_dr_in = i_Fr_dr_in2 + NV, &
         i_Fr_dr_00 = i_Fr_dr_in + NV, &
         i_Fr_dr_out = i_Fr_dr_00 + NV, &
         
         i_Fr_dT_00 = i_Fr_dT_in + NV, &
         i_Fr_dT_out = i_Fr_dT_00 + NV, &
         i_Fr_dT_out2 = i_Fr_dT_out + NV, &
         
         i_Fr_dw_00 = i_Fr_dw_in + NV, &
         i_Fr_dw_out = i_Fr_dw_00 + NV, &
         i_Fr_dw_out2 = i_Fr_dw_out + NV, &
         
         i_Fr_der_00 = i_Fr_der_in + NV, &
         i_Fr_der_out = i_Fr_der_00 + NV, &
         i_Fr_der_out2 = i_Fr_der_out + NV
      
      integer :: iter, min_k_for_turbulent_flux
      
      
      contains


      subroutine HYD(s,ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: DXXT, DXXC, DXXE, DXXL, PREC2, DXH, EZH, rel_err_energy, &
            T_surf, P_surf, kap_surf, kap_guess, Teff_atm, &
            R_center_start, dt_max, total
         integer :: &
            i_min, i_max, num_tries, max_retries, max_iters, i, IT, k, nz, &
            kT_max, kW_max, kE_max, kL_max, iter_for_dfridr, test_partials_k
         integer(8) :: time0
         logical :: converged, dbg_msg, trace
         include 'formats'
         
         ierr = 0
         dbg_msg = s% report_solver_progress
         trace = s% trace_evolve
         iter = 0
         iter_for_dfridr = - 1
         test_partials_k = s% solver_test_partials_k
         s% solver_test_partials_k = 0
         if (s% model_number == s% solver_test_partials_call_number &
              .and. s% solver_test_partials_dx_0 > 0 &
              .and. s% solver_test_partials_iter_number > 0 &
              .and. test_partials_k > 0) then
            iter_for_dfridr = s% solver_test_partials_iter_number   
            s% solver_test_partials_var = 0
            s% solver_test_partials_val = 0
            s% solver_test_partials_dval_dx = 0
         end if   
         nz = s% nz
         i_min = 1
         i_max = nz
         PREC2 = s% RSP_tol_max_corr
         DXH = 0.3d0
         max_retries = s% RSP_max_retries_per_step
         max_iters = s% RSP_max_iters_per_try
         R_center_start = s% R_center
         call save_start_vars(s)
         do k=1,s% nz
            s% L(k) = s% Fr(k)*4*pi*s% r(k)**2 + s% Lc(k) + s% Lt(k)
            s% L_start(k) = s% L(k)
         end do
         call set_f_Edd(s,ierr)
         if (ierr /= 0) return
         P_surf = get_Psurf(s,ierr)
         if (ierr /= 0) return
         if (ALFAT > 0d0) then
            call set_min_k_for_turbulent_flux
         else
            min_k_for_turbulent_flux = 0
         end if
         
         if (s% v_center > s% v(nz) .and. s% v_center > 0d0) then ! compressing innermost cell
            dt_max = 1d-2*(s% r(nz) - s% r_center)/(s% v_center - s% v(nz))
            if (s% dt > dt_max) then
               !write(*,2) 'reduce dt in HYD', s% model_number, s% dt, dt_max
               !write(*,2) 's% r(nz) - s% r_center', s% model_number, s% r(nz) - s% r_center
               !write(*,2) 's% v_center - s% v(nz)', s% model_number, s% v_center - s% v(nz)
               s% dt = dt_max
               if (call_is_bad) then
                  if (is_bad(s% dt)) then
                     write(*,1) 'dt', s% dt
                     stop 'HYD compressing innermost cell'
                  end if
               end if
               if (s% RSP_report_limit_dt) &
                  write(*,4) 'limit dt to max_dt set by compressing innermost cell', s% model_number
            end if
         end if
         if (s% dt < s% force_timestep_min .and. s% force_timestep_min > 0) &
            s% dt = s% force_timestep_min
         if (s% force_timestep > 0) s% dt = s% force_timestep

         retry_loop: do num_tries = 1, max_retries+1

            converged = .false.
            call set_1st_iter_R_using_v_start(s)
            s% R_center = s% R_center + s% dt*s% v_center
            
            iter_loop: do iter = 1, max_iters  
               s% solver_iter = iter
               if (iter == iter_for_dfridr) then
                  s% solver_test_partials_k = test_partials_k
               end if
               call eval_vars(s,iter,i_min,i_max,ierr)
               if (ierr /= 0) return
               call eval_eqns(s,P_surf)
               if (iter == iter_for_dfridr) call check_partial()
               if (converged) then
                  s% num_solver_iterations = iter - 1
                  s% solver_test_partials_k = test_partials_k
                  return      
               end if
               if (s% doing_timing) call start_time(s, time0, total)
               call solve_for_corrections(s,iter)
               call apply_corrections(s, &
                  DXH,DXXT,DXXC,DXXE,DXXL,EZH,kT_max,kW_max,kE_max,kL_max)
               if (s% doing_timing) call update_time(s, time0, total, s% time_solver_matrix)
               if (dbg_msg) call write_msg
               if (iter == 1) cycle iter_loop
               converged = (abs(DXXT) < PREC2 .and. abs(DXXC) < PREC2)
            end do iter_loop      
            call doing_retry
         end do retry_loop
      
         write(*,*) ' NO CONVERGENCE IN HYD, TIME STEP: ',s% model_number
         stop
      
         contains
         
         subroutine set_min_k_for_turbulent_flux
            real(dp) :: tau, rmid
            integer :: k
            tau = 0
            min_k_for_turbulent_flux = nz
            do k = 1, nz-1
               rmid = 0.5d0*(s% r(k) + s% r(k+1))
               tau = tau + s% dm(k)*s% opacity(k)/(4*pi*rmid**2)
               if (tau >= s% RSP_min_tau_for_turbulent_flux) then
                  min_k_for_turbulent_flux = k
                  exit
               end if
            end do
         end subroutine set_min_k_for_turbulent_flux
               
         subroutine doing_retry
            integer :: i
            include 'formats'
            if (num_tries  ==  max_retries+1) then
               write(*,*) 'NO CONVERGENCE IN HYD, TIME STEP num tries, max allowed', &
                  s% model_number, num_tries, max_retries+1
               stop 'RSP: step num_retries = RSP_max_retries_per_step'
            end if  
            call restore_start_vars(s)
            s% R_center = R_center_start
            s% dt = s% dt/2.d0   
            s% num_retries = s% num_retries + 1
            if (s% max_number_retries < 0) return
            if (s% num_retries > s% max_number_retries) then
               write(*,3) 'model max_number_retries', s% model_number, s% max_number_retries
               stop 'RSP: num_retries > max_number_retries'
            end if
         end subroutine doing_retry
         
         subroutine write_msg
            integer :: i, k
            include 'formats'
            !if (EZH < 1d0) write(*,3) 'undercorrection factor', s% model_number, iter, EZH
            write(*,'(i6, 2x, i3, 4(4x, a, 1x, i4, 1x, 1pe11.4, 1x, 1pe11.4))') &
               s% model_number, iter, &
               'T', kT_max, DXXT, s% T(max(1,kT_max)), &
               'w', kW_max, DXXC, max(1d-99,s% RSP_w(max(1,kW_max))), &
               'erad', kE_max, DXXE, s% erad(max(1,kE_max)), &
               'Fr', kL_max, DXXL, s% Fr(max(1,kL_max))
         end subroutine write_msg
      
         subroutine check_partial
            integer :: i_var
            real(dp) :: dvardx_0, dx_0, dvardx, xdum, err
            include 'formats'
            i_var = s% solver_test_partials_var 
            dvardx_0 = s% solver_test_partials_dval_dx ! analytic partial
            if (i_var <= 0) then
               write(*,2) 'need to set test_partials_var', i_var
               stop 'check_partial'
            end if            
            dx_0 = get1_val(i_var, s% solver_test_partials_k)
            dx_0 = s% solver_test_partials_dx_0*max(1d-99, abs(dx_0))
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
            write(*,1) 'analytic numeric err rel_diff',dvardx_0,dvardx,err,xdum
            !write(*,*)
            stop 'check_partial'            
         end subroutine check_partial         
         
         real(dp) function get1_val(i_var,k) result(val)
            integer, intent(in) :: i_var, k
            include 'formats'
            if (i_var == i_var_w) then
               val = s% RSP_w(k)
            else if (i_var == i_var_R) then
               val = s% r(k)
            else if (i_var == i_var_T) then
               val = s% T(k)
            else if (i_var == i_var_er) then
               val = s% erad(k)
            else if (i_var == i_var_Fr) then
               val = s% Fr(k)
            else if (i_var == i_var_Vol) then
               val = s% Vol(k)
            else 
               write(*,2) 'bad value for solver_test_partials_var', i_var
               stop 'solver_test_partials'
            end if
         end function get1_val
         
         subroutine store1_val(i_var, k, val) 
            integer, intent(in) :: i_var, k
            real(dp), intent(in) :: val
            include 'formats'
            if (i_var == i_var_w) then
               s% RSP_w(k) = val
            else if (i_var == i_var_R) then
               s% r(k) = val
               s% v(k) = 2.d0*(s% r(k) - s% r_start(k))/s% dt - s% v_start(k)
               ! partials wrt R assume v changes along with R,
               ! so need to do it here.
            else if (i_var == i_var_T) then
               s% T(k) = val
            else if (i_var == i_var_er) then
               s% erad(k) = val
            else if (i_var == i_var_Fr) then
               s% Fr(k) = val
            else if (i_var == i_var_Vol) then
               s% Vol(k) = val
            else 
               write(*,2) 'bad value for solver_test_partials_var', i_var
               stop 'solver_test_partials'
            end if
         end subroutine store1_val

         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: i_var, k, ierr
            real(dp) :: save1
            include 'formats'
            i_var = s% solver_test_partials_var
            k = s% solver_test_partials_k
            save1 = get1_val(i_var, k)
            call store1_val(i_var, k, save1 + delta_x)
            call eval_vars(s,0,i_min,i_max,ierr)
            if (ierr /= 0) stop 'failed in eval_vars'
            call eval_eqns(s,P_surf)
            val = s% solver_test_partials_val
            !write(*,2) 'dfridr val', k, val
            call store1_val(i_var, k, save1)
         end function dfridr_func

         real(dp) function dfridr(hx,err)
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: errt,fac,hh,a(ntab,ntab),f1,f2
            real(dp), parameter :: con2 = 2d0, con = sqrt(con2), big = 1d199, safe = 2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            f1 = dfridr_func(hh)
            !write(*,2) 'f1', 1, f1, save_dx(i_var,k) + hh
            f2 = dfridr_func( - hh)
            !write(*,2) 'f2', 1, f2, save_dx(i_var,k) - hh
            a(1,1) = (f1 - f2)/(2d0*hh)
            !write(*,2) 'dfdx', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i = 2,ntab
               hh = hh/con
               f1 = dfridr_func(hh)
               !write(*,2) 'f1', i, f1
               f2 = dfridr_func( - hh)
               !write(*,2) 'f2', i, f2
               a(1,i) = (f1 - f2)/(2d0*hh)
               !write(*,2) 'dfdx', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j = 2,i
                  a(j,i) = (a(j - 1,i)*fac - a(j - 1,i - 1))/(fac - 1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i) - a(j - 1,i)),abs(a(j,i) - a(j - 1,i - 1)))
                  !write(*,2) 'errt, err', j, errt, err
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i - 1,i - 1)) >= safe*err) then
                  !write(*,1) 'higher order is worse'
                  write(*,*)
                  return
               end if
            end do
         end function dfridr
      
      end subroutine HYD
      
      
      real(dp) function get_Psurf(s,ierr) result(P_surf)
         use rsp_eval_eos_and_kap, only: get_surf_P_T_kap
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: kap_guess, T_surf, kap_surf, Teff_atm
         include 'formats'
         ierr = 0
         if (s% RSP_fixed_Psurf) then
            P_surf = Psurf_from_atm
         else if (s% RSP_use_atm_grey_with_kap_for_Psurf) then
            ierr = 0
            kap_guess = 1d-2
            call get_surf_P_T_kap(s, &
               s% M(1), s% r(1), s% L(1), &
               (2d0/3d0)*s% tau_factor, kap_guess, &
               T_surf, P_surf, kap_surf, Teff_atm, ierr)
            if (ierr/= 0) then
               write(*,*) 'ierr from get_surf_P_T_kap'
               return
            end if
         else if (s% RSP_use_Prad_for_Psurf) then
            P_surf = crad*s% T(1)**4/3d0
         else if (s% RSP_Psurf >= 0d0) then
            P_surf = s% RSP_Psurf
         else
            P_surf = 0d0
         end if
      end function get_Psurf


      subroutine calculate_energies(s,total_radiation)
         use star_utils, only: cell_specific_KE, cell_specific_PE
         type (star_info), pointer :: s
         real(dp), intent(out) :: total_radiation
         integer :: i, k
         real(dp) :: d_dlnR00, d_dlnRp1, d_dv00, d_dvp1
         include 'formats'
         EGRV = 0.d0
         ETHE = 0.d0
         EKIN = 0.d0
         ECON = 0.d0
         do i=1,NZN
            k = NZN+1-i
            ETHE = ETHE + (s% egas(k)+s% erad(k))*s% dm(k)
            ECON = ECON + s% RSP_w(k)**2*s% dm(k)
            EKIN = EKIN + cell_specific_KE(s,k,d_dv00,d_dvp1)*s% dm(k)
            EGRV = EGRV + cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)*s% dm(k)           
            
            !EKIN = EKIN + 0.5d0*s% v(k)**2*s% dm_bar(k)               
            !if (k < NZN) then
            !   EGRV = EGRV - s% cgrav(k) * (s%m(k)-0.5d0*s%dm(k))*s%dm(k)/(0.5d0*(s%r(k)+s%r(k+1)))
            !else
            !   EGRV = EGRV - s% cgrav(k) * (s%m(k)-0.5d0*s%dm(k))*s%dm(k)/(0.5d0*(s%r(k)+s%r_center))
            !end if
            
         enddo
         if (s% RSP_hydro_only) then
            total_radiation = 0d0
         else
            total_radiation = s% dt*(s% L_center - (WTR*s% L(1) + WTR1*s% L_start(1)))
         end if
         ETOT = ETHE+EKIN+ECON+EGRV
         EDE_start = EDE_start-E0
         if (EKIN > EKMAX) EKMAX = EKIN
         if (EKIN < EKMIN) EKMIN = EKIN
         !write(*,1) 'ETHE', ETHE
         !write(*,1) 'ECON', ECON
         !write(*,1) 'EKIN', EKIN
         !write(*,1) 'EGRV', EGRV
         !write(*,1) 'ETOT', ETOT
      end subroutine calculate_energies


      subroutine init_HYD(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: i_min, i_max
         include 'formats'
         i_min = 1
         i_max = s% nz
         ! setup so can save start vals in 1st call to HYD
         s% f_Edd(1:NZN) = f_Edd_isotropic ! fake it for 1st call on eval_vars
         call eval_vars(s,0,i_min,i_max,ierr)
         if (ierr /= 0) return
         call set_f_Edd(s,ierr) ! needs opacities
         if (ierr /= 0) return
         call save_start_vars(s) ! needed by eval_eqns
         call eval_eqns(s,0d0)
      end subroutine init_HYD


      subroutine set_1st_iter_R_using_v_start(s)
         type (star_info), pointer :: s
         real(dp) :: EH1,EHJT
         integer :: k
         include 'formats'         
         EHJT = 1.d0
         do k = 1,NZN
            if (k /= NZN) then
               EH1 = - (s% v_start(k) - s% v_start(k+1))*s% dt/ &
                  (s% r_start(k) - s% r_start(k+1))/0.8d0
            else
               EH1 = - (s% v_start(k) - s% v_center)*s% dt/ &
                  (s% r_start(k) - s% R_center)/0.8d0
            end if
            EHJT = 1.d0/max(1.d0/EHJT,EH1,1d0)
         end do
         do k = 1,NZN
            s% r(k) = s% r_start(k) + EHJT*s% dt*s% v_start(k) ! initial guess for R
            s% v(k) = (2d0*EHJT - 1d0)*s% v_start(k) ! initial guess for v
         end do
      end subroutine set_1st_iter_R_using_v_start


      subroutine rsp_set_Teff(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: tau00, taup1, dtau, tau_phot, Tface_0, Tface_1
         include 'formats'
         tau_phot = 2d0/3d0
         tau00 = s% tau_factor*s% tau_base
         s% Teff = s% T(1)
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(2*pi*(s% r(k)**2 + s% r(k+1)**2))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               if (k == 1) then
                  Tface_0 = s% T(k)
               else
                  Tface_0 = 0.5d0*(s% T(k) + s% T(k-1))
               end if
               Tface_1 = 0.5d0*(s% T(k) + s% T(k+1))
               s% Teff = Tface_0 + (Tface_1 - Tface_0)*(tau_phot - tau00)/dtau
               exit
            end if
            tau00 = taup1
         end do
      end subroutine rsp_set_Teff
      
      
      subroutine set_f_Edd(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: lim_f_Edd, lim_g_Edd
         integer :: k, j
         include 'formats'
         ierr = 0
         s% g_Edd = 0.5d0
         s% f_Edd(1:NZN) = f_Edd_isotropic
      end subroutine set_f_Edd
      
      
      subroutine save_start_vars(s)
         type (star_info), pointer :: s
         integer :: I, k
         do I = 1,NZN
            k = NZN+1-i
            s% T_start(k) = s% T(k)
            s% r_start(k) = s% r(k)
            s% Pgas_start(k) = s% Pgas(k)
            s% Prad_start(k) = s% Prad(k)
            s% Pvsc_start(k) = s% Pvsc(k)
            s% Vol_start(k) = s% Vol(k)
            s% csound_start(k) = s% csound(k)
            s% opacity_start(k) = s% opacity(k)
            s% egas_start(k) = s% egas(k)
            s% erad_start(k) = s% erad(k)
            s% RSP_w_start(k) = s% RSP_w(k)
            s% Ptrb_start(k) = s% Ptrb(k)
            s% Chi_start(k) = s% Chi(k)
            s% v_start(k) = s% v(k)
            s% Fr_start(k) = s% Fr(k)
            s% Lc_start(k) = s% Lc(k)
            s% Lt_start(k) = s% Lt(k)
            s% COUPL_start(k) = s% COUPL(k)
         end do
      end subroutine save_start_vars
      
      
      subroutine restore_start_vars(s)
         type (star_info), pointer :: s
         integer :: I, k
         do I = 1,NZN
            k = NZN+1-i
            s% T(k) = s% T_start(k)
            s% r(k) = s% r_start(k)
            s% Pgas(k) = s% Pgas_start(k)
            s% Prad(k) = s% Prad_start(k)
            s% Pvsc(k) = s% Pvsc_start(k)
            s% Vol(k) = s% Vol_start(k)
            s% csound(k) = s% csound_start(k)
            s% opacity(k) = s% opacity_start(k)
            s% egas(k) = s% egas_start(k)
            s% erad(k) = s% erad_start(k)
            s% RSP_w(k) = s% RSP_w_start(k)
            s% Ptrb(k) = s% Ptrb_start(k)
            s% Chi(k) = s% Chi_start(k)
            s% v(k) = s% v_start(k)
            s% L(k) = s% L_start(k)
            s% Fr(k) = s% Fr_start(k)
            s% Lc(k) = s% Lc_start(k)
            s% Lt(k) = s% Lt_start(k)
            s% COUPL(k) = s% COUPL_start(k)
         end do            
      end subroutine restore_start_vars
      
      
      subroutine eval_vars(s,iter,i_min,i_max,ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         integer, intent(in) :: iter,i_min,i_max
         integer, intent(out) :: ierr
         integer :: i, k, op_err
         integer(8) :: time0
         real(dp) :: total
         include 'formats'
         ierr = 0
         if (s% doing_timing) call start_time(s, time0, total)
         !$OMP PARALLEL DO PRIVATE(I,op_err) SCHEDULE(dynamic,2)
         do i = i_min,i_max
            call do1_specific_volume(s,i)
            call do1_eos_and_kap(s,i,op_err)
            if (op_err /= 0) ierr = op_err
            call calc_Prad(s,i)
         end do
         !$OMP END PARALLEL DO
         if (s% doing_timing) call update_time(s, time0, total, s% time_eos)
         if (ierr /= 0) return
         !$OMP PARALLEL DO PRIVATE(I) SCHEDULE(dynamic,2)
         do i = i_min,i_max
            call calc_Hp_face(s,i)
            call calc_Y_face(s,i)
            call calc_PII_face(s,i)
            call calc_Pvsc(s,i)
         end do
         !$OMP END PARALLEL DO
         if (iter == 1) then
            do i = i_min,i_max
               call check_omega(s,i)
            end do
         end if
         !$OMP PARALLEL DO PRIVATE(I) SCHEDULE(dynamic,2)
         do i = i_min,i_max
            call calc_Pturb(s,i)
            call calc_Chi(s,i)
            call calc_Eq(s,i)
            call calc_source_sink(s,i)
         end do
         !$OMP END PARALLEL DO
         call zero_boundaries(s)               
      end subroutine eval_vars
      
      
      subroutine eval_eqns(s,P_surf)
         type (star_info), pointer :: s
         real(dp), intent(in) :: P_surf
         integer :: i, k
         include 'formats'
         !$OMP PARALLEL DO PRIVATE(I) SCHEDULE(dynamic,2)
         do I = 1,NZN
            call calc_equations(s, I, P_surf)
         end do
         !$OMP END PARALLEL DO
      end subroutine eval_eqns


      subroutine do1_specific_volume(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: k
         real(dp) :: T1
         real(qp) :: q1, q2, q3, q4
         include 'formats'
         k = NZN+1-i
         T1 = P43/s% dm(k)
         if (I /= 1) then
            q1 = T1
            q2 = s% r(k)
            q3 = s% r(k+1)
            q4 = q1*(q2**3 - q3**3)
            s% Vol(k) = dble(q4)
            dVol_dr_in(I) = - 3.0d0*T1*s% r(k+1)**2
         else
            s% Vol(k) = T1*(s% r(k)**3 - s% R_center**3)
            dVol_dr_in(I) = 0d0
         end if
         if (s% Vol(k) <= 0d0) then
            write(*,2) 'bad Vol', k, s% Vol(k)
            stop 'do1_specific_volume'
         end if
         dVol_dr_00(I) = 3.d0*T1*s% r(k)**2
      end subroutine do1_specific_volume
      
      
      subroutine solve_for_corrections(s,iter)
         type (star_info), pointer :: s
         integer, intent(in) :: iter
         integer :: i, j, info, IR, IT, IW, IE, IL, N
         logical :: okay
         include 'formats'
         
         !if (s% model_number >= s% max_model_number-1 .and. iter == 1) then
         !if (.false. .and. s% model_number == s% max_model_number) then !.and. iter == 1) then
         if (.false.) then
            ! this can be useful for finding restart bugs.
            write(*,*) 'solve_for_corrections'
            okay = .true.
            do i = 1,NZN
               IR = i_var_R + NV*(i-1)
               IT = i_var_T + NV*(i-1)
               IW = i_var_w + NV*(i-1)
               IE = i_var_er + NV*(i-1)
               IL = i_var_Fr + NV*(i-1)
               
               if (is_bad(HR(IR))) then
                  write(*,4) 'HR(IR)', iter, i, IR, HR(IR)
                  okay = .false.
               end if
               if (is_bad(HR(IT)))  then
                  write(*,4) 'HR(IT)', iter, i, IT, HR(IT)
                  okay = .false.
               end if
               if (is_bad(HR(IW)))  then
                  write(*,4) 'HR(IW)', iter, i, IW, HR(IW)
                  okay = .false.
               end if
               if (is_bad(HR(IL)))  then
                  write(*,4) 'HR(IL)', iter, i, IL, HR(IL)
                  okay = .false.
               end if
               
               if (.not. okay) stop 'solve_for_corrections'
               
               do j = 1,LD_HD
                  if (is_bad(HD(j,IR))) then
                     write(*,5) 'HD(j,IR)', iter, i, j, IR, HD(j,IR)
                     okay = .false.
                  end if
                  if (is_bad(HD(j,IT))) then
                     write(*,5) 'HD(j,IT)', iter, i, j, IT, HD(j,IT)
                     okay = .false.
                  end if
                  if (is_bad(HD(j,IW))) then
                     write(*,5) 'HD(j,IW)', iter, i, j, IW, HD(j,IW)
                     okay = .false.
                  end if
                  if (is_bad(HD(j,IE))) then
                     write(*,5) 'HD(j,IE)', iter, i, j, IE, HD(j,IE)
                     okay = .false.
                  end if
                  if (is_bad(HD(j,IL))) then
                     write(*,5) 'HD(j,IL)', iter, i, j, IL, HD(j,IL)
                     okay = .false.
                  end if
               
                  if (.not. okay) stop 'solve_for_corrections'
                  
               end do

               cycle

               write(*,3) 'HR(IR)', iter, i, HR(IR)
               write(*,3) 'HR(IT)', iter, i, HR(IT)
               write(*,3) 'HR(IW)', iter, i, HR(IW)
               write(*,3) 'HR(IL)', iter, i, HR(IL)
               
               do j = 1,13
                  write(*,5) 'HD(j,IR)', iter, j, IR, i, HD(j,IR)
               end do
               
               do j = 1,13
                  write(*,5) 'HD(j,IT)', iter, j, IT, i, HD(j,IT)
               end do
               
               do j = 1,13
                  write(*,5) 'HD(j,IW)', iter, j, IW, i, HD(j,IW)
               end do
               
               write(*,3) 'HR(IE)', iter, i, HR(IE)
               do j = 1,13
                  write(*,5) 'HD(j,IE)', iter, j, IE, i, HD(j,IE)
               end do
               
               do j = 1,13
                  write(*,5) 'HD(j,IL)', iter, j, IL, i, HD(j,IL)
               end do
               
            end do
            !stop 'solve_for_corrections'
         end if
         
         N = NV*NZN+1
         
         if (.false.) then ! check HR and HD for NaN's
            do I = 1,N
               if (is_bad(HR(I))) then
                  write(*,3) 'HR(I)', iter, I, HR(I)
                  stop 'solve_for_corrections'
               end if
            end do
            do I = 1,N
               do j = 1,LD_HD
                  if (is_bad(HD(j,i))) then
                     write(*,4) 'HD(j,i)', iter, j, i, HD(j,i)
                     stop 'solve_for_corrections'
                  end if
               end do
            end do
         end if
     
         do J = 1,2*NV     ! translate hd into band storage of LAPACK
            do I = 1,N
               ABB(J,I) = 0.0d0
            end do
         end do      
         do J = 1,2*NV      
            do I = 1,N - J
               ABB(LD_HD - J,I + J) = HD(HD_DIAG + J,I) ! upper diagonals
            end do  
         end do
         do J = 1,2*NV     
            do I = 1,N - J
               ABB(LD_HD + J,I) = HD(HD_DIAG - J,I + J) ! lower diagonals
            end do  
         end do
         do I = 1,N
            ABB(LD_HD,I) = HD(HD_DIAG,I)
         end do  
      
         call DGBTRF(N,N,2*NV,2*NV,ABB,LD_ABB,IPVT,INFO)
         if (INFO/= 0) then
            write(*,*) 'hyd: LAPACK/dgbtrf problem',INFO
            stop
         end if
      
         call DGBTRS('n',N,2*NV,2*NV,1,ABB,LD_ABB,IPVT,HR,N,INFO)
         if (INFO/= 0) then
            write(*,*) 'hyd: LAPACK/dgbtrs problem',INFO
            stop
         end if             

         do I = 1,N
            DX(I) = HR(I)
            if (call_is_bad) then
               if (is_bad(DX(I))) then
                  write(*,2) 'DX(I)', I, DX(I)
                  stop 'solve_for_corrections'
               end if
            end if
         end do

      end subroutine solve_for_corrections

      
      subroutine apply_corrections(s, &
            DXH, XXT, XXC, XXE, XXL, EZH, &
            kT_max, kW_max, kE_max, kL_max)
         use star_utils, only: rand
         type (star_info), pointer :: s
         real(dp), intent(in) :: DXH
         real(dp), intent(out) :: XXT, XXC, XXE, XXL, EZH
         integer, intent(out) :: kT_max,kW_max,kE_max,kL_max
         integer :: i, k, IR, IT, IW, IU, IE, IL, kEZH, &
            iTM, kTM, iRM, kRM, iEM, kEM, iCM, kCM, iLM, kLM
         real(dp) :: XXR, DXXT, DXXC, DXXE, DXXL, DXRM, DXR, &
            EZH1, XXTM, XXCM, XXRM, XXEM, XXLM, DXKT, DXKC, DXKE, DXKL
         include 'formats'
         EZH = 1.0d0; kEZH = 0
         XXTM = 0d0; kTM = 0; iTM = 0
         XXRM = 0d0; kRM = 0; iRM = 0
         XXEM = 0d0; kEM = 0; iEM = 0
         XXCM = 0d0; kCM = 0; iCM = 0
         XXLM = 0d0; kLM = 0; iLM = 0
         do I = 1,NZN
            k = NZN+1-i
            IR = i_var_R + NV*(i-1)
            IT = i_var_T + NV*(i-1)
            IW = i_var_w + NV*(i-1)
            IE = i_var_er + NV*(i-1)
            IL = i_var_Fr + NV*(i-1)
            if (s% RSP_w(k) > (1.d+2)*EFL0) then
               XXC = abs(DX(IW)/s% RSP_w(k))/DXH
               if (XXC > XXCM) then
                  XXCM = XXC; kCM = k; iCM = IW
               end if
            else
               XXC = 0d0
            end if      
            if (i == 1) then
               DXR = -DX(IR)
               XXR = (DXR/(s% r(k) - s% r_center))/DXH
               if (XXR > XXRM) then
                  DXRM = DXR/(s% r(k) - s% r_center)
                  XXRM = XXR; kRM = k; iRM = IR
               end if      
            else      
               DXR = DX(IR-NV) - DX(IR)
               XXR = (DXR/(s% r(k) - s% r(k+1)))/DXH
               if (XXR > XXRM) then
                  DXRM = DXR/(s% r(k) - s% r(k+1))
                  XXRM = XXR; kRM = k; iRM = IR
               end if      
            end if      
            XXE = abs(DX(IE)/s% erad(k))/DXH
            if (XXE > XXEM) then
               XXEM = XXE; kEM = k; iEM = IE
            end if            
            XXL = abs(DX(IL)/s% Fr(k))/DXH
            if (XXL > XXLM) then
               XXLM = XXL; kLM = k; iLM = IL
            end if            
            XXT = abs(DX(IT)/s% T(k))/DXH
            if (XXT > XXTM) then
               XXTM = XXT; kTM = k; iTM = IT
            end if    
            EZH1 = EZH       
            EZH = 1.d0/max(1.d0/EZH,XXR,XXT,XXC,XXE)
            if (EZH1 /= EZH) kEZH = k
         end do
         
         if (EZH < 1d0 .and. s% RSP_report_undercorrections) then
            write(*,'(i6, 2x, i3, 6(4x, a, 1x, i4, 1x, 1pe10.3))') &
               s% model_number, iter, &
               'EZH', kEZH, EZH, &
               'r', kRM, DXRM, &
               'T', kTM, DX(iTM)/s% T(kTM), &
               'w', kCM, DX(iCM)/max(1d-99,s% RSP_w(kCM)), &
               'erad', kEM, DX(iEM)/s% erad(kEM), &
               'Fr', kLM, DX(iLM)/s% Fr(kLM)
         end if

         DXXT = 0.0d0
         DXXC = 0.0d0
         DXXE = 0.0d0
         DXXL = 0.0d0
         XXT = 0.0d0
         XXC = 0.0d0
         XXE = 0.0d0
         XXL = 0.0d0
         kT_max = 0
         kW_max = 0
         kE_max = 0
         kL_max = 0

         do I = 1,NZN
            k = NZN+1-i
            IT = i_var_T + NV*(i-1)
            IR = i_var_R + NV*(i-1)
            IW = i_var_w + NV*(i-1)
            IE = i_var_er + NV*(i-1)
            IL = i_var_Fr + NV*(i-1)
            s% T(k) = s% T(k) + EZH*DX(IT)
            s% erad(k) = s% erad(k) + EZH*DX(IE)
            if (I > IBOTOM .and. I < NZN)then
               if ((s% RSP_w(k) + EZH*DX(IW)) <= 0d0)then
                  s% RSP_w(k) = EFL0*rand(s)*1d-6 ! RSP NEEDS THIS to give seed for SOURCE
               else
                  s% RSP_w(k) = s% RSP_w(k) + EZH*DX(IW)
               end if
            end if
            s% v(k) = s% v(k) - &
               (s% v(k) + s% v_start(k)) + &
                2.d0/s% dt*(EZH*DX(IR) + (s% r(k) - s% r_start(k)))
            s% r(k) = s% r(k) + EZH*DX(IR)
            s% Fr(k) = s% Fr(k) + EZH*DX(IL)
            s% Lr(k) = s% Fr(k)*4d0*pi*s% r(k)**2
            DXKT = DXXT
            DXKC = DXXC
            DXKE = DXXE
            DXKL = DXXL
            DXXT = max(DXXT,abs(DX(IT)/s% T(k)))
            DXXE = max(DXXE,abs(DX(IE)/s% erad(k)))
            DXXL = max(DXXL,abs(DX(IL)/s% Fr(k)))
            if (s% RSP_w(k) > (1.d-2)*EFL0) &
               DXXC = max(DXXC,abs(DX(IW)/s% RSP_w(k)))
            if (DXXC > DXKC) then
               kW_max = k; XXC = DX(IW)/max(1d-99,s% RSP_w(k))
            end if
            if (DXXT > DXKT) then
               kT_max = k; XXT = DX(IT)/s% T(k)
            end if
            if (DXXE > DXKE) then
               kE_max = k; XXE = DX(IE)/s% erad(k)
            end if
            if (DXXL > DXKL) then
               kL_max = k; XXL = DX(IL)/s% Fr(k)
            end if
         end do
      
      end subroutine apply_corrections


      subroutine do1_eos_and_kap(s,i,ierr)
         use rsp_eval_eos_and_kap, only : eval_mesa_eos_and_kap
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer, intent(out) :: ierr
         real(dp) :: Prad, d_Pr_dT, erad, d_erad_dVol, d_erad_dT
         integer :: k
         include 'formats'
         k = NZN + 1 - i
         call eval_mesa_eos_and_kap(&
            s, k, s% T(k), s% Vol(k), &
            s% Pgas(k), d_Pg_dVol(I), d_Pg_dT(I), &
            Prad, d_Pr_dT, &
            s% egas(k), d_egas_dVol(I), d_egas_dT(I), &
            erad, d_erad_dVol, d_erad_dT, &
            s% csound(k), s% Cp(k), dCp_dVol(I), dCp_dT(I), &
            s% QQ(k), dQQ_dVol(I), dQQ_dT(I), &
            s% opacity(k), dK_dVol(I), dK_dT(I),ierr)   
         if (ierr /= 0) return         
         d_Pg_dr_00(I) = d_Pg_dVol(I)*dVol_dr_00(I)
         d_Pg_dr_in(I) = d_Pg_dVol(I)*dVol_dr_in(I)
         d_egas_dr_00(I) = d_egas_dVol(I)*dVol_dr_00(I)
         d_egas_dr_in(I) = d_egas_dVol(I)*dVol_dr_in(I)
         dCp_dr_00(I) = dCp_dVol(I)*dVol_dr_00(I)
         dCp_dr_in(I) = dCp_dVol(I)*dVol_dr_in(I)
         dQQ_dr_00(I) = dQQ_dVol(I)*dVol_dr_00(I)
         dQQ_dr_in(I) = dQQ_dVol(I)*dVol_dr_in(I)
         dK_dr_00(I) = dK_dVol(I)*dVol_dr_00(I)
         dK_dr_in(I) = dK_dVol(I)*dVol_dr_in(I)
         if (call_is_bad) then
            if (is_bad(s% opacity(k))) then
!$OMP critical
               write(*,2) 's% opacity(k)', k, s% opacity(k)
               stop 'do1_eos_and_kap'
!$OMP end critical
            end if
         end if
      end subroutine do1_eos_and_kap
      
      
      subroutine calc_Prad(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: k
         real(dp) :: V
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         V = s% Vol(k)
         s% Prad(k) = s% f_Edd(k)*s% erad(k)/V
         d_Pr_der(i) = s% f_Edd(k)/V
         d_Pr_dVol(i) = -s% Prad(k)/V
         d_Pr_dr_00(i) = d_Pr_dVol(i)*dVol_dr_00(i)
         d_Pr_dr_in(i) = d_Pr_dVol(i)*dVol_dr_in(i)

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (test_partials) then
            s% solver_test_partials_val = s% Prad(k)
            s% solver_test_partials_var = i_var_er
            s% solver_test_partials_dval_dx = d_Pr_der(i)
            write(*,*) 'calc_Prad', s% solver_test_partials_var
            write(*,2) 'erad Prad f_Edd', k, s% erad(k), s% Prad(k), s% f_Edd(k)
         end if
         
      end subroutine calc_Prad


      subroutine calc_Hp_face(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: POM
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (I < NZN) then
         
            POM = (s% r(k)**2)/(2.d0*s% cgrav(k)*s% m(k))
            s% Hp_face(k) = POM*( &
               (s% Pgas(k) + s% Prad(k))*s% Vol(k) &
             + (s% Pgas(k-1) + s% Prad(k-1))*s% Vol(k-1))
            
            dHp_dVol_00(I) = POM*( &
               + s% Vol(k)*(d_Pg_dVol(i) + d_Pr_dVol(i)) &
               + s% Pgas(k) + s% Prad(k))
            dHp_dVol_out(I) = POM*( &
               + s% Vol(k-1)*(d_Pg_dVol(i+1) + d_Pr_dVol(i+1)) &
               + s% Pgas(k-1) + s% Prad(k-1))
               
            dHp_dr_in(I) = POM*( &
               (s% Pgas(k) + s% Prad(k))*dVol_dr_in(I) &
             + (d_Pg_dr_in(I) + d_Pr_dr_in(I))*s% Vol(k)) !             
            dHp_dr_00(I) = 2.d0*s% Hp_face(k)/s% r(k) + POM*( &
                (s% Pgas(k) + s% Prad(k))*dVol_dr_00(I) &
              + (d_Pg_dr_00(I) + d_Pr_dr_00(I))*s% Vol(k) &
              + (s% Pgas(k-1) + s% Prad(k-1))*dVol_dr_in(i+1) &
              + (d_Pg_dr_in(i+1) + d_Pr_dr_in(i+1))*s% Vol(k-1)) !                   
            dHp_dr_out(I) = POM*( &
               (s% Pgas(k-1) + s% Prad(k-1))*dVol_dr_00(i+1) &
             + (d_Pg_dr_00(i+1) + d_Pr_dr_00(i+1))*s% Vol(k-1)) ! 
             
            dHp_dT_00(I) = POM*s% Vol(k)*d_Pg_dT(I) ! 
            dHp_dT_out(I) = POM*s% Vol(k-1)*d_Pg_dT(I+1) ! 
              
            dHp_der_00(I) = POM*s% Vol(k)*d_Pr_der(I) ! 
            dHp_der_out(I) = POM*s% Vol(k-1)*d_Pr_der(I+1) ! 

         else ! surface
         
            POM = (s% r(k)**2)/(s% cgrav(k)*s% M(k))
            s% Hp_face(k) = POM*(s% Pgas(k) + s% Prad(k))*s% Vol(k)
            
            dHp_dVol_00(I) = POM*( &
               + s% Vol(k)*(d_Pg_dVol(i) + d_Pr_dVol(i)) &
               + s% Pgas(k) + s% Prad(k))
            dHp_dVol_out(I) = 0d0
            
            dHp_dr_in(i) = POM*( &
               (s% Pgas(k) + s% Prad(k))*dVol_dr_in(i) &
             + (d_Pg_dr_in(i) + d_Pr_dr_in(i))*s% Vol(k))  
            dHp_dr_00(i) = 2.d0*s% Hp_face(k)/s% r(k) + POM*( &
               (s% Pgas(k) + s% Prad(k))*dVol_dr_00(i) &
             + (d_Pg_dr_00(i) + d_Pr_dr_00(i))*s% Vol(k))
            dHp_dr_out(i) = 0.d0
            dHp_dT_00(i) = POM*s% Vol(k)*d_Pg_dT(i)
            dHp_dT_out(i) = 0.d0
            dHp_der_00(i) = POM*s% Vol(k)*d_Pr_der(i)
            dHp_der_out(i) = 0.d0
            
         end if

         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% Hp_face(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dHp_dVol_out(I)
            write(*,*) 'calc_Hp_face', s% solver_test_partials_var
         end if

      end subroutine calc_Hp_face
      
      
      subroutine calc_Y_face(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: POM, POM2, &
            Y1, d_Y1_dr_00, d_Y1_dr_in, d_Y1_dr_out, &
            d_Y1_dVol_00, d_Y1_dVol_out, &
            d_Y1_dT_00, d_Y1_dT_out, &
            d_Y1_der_00, d_Y1_der_out, &
            Y2, d_Y2_dr_00, d_Y2_dr_in, d_Y2_dr_out, &
            d_Y2_dVol_00, d_Y2_dVol_out, &
            d_Y2_dT_00, d_Y2_dT_out, &
            d_Y2_der_00, d_Y2_der_out
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (i < NZN .and. ALFA /= 0d0) then
            POM = 0.5d0*(s% QQ(k)/s% Cp(k) + s% QQ(k-1)/s% Cp(k-1))
            POM2 = 0.5d0*((s% Pgas(k-1) + s% Prad(k-1)) - (s% Pgas(k) + s% Prad(k)))
            Y1 = &
               0.5d0*(s% QQ(k)/s% Cp(k)+ s% QQ(k-1)/s% Cp(k-1))* &
                     ((s% Pgas(k-1) + s% Prad(k-1)) - (s% Pgas(k) + s% Prad(k))) &
               - (s% lnT(k-1) - s% lnT(k))
            
            d_Y1_dVol_00 = &
               - POM*(d_Pg_dVol(i) + d_Pr_dVol(i)) &
               + POM2*(dQQ_dVol(i) - s% QQ(k)*dCp_dVol(i)/s% Cp(k))/s% Cp(k)
            d_Y1_dVol_out = &
               + POM*(d_Pg_dVol(i+1) + d_Pr_dVol(i+1)) &
               + POM2*(dQQ_dVol(i+1) - s% QQ(k-1)*dCp_dVol(i+1)/s% Cp(k-1))/s% Cp(k-1)
               
            d_Y1_dr_00 = &
               POM2*( &
                  (dQQ_dr_00(I) - s% QQ(k)/s% Cp(k)*dCp_dr_00(I))/s% Cp(k)&
                 +(dQQ_dr_in(i+1) - s% QQ(k-1)/s% Cp(k-1)*dCp_dr_in(i+1))/s% Cp(k-1)) &
              + POM*(d_Pg_dr_in(i+1) - d_Pg_dr_00(I) + d_Pr_dr_in(i+1) - d_Pr_dr_00(I))
                    
            d_Y1_dr_in = &
                 POM2*(dQQ_dr_in(I) - s% QQ(k)/s% Cp(k)*dCp_dr_in(I))/s% Cp(k) &
               + POM*(- d_Pg_dr_in(I) - d_Pr_dr_in(I))
               
            d_Y1_dr_out = &
                 POM2*(dQQ_dr_00(i+1) - s% QQ(k-1)/s% Cp(k-1)*dCp_dr_00(i+1))/s% Cp(k-1) &
               + POM*(d_Pg_dr_00(i+1) + d_Pr_dr_00(i+1))
               
            d_Y1_dT_00 = &
                 POM2*(dQQ_dT(I) - s% QQ(k)/s% Cp(k)*dCp_dT(I))/s% Cp(k) &
               - POM*d_Pg_dT(I) &
               + 1.d0/s% T(k)
               
            d_Y1_dT_out = &
                 POM2*(dQQ_dT(i+1) - s% QQ(k-1)/s% Cp(k-1)*dCp_dT(i+1))/s% Cp(k-1) &
               + POM*d_Pg_dT(I+1) &
               - 1.d0/s% T(k-1)
               
            d_Y1_der_00 = -POM*d_Pr_der(I)
               
            d_Y1_der_out = POM*d_Pr_der(I+1)

            POM = 2.d0/(s% Vol(k) + s% Vol(k-1))
            POM2 = 8.d0*PI*(s% r(k)**2)/s% dm_bar(k)*s% Hp_face(k)
            Y2 = 4.d0*PI*(s% r(k)**2)*s% Hp_face(k)*POM/s% dm_bar(k)
            
            d_Y2_dVol_00 = &
               + Y2/s% Hp_face(k)*dHp_dVol_00(i) & 
               - POM2/(s% Vol(k) + s% Vol(k-1))**2
            d_Y2_dVol_out = &
               + Y2/s% Hp_face(k)*dHp_dVol_out(i) & 
               - POM2/(s% Vol(k) + s% Vol(k-1))**2
            
            d_Y2_dr_00 = 2.d0*Y2/s% r(k) &
               + Y2/s% Hp_face(k)*dHp_dr_00(I) &
               - POM2/(s% Vol(k) + s% Vol(k-1))**2*(dVol_dr_00(I) + dVol_dr_in(i+1))
               
            d_Y2_dr_in = - POM2/(s% Vol(k) &
               + s% Vol(k-1))**2*dVol_dr_in(I) &
               + Y2/s% Hp_face(k)*dHp_dr_in(I)
               
            d_Y2_dr_out = - POM2/(s% Vol(k) &
               + s% Vol(k-1))**2*dVol_dr_00(i+1) &
               + Y2/s% Hp_face(k)*dHp_dr_out(I)
               
            d_Y2_dT_00 = Y2/s% Hp_face(k)*dHp_dT_00(I)
            
            d_Y2_dT_out = Y2/s% Hp_face(k)*dHp_dT_out(I)
               
            d_Y2_der_00 = Y2/s% Hp_face(k)*dHp_der_00(I)
            
            d_Y2_der_out = Y2/s% Hp_face(k)*dHp_der_out(I)

            s% Y_face(k) = Y1*Y2
            
            if (k==-109) write(*,3) 'Y_face Y1 Y2', k, s% solver_iter, &
               s% Y_face(k), Y1, Y2

            dY_dr_00(I) = Y1*d_Y2_dr_00 + Y2*d_Y1_dr_00 ! 
            dY_dr_in(I) = Y1*d_Y2_dr_in + Y2*d_Y1_dr_in ! 
            dY_dr_out(I) = Y1*d_Y2_dr_out + Y2*d_Y1_dr_out ! 
            dY_dVol_00(I) = Y1*d_Y2_dVol_00 + Y2*d_Y1_dVol_00 ! 
            dY_dVol_out(I) = Y1*d_Y2_dVol_out + Y2*d_Y1_dVol_out ! 
            dY_dT_00(I) = Y1*d_Y2_dT_00 + Y2*d_Y1_dT_00 ! 
            dY_dT_out(I) = Y1*d_Y2_dT_out + Y2*d_Y1_dT_out ! 
            dY_der_00(I) = Y1*d_Y2_der_00 + Y2*d_Y1_der_00 ! 
            dY_der_out(I) = Y1*d_Y2_der_out + Y2*d_Y1_der_out ! 
            
            if (call_is_bad) then
               if (is_bad(s% Y_face(k))) then
                  !$OMP critical
                  write(*,2) 's% Y_face(k)', k, s% Y_face(k)
                  write(*,2) 'Y1', k, Y1
                  write(*,2) 'Y2', k, Y2
                  write(*,2) 's% QQ(k)', k, s% QQ(k)
                  write(*,2) 's% Cp(k)', k, s% Cp(k)
                  write(*,2) 's% Pgas(k)', k, s% Pgas(k)
                  write(*,2) 's% Prad(k)', k, s% Prad(k)
                  write(*,2) 's% T(k)', k, s% T(k)
                  stop 'calc_Y_face'
                  !$OMP end critical
               end if
            end if
         
         else
            s% Y_face(k) = 0
            dY_dr_00(I) = 0
            dY_dr_in(I) = 0
            dY_dr_out(I) = 0
            dY_dVol_00(I) = 0
            dY_dVol_out(I) = 0
            dY_dT_00(I) = 0
            dY_dT_out(I) = 0
            dY_der_00(I) = 0
            dY_der_out(I) = 0
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% Y_face(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dY_dVol_00(I)
            write(*,*) 'calc_Y_face', s% solver_test_partials_var
         end if
      
      end subroutine calc_Y_face

      
      subroutine calc_PII_face(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I
         real(dp) :: POM, POM2
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (k == 1 .or. k == s% nz .or. ALFA == 0d0) then
            s% PII(k) = 0
            dPII_dr_00(I) = 0
            dPII_dr_in(I) = 0 
            dPII_dr_out(I) = 0
            dPII_dVol_00(I) = 0
            dPII_dVol_out(I) = 0
            dPII_dT_00(I) = 0
            dPII_dT_out(I) = 0
            dPII_der_00(I) = 0
            dPII_der_out(I) = 0
         else         
            POM = ALFAS*ALFA
            POM2 = 0.5d0*(s% Cp(k) + s% Cp(k-1))
            s% PII(k) = POM*POM2*s% Y_face(k)

            dPII_dVol_00(I) = &
               POM*(POM2*dY_dVol_00(I) + s% Y_face(k)*0.5d0*dCp_dVol(I))
            dPII_dVol_out(I) = &
               POM*(POM2*dY_dVol_out(I) + s% Y_face(k)*0.5d0*dCp_dVol(I+1))
            dPII_dr_in(I) = & ! 
               POM*(POM2*dY_dr_in(I) + s% Y_face(k)*0.5d0*dCp_dr_in(I))
            dPII_dr_00(I) = & ! 
               POM*(POM2*dY_dr_00(I) + s% Y_face(k)*0.5d0*(dCp_dr_00(I) + dCp_dr_in(i+1)))
            dPII_dr_out(I) = & ! 
               POM*(POM2*dY_dr_out(I) + s% Y_face(k)*0.5d0*dCp_dr_00(i+1))
            dPII_dT_00(I) = & ! 
               POM*(POM2*dY_dT_00(I) + s% Y_face(k)*0.5d0*dCp_dT(I))
            dPII_dT_out(I) = & ! 
               POM*(POM2*dY_dT_out(I) + s% Y_face(k)*0.5d0*dCp_dT(i+1))
            dPII_der_00(I) = POM*POM2*dY_der_00(I) ! 
            dPII_der_out(I) = POM*POM2*dY_der_out(I) ! 
         end if

         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% PII(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dPII_dVol_out(I)
            write(*,*) 'calc_PII_face', s% solver_test_partials_var
         end if
      
      end subroutine calc_PII_face
      
      
      subroutine calc_Pvsc(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: &
            dv, dv1, P, dP_dT, dP_der, dP_dr_in, dP_dr_00, V, sqrt_PV, &
            d_PV_dT, d_PV_der, d_PV_dr_in, d_PV_dr_00, &
            d_sqrt_PV_dT, d_sqrt_PV_der, d_sqrt_PV_dr_in, d_sqrt_PV_dr_00, &
            d_dv_dT, d_dv_der, d_dv_dr_in, d_dv_dr_00, &
            dP_dVol, d_PV_dVol, d_sqrt_PV_dVol, d_dv_dVol
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (CQ == 0d0) then
            s% Pvsc(k) = 0d0
            d_Pvsc_dVol(i) = 0d0
            d_Pvsc_dT(i) = 0d0
            d_Pvsc_der(i) = 0d0
            d_Pvsc_dr_in(i) = 0d0
            d_Pvsc_dr_00(i) = 0d0
            return
         end if
         if (I > 1) then
            dv = s% v(k+1) - s% v(k)
         else
            dv = s% v_center - s% v(NZN)
         end if
         P = s% Pgas(k) + s% Prad(k)
         dP_dT = d_Pg_dT(I)
         dP_dVol = d_Pg_dVol(I) + d_Pr_dVol(I)
         dP_der = d_Pr_der(I)
         dP_dr_in = d_Pg_dr_in(I) + d_Pr_dr_in(I)
         dP_dr_00 = d_Pg_dr_00(I) + d_Pr_dr_00(I)
         V = s% Vol(k)
         sqrt_PV = sqrt(P*V)
         if (dv <= ZSH*sqrt_PV) then
            s% Pvsc(k) = 0d0
            d_Pvsc_dVol(i) = 0d0
            d_Pvsc_dT(i) = 0d0
            d_Pvsc_der(i) = 0d0
            d_Pvsc_dr_in(i) = 0d0
            d_Pvsc_dr_00(i) = 0d0
            return
         end if

         d_PV_dT = dP_dT*V
         d_PV_dVol = P + dP_dVol*V
         d_PV_der = dP_der*V
         d_PV_dr_in = dP_dr_in*V + P*dVol_dr_in(I)
         d_PV_dr_00 = dP_dr_00*V + P*dVol_dr_00(I)
         
         d_sqrt_PV_dT = 0.5d0*d_PV_dT/sqrt_PV
         d_sqrt_PV_dVol = 0.5d0*d_PV_dVol/sqrt_PV
         d_sqrt_PV_der = 0.5d0*d_PV_der/sqrt_PV
         d_sqrt_PV_dr_in = 0.5d0*d_PV_dr_in/sqrt_PV
         d_sqrt_PV_dr_00 = 0.5d0*d_PV_dr_00/sqrt_PV
         
         dv1 = dv
         dv = dv - ZSH*sqrt_PV
         
         d_dv_dT = - ZSH*d_sqrt_PV_dT
         d_dv_dVol = - ZSH*d_sqrt_PV_dVol
         d_dv_der = - ZSH*d_sqrt_PV_der
         
         d_dv_dr_in = 2d0/s% dt - ZSH*d_sqrt_PV_dr_in ! not used if I == 1
         d_dv_dr_00 = -2d0/s% dt - ZSH*d_sqrt_PV_dr_00
         
         ! Pvsc = CQ*P*(dv1/sqrt_PV - cut)^2  eqn 3.60
         !     = CQ*P*((dv1 - cut*sqrt_PV)/sqrt_PV)^2
         !     = CQ*P/(P*V)*dv^2
         !     = CQ/V*dv^2
         
         s% Pvsc(k) = CQ/V*dv**2
         d_Pvsc_dVol(i) = -s% Pvsc(k)/V + 2d0*d_dv_dVol*CQ/V*dv
         d_Pvsc_dT(i) = CQ/V*2d0*dv*d_dv_dT
         d_Pvsc_der(i) = CQ/V*2d0*dv*d_dv_der
         d_Pvsc_dr_in(i) = CQ/V*2d0*dv*d_dv_dr_in - CQ*dv**2*dVol_dr_in(I)/V**2
         d_Pvsc_dr_00(i) = CQ/V*2d0*dv*d_dv_dr_00 - CQ*dv**2*dVol_dr_00(I)/V**2

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% Pvsc(k)
            s% solver_test_partials_var = i_var_T
            s% solver_test_partials_dval_dx = d_Pvsc_dT(i)
            write(*,*) 'calc_Pvsc', s% solver_test_partials_var
         end if
      end subroutine calc_Pvsc

      
      subroutine check_omega(s,i) ! needs cleanup
         type (star_info), pointer :: s
         integer, intent(in) :: i   
         real(dp) :: SOURS, DAMPS, DAMPRS, DELTA, SOL, POM, POM2, POM3
         integer :: k
         if (I > IBOTOM .and. I < NZN .and. ALFA /= 0d0) then
         !     JAK OKRESLIC OMEGA DLA PIERWSZEJ ITERACJI
            k = NZN+1-i
            if (s% RSP_w(k) > EFL0) return
            POM = (s% PII(k)/s% Hp_face(k) + s% PII(k+1)/s% Hp_face(k+1))*0.5d0
            POM2 = s% T(k)*(s% Pgas(k) + s% Prad(k))*s% QQ(k)/s% Cp(k)
            SOURS = POM*POM2
            DAMPS = (CEDE/ALFA)/((s% Hp_face(k) + s% Hp_face(k+1))*0.5d0)
            POM3 = (GAMMAR**2)/(ALFA**2)*4.d0*SIG
            POM2 = (s% T(k)**3)*(s% Vol(k)**2)/(s% Cp(k)*s% opacity(k))       
            DAMPRS = POM3*POM2/((s% Hp_face(k)**2 + s% Hp_face(k+1)**2)*0.5d0)
            DELTA = DAMPRS**2 + 4.d0*DAMPS*SOURS
            if (DELTA >= 0.d0) SOL = ( - DAMPRS + sqrt(DELTA))/(2.d0*DAMPS)
            if (DELTA < 0.d0) SOL = - 99.99d0
            if (SOL >= 0.d0) SOL = SOL**2
            if (SOL > 0.d0) s% RSP_w(k) = sqrt(SOL)
         end if
      end subroutine check_omega
      
      
      subroutine calc_Pturb(s,i) ! TURBULENT PRESSURE (ZONE)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: TEM1, Vol
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (ALFA == 0d0 .or. ALFAP == 0d0 .or. &
             I <= IBOTOM .or. I >= NZN) then
            s% Ptrb(k) = 0.d0
            dPtrb_dVol_00(I) = 0.d0
            dPtrb_dw_00(I) = 0.d0
            dPtrb_dr_00(I) = 0.d0
            dPtrb_dr_in(I) = 0.d0
         else
            Vol = s% Vol(k)
            s% Ptrb(k) = ALFAP*s% RSP_w(k)**2/Vol
            dPtrb_dVol_00(I) = -s% Ptrb(k)/Vol
            dPtrb_dw_00(I) = 2.d0*ALFAP*s% RSP_w(k)/Vol
            TEM1 = - ALFAP*s% RSP_w(k)**2/Vol**2
            dPtrb_dr_00(I) = TEM1*dVol_dr_00(I)
            dPtrb_dr_in(I) = TEM1*dVol_dr_in(I)
         end if
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = 0 ! residual
            s% solver_test_partials_var = i_var_R
            s% solver_test_partials_dval_dx = 0 ! d_residual_dr_00
            write(*,*) 'calc_Pturb', s% solver_test_partials_var
         end if
      end subroutine calc_Pturb
      
      
      subroutine calc_Chi(s,i) ! eddy viscosity (Kuhfuss 1986)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: POM, POM1, POM2, POM3, POM4, &
            POMT1, POMT2, POMT3, POMT4, Vol
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (ALFA == 0d0 .or. I <= IBOTOM .or. I >= NZN) then
            s% Chi(k) = 0      
            dChi_dT_out(I) = 0  
            dChi_dT_00(I) = 0
            dChi_dT_in(I) = 0
            dChi_der_out(I) = 0  
            dChi_der_00(I) = 0
            dChi_der_in(I) = 0
            dChi_dw_00(I) = 0            
            dChi_dr_out(I) = 0
            dChi_dr_00(I) = 0
            dChi_dr_in(I) = 0
            dChi_dr_in2(I) = 0
            !s% profile_extra(k,2) = 0
            !s% profile_extra(k,3) = 0
            !s% profile_extra(k,4) = 0
         else
            POM = (16.d0/3.d0)*PI*ALFA*ALFAM/s% dm(k)  
            Vol = s% Vol(k)
            POM1 = s% RSP_w(k)/Vol**2
            POM2 = 0.5d0*(s% r(k)**6 + s% r(k+1)**6)
            POM4 = 0.5d0*(s% Hp_face(k) + s% Hp_face(k+1))
            POM3 = s% v(k)/s% r(k) - s% v(k+1)/s% r(k+1)
            POMT3 = POM*POM1*POM2*POM4
            POMT1 = POM*POM2*POM3*POM4
            POMT2 = POM*POM1*POM3*POM4
            POMT4 = POM*POM1*POM2*POM3

            s% Chi(k) = POMT1*POM1

            if (call_is_bad) then
               if (is_bad(s% Chi(k))) then
                  !$OMP critical
                  write(*,2) 'POM', k, POM
                  write(*,2) 'POM1', k, POM1
                  write(*,2) 'POM2', k, POM2
                  write(*,2) 'POM3', k, POM3
                  write(*,2) 'POM4', k, POM4
                  write(*,2) 's% RSP_w(k)', k, s% RSP_w(k)
                  write(*,2) 's% Volk)', k, s% Vol(k)
                  write(*,2) 's% r(k)', k, s% r(k)
                  write(*,2) 's% r(k+1)', k+1, s% r(k+1)
                  write(*,2) 's% v(k)', k, s% v(k)
                  write(*,2) 's% v(k+1)', k+1, s% v(k+1)
                  write(*,2) 's% Vol(k+1)', k+1, s% Vol(k+1)
                  write(*,2) 's% rho(k+1)', k+1, s% rho(k+1)
                  stop 'calc_Chi'
                  !$OMP end critical
               end if
            end if
         
            dChi_dVol_out(I) = &
                 POMT4*0.5d0*dHp_dVol_out(I)
            dChi_dVol_00(I) = &
               + POMT4*0.5d0*(dHp_dVol_00(I) + dHp_dVol_out(i-1)) &
               - 2d0*s% Chi(k)/Vol
            dChi_dVol_in(I) = &
                 POMT4*0.5d0*dHp_dVol_00(i-1)

            dChi_dT_out(I) = POMT4*0.5d0*dHp_dT_out(I) ! 
            dChi_dT_00(I) = POMT4*0.5d0*(dHp_dT_00(I) + dHp_dT_out(i-1)) ! 
            dChi_dT_in(I) = POMT4*0.5d0*dHp_dT_00(i-1) ! 

            dChi_der_out(I) = POMT4*0.5d0*(dHp_der_out(I)) ! 
            dChi_der_00(I) = POMT4*0.5d0*(dHp_der_00(I) + dHp_der_out(i-1)) ! 
            dChi_der_in(I) = POMT4*0.5d0*(dHp_der_00(i-1)) ! 

            dChi_dw_00(I) = POMT1/Vol**2 ! 
            
            dChi_dr_out(I) = POMT4*0.5d0*(dHp_dr_out(I)) ! 
            dChi_dr_00(I) = &
               - 2.d0*s% Chi(k)/Vol*dVol_dr_00(I) & ! 
               + POMT2*3.d0*s% r(k)**5 &
               + POMT3*(2.d0/s% dt/s% r(k) - s% v(k)/s% r(k)**2) &
               + POMT4*0.5d0*(dHp_dr_00(I) + dHp_dr_out(i-1))
            dChi_dr_in(I) = &
               - 2.d0*s% Chi(k)/Vol*dVol_dr_in(I) & ! 
               + POMT2*3.d0*s% r(k+1)**5 &
               - POMT3*(2.d0/s% dt/s% r(k+1) - s% v(k+1)/s% r(k+1)**2) &
               + POMT4*0.5d0*(dHp_dr_in(I) + dHp_dr_00(i-1))
            dChi_dr_in2(I) = POMT4*0.5d0*dHp_dr_in(i-1) ! 
            
            if (call_is_bad) then
               if (is_bad(dChi_dr_in(I))) then
                  !$OMP critical
                  write(*,2) 's% Chi(k)', k, s% Chi(k)
                  write(*,2) 'Vol', k, Vol
                  write(*,2) 'POMT2', k, POMT2
                  write(*,2) 'POMT3', k, POMT3
                  write(*,2) 'POMT4', k, POMT4
                  write(*,2) 'dVol_dr_in(I)', I, dVol_dr_in(I)
                  write(*,2) 'dChi_dr_in(I)', I, dChi_dr_in(I)
                  write(*,2) 'dHp_dr_00(i-1)', I-1, dHp_dr_00(i-1)
                  write(*,2) 'dHp_dr_in(I)', I, dHp_dr_in(I)
                  write(*,2) 's% r(k+1)', k+1, s% r(k+1)
                  write(*,2) 's% v(k+1)', k+1, s% v(k+1)
                  stop 'calc_Chi'
                  !$OMP end critical
               end if
            end if
            
         end if
            
            !if (k==194) then
            !   write(*,2) 'RSP Chi', k, s% Chi(k)
            !end if
         
         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% Chi(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dChi_dVol_out(I)
            write(*,*) 'calc_Chi', s% solver_test_partials_var
         end if
      end subroutine calc_Chi


      subroutine calc_Eq(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: RRI, RRM, UUI, UUM, POM, POM2
         integer :: k
         logical :: test_partials, smoothed
         include 'formats'
         k = NZN+1-i
         if (ALFA == 0d0 .or. I <= IBOTOM .or. I >= NZN) then
            s% Eq(k) = 0
            dEq_dr_out(I) = 0
            dEq_dr_00(I) = 0
            dEq_dr_in(I) = 0
            dEq_dr_in2(I) = 0
            dEq_dVol_out(I) = 0
            dEq_dVol_00(I) = 0
            dEq_dVol_in(I) = 0
            dEq_dT_out(I) = 0
            dEq_dT_00(I) = 0
            dEq_dT_in(I) = 0
            dEq_der_out(I) = 0
            dEq_der_00(I) = 0
            dEq_der_in(I) = 0
            dEq_dw_00(I) = 0
         else
         
            RRI = 0.5d0*(s% r(k) + s% r_start(k))
            RRM = 0.5d0*(s% r(k+1) + s% r_start(k+1))
            UUI = 0.5d0*(s% v(k) + s% v_start(k))
            UUM = 0.5d0*(s% v(k+1) + s% v_start(k+1))
            
            POM = P4/s% dm(k)*(UUI/RRI - UUM/RRM)
            POM2 = P4/s% dm(k)*(THETAU*s% Chi(k) + THETAU1*s% Chi_start(k))
         
            s% Eq(k) = POM*(THETAU*s% Chi(k) + THETAU1*s% Chi_start(k))
      
            POM = POM*THETAU

            dEq_dVol_out(I) = POM*THETAU*dChi_dVol_out(i)
            dEq_dVol_00(I) = POM*THETAU*dChi_dVol_00(i)
            dEq_dVol_in(I) = POM*THETAU*dChi_dVol_in(i)
            
            dEq_dr_out(I) = POM*dChi_dr_out(I) ! 
            dEq_dr_00(I) = POM*dChi_dr_00(I) & ! 
               + POM2*(1.d0/(RRI*s% dt) - 0.5d0*UUI/(RRI**2))
            dEq_dr_in(I) = POM*dChi_dr_in(I) & ! 
               - POM2*(1.d0/(RRM*s% dt) - 0.5d0*UUM/(RRM**2))
            dEq_dr_in2(I) = POM*dChi_dr_in2(I) ! 
            dEq_dT_out(I) = POM*dChi_dT_out(I) ! 
            dEq_dT_00(I) = POM*dChi_dT_00(I) ! 
            dEq_dT_in(I) = POM*dChi_dT_in(I) ! 
            dEq_der_out(I) = POM*dChi_der_out(I) ! 
            dEq_der_00(I) = POM*dChi_der_00(I) ! 
            dEq_der_in(I) = POM*dChi_der_in(I) ! 
            dEq_dw_00(I) = POM*dChi_dw_00(I) ! 
         end if
         
         !test_partials = (k+1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% Eq(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dEq_dVol_in(I)
            write(*,*) 'calc_Eq', s% solver_test_partials_var
         end if
      end subroutine calc_Eq
      
      
      subroutine calc_source_sink(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: I      
         real(dp) :: POM, POM2, POM3, POM4, TEM1, TEMI, TEMM, &
            dsrc_dr_in2, dsrc_dr_in, dsrc_dr_00, dsrc_dr_out, &
            dsrc_dVol_out, dsrc_dVol_00, dsrc_dVol_in, &
            dsrc_dT_out, dsrc_dT_00, dsrc_dT_in, &
            dsrc_der_out, dsrc_der_00, dsrc_der_in, &
            dsrc_dw_00, &
            d_damp_dr_in2, d_damp_dr_in, d_damp_dr_00, d_damp_dr_out, &
            d_damp_dVol_out, d_damp_dVol_00, d_damp_dVol_in, &
            d_damp_dT_out, d_damp_dT_00, d_damp_dT_in, &
            d_damp_der_out, d_damp_der_00, d_damp_der_in, &
            d_damp_dw_00, &
            d_dampR_dr_in2, d_dampR_dr_in, d_dampR_dr_00, d_dampR_dr_out, &
            d_dampR_dVol_out, d_dampR_dVol_00, d_dampR_dVol_in, &
            d_dampR_dT_out, d_dampR_dT_00, d_dampR_dT_in, &
            d_dampR_der_out, d_dampR_der_00, d_dampR_der_in, &
            d_dampR_dw_00, d_POM_dVol_00, d_POM_dVol_in, &
            d_QQ_div_Cp_d_Vol_00, d_POM2_dVol_00, QQ_div_Cp, &
            POM2a, POM2b, d_POM2a_dVol_00, d_POM2b_dVol_00, &
            d_POM4_dVol_00
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (ALFA == 0d0 .or. I <= IBOTOM .or. I >= NZN) then
            s% SOURCE(k) = 0
            s% DAMP(k) = 0
            s% DAMPR(k) = 0
            s% COUPL(k) = 0            
         else
            ! SOURCE TERM
            POM = 0.5d0*(s% PII(k)/s% Hp_face(k) + s% PII(k+1)/s% Hp_face(k+1))
            QQ_div_Cp = s% QQ(k)/s% Cp(k)
            POM2 = s% T(k)*(s% Pgas(k) + s% Prad(k))*QQ_div_Cp
            POM3 = s% RSP_w(k)            
            s% SOURCE(k) = POM*POM2*POM3
         
            ! P*QQ/Cp = grad_ad
            if (k==-109) write(*,3) 'w grada PII_00 PII_p1 SOURCE', k, s% solver_iter, &
               s% RSP_w(k), (s% Pgas(k) + s% Prad(k))*QQ_div_Cp, s% PII(k), &
               s% PII(k+1), s% SOURCE(k)
      
            TEM1 = POM2*POM3*0.5d0
            TEMI = - s% PII(k)/s% Hp_face(k)**2
            TEMM = - s% PII(k+1)/s% Hp_face(k+1)**2

            d_POM_dVol_00 =  &
               0.5d0*(dPII_dVol_00(I) - s% PII(k)/s% Hp_face(k)*dHp_dVol_00(I))/s% Hp_face(k) + &
               0.5d0*(dPII_dVol_out(I-1) - s% PII(k+1)/s% Hp_face(k+1)*dHp_dVol_out(I-1))/s% Hp_face(k+1)
            d_POM_dVol_in = 0.5d0*( &
               dPII_dVol_00(I-1) - s% PII(k+1)/s% Hp_face(k+1)*dHp_dVol_00(I-1))/s% Hp_face(k+1)
            d_QQ_div_Cp_d_Vol_00 = (dQQ_dVol(i) - dCp_dVol(i)*s% QQ(k)/s% Cp(k))/s% Cp(k)
            d_POM2_dVol_00 = s% T(k)*( &
                  (d_Pg_dVol(i) + d_Pr_dVol(i))*QQ_div_Cp + &
                  (s% Pgas(k) + s% Prad(k))*d_QQ_div_Cp_d_Vol_00)
            dsrc_dVol_out = TEM1*( & ! ok
                  dPII_dVol_out(I)/s% Hp_face(k) &
                + TEMI*dHp_dVol_out(I))
            dsrc_dVol_00 = POM3*(d_POM_dVol_00*POM2 + POM*d_POM2_dVol_00)
            dsrc_dVol_in = TEM1*( & ! ok
                  dPII_dVol_00(i-1)/s% Hp_face(k+1) &
                + TEMM*dHp_dVol_00(i-1))

            dsrc_dr_out = TEM1*(dPII_dr_out(I)/s% Hp_face(k) &
               + TEMI*dHp_dr_out(I))
            dsrc_dr_00 = TEM1*(dPII_dr_00(I)/s% Hp_face(k) + dPII_dr_out(i-1)/s% Hp_face(k+1) &
               + TEMI*dHp_dr_00(I) + TEMM*dHp_dr_out(i-1))
            dsrc_dr_in = TEM1*(dPII_dr_in(I)/s% Hp_face(k) + dPII_dr_00(i-1)/s% Hp_face(k+1) &
               + TEMI*dHp_dr_in(I) + TEMM*dHp_dr_00(i-1))
            dsrc_dr_in2 = TEM1*(dPII_dr_in(i-1)/s% Hp_face(k+1) &
               + TEMM*dHp_dr_in(i-1))

            dsrc_dT_out = TEM1*( &
                  dPII_dT_out(I)/s% Hp_face(k) &
                + TEMI*dHp_dT_out(I))
            dsrc_dT_00 = TEM1*( &
                  dPII_dT_00(I)/s% Hp_face(k) &
                + dPII_dT_out(i-1)/s% Hp_face(k+1) &
                + TEMI*dHp_dT_00(I) &
                + TEMM*dHp_dT_out(i-1))
            dsrc_dT_in = TEM1*( &
                  dPII_dT_00(i-1)/s% Hp_face(k+1) &
                + TEMM*dHp_dT_00(i-1))

            dsrc_der_out = TEM1*( &
                  dPII_der_out(I)/s% Hp_face(k) &
                + TEMI*dHp_der_out(I))
            dsrc_der_00 = TEM1*( &
                  dPII_der_00(I)/s% Hp_face(k) &
                + dPII_der_out(i-1)/s% Hp_face(k+1) &
                + TEMI*dHp_der_00(I) &
                + TEMM*dHp_der_out(i-1))
            dsrc_der_in = TEM1*( &
                  dPII_der_00(i-1)/s% Hp_face(k+1) &
                + TEMM*dHp_der_00(i-1))

            dsrc_dw_00 = POM*POM2

            POM = POM*POM3

            dsrc_dT_00 = dsrc_dT_00 + POM/s% Cp(k)*( &
                 (s% Pgas(k) + s% Prad(k))*s% QQ(k) &
               + s% T(k)*s% QQ(k)*d_Pg_dT(I) &
               + s% T(k)*(s% Pgas(k) + s% Prad(k))*dQQ_dT(I) &
               - s% T(k)*(s% Pgas(k) + s% Prad(k))*s% QQ(k)/s% Cp(k)*dCp_dT(I))
               
            dsrc_der_00 = dsrc_der_00 &
               + POM/s% Cp(k)*s% T(k)*s% QQ(k)*d_Pr_der(I)

            dsrc_dr_00 = dsrc_dr_00 + POM*s% T(k)/s% Cp(k)*( &
                 s% QQ(k)*(d_Pg_dr_00(I) + d_Pr_dr_00(I)) &
               + (s% Pgas(k) + s% Prad(k))*dQQ_dr_00(I) &
               - (s% Pgas(k) + s% Prad(k))*s% QQ(k)/s% Cp(k)*dCp_dr_00(I))
            dsrc_dr_in = dsrc_dr_in + POM*s% T(k)/s% Cp(k)*( &
                 s% QQ(k)*(d_Pg_dr_in(I) + d_Pr_dr_in(I)) &
               + (s% Pgas(k) + s% Prad(k))*dQQ_dr_in(I) &
               - (s% Pgas(k) + s% Prad(k))*s% QQ(k)/s% Cp(k)*dCp_dr_in(I))

            ! DAMP TERM
            POM = (CEDE/ALFA)*(s% RSP_w(k)**3 - EFL0**3)
            POM2 = 0.5d0*(s% Hp_face(k) + s% Hp_face(k+1))
            s% DAMP(k) = POM/POM2
      
            TEM1 = - 0.5d0*POM/POM2**2
            
            d_damp_dVol_out = TEM1*dHp_dVol_out(I)
            d_damp_dVol_00 = TEM1*(dHp_dVol_00(I) + dHp_dVol_out(i-1))
            d_damp_dVol_in = TEM1*dHp_dVol_00(i-1)
            d_damp_dr_out = TEM1*dHp_dr_out(I)
            d_damp_dr_00 = TEM1*(dHp_dr_00(I) + dHp_dr_out(i-1))
            d_damp_dr_in = TEM1*(dHp_dr_in(I) + dHp_dr_00(i-1))
            d_damp_dr_in2 = TEM1*(dHp_dr_in(i-1))
            d_damp_dT_00 = TEM1*(dHp_dT_00(I) + dHp_dT_out(i-1))
            d_damp_dT_out = TEM1*dHp_dT_out(I)
            d_damp_dT_in = TEM1*dHp_dT_00(i-1)
            d_damp_der_00 = TEM1*(dHp_der_00(I) + dHp_der_out(i-1))
            d_damp_der_out = TEM1*dHp_der_out(I)
            d_damp_der_in = TEM1*dHp_der_00(i-1)
            d_damp_dw_00 = 3.0d0*(CEDE/ALFA)/POM2*s% RSP_w(k)**2
            
            ! RADIATIVE DAMP TERM
            if (GAMMAR == 0.d0)then
               s% DAMPR(k) = 0.d0
               d_dampR_dr_out = 0.d0
               d_dampR_dr_00 = 0.d0
               d_dampR_dr_in = 0.d0
               d_dampR_dr_in2 = 0.d0
               d_dampR_dVol_out = 0.d0
               d_dampR_dVol_in = 0.d0
               d_dampR_dVol_00 = 0.d0
               d_dampR_dT_out = 0.d0
               d_dampR_dT_in = 0.d0
               d_dampR_dT_00 = 0.d0
               d_dampR_der_out = 0.d0
               d_dampR_der_in = 0.d0
               d_dampR_der_00 = 0.d0
               d_dampR_dw_00 = 0.d0
            else   
               POM = (GAMMAR**2)/(ALFA**2)*4.d0*SIG
               POM2a = s% T(k)**3*s% Vol(k)**2
               POM2b = 1d0/(s% Cp(k)*s% opacity(k))
               POM2 = POM2a*POM2b
               POM3 = s% RSP_w(k)**2
               POM4 = 0.5d0*(s% Hp_face(k)**2 + s% Hp_face(k+1)**2)
               s% DAMPR(k) = POM*POM2*POM3/POM4
      
               TEM1 = - s% DAMPR(k)/POM4

               d_POM2a_dVol_00 = 2d0*s% T(k)**3*s% Vol(k)
               d_POM2b_dVol_00 = &
                  -POM2b*(dCp_dVol(i)/s% Cp(k) + dK_dVol(i)/s% opacity(k))
               d_POM2_dVol_00 = d_POM2a_dVol_00*POM2b + POM2a*d_POM2b_dVol_00
               d_POM4_dVol_00 = &
                  s% Hp_face(k)*dHp_dVol_00(I) + &
                  s% Hp_face(k+1)*dHp_dVol_out(I-1)
               
               d_dampR_dVol_out = TEM1*s% Hp_face(k)*dHp_dVol_out(I)
               d_dampR_dVol_00 = POM*POM3*( &
                  d_POM2_dVol_00 - POM2*d_POM4_dVol_00/POM4)/POM4
               
               d_dampR_dVol_in = TEM1*s% Hp_face(k+1)*dHp_dVol_00(i-1)
               
               d_dampR_dr_out = TEM1*s% Hp_face(k)*dHp_dr_out(I)
               d_dampR_dr_00 = TEM1*(s% Hp_face(k)*dHp_dr_00(I) &
                  + s% Hp_face(k+1)*dHp_dr_out(i-1))
               d_dampR_dr_in = TEM1*(s% Hp_face(k)*dHp_dr_in(I) &
                  + s% Hp_face(k+1)*dHp_dr_00(i-1))
               d_dampR_dr_in2 = TEM1*s% Hp_face(k+1)*dHp_dr_in(i-1)
               
               d_dampR_dT_out = TEM1*s% Hp_face(k)*dHp_dT_out(I)
               d_dampR_dT_00 = TEM1*(s% Hp_face(k)*dHp_dT_00(I) &
                  + s% Hp_face(k+1)*dHp_dT_out(i-1))
               d_dampR_dT_in = TEM1*s% Hp_face(k+1)*dHp_dT_00(i-1)
               
               d_dampR_der_out = TEM1*s% Hp_face(k)*dHp_der_out(I)
               d_dampR_der_00 = TEM1*(s% Hp_face(k)*dHp_der_00(I) &
                  + s% Hp_face(k+1)*dHp_der_out(i-1))
               d_dampR_der_in = TEM1*s% Hp_face(k+1)*dHp_der_00(i-1)

               d_dampR_dw_00 = POM*POM2/POM4*2.d0*s% RSP_w(k)

               TEM1 = POM*POM3/POM4
               d_dampR_dr_00 = d_dampR_dr_00 &
                  + TEM1*s% T(k)**3*(2.d0*s% Vol(k)*dVol_dr_00(I) &
                  - s% Vol(k)**2*(1.d0/s% Cp(k)*dCp_dr_00(I) &
                  + 1.d0/s% opacity(k)*dK_dr_00(I))) &
                  /(s% Cp(k)*s% opacity(k))              

               d_dampR_dr_in = d_dampR_dr_in &
                  + TEM1*s% T(k)**3*(&
                       2.d0*s% Vol(k)*dVol_dr_in(I) &
                     - s% Vol(k)**2* &
                        (dCp_dr_in(I)/s% Cp(k) + dK_dr_in(I)/s% opacity(k))) &
                  /(s% Cp(k)*s% opacity(k))

               d_dampR_dT_00 = d_dampR_dT_00 &
                  + TEM1*s% Vol(k)**2*(3.d0*s% T(k)**2 &
                  - s% T(k)**3*(1.d0/s% Cp(k)*dCp_dT(I) &
                  + 1.d0/s% opacity(k)*dK_dT(I))) &
                  /(s% Cp(k)*s% opacity(k))

            end if
 
            s% COUPL(k) = s% SOURCE(k) - s% DAMP(k) - s% DAMPR(k)
            dC_dr_00(I) = dsrc_dr_00 - d_damp_dr_00 - d_dampR_dr_00 ! 
            dC_dr_out(I) = dsrc_dr_out - d_damp_dr_out - d_dampR_dr_out ! 
            dC_dr_in(I) = dsrc_dr_in - d_damp_dr_in - d_dampR_dr_in ! 
            dC_dr_in2(I) = dsrc_dr_in2 - d_damp_dr_in2 - d_dampR_dr_in2 ! 
            dC_dVol_00(I) = dsrc_dVol_00 - d_damp_dVol_00 - d_dampR_dVol_00 ! 
            dC_dVol_out(I) = dsrc_dVol_out - d_damp_dVol_out - d_dampR_dVol_out ! 
            dC_dVol_in(I) = dsrc_dVol_in - d_damp_dVol_in - d_dampR_dVol_in ! 
            dC_dT_00(I) = dsrc_dT_00 - d_damp_dT_00 - d_dampR_dT_00 ! 
            dC_dT_out(I) = dsrc_dT_out - d_damp_dT_out - d_dampR_dT_out ! 
            dC_dT_in(I) = dsrc_dT_in - d_damp_dT_in - d_dampR_dT_in ! 
            dC_der_00(I) = dsrc_der_00 - d_damp_der_00 - d_dampR_der_00 ! 
            dC_der_out(I) = dsrc_der_out - d_damp_der_out - d_dampR_der_out ! 
            dC_der_in(I) = dsrc_der_in - d_damp_der_in - d_dampR_der_in ! 
            dC_dw_00(I) = dsrc_dw_00 - d_damp_dw_00 - d_dampR_dw_00 ! 

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = s% COUPL(k)
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dC_dVol_00(I)
            write(*,*) 'calc_source_sink', s% solver_test_partials_var
         end if

         end if

      end subroutine calc_source_sink


      subroutine zero_boundaries(s)
         type (star_info), pointer :: s
         integer :: I, k
         do I = 1,IBOTOM
            k = NZN+1-i
            s% Eq(k) = 0.d0
            dEq_dr_in2(I) = 0.d0
            dEq_dr_in(I) = 0.d0
            dEq_dr_00(I) = 0.d0
            dEq_dr_out(I) = 0.d0
            dEq_dT_in(I) = 0.d0
            dEq_dT_00(I) = 0.d0
            dEq_dT_out(I) = 0.d0
            dEq_dw_00(I) = 0.d0! - 
            s% Chi(k) = 0.d0
            dChi_dr_in2(I) = 0.d0
            dChi_dr_in(I) = 0.d0
            dChi_dr_00(I) = 0.d0
            dChi_dr_out(I) = 0.d0
            dChi_dT_in(I) = 0.d0
            dChi_dT_00(I) = 0.d0
            dChi_dT_out(I) = 0.d0
            dChi_dw_00(I) = 0.d0! - 
            s% COUPL(k) = 0.d0
            dC_dr_00(I) = 0.d0
            dC_dr_out(I) = 0.d0
            dC_dr_in(I) = 0.d0
            dC_dr_in2(I) = 0.d0
            dC_dT_in(I) = 0.d0
            dC_dT_00(I) = 0.d0
            dC_dT_out(I) = 0.d0
            dC_dw_00(I) = 0.d0! - 
            s% Ptrb(k) = 0.d0
            dPtrb_dr_00(I) = 0.d0 
            dPtrb_dr_in(I) = 0.d0
            dPtrb_dw_00(I) = 0.d0! - 
         end do
         do I = NZN,NZN
            k = NZN+1-i
            s% Eq(k) = 0.d0
            dEq_dr_in2(I) = 0.d0
            dEq_dr_in(I) = 0.d0
            dEq_dr_00(I) = 0.d0
            dEq_dr_out(I) = 0.d0
            dEq_dT_in(I) = 0.d0
            dEq_dT_00(I) = 0.d0
            dEq_dT_out(I) = 0.d0
            dEq_dw_00(I) = 0.d0! - 
            s% Chi(k) = 0.d0
            dChi_dr_in2(I) = 0.d0
            dChi_dr_in(I) = 0.d0
            dChi_dr_00(I) = 0.d0
            dChi_dr_out(I) = 0.d0
            dChi_dT_in(I) = 0.d0
            dChi_dT_00(I) = 0.d0
            dChi_dT_out(I) = 0.d0
            dChi_dw_00(I) = 0.d0! - 
            s% COUPL(k) = 0.d0
            dC_dr_00(I) = 0.d0
            dC_dr_out(I) = 0.d0
            dC_dr_in(I) = 0.d0
            dC_dr_in2(I) = 0.d0
            dC_dT_in(I) = 0.d0
            dC_dT_00(I) = 0.d0
            dC_dT_out(I) = 0.d0
            dC_dw_00(I) = 0.d0! - 
            s% Ptrb(k) = 0.d0
            dPtrb_dr_00(I) = 0.d0
            dPtrb_dr_in(I) = 0.d0
            dPtrb_dw_00(I) = 0.d0! - 
         end do
      end subroutine zero_boundaries

         
      subroutine calc_Lt(s,i,Lt_00, &
            dLt_dr_00, dLt_dr_in, dLt_dr_out, &
            dLt_dVol_00, dLt_dVol_out, &
            dLt_dT_00, dLt_dT_out, &
            dLt_der_00, dLt_der_out, &
            dLt_dw_00, dLt_dw_out)
         type (star_info), pointer :: s
         integer, intent(in) :: I         
         real(dp), intent(out) :: &
            Lt_00, &
            dLt_dr_00, dLt_dr_in, dLt_dr_out, &
            dLt_dVol_00, dLt_dVol_out, &
            dLt_dT_00, dLt_dT_out, &
            dLt_der_00, dLt_der_out, &
            dLt_dw_00, dLt_dw_out
         real(dp) :: POM, POM2, POM3, TEM1, TEM2, &
            d_POM2_dVol_00, d_POM2_dVol_out, rho2_face, &
            d_rho2_face_dVol_00, d_rho2_face_dVol_out
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (I <= IBOTOM .or. I == NZN .or. ALFA == 0d0 .or. &
             ALFAT == 0.d0 .or. k < min_k_for_turbulent_flux) then
            Lt_00 = 0.d0
            dLt_dr_00 = 0.d0
            dLt_dr_in = 0.d0
            dLt_dr_out = 0.d0
            dLt_dVol_00 = 0.d0 
            dLt_dVol_out = 0.d0
            dLt_dT_00 = 0.d0 
            dLt_dT_out = 0.d0
            dLt_der_00 = 0.d0 
            dLt_der_out = 0.d0
            dLt_dw_00 = 0.d0
            dLt_dw_out = 0.d0
         else
            POM3 = (s% RSP_w(k-1)**3 - s% RSP_w(k)**3)/s% dm_bar(k)
            POM = - 2.d0/3.d0*ALFA*ALFAT*(P4*(s% r(k)**2))**2
            rho2_face = 0.5d0*(1.d0/s% Vol(k)**2 + 1.d0/s% Vol(k-1)**2)
            POM2 = s% Hp_face(k)*rho2_face
            Lt_00 = POM*POM2*POM3       
                 
            TEM1 = Lt_00/s% Hp_face(k)
            TEM2 = Lt_00/POM2*s% Hp_face(k)
            
            d_POM2_dVol_00 = dHp_dVol_00(i)*rho2_face - s% Hp_face(k)/s% Vol(k)**3
            d_POM2_dVol_out = dHp_dVol_out(i)*rho2_face - s% Hp_face(k)/s% Vol(k-1)**3
            dLt_dVol_00 = POM*POM3*d_POM2_dVol_00
            dLt_dVol_out = POM*POM3*d_POM2_dVol_out
            
            dLt_dr_00 = 4.d0*Lt_00/s% r(k) & ! 
               - TEM2/s% Vol(k)**3*dVol_dr_00(I) &
               - TEM2/s% Vol(k-1)**3*dVol_dr_in(i+1) &
               + TEM1*dHp_dr_00(I)
            dLt_dr_in = &
               - TEM2/s% Vol(k)**3*dVol_dr_in(I) & ! 
               + TEM1*dHp_dr_in(I)
            dLt_dr_out = &
               - TEM2/s% Vol(k-1)**3*dVol_dr_00(i+1) & ! 
               + TEM1*dHp_dr_out(I)
            dLt_dT_00 = TEM1*dHp_dT_00(I) ! 
            dLt_dT_out = TEM1*dHp_dT_out(I) ! 
            dLt_der_00 = TEM1*dHp_der_00(I) ! 
            dLt_der_out = TEM1*dHp_der_out(I) ! 
            TEM1 = POM*POM2*3.0d0/s% dm_bar(k) 
            dLt_dw_00 = -TEM1*s% RSP_w(k)**2 ! 
            dLt_dw_out = TEM1*s% RSP_w(k-1)**2 ! 
         end if
         if (i > 0) s% Lt(k) = Lt_00
         
         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = Lt_00
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dLt_dVol_out
            write(*,*) 'calc_Lt', s% solver_test_partials_var
         end if
      end subroutine calc_Lt
      
      
      subroutine calc_Lc(s,i,Lc_00, &
            dLc_dr_in, dLc_dr_00, dLc_dr_out, &
            dLc_dVol_00, dLc_dVol_out, &
            dLc_dT_00, dLc_dT_out, &
            dLc_der_00, dLc_der_out, &
            dLc_dw_00, dLc_dw_out)
         type (star_info), pointer :: s
         integer, intent(in) :: I         
         real(dp), intent(out) :: &
            Lc_00, dLc_dr_in, dLc_dr_00, dLc_dr_out, &
            dLc_dVol_00, dLc_dVol_out, &
            dLc_dT_00, dLc_dT_out, &
            dLc_der_00, dLc_der_out, &
            dLc_dw_00, dLc_dw_out
         real(dp) :: POM,POM2,POM3
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (I <= IBOTOM .or. I == NZN .or. ALFA == 0d0)then 
            Lc_00 = 0.d0 !DODANE SMOLEC
            dLc_dr_in = 0.d0
            dLc_dr_00 = 0.d0
            dLc_dr_out = 0.d0
            dLc_dVol_00 = 0.d0
            dLc_dVol_out = 0.d0
            dLc_dT_00 = 0.d0
            dLc_dT_out = 0.d0
            dLc_der_00 = 0.d0
            dLc_der_out = 0.d0
            dLc_dw_00 = 0.d0
            dLc_dw_out = 0.d0! - 
         else if (s% RSP_w(k) < EFL0*1d-8)then
            Lc_00 = 0.d0
            dLc_dr_00 = 0.d0
            dLc_dr_in = 0.d0
            dLc_dr_out = 0.d0
            dLc_dVol_00 = 0.d0
            dLc_dVol_out = 0.d0
            dLc_dT_00 = 0.d0
            dLc_dT_out = 0.d0
            dLc_der_00 = 0.d0
            dLc_der_out = 0.d0
            dLc_dw_00 = 0.d0
            dLc_dw_out = 0.d0
         else
            
            POM3 = 0.5d0*(s% RSP_w(k) + s% RSP_w(k-1))

            POM = P4*(s% r(k)**2)*(ALFAC/ALFAS)* &
               0.5d0*(s% T(k)/s% Vol(k) + s% T(k-1)/s% Vol(k-1))
            Lc_00 = POM*s% PII(k)*POM3 
            
            dLc_dw_00 = POM*s% PII(k)*0.5d0 ! 
            if (I >= NZN - 1) then
               dLc_dw_out = 0.d0
            else
               dLc_dw_out = POM*s% PII(k)*0.5d0 ! 
            end if
            
            POM2 = P4*(s% r(k)**2)*s% PII(k)*POM3*0.5d0*(ALFAC/ALFAS)
            POM = POM*POM3
            
            dLc_dr_00 = & ! 
                 POM*dPII_dr_00(I) &
               + 2.d0*Lc_00/s% r(k) &
               - POM2*(s% T(k)/(s% Vol(k)**2)*dVol_dr_00(I) + &
                        s% T(k-1)/(s% Vol(k-1)**2)*dVol_dr_in(i+1))
            dLc_dr_in = &
                 POM*dPII_dr_in(I) &
               - POM2*s% T(k)/(s% Vol(k)**2)*dVol_dr_in(I) ! 
            dLc_dr_out = &
                 POM*dPII_dr_out(I) &
               - POM2*s% T(k-1)/(s% Vol(k-1)**2)*dVol_dr_00(i+1) ! 

            dLc_dVol_00 = &
                 POM*dPII_dVol_00(I) &
               - POM2*s% T(k)/(s% Vol(k)**2) ! 
            dLc_dVol_out = &
                 POM*dPII_dVol_out(I) &
               - POM2*s% T(k-1)/(s% Vol(k-1)**2) ! 
               
            dLc_dT_00 = POM*dPII_dT_00(I) + POM2/s% Vol(k) ! 
            dLc_dT_out = POM*dPII_dT_out(I) + POM2/s% Vol(k-1) ! 
               
            dLc_der_00 = POM*dPII_der_00(I) ! 
            dLc_der_out = POM*dPII_der_out(I) ! 
            
         end if
         
         if (i > 0) s% Lc(k) = Lc_00

         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = Lc_00
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dLc_dVol_out
            write(*,*) 'calc_Lc', s% solver_test_partials_var
         end if
      end subroutine calc_Lc
      
      ! in diffusion limit, radiative flux equation reduces to Fr calculated from d_erad_dm as below.
      ! note that can have nonequilibrium diffusion regime with different T for gas and photons.
      ! this happens when absorption mean opacity is different than Planck mean opacity.
      subroutine calc_Fr(s, i, Fr_00, & !rs Stellingwerf 1975, Appendix A 
            dFr_dr_out, dFr_dr_00, dFr_dr_in, &
            dFr_dVol_out, dFr_dVol_00, &
            dFr_dT_out, dFr_dT_00, &
            dFr_der_out, dFr_der_00)
         type (star_info), pointer :: s
         integer, intent(in) :: I         
         real(dp), intent(out) :: &
            Fr_00, dFr_dr_out, dFr_dr_00, dFr_dr_in, &
            dFr_dVol_out, dFr_dVol_00, &
            dFr_dT_out, dFr_dT_00, &
            dFr_der_out, dFr_der_00
         real(dp) :: &
            Vol_00, W_00, d_W_00_dr_in, d_W_00_dr_00, d_W_00_dVol_00, d_W_00_der_00, &
            Vol_out, W_out, d_W_out_dr_00, d_W_out_dr_out, d_W_out_dVol_out, d_W_out_der_out, &
            Prad_factor, Fr2a, d_Fr2a_dW_00, d_Fr2a_dW_out, &
            Fr2b, d_Fr2b_dW_00, d_Fr2b_dW_out, &
            BW, kap_00, kap_out, BK, Fr1, Fr2, Fr3, &
            d_Fr_dK_00, d_Fr_dK_out, d_Fr_dW_00, d_Fr_dW_out
         integer :: k
         logical :: test_partials
         include 'formats'
         
         k = NZN+1-i
         
         if (i < 1) then
            if (s% RSP_hydro_only) then
               Fr_00 = 0d0
            else
               Fr_00 = s% L_center
            end if
            Fr_00 = Fr_00/(4d0*pi*s% r_center**2)
            dFr_dr_out = 0
            dFr_dr_00 = 0 
            dFr_dr_in = 0
            dFr_dVol_out = 0
            dFr_dVol_00 = 0
            dFr_dT_out = 0
            dFr_dT_00 = 0
            dFr_der_out = 0
            dFr_der_00 = 0
            return
         end if
         
         Prad_factor = 3d0/crad ! 3d0 to cancel the 1/3d0 factor in CL below
         W_00 = Prad_factor*s% Prad(k) ! replaces s% T(k)**4
         d_W_00_dVol_00 = Prad_factor*d_Pr_dVol(i)
         d_W_00_dr_in = Prad_factor*d_Pr_dr_in(i)
         d_W_00_dr_00 = Prad_factor*d_Pr_dr_00(i)
         d_W_00_der_00 = Prad_factor*d_Pr_der(i) 
         
         if (k == 1) then ! surface
            Fr1 = s% g_Edd*4d0*SIG
            Fr_00 = Fr1*W_00 ! s% T(k)**4 => W_00
            dFr_dr_out = 0
            dFr_dr_in = Fr1*d_W_00_dr_in
            dFr_dr_00 = Fr1*d_W_00_dr_00
            dFr_dVol_out = 0
            dFr_dVol_00 = Fr1*d_W_00_dVol_00
            dFr_dT_out = 0
            dFr_dT_00 = 0
            dFr_der_out = 0
            dFr_der_00 = Fr1*d_W_00_der_00
            return
         end if
         
         W_out = Prad_factor*s% Prad(k-1) ! replaces s% T(k-1)**4
         d_W_out_dVol_out = Prad_factor*d_Pr_dVol(i+1)
         d_W_out_dr_00 = Prad_factor*d_Pr_dr_in(i+1)
         d_W_out_dr_out = Prad_factor*d_Pr_dr_00(i+1)
         d_W_out_der_out = Prad_factor*d_Pr_der(i+1)
         
         BW = log(W_out/W_00) 
         if (abs(BW) < 1d-30) then
            Fr_00 = 0
            dFr_dr_out = 0
            dFr_dr_in = 0
            dFr_dr_00 = 0
            dFr_dVol_out = 0
            dFr_dVol_00 = 0
            dFr_dT_out = 0
            dFr_dT_00 = 0
            dFr_der_out = 0
            dFr_der_00 = 0
            return
         end if
         
         kap_00 = s% opacity(k)
         kap_out = s% opacity(k-1)
         BK = log(kap_out/kap_00)

         Fr1 = -CL*s% r(k)**2/(4d0*pi*s% dm_bar(k))   ! CL = 4d0*(4d0*PI)**2*SIG/3d0
         
         Fr2a = W_out/kap_out - W_00/kap_00
         d_Fr2a_dW_00 = -1d0/kap_00
         d_Fr2a_dW_out = 1d0/kap_out
         
         Fr2b = 1d0 - BK/BW
         d_Fr2b_dW_00 = -BK/BW**2/W_00
         d_Fr2b_dW_out = BK/BW**2/W_out
         
         Fr2 = Fr2a/Fr2b
         Fr_00 = Fr1*Fr2
         d_Fr_dW_00 = Fr1*(d_Fr2a_dW_00 - Fr2*d_Fr2b_dW_00)/Fr2b
         d_Fr_dW_out = Fr1*(d_Fr2a_dW_out - Fr2*d_Fr2b_dW_out)/Fr2b

         Fr3 = Fr1/(BW - BK)
         d_Fr_dK_00 = (Fr3/kap_00)*(W_00*BW/kap_00 - Fr2)
         d_Fr_dK_out = -(Fr3/kap_out)*(W_out*BW/kap_out - Fr2)
         
         dFr_dr_in = & ! 
            + d_Fr_dK_00*dK_dr_in(i) &
            + d_Fr_dW_00*d_W_00_dr_in
         dFr_dr_00 = 2d0*Fr_00/s% r(k) & ! 
            + d_Fr_dK_00*dK_dr_00(i) &
            + d_Fr_dK_out*dK_dr_in(i+1) &
            + d_Fr_dW_00*d_W_00_dr_00 &
            + d_Fr_dW_out*d_W_out_dr_00
         dFr_dr_out = & ! 
            + d_Fr_dK_out*dK_dr_00(i+1) &
            + d_Fr_dW_out*d_W_out_dr_out
         
         dFr_dVol_00 = &
            + d_Fr_dK_00*dK_dVol(i) &
            + d_Fr_dW_00*d_W_00_dVol_00
         dFr_dVol_out = &
            + d_Fr_dK_out*dK_dVol(i+1) &
            + d_Fr_dW_out*d_W_out_dVol_out

         dFr_dT_out = & ! 
            + d_Fr_dK_out*dK_dT(i+1)
         dFr_dT_00 = & ! 
            + d_Fr_dK_00*dK_dT(i)
            
         dFr_der_out = & ! 
            + d_Fr_dW_out*d_W_out_der_out
         dFr_der_00 = & ! 
            + d_Fr_dW_00*d_W_00_der_00
         
         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = Fr_00
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = dFr_dVol_out
            write(*,*) 'calc_Fr', s% solver_test_partials_var
         end if
         
         if (call_is_bad) then
            if (is_bad(Fr_00)) then
   !$OMP critical
            write(*,2) 'Fr_00', k, Fr_00
            write(*,2) 'Fr1', k, Fr1
            write(*,2) 'Fr2', k, Fr2
            write(*,2) 'Fr2a', k, Fr2a
            write(*,2) 'Fr2b', k, Fr2b
            write(*,2) 'W_00', k, W_00
            write(*,2) 'W_out', k, W_out
            write(*,2) 'kap_00', k, kap_00
            write(*,2) 'kap_out', k, kap_out
            write(*,2) 'r(k)', k, s% r(k)
            write(*,2) 'dm_bar(k)', k, s% dm_bar(k)
            write(*,2) 'erad(k)', k, s% erad(k)
            write(*,2) 'erad(k-1)', k-1, s% erad(k-1)
            write(*,2) 'BK', k, BK
            write(*,2) 'BW', k, BW
            write(*,2) 'nz', s% nz
            stop 'calc_Fr'
   !$OMP end critical
            end if
         end if
         
      end subroutine calc_Fr


      subroutine rsp_calc_XP(s, P_surf, i, with_Prad, & ! time weighted combined pressure
            XP, d_XP_dVol_00, d_XP_dT_00, d_XP_der_00, &
            d_XP_dw_00, d_XP_dr_in, d_XP_dr_00)
         type (star_info), pointer :: s
         real(dp), intent(in) :: P_surf
         integer, intent(in) :: i
         logical, intent(in) :: with_Prad
         real(dp), intent(out) :: &
            XP, d_XP_dVol_00, d_XP_dT_00, d_XP_der_00, &
            d_XP_dw_00, d_XP_dr_in, d_XP_dr_00
         real(dp) :: T_surf, Prad_factor
         logical :: test_partials
         integer :: k
         include 'formats'
         
         k = NZN+1-i
         if (k == 0) then ! pressure outside of surface
            if (s% RSP_use_atm_grey_with_kap_for_Psurf) then
               XP = P_surf
            else if (s% RSP_use_Prad_for_Psurf) then
               T_surf = s% T_start(1)
               XP = crad*T_surf**4/3d0
            else
               XP = 0.0d0
            end if
            d_XP_dVol_00 = 0d0
            d_XP_dT_00 = 0d0
            d_XP_der_00 = 0d0
            d_XP_dw_00 = 0d0
            d_XP_dr_in = 0d0
            d_XP_dr_00 = 0d0
            return
         end if
         if (with_Prad) then
            Prad_factor = 1d0
         else
            Prad_factor = 0d0
         end if
         
         XP = THETA*(s% Pgas(k) + Prad_factor*s% Prad(k)) &
            + THETA1*(s% Pgas_start(k) + Prad_factor*s% Prad_start(k)) &
            + THETAQ*s% Pvsc(k) + THETAQ1*s% Pvsc_start(k) &
            + THETAT*s% Ptrb(k) + THETAT1*s% Ptrb_start(k)
         d_XP_dVol_00 = &
              THETA*(d_Pg_dVol(i) + Prad_factor*d_Pr_dVol(i)) &
            + THETAQ*d_Pvsc_dVol(i) &
            + THETAT*dPtrb_dVol_00(i)
         d_XP_dT_00 = THETA*d_Pg_dT(i) + THETAQ*d_Pvsc_dT(i) ! 
         d_XP_der_00 = THETA*Prad_factor*d_Pr_der(i) ! 
         d_XP_dw_00 = THETAT*dPtrb_dw_00(i) ! 
         d_XP_dr_in = & ! 
              THETA*(d_Pg_dr_in(i) + Prad_factor*d_Pr_dr_in(i)) &
            + THETAQ*d_Pvsc_dr_in(i) &
            + THETAT*dPtrb_dr_in(i)
         d_XP_dr_00 = & ! 
              THETA*(d_Pg_dr_00(i) + Prad_factor*d_Pr_dr_00(i)) &
            + THETAQ*d_Pvsc_dr_00(i) &
            + THETAT*dPtrb_dr_00(i)
            
         if (call_is_bad) then
            if (is_bad(XP)) then
   !$OMP critical
               write(*,2) 'XP', k, XP
               write(*,2) 's% Pgas(k)', k, s% Pgas(k)
               write(*,2) 's% Prad(k)', k, s% Prad(k)
               write(*,2) 's% Pvsc(k)', k, s% Pvsc(k)
               write(*,2) 's% Ptrb(k)', k, s% Ptrb(k)
               write(*,2) 'THETA', k, THETA
               write(*,2) 'THETAQ', k, THETAQ
               write(*,2) 'THETAT', k, THETAT
               stop 'rsp_calc_XP'
   !$OMP end critical
            end if
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = XP
            s% solver_test_partials_var = i_var_Vol
            s% solver_test_partials_dval_dx = d_XP_dVol_00
            write(*,*) 'rsp_calc_XP', s% solver_test_partials_var
         end if        
      end subroutine rsp_calc_XP
      
      
      subroutine calc_equations(s,i,P_surf)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: P_surf
         call calc_face_equations(s,i,P_surf)
         call calc_cell_equations(s,i,P_surf,.true.,.true.,.true.)
      end subroutine calc_equations
      
      
      subroutine calc_face_equations(s,i,P_surf)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: P_surf
         call acceleration_eqn(s,i,P_surf)
         call Fr_eqn(s,i)         
      end subroutine calc_face_equations
      
      
      subroutine calc_cell_equations( &
            s, i, P_surf, do_etot, do_eturb, do_erad)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: P_surf
         logical, intent(in) :: do_etot, do_eturb, do_erad
         integer :: k
         real(dp) :: &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00                        

         include 'formats'
         
         ! HR = -residual
         ! partials of residual go in HD            

         k = NZN+1-i

         call get_Lt 
         
         if (do_etot) call total_energy_eqn(s, i, P_surf, &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00)                    
            
         if (do_eturb) call turbulent_energy_eqn(s, i, &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00)                    
         
         if (do_erad) call erad_eqn(s,i)
            
         if (I == 1) call inner_boundary_eqn
         
         contains
         
         subroutine inner_boundary_eqn 
            HR(1) = 0.d0              
            HD(1:LD_HD,1) = 0.d0            
            HD(i_r_dr_00,1) = 1.d0              
         end subroutine inner_boundary_eqn          
                  
         subroutine get_Lt
            integer :: k
            k = NZN+1-i
            call calc_Lt(s,i - 1,Lt_in, &
               dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
               dLt_in_dVol_in, dLt_in_dVol_00, &
               dLt_in_dT_in, dLt_in_dT_00, &
               dLt_in_der_in, dLt_in_der_00, &
               dLt_in_dw_in, dLt_in_dw_00)        
            call calc_Lt(s,i,Lt_00, &
               dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
               dLt_00_dVol_00, dLt_00_dVol_out, &
               dLt_00_dT_00, dLt_00_dT_out, &
               dLt_00_der_00, dLt_00_der_out, &
               dLt_00_dw_00, dLt_00_dw_out)
            if (I == NZN) then
               Lt_00_start = 0.d0   
            else
               Lt_00_start = s% Lt_start(k)
            end if
            if (i == 1) then
               Lt_in_start = 0
            else
               Lt_in_start = s% Lt_start(k+1)
            end if         
         end subroutine get_Lt   

      end subroutine calc_cell_equations
      
      
      subroutine acceleration_eqn(s, i, P_surf)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: P_surf
         integer :: IR, k
         real(dp) :: dt, dm_bar, residual, area, d_area_dr_00, &
            grav, dXP_dm, Uq1, A_dm, R_00, dv_dr, dvdt_factor, &
            XP_00, dXP_00_dVol_00, dXP_00_dT_00, dXP_00_der_00, &
            dXP_00_dw_00, dXP_00_dr_in, dXP_00_dr_00, &
            XP_out, dXP_out_dVol_out, dXP_out_dT_out, dXP_out_der_out, &
            dXP_out_dw_out, dXP_out_dr_00, dXP_out_dr_out, &
            Chi_00, Chi_out, d_Chi_out_dVol_00, &
            d_Chi_out_dr_in, d_Chi_out_dr_00, d_Chi_out_dr_out, d_Chi_out_dr_out2, &
            d_Chi_out_dT_00, d_Chi_out_dT_out, d_Chi_out_dT_out2, &
            d_Chi_out_der_00, d_Chi_out_der_out, d_Chi_out_der_out2, &
            d_Chi_out_dw_out, &
            d_Chi_00_dr_in2, d_Chi_00_dr_in, d_Chi_00_dr_00, d_Chi_00_dr_out, &
            d_Chi_00_dT_in, d_Chi_00_dT_00, d_Chi_00_dT_out, &
            d_Chi_00_der_in, d_Chi_00_der_00, d_Chi_00_der_out, &
            d_Chi_00_dw_00, d_Chi_00_dVol_00, &
            d_Uq1_dr_00, d_Uq_dVol_00, d_Uq_dr_in2, d_Uq_dw_00, d_Uq_dw_out, &
            d_Uq_dT_in, d_Uq_dT_00, d_Uq_dT_out, d_Uq_dT_out2, &
            d_Uq_der_in, d_Uq_der_00, d_Uq_der_out, d_Uq_der_out2, &
            d_Uq_dr_in, d_Uq_dr_00, d_Uq_dr_out, d_Uq_dr_out2, &
            kap_face, d_kap_dVol_00, d_kap_dT_00, d_kap_dT_out, &
            d_kap_dr_in, d_kap_dr_00, d_kap_dr_out, &
            Fr_tw, Fr_term, d_Fr_term_dFr_00, d_Fr_term_dVol_00, d_Fr_term_dT_00, d_Fr_term_dT_out, &
            d_Fr_term_dr_in, d_Fr_term_dr_00, d_Fr_term_dr_out
         logical :: test_partials, use_Prad

         include 'formats'
         
         use_Prad = .true. ! s% RSP_accel_eqn_use_Prad_instead_of_Fr_term
         
         dvdt_factor = 1d0

         k = NZN+1-i
         IR = i_var_R + NV*(i-1)
         HD(1:LD_HD,IR) = 0.d0   
         
         dt = s% dt
         if (s% use_compression_outer_BC .and. I == NZN) then
            stop 'no rsp support for use_compression_outer_BC'
         end if
         
         ! XP doesn't include Prad for acceleration equation
         ! instead introduce term using Fr
         
         call rsp_calc_XP(s, P_surf, i+1, use_Prad, &
            XP_out, dXP_out_dVol_out, dXP_out_dT_out, dXP_out_der_out, &
            dXP_out_dw_out, dXP_out_dr_00, dXP_out_dr_out)
            
         call rsp_calc_XP(s, P_surf, i, use_Prad, &
            XP_00, dXP_00_dVol_00, dXP_00_dT_00, dXP_00_der_00, &
            dXP_00_dw_00, dXP_00_dr_in, dXP_00_dr_00)

         area = P43*(s% r(k)**2 + s% r(k)*s% r_start(k) + s% r_start(k)**2)
         d_area_dr_00 = P43*(2d0*s% r(k) + s% r_start(k))
         dm_bar = s% dm_bar(k)
         A_dm = area/dm_bar
         dv_dr = 2d0/dt
         dXP_dm = (XP_out - XP_00)/dm_bar
         grav = - s% cgrav(k)*s% m(k)/(s% r(k)*s% r_start(k))
        
         R_00 = 0.5d0*(s% r(k) + s% r_start(k))
         Uq1 = P4/(dm_bar*R_00)
         d_Uq1_dr_00 = -Uq1*0.5d0/R_00

         Chi_00 = THETAU*s% Chi(k) + THETAU1*s% Chi_start(k)
         d_Chi_00_dVol_00 = THETAU*dChi_dVol_00(i)
         d_Chi_00_dT_in = THETAU*dChi_dT_in(i)
         d_Chi_00_dT_00 = THETAU*dChi_dT_00(i)
         d_Chi_00_dT_out = THETAU*dChi_dT_out(i)
         d_Chi_00_der_in = THETAU*dChi_der_in(i)
         d_Chi_00_der_00 = THETAU*dChi_der_00(i)
         d_Chi_00_der_out = THETAU*dChi_der_out(i)
         d_Chi_00_dw_00 = THETAU*dChi_dw_00(i)
         d_Chi_00_dr_in2 = THETAU*dChi_dr_in2(i)
         d_Chi_00_dr_in = THETAU*dChi_dr_in(i)
         d_Chi_00_dr_00 = THETAU*dChi_dr_00(i)
         d_Chi_00_dr_out = THETAU*dChi_dr_out(i)   
         if (I == NZN) then
            Chi_out = 0d0
            d_Chi_out_dVol_00 = 0d0
            d_Chi_out_dT_00 = 0d0
            d_Chi_out_dT_out = 0d0
            d_Chi_out_dT_out2 = 0d0
            d_Chi_out_der_00 = 0d0
            d_Chi_out_der_out = 0d0
            d_Chi_out_der_out2 = 0d0
            d_Chi_out_dw_out = 0d0
            d_Chi_out_dr_in = 0d0
            d_Chi_out_dr_00 = 0d0
            d_Chi_out_dr_out = 0d0
            d_Chi_out_dr_out2 = 0d0
         else
            Chi_out = THETAU*s% Chi(k-1) + THETAU1*s% Chi_start(k-1)         
            d_Chi_out_dVol_00 = THETAU*dChi_dVol_in(i+1)
            d_Chi_out_dT_00 = THETAU*dChi_dT_in(i+1)
            d_Chi_out_dT_out = THETAU*dChi_dT_00(i+1)
            d_Chi_out_dT_out2 = THETAU*dChi_dT_out(i+1)
            d_Chi_out_der_00 = THETAU*dChi_der_in(i+1)
            d_Chi_out_der_out = THETAU*dChi_der_00(i+1)
            d_Chi_out_der_out2 = THETAU*dChi_der_out(i+1)
            d_Chi_out_dw_out = THETAU*dChi_dw_00(i+1)
            d_Chi_out_dr_in = THETAU*dChi_dr_in2(i+1)
            d_Chi_out_dr_00 = THETAU*dChi_dr_in(i+1)
            d_Chi_out_dr_out = THETAU*dChi_dr_00(i+1)
            d_Chi_out_dr_out2 = THETAU*dChi_dr_out(i+1)   
         end if                          

         s% Uq(k) = Uq1*(Chi_out - Chi_00)         
         d_Uq_dVol_00 = Uq1*(d_Chi_out_dVol_00 - d_Chi_00_dVol_00)
         d_Uq_dT_in = -Uq1*d_Chi_00_dT_in
         d_Uq_dT_00 = Uq1*(d_Chi_out_dT_00 - d_Chi_00_dT_00)
         d_Uq_dT_out = Uq1*(d_Chi_out_dT_out - d_Chi_00_dT_out)
         d_Uq_dT_out2 = Uq1*d_Chi_out_dT_out2
         d_Uq_der_in = -Uq1*d_Chi_00_der_in
         d_Uq_der_00 = Uq1*(d_Chi_out_der_00 - d_Chi_00_der_00)
         d_Uq_der_out = Uq1*(d_Chi_out_der_out - d_Chi_00_der_out)
         d_Uq_der_out2 = Uq1*d_Chi_out_der_out2         
         d_Uq_dw_00 = -Uq1*d_Chi_00_dw_00
         d_Uq_dw_out = Uq1*d_Chi_out_dw_out         
         d_Uq_dr_in2 = -Uq1*d_Chi_00_dr_in2
         d_Uq_dr_in = Uq1*(d_Chi_out_dr_in - d_Chi_00_dr_in)
         d_Uq_dr_00 = &
              Uq1*(d_Chi_out_dr_00 - d_Chi_00_dr_00) &
            + d_Uq1_dr_00*(Chi_out - Chi_00)
         d_Uq_dr_out = Uq1*(d_Chi_out_dr_out - d_Chi_00_dr_out)
         d_Uq_dr_out2 = Uq1*d_Chi_out_dr_out2

         if (use_Prad) then
            Fr_term = 0
            d_Fr_term_dFr_00 = 0
            d_Fr_term_dVol_00 = 0
            d_Fr_term_dT_00 = 0
            d_Fr_term_dT_out = 0
            d_Fr_term_dr_in = 0
            d_Fr_term_dr_00 = 0
            d_Fr_term_dr_out = 0
         else ! include radiative force, Fr*kap_face/clight
            if (k == 1) then
               kap_face = s% opacity(k)
               d_kap_dVol_00 = dK_dVol(i)
               d_kap_dT_00 = dK_dT(i)
               d_kap_dT_out = 0d0
               d_kap_dr_in = dK_dr_in(i)
               d_kap_dr_00 = dK_dr_00(i)
               d_kap_dr_out = 0d0
            else
               kap_face = 0.5d0*(s% opacity(k) + s% opacity(k-1))
               d_kap_dVol_00 = 0.5d0*dK_dVol(i)
               d_kap_dT_00 = 0.5d0*dK_dT(i)
               d_kap_dT_out = 0.5d0*dK_dT(i+1)
               d_kap_dr_in = 0.5d0*dK_dr_in(i)
               d_kap_dr_00 = 0.5d0*(dK_dr_00(i) + dK_dr_in(i+1))
               d_kap_dr_out = 0.5d0*dK_dr_00(i+1)
            end if
            Fr_tw = WTR*s% Fr(k) + WTR1*s% Fr_start(k)
            Fr_term = Fr_tw*kap_face/clight
            d_Fr_term_dFr_00 = WTR*kap_face/clight
            d_Fr_term_dVol_00 = Fr_term*d_kap_dVol_00/kap_face
            d_Fr_term_dT_00 = Fr_term*d_kap_dT_00/kap_face
            d_Fr_term_dT_out = Fr_term*d_kap_dT_out/kap_face
            d_Fr_term_dr_in = Fr_term*d_kap_dr_in/kap_face
            d_Fr_term_dr_00 = Fr_term*d_kap_dr_00/kap_face
            d_Fr_term_dr_out = Fr_term*d_kap_dr_out/kap_face
         end if
         
         residual = &
            dvdt_factor*(s% v(k) - s% v_start(k))/dt &
           + area*dXP_dm - grav - s% Uq(k) - Fr_term
         HR(IR) = -residual
         
         !s% xtra1_array(k) = s% Pgas(k) + s% Prad(k)
         !s% xtra2_array(k) = s% Vol(k)
         !s% xtra3_array(k) = s% T(k)
         !s% xtra4_array(k) = s% v(k)
         !s% xtra5_array(k) = s% RSP_w(k)**2
         !s% xtra6_array(k) = s% r(k)
            
         HD(i_r_dFr_00,IR) = - d_Fr_term_dFr_00

         HD(i_r_dT_in,IR) = &
            - d_Uq_dT_in
         HD(i_r_dT_00,IR) = &
            + A_dm*(-dXP_00_dT_00) & 
            - d_Uq_dT_00 &
            - d_Fr_term_dT_00
         HD(i_r_dT_out,IR) = &
            + A_dm*dXP_out_dT_out & 
            - d_Uq_dT_out &
            - d_Fr_term_dT_out
         HD(i_r_dT_out2,IR) = & 
            - d_Uq_dT_out2
                  
         HD(i_r_dr_in2,IR) = - d_Uq_dr_in2
         HD(i_r_dr_in,IR) = & ! 
            + A_dm*(-dXP_00_dr_in) &
            - d_Uq_dr_in &
            - d_Fr_term_dr_in
         HD(i_r_dr_00,IR) = &
            dvdt_factor*dv_dr/dt &
            + d_area_dr_00*dXP_dm &
            + A_dm*(dXP_out_dr_00 - dXP_00_dr_00) &
            + grav/s% r(k) &
            - d_Uq_dr_00 &
            - d_Fr_term_dr_00
         HD(i_r_dr_00,IR) = &
            HD(i_r_dr_00,IR) + 0.5d0*Uq1/R_00*(Chi_out - Chi_00)
         HD(i_r_dr_out,IR) = &
            + A_dm*dXP_out_dr_out &
            - d_Uq_dr_out &
            - d_Fr_term_dr_out
         HD(i_r_dr_out2,IR) = & ! 
            - d_Uq_dr_out2

         HD(i_r_der_in,IR) = &
            - d_Uq_der_in
         HD(i_r_der_00,IR) = &
            + A_dm*(-dXP_00_der_00) & 
            - d_Uq_der_00
         HD(i_r_der_out,IR) = &
            + A_dm*dXP_out_der_out & 
            - d_Uq_der_out
         HD(i_r_der_out2,IR) = & 
            - d_Uq_der_out2

         HD(i_r_dw_in,IR) = 0.d0      
         if (I <= IBOTOM .or. I == NZN) then
            HD(i_r_dw_00,IR) = 0.d0   
         else     
            HD(i_r_dw_00,IR) = &
               + A_dm*(-dXP_00_dw_00) & 
               - d_Uq_dw_00
         end if
         if (I <= IBOTOM - 1 .or. I >= NZN - 1) then
            HD(i_r_dw_out,IR) = 0.d0  
         else
            HD(i_r_dw_out,IR) = &
               + A_dm*dXP_out_dw_out &
               - d_Uq_dw_out
         end if
         HD(i_r_dw_out2,IR) = 0.0d0   
         
         !call check_is_bad
         

         !HD(i_r_dr_in2,IR) !  
         !HD(i_r_dr_in,IR) ! ok
         !HD(i_r_dr_00,IR) !  ok
         !HD(i_r_dr_out,IR) ! ok
         !HD(i_r_dr_out2,IR) ! 
         
         !HD(i_r_dT_in,IR) ! 0
         !HD(i_r_dT_00,IR) ! ok
         !HD(i_r_dT_out,IR) ! need 1d-3 like for er
         !HD(i_r_dT_out2,IR) ! 0
         
         ! NOTE: may need solver_test_partials_dx_0 = 1d-3 for er
         !HD(i_r_der_in,IR) ! 
         !HD(i_r_der_00,IR) ! 0
         !HD(i_r_der_out,IR) ! 0
         !HD(i_r_der_out2,IR) ! 
         
         !HD(i_r_dw_00,IR) ! 
         !HD(i_r_dw_out,IR) ! 

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_r
            s% solver_test_partials_dval_dx = i_r_dr_00
            write(*,*) 'acceleration_eqn', s% solver_test_partials_var
         end if    
         
         contains
         
         subroutine check_is_bad
            include 'formats'
            if (is_bad(residual)) then
            !$OMP critical
               write(*,2) 'residual', k, residual
               write(*,2) 's% v(k)', k, s% v(k)
               write(*,2) 's% v_start(k)', k, s% v_start(k)
               write(*,2) 'area', k, area
               write(*,2) 'dXP_dm', k, dXP_dm
               write(*,2) 'XP_out', k, XP_out
               write(*,2) 'XP_00', k, XP_00
               write(*,2) 's% dm_bar(k)', k, s% dm_bar(k)
               write(*,2) 'grav', k, grav
               write(*,2) 's% Uq(k)', k, s% Uq(k)
               write(*,2) 'Fr_term', k, Fr_term
               write(*,2) 'dt', k, dt
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dr_in2,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dr_in2,IR)', k, HD(i_r_dr_in2,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dr_in,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dr_in,IR)', k, HD(i_r_dr_in,IR)
               write(*,2) 'dXP_00_dr_in', k, dXP_00_dr_in
               write(*,2) 'd_Uq_dr_in', k, d_Uq_dr_in
               write(*,2) 'd_Fr_term_dr_in', k, d_Fr_term_dr_in
               write(*,2) 'd_Chi_out_dr_in', k, d_Chi_out_dr_in
               write(*,2) 'd_Chi_00_dr_in', k, d_Chi_00_dr_in
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dr_00,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dr_00,IR)', k, HD(i_r_dr_00,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dr_out,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dr_out,IR)', k, HD(i_r_dr_out,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dr_out2,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dr_out2,IR)', k, HD(i_r_dr_out2,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dT_in,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dT_in,IR)', k, HD(i_r_dT_in,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dT_00,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dT_00,IR)', k, HD(i_r_dT_00,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dT_out,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dT_out,IR)', k, HD(i_r_dT_out,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dT_out2,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dT_out2,IR)', k, HD(i_r_dT_out2,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_der_in,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_der_in,IR)', k, HD(i_r_der_in,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_der_00,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_der_00,IR)', k, HD(i_r_der_00,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_der_out,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_der_out,IR)', k, HD(i_r_der_out,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_der_out2,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_der_out2,IR)', k, HD(i_r_der_out2,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dw_00,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dw_00,IR)', k, HD(i_r_dw_00,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if

            if (is_bad(HD(i_r_dw_out,IR))) then
            !$OMP critical
               write(*,2) 'HD(i_r_dw_out,IR)', k, HD(i_r_dw_out,IR)
               stop 'acceleration_eqn'
            !$OMP end critical
            end if
         end subroutine check_is_bad

      end subroutine acceleration_eqn
               
      
      subroutine total_energy_eqn(s, i, P_surf, &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00)                    
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: P_surf, &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00             
         integer :: IT, k
         real(dp) :: dt, dm, residual, erad, erad_tw, DV, dt_div_dm, &
            area_00, area_00_start, area_in, area_in_start, &
            L_00, Lr_00, Lr_00_start, d_Lr_00_dFr_00, d_Lr_00_dr_00, &
            L_in, Lr_in, Lr_in_start, d_Lr_in_dFr_in, d_Lr_in_dr_in, &
            Lc_00, Lc_00_start, Lc_in, Lc_in_start, &
            dLc_in_dr_in, dLc_in_dr_in2, dLc_in_dr_00, &
            dLc_in_dVol_in, dLc_in_dVol_00, &
            dLc_in_dT_in, dLc_in_dT_00, &
            dLc_in_der_in, dLc_in_der_00, &
            dLc_in_dw_in, dLc_in_dw_00, &    
            dLc_00_dr_00, dLc_00_dr_in, dLc_00_dr_out, &
            dLc_00_dVol_00, dLc_00_dVol_out, &
            dLc_00_dT_00, dLc_00_dT_out, &
            dLc_00_der_00, dLc_00_der_out, &
            dLc_00_dw_00, dLc_00_dw_out, &        
            XP_00, dXP_00_dr_00, dXP_00_dr_in, &
            dXP_00_dVol_00, dXP_00_dT_00, dXP_00_der_00, dXP_00_dw_00, &
            u_div_r, d_u_div_r_dr_00, d_u_div_r_dr_in, u_div_r_factor
         logical :: test_partials

         include 'formats'

         k = NZN+1-i

         IT = i_var_T + NV*(i-1)
         HD(1:LD_HD,IT) = 0.d0

         call calc_Lc(s, i-1, Lc_in, &
            dLc_in_dr_in2, dLc_in_dr_in, dLc_in_dr_00, &
            dLc_in_dVol_in, dLc_in_dVol_00, &
            dLc_in_dT_in, dLc_in_dT_00, &
            dLc_in_der_in, dLc_in_der_00, &
            dLc_in_dw_in, dLc_in_dw_00)    
         call calc_Lc(s, i, Lc_00, &
            dLc_00_dr_in, dLc_00_dr_00, dLc_00_dr_out, &
            dLc_00_dVol_00, dLc_00_dVol_out, &
            dLc_00_dT_00, dLc_00_dT_out, &
            dLc_00_der_00, dLc_00_der_out, &
            dLc_00_dw_00, dLc_00_dw_out)        

         if (I == NZN) then
            Lc_00_start = 0.d0    
         else
            Lc_00_start = s% Lc_start(k)
         end if
         if (i == 1) then
            Lc_in_start = 0
         else
            Lc_in_start = s% Lc_start(k+1)  
         end if         
         
         area_00 = 4d0*pi*s% r(k)**2
         Lr_00 = s% Fr(k)*area_00
         d_Lr_00_dFr_00 = area_00    
         d_Lr_00_dr_00 = 2d0*Lr_00/s% r(k)
         area_00_start = 4d0*pi*s% r_start(k)**2
         Lr_00_start = s% Fr_start(k)*area_00_start
         if (i == 1) then
            if (s% RSP_hydro_only) then
               Lr_in = 0d0
            else
               Lr_in = s% L_center
            end if
            d_Lr_in_dFr_in = 0d0      
            d_Lr_in_dr_in = 0d0      
            Lr_in_start = Lr_in 
         else
            area_in = 4d0*pi*s% r(k+1)**2
            Lr_in = s% Fr(k+1)*area_in
            d_Lr_in_dFr_in = area_in       
            d_Lr_in_dr_in = 2d0*Lr_in/s% r(k+1)
            area_in_start = 4d0*pi*s% r_start(k+1)**2
            Lr_in_start = s% Fr_start(k+1)*area_in_start
         end if         
         
         L_00 = &
            WTR*Lr_00 + WTC*Lc_00 + WTT*Lt_00 + &
            WTR1*Lr_00_start + WTC1*Lc_00_start + WTT1*Lt_00_start
         L_in = &
            WTR*Lr_in + WTC*Lc_in + WTT*Lt_in + &
            WTR1*Lr_in_start + WTC1*Lc_in_start + WTT1*Lt_in_start
         
         dt = s% dt
         dm = s% dm(k)
         dt_div_dm = dt/dm
         DV = s% Vol(k) - s% Vol_start(k)
         
         if (s% f_Edd(k) == f_Edd_isotropic .or. k == NZN) then
            u_div_r = 0d0
            d_u_div_r_dr_00 = 0d0
            d_u_div_r_dr_in = 0d0
            u_div_r_factor = 0d0
         else
            u_div_r = 0.5d0*(s% v(k)/s% r(k) + s% v(k+1)/s% r(k+1))
            d_u_div_r_dr_00 = 0.5d0*(2d0/dt - s% v(k)/s% r(k))/s% r(k)
            d_u_div_r_dr_in = 0.5d0*(2d0/dt - s% v(k+1)/s% r(k+1))/s% r(k+1)
            u_div_r_factor = dt*(1d0 - 3d0*s% f_Edd(k))
         end if
         
         erad = s% erad(k)
         erad_tw = THETAE*erad + THETAE1*s% erad_start(k)
         
         call rsp_calc_XP(s, P_surf, i, .true., &
            XP_00, dXP_00_dVol_00, dXP_00_dT_00, dXP_00_der_00, &
            dXP_00_dw_00, dXP_00_dr_in, dXP_00_dr_00)

         residual = &
              s% egas(k) - s% egas_start(k) &
            + erad - s% erad_start(k) &
            + s% RSP_w(k)**2 - s% RSP_w_start(k)**2 &
            + dt_div_dm*(L_00 - L_in) &
            + XP_00*DV &
            + erad_tw*u_div_r_factor*u_div_r &
            - dt*s% Eq(k)

         s% ergs_error(k) = s% dm(k)*residual

         HR(IT) = -residual

         HD(i_T_dFr_in,IT) = -dt_div_dm*WTR*d_Lr_in_dFr_in
         HD(i_T_dFr_00,IT) =  dt_div_dm*WTR*d_Lr_00_dFr_00

         HD(i_T_dr_in2,IT) = &
            - dt_div_dm*WTC*dLc_in_dr_in2 &
            - dt_div_dm*WTT*dLt_in_dr_in2 &
            - dt*dEq_dr_in2(I)
         HD(i_T_dr_in,IT) = &
            + d_egas_dr_in(i) & 
            - dt_div_dm*WTR*d_Lr_in_dr_in &
            + dt_div_dm*WTC*(dLc_00_dr_in - dLc_in_dr_in) &
            + dt_div_dm*WTT*(dLt_00_dr_in - dLt_in_dr_in) &
            + dVol_dr_in(I)*XP_00 & 
            + DV*dXP_00_dr_in &
            + erad_tw*u_div_r_factor*d_u_div_r_dr_in &
            - dt*dEq_dr_in(I)             
         HD(i_T_dr_00,IT) = &
            + d_egas_dr_00(i) & 
            + dt_div_dm*WTR*d_Lr_00_dr_00 &
            + dt_div_dm*WTC*(dLc_00_dr_00 - dLc_in_dr_00) &
            + dt_div_dm*WTT*(dLt_00_dr_00 - dLt_in_dr_00) &
            + dVol_dr_00(I)*XP_00 & 
            + DV*dXP_00_dr_00 &
            + erad_tw*u_div_r_factor*d_u_div_r_dr_00 &
            - dt*dEq_dr_00(I)
         HD(i_T_dr_out,IT) = &
          + dt_div_dm*WTC*dLc_00_dr_out &
          + dt_div_dm*WTT*dLt_00_dr_out &
          - dt*dEq_dr_out(I)

         HD(i_T_dT_in,IT) = &
            - dt_div_dm*WTC*dLc_in_dT_in &
            - dt_div_dm*WTT*dLt_in_dT_in &
            - dt*dEq_dT_in(I)  
         HD(i_T_dT_00,IT) = &
              d_egas_dT(i) & 
            + DV*dXP_00_dT_00 &     
            + dt_div_dm*WTC*(dLc_00_dT_00 - dLc_in_dT_00) &
            + dt_div_dm*WTT*(dLt_00_dT_00 - dLt_in_dT_00) &
            - dt*dEq_dT_00(I)
         HD(i_T_dT_out,IT) = & ! 
            + dt_div_dm*WTC*dLc_00_dT_out &
            + dt_div_dm*WTT*dLt_00_dT_out &
            - dt*dEq_dT_out(I)                    

         HD(i_T_der_in,IT) = &
            - dt_div_dm*WTC*dLc_in_der_in &
            - dt_div_dm*WTT*dLt_in_der_in &
            - dt*dEq_der_in(I)  
         HD(i_T_der_00,IT) = &
              1d0 & 
            + DV*dXP_00_der_00 &     
            + dt_div_dm*WTC*(dLc_00_der_00 - dLc_in_der_00) &
            + dt_div_dm*WTT*(dLt_00_der_00 - dLt_in_der_00) &
            + THETAE*u_div_r_factor*u_div_r &
            - dt*dEq_der_00(I)
         HD(i_T_der_out,IT) = & ! 
            + dt_div_dm*WTC*dLc_00_der_out &
            + dt_div_dm*WTT*dLt_00_der_out &
            - dt*dEq_der_out(I)                    
    
         if (I <= IBOTOM + 1) then
            HD(i_T_dw_in,IT) = 0.d0   
         else
            HD(i_T_dw_in,IT) = & 
               - dt_div_dm*WTC*dLc_in_dw_in &
               - dt_div_dm*WTT*dLt_in_dw_in      
         end if         
         if (I <= IBOTOM .or. I == NZN) then
            HD(i_T_dw_00,IT) = 0.d0  
         else
            HD(i_T_dw_00,IT) = &
               2.d0*s% RSP_w(k) & 
               + dt_div_dm*WTC*(dLc_00_dw_00 - dLc_in_dw_00) &
               + dt_div_dm*WTT*(dLt_00_dw_00 - dLt_in_dw_00) &
               + DV*dXP_00_dw_00 &
               - dt*dEq_dw_00(I) 
         end if                   
         if (I <= IBOTOM - 1 .or. I >= NZN - 1) then
            HD(i_T_dw_out,IT) = 0.d0
         else
            HD(i_T_dw_out,IT) = &
               dt_div_dm*WTT*dLt_00_dw_out & 
               + dt_div_dm*WTC*dLc_00_dw_out
         end if
         
         if (i == -6 .and. s% model_number == s% max_model_number) then
            write(*,5) 'HD(i_T_dw_00,IT)', k, i, iter, s% model_number, HD(i_T_dw_00,IT)
            write(*,5) 's% RSP_w(k)', k, i, iter, s% model_number, s% RSP_w(k)
            write(*,5) 'dt_div_dm', k, i, iter, s% model_number, dt_div_dm
            write(*,5) 'WTC', k, i, iter, s% model_number, WTC
            write(*,5) 'WTT', k, i, iter, s% model_number, WTT
            write(*,5) 'DV', k, i, iter, s% model_number, DV
            write(*,5) 'dt', k, i, iter, s% model_number, dt
            write(*,5) 'dLc_00_dw_00', k, i, iter, s% model_number, dLc_00_dw_00
            write(*,5) 'dLc_in_dw_00', k, i, iter, s% model_number, dLc_in_dw_00
            write(*,5) 'dLt_00_dw_00', k, i, iter, s% model_number, dLt_00_dw_00
            write(*,5) 'dLt_in_dw_00', k, i, iter, s% model_number, dLt_in_dw_00
            write(*,5) 'dXP_00_dw_00', k, i, iter, s% model_number, dXP_00_dw_00
            write(*,5) 'dEq_dw_00(I)', k, i, iter, s% model_number, dEq_dw_00(I)
         endif

         !HD(i_T_dr_in2,IT) ! 
         !HD(i_T_dr_in,IT) ! ok
         !HD(i_T_dr_00,IT) ! ok
         !HD(i_T_dr_out,IT) ! 

         !HD(i_T_dT_in,IT) ! 0
         !HD(i_T_dT_00,IT) ! ok
         !HD(i_T_dT_out,IT) ! 0

         !HD(i_T_der_in,IT) ! 
         !HD(i_T_der_00,IT) ! ok
         !HD(i_T_der_out,IT) ! 
    
         !HD(i_T_dw_in,IT) ! 
         !HD(i_T_dw_00,IT) ! 
         !HD(i_T_dw_out,IT) ! 

         !test_partials = (k+1 == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_r
            s% solver_test_partials_dval_dx = HD(i_T_dr_00,IT)
            write(*,*) 'total_energy_eqn', s% solver_test_partials_var
         end if        

      end subroutine total_energy_eqn
      
      
      subroutine turbulent_energy_eqn(s, i, &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00)                    
         type (star_info), pointer :: s
         integer, intent(in) :: i
         real(dp), intent(in) :: &
            Lt_00, Lt_00_start, Lt_in, Lt_in_start, &
            dLt_00_dr_00, dLt_00_dr_in, dLt_00_dr_out, &
            dLt_00_dVol_00, dLt_00_dVol_out, &
            dLt_00_dT_00, dLt_00_dT_out, &
            dLt_00_der_00, dLt_00_der_out, &
            dLt_00_dw_00, dLt_00_dw_out, &
            dLt_in_dr_in, dLt_in_dr_in2, dLt_in_dr_00, &
            dLt_in_dVol_in, dLt_in_dVol_00, &
            dLt_in_dT_in, dLt_in_dT_00, &
            dLt_in_der_in, dLt_in_der_00, &
            dLt_in_dw_in, dLt_in_dw_00
         integer :: IW, k
         real(dp) :: dt, dm, residual, Ptrb_tw, DV, dt_div_dm, L_00, L_in
         logical :: test_partials

         include 'formats'

         k = NZN+1-i
      
         IW = i_var_w + NV*(i-1) 
         HD(1:LD_HD,IW) = 0.d0            

         if (ALFA == 0d0 .or. I <= IBOTOM .or. I == NZN) then         
            HD(i_w_dw_00,IW) = 1.d0
            HR(IW) = 0.d0
            return
         end if

         L_00 = WTT*Lt_00 + WTT1*Lt_00_start
         L_in = WTT*Lt_in + WTT1*Lt_in_start
         
         dt = s% dt
         dm = s% dm(k)
         dt_div_dm = dt/dm
         DV = s% Vol(k) - s% Vol_start(k)
         Ptrb_tw = THETAT*s% Ptrb(k) + THETAT1*s% Ptrb_start(k)

         residual = &
            (s% RSP_w(k)**2 - s% RSP_w_start(k)**2) &
            + dt_div_dm*(L_00 - L_in) &
            + DV*Ptrb_tw &
            - dt*(GAM*s% COUPL(k) + GAM1*s% COUPL_start(k) + s% Eq(k))
         HR(IW) = -residual
         
         if (k==-109) then
            write(*,3) 'RSP w dEt PdV dtC dtEq', k, iter, &
               s% RSP_w(k), s% RSP_w(k)**2 - s% RSP_w_start(k)**2, DV*Ptrb_tw, &
               dt*(GAM*s% COUPL(k) + GAM1*s% COUPL_start(k)), dt*s% Eq(k)
            !write(*,2) 'RSP w COUPL SOURCE DAMP DAMPR', k, &
            !   s% RSP_w(k), s% COUPL(k), s% SOURCE(k), s% DAMP(k), s% DAMPR(k)
            !write(*,2) 'RSP w SOURCE PII/Hp P*QQ_div_Cp P T', k, &
            !   s% RSP_w(k), s% SOURCE(k), &
            !   0.5d0*(s% PII(k)/s% Hp_face(k) + s% PII(k+1)/s% Hp_face(k+1)), &
            !   (s% Pgas(k) + s% Prad(k))*s% QQ(k)/s% Cp(k), s% Pgas(k) + s% Prad(k), s% T(k)
            !write(*,2) 'RSP PII_00 PII_p1 Hp_00 Hp_p1', k, &
            !   s% PII(k), s% PII(k+1), s% Hp_face(k), s% Hp_face(k+1)
         end if
         
         HD(i_w_dw_in2,IW) = 0.d0          
         HD(i_w_dw_in,IW) = - dt_div_dm*WTT*dLt_in_dw_in  
         HD(i_w_dw_00,IW) = &
            2.d0*s% RSP_w(k) & 
            - dt*GAM*dC_dw_00(I) &
            - dt*dEq_dw_00(I) &
            + dt_div_dm*WTT*(dLt_00_dw_00 - dLt_in_dw_00) &
            + DV*THETAT*dPtrb_dw_00(I)
         HD(i_w_dw_out,IW) = dt_div_dm*WTT*dLt_00_dw_out  
         HD(i_w_dw_out2,IW) = 0.d0           

         HD(i_w_dr_in2,IW) = &
            - dt*GAM*dC_dr_in2(I) & 
            - dt_div_dm*WTT*dLt_in_dr_in2 &
            - dt*dEq_dr_in2(I)                      
         HD(i_w_dr_in,IW) = &
            - dt*GAM*dC_dr_in(I) &  
            - dt*dEq_dr_in(I) &
            + dt_div_dm*WTT*(dLt_00_dr_in - dLt_in_dr_in) &
            + DV*THETAT*dPtrb_dr_in(I) &
            + (THETAT*s% Ptrb(k) + THETAT1*s% Ptrb_start(k))*dVol_dr_in(I)
         HD(i_w_dr_00,IW) = &
            - dt*GAM*dC_dr_00(I) & 
            - dt*dEq_dr_00(I) &
            + dt_div_dm*WTT*(dLt_00_dr_00 - dLt_in_dr_00) &
            + DV*THETAT*dPtrb_dr_00(I) &
            + (THETAT*s% Ptrb(k) + THETAT1*s% Ptrb_start(k))*dVol_dr_00(I)
         HD(i_w_dr_out,IW) = &
            - dt*GAM*dC_dr_out(I) &   
            + dt_div_dm*WTT*dLt_00_dr_out &
            - dt*dEq_dr_out(I)

         if (I <= IBOTOM + 1) then
            HD(i_w_dT_in,IW) = 0.d0
         else
            HD(i_w_dT_in,IW) = &
               - dt_div_dm*WTT*dLt_in_dT_in &  
               - dt*dEq_dT_in(I) &
               - dt*GAM*dC_dT_in(I)  
         end if
         HD(i_w_dT_00,IW) = &
            - dt*GAM*dC_dT_00(I) & 
            - dt*dEq_dT_00(I) &
            + dt_div_dm*WTT*(dLt_00_dT_00 - dLt_in_dT_00)
         if (I <= IBOTOM - 1 .or. I >= NZN - 1) then 
            HD(i_w_dT_out,IW) = 0.d0
         else
            HD(i_w_dT_out,IW) = &
               - dt*GAM*dC_dT_out(I) &  
               + dt_div_dm*WTT*dLt_00_dT_out &
               - dt*dEq_dT_out(I)
         end if

         if (I <= IBOTOM + 1) then
            HD(i_w_der_in,IW) = 0.d0
         else
            HD(i_w_der_in,IW) = &
               - dt_div_dm*WTT*dLt_in_der_in &  
               - dt*dEq_der_in(I) &
               - dt*GAM*dC_der_in(I)  
         end if
         HD(i_w_der_00,IW) = &
            - dt*GAM*dC_der_00(I) & 
            - dt*dEq_der_00(I) &
            + dt_div_dm*WTT*(dLt_00_der_00 - dLt_in_der_00)
         if (I <= IBOTOM - 1 .or. I >= NZN - 1) then 
            HD(i_w_der_out,IW) = 0.d0
         else
            HD(i_w_der_out,IW) = &
               - dt*GAM*dC_der_out(I) &  
               + dt_div_dm*WTT*dLt_00_der_out &
               - dt*dEq_der_out(I)
         end if

         !HD(i_w_dw_in,IW) ! 
         !HD(i_w_dw_00,IW) ! 
         !HD(i_w_dw_out,IW) ! 
         
         !HD(i_w_dr_in2,IW) !           
         !HD(i_w_dr_in,IW) ! 
         !HD(i_w_dr_00,IW) ! 
         !HD(i_w_dr_out,IW) ! 
         
         !HD(i_w_dT_in,IW) ! 
         !HD(i_w_dT_00,IW) ! 
         !HD(i_w_dT_out,IW) ! 
         
         !HD(i_w_der_in,IW) ! 
         !HD(i_w_der_00,IW) ! 
         !HD(i_w_der_out,IW) ! 
         
         !test_partials = (k+1 == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_r
            s% solver_test_partials_dval_dx = HD(i_w_dr_00,IW)
            write(*,*) 'turbulent_energy_eqn', s% solver_test_partials_var
         end if        
      
      end subroutine turbulent_energy_eqn

      
      subroutine erad_eqn(s, i)  
         use const_def, only: crad, clight
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: IE, k
         real(dp) :: time_scale, residual, XP, d_XP_der, &
            d_XP_dVol_00, d_XP_dr_00, d_XP_dr_in, &
            area_00, area_00_start, area_in, area_in_start, &
            Lr_00, Lr_00_start, d_Lr_00_dFr_00, d_Lr_00_dr_00, &
            Lr_in, Lr_in_start, d_Lr_in_dFr_in, d_Lr_in_dr_in, &
            L_00, L_in, dt, dm, dt_div_dm, DV, opacity, &
            erad, erad_tw, Vol, T, COUPL_factor, COUPL, COUPL1, &
            d_COUPL_dVol, d_COUPL_dT, d_COUPL_der, d_COUPL_dr_in, d_COUPL_dr_00, &
            heat_exchange_time_scale, kE_fac, kP_fac, dt_term, &
            d_dt_term_dT, d_dt_term_dr_in, d_dt_term_dr_00, &
            u_div_r, d_u_div_r_dr_00, d_u_div_r_dr_in, u_div_r_factor
         logical :: test_partials
         include 'formats'
         
         call T_form_of_erad_eqn(s, i)

      end subroutine erad_eqn
      
      
      subroutine Fr_eqn(s, i)
         use const_def, only: clight
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: IL, k
         real(dp) :: residual, &
            XP_00, d_XP_00_dVol_00, d_XP_00_der_00, d_XP_00_dr_00, d_XP_00_dr_in, &
            XP_out, d_XP_out_der_out, d_XP_out_dr_out, d_XP_out_dr_00, &
            u_div_r, d_u_div_r_dr_00, &
            d_lnV, d_dlnV_dVol_00, d_dlnV_dr_in, d_dlnV_dr_00, d_dlnV_dr_out, &
            Vol_out, d_Vol_out_dr_00, d_Vol_out_dr_out, &
            rho_out, d_rho_out_dr_out, d_rho_out_dr_00, &
            Vol_00, d_Vol_00_dr_in, d_Vol_00_dr_00, &
            rho_00, d_rho_00_dVol_00, d_rho_00_dr_00, d_rho_00_dr_in, &
            rho_face, d_rho_dVol_00, d_rho_dr_in, d_rho_dr_00, d_rho_dr_out, &
            r, dt, dm_bar, dt_div_dm_bar, d_dXP_factor_dr_in, d_dXP_factor_dr_00, &
            d_dXP_factor_dr_out, Fr, Fr_start, Fr_tw, &
            kap_00, kap_out, kap_rho_face, rho_face_start, d_dlnV_drho_face, &
            d_kap_00_dVol_00, d_kap_00_dT_00, d_kap_00_dr_00, d_kap_00_dr_in, &
            d_kap_out_dT_out, d_kap_out_dr_00, d_kap_out_dr_out, &
            area, d_area_dr_00, d_kaprho_dVol_00, d_kaprho_dT_00, d_kaprho_dT_out, &
            d_kaprho_dr_in, d_kaprho_dr_00, d_kaprho_dr_out, &
            dXP, dXP_factor, Fr_factor, d_dXP_factor_dVol_00, &
            d_Fr_factor_dVol_00, d_Fr_factor_dT_00, d_Fr_factor_dT_out, &
            d_Fr_factor_dr_in, d_Fr_factor_dr_00, d_Fr_factor_dr_out, &
            erad_00, erad_out, erad_tw, erad_factor1, erad_factor, &
            d_erad_factor_dVol_00, d_erad_factor_dr_in, &
            d_erad_factor_dr_00, d_erad_factor_dr_out
         logical :: test_partials, okay

         include 'formats'
         
         if (s% RSP_hydro_only) then
            k = NZN+1-i
            IL = i_var_Fr + NV*(i-1)
            HD(1:LD_HD,IL) = 0d0        
            residual = - s% Fr(k)  ! want Fr = 0d0
            HR(IL) = -residual
            HD(i_Fr_dFr_00,IL) = -1d0
            return
         end if
         
         call T_form_of_Fr_eqn(s,i)

      end subroutine Fr_eqn
      
      
      subroutine d_Prad_dm_Fr_eqn(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: IL, k
         real(dp) :: residual, Fr_00, &
            dFr_dr_out, dFr_dr_00, dFr_dr_in, &
            dFr_dVol_out, dFr_dVol_00, &
            dFr_dT_out, dFr_dT_00, &
            dFr_der_out, dFr_der_00

         logical :: test_partials

         include 'formats'

         k = NZN+1-i
         IL = i_var_Fr + NV*(i-1)
         HD(1:LD_HD,IL) = 0d0        

         call calc_Fr(s, i, Fr_00, &
            dFr_dr_out, dFr_dr_00, dFr_dr_in, &
            dFr_dVol_out, dFr_dVol_00, &
            dFr_dT_out, dFr_dT_00, &
            dFr_der_out, dFr_der_00)
             
         residual = Fr_00 - s% Fr(k) 
         HR(IL) = -residual
                   
         HD(i_Fr_der_00,IL) = dFr_der_00
         HD(i_Fr_der_out,IL) = dFr_der_out
         HD(i_Fr_dr_in,IL) = dFr_dr_in ! 
         HD(i_Fr_dr_00,IL) = dFr_dr_00 ! 
         HD(i_Fr_dr_out,IL) = dFr_dr_out ! 
         HD(i_Fr_dT_00,IL) = dFr_dT_00
         HD(i_Fr_dT_out,IL) = dFr_dT_out
         HD(i_Fr_dFr_00,IL) = -1d0
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_er
            s% solver_test_partials_dval_dx = dFr_der_00
            write(*,*) 'd_Prad_dm_Fr_eqn', s% solver_test_partials_var
         end if        
            
      end subroutine d_Prad_dm_Fr_eqn

      
      subroutine T_form_of_erad_eqn(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: IE, k
         real(dp) :: T, V, erad_expected, residual
         logical :: test_partials

         include 'formats'

         k = NZN+1-i
         IE = i_var_er + NV*(i-1) 
         HD(1:LD_HD,IE) = 0.d0            
         
         T = s% T(k)
         V = s% Vol(k)
         erad_expected = crad*T**4*V ! ergs/gm
         residual = erad_expected - s% erad(k)
         
         HR(IE) = -residual
         
         HD(i_er_der_00,IE) = -1d0
         HD(i_er_dT_00,IE) = 4d0*crad*T**3*V
         HD(i_er_dr_00,IE) = crad*T**4*dVol_dr_00(i)
         HD(i_er_dr_in,IE) = crad*T**4*dVol_dr_in(i)
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_T
            s% solver_test_partials_dval_dx = HD(i_er_dT_00,IE)
            write(*,*) 'T_form_of_erad_eqn', s% solver_test_partials_var
         end if        
      end subroutine T_form_of_erad_eqn
      
      
      subroutine T_form_of_calc_Fr(s, i, Fr_00, & !rs Stellingwerf 1975, Appendix A 
            dFr_dr_out, dFr_dr_in, dFr_dr_00, &
            dFr_dT_out, dFr_dT_00, dFr_dVol_00)
         type (star_info), pointer :: s
         integer, intent(in) :: I         
         real(dp), intent(out) :: &
            Fr_00, dFr_dr_out, dFr_dr_in, &
            dFr_dr_00, dFr_dT_out, dFr_dT_00, dFr_dVol_00
         real(dp) :: dFr_dK_00, dFr_dK_out, W, WP, BW, BK, T1, T2, T3
         integer :: k
         logical :: test_partials
         include 'formats'
         k = NZN+1-i
         if (i < 1) then
            if (s% RSP_hydro_only) then
               Fr_00 = 0d0
            else
               Fr_00 = s% L_center
            end if
            Fr_00 = Fr_00/(4d0*pi*s% r_center**2)
            dFr_dr_00 = 0 
            dFr_dT_00 = 0
            dFr_dK_00 = 0
            dFr_dK_out = 0
            dFr_dr_out = 0
            dFr_dr_in = 0
            dFr_dT_out = 0
         else if (i == NZN) then
            Fr_00 = 2d0*SIG*s% T(k)**4 !EDDI
            dFr_dT_00 = 4.d0*Fr_00/s% T(k)   
            dFr_dK_00 = 0
            dFr_dK_out = 0
            dFr_dr_out = 0
            dFr_dr_in = 0
            dFr_dr_00 = 0
            dFr_dT_out = 0
         else
            Fr_00 = 0d0
            dFr_dr_00 = 0d0 
            dFr_dT_00 = 0d0
            dFr_dK_00 = 0d0
            dFr_dK_out = 0d0
            dFr_dr_out = 0d0
            dFr_dr_in = 0d0
            dFr_dT_out = 0d0
            W = s% T(k)**4
            WP = s% T(k-1)**4
            BW = 4d0*(s% lnT(k-1) - s% lnT(k)) ! log(WP/W)
            if (abs(BW) < 1d-20) return
            BK = log(s% opacity(k-1)/s% opacity(k))
            if (abs(1.d0 - BK/BW) < 1d-15 .or. abs(BW - BK) < 1d-15) return
            T1 = - CL*s% r(k)**2/(4d0*pi*s% dm_bar(k))   ! CL = 4d0*(4d0*PI)**2*SIG/3d0
            T2 = (WP/s% opacity(k-1) - W/s% opacity(k))/(1.d0 - BK/BW)
            Fr_00 = T1*T2
            T3 = T1/(BW - BK)
            !rs radiative luminosity derivatives
            dFr_dK_00 = (T3/s% opacity(k))*(W*BW/s% opacity(k) - T2)
            dFr_dK_out = -(T3/s% opacity(k-1))*(WP*BW/s% opacity(k-1) - T2)            
            dFr_dr_out = dFr_dK_out*dK_dr_00(i+1) ! 
            dFr_dr_in = dFr_dK_00*dK_dr_in(I) ! 
            dFr_dr_00 = 2d0*Fr_00/s% r(k) & ! 
               + dFr_dK_00*dK_dr_00(I) &
               + dFr_dK_out*dK_dr_in(i+1)
            dFr_dT_out = & ! 
                 4.d0*(T3/s% T(k-1))*(WP*BW/s% opacity(k-1) &
               - T2*BK/BW) + dFr_dK_out*dK_dT(i+1)
            dFr_dT_00 = & ! 
               - 4.d0*(T3/s% T(k))*(W*BW/s% opacity(k) - T2*BK/BW) &
               + dFr_dK_00*dK_dT(I)
            
            if (call_is_bad) then
               if (is_bad(dFr_dT_out + dFr_dT_00)) then
                  write(*,3) 'dFr_dT_out', i, k, dFr_dT_out
                  write(*,3) 'dFr_dT_00', i, k, dFr_dT_00
                  write(*,3) 's% T(k-1)', i, k, s% T(k-1)
                  write(*,3) 's% T(k)', i, k, s% T(k)
                  write(*,3) 's% opacity(k-1)', i, k, s% opacity(k-1)
                  write(*,3) 's% opacity(k)', i, k, s% opacity(k)
                  write(*,3) 'dK_dT(I)', i, k, dK_dT(I)
                  write(*,3) 'dK_dT(i+1)', i, k, dK_dT(i+1)
                  write(*,3) 'T2', i, k, T2
                  write(*,3) 'T3', i, k, T3
                  write(*,3) 'BK', i, k, BK
                  write(*,3) 'BW', i, k, BW
                  write(*,3) '1.d0 - BK/BW', i, k, 1.d0 - BK/BW
                  write(*,3) 'W', i, k, W
                  write(*,3) 'WP', i, k, WP
                  stop 'T_form_of_calc_Fr'
               end if
            end if
            
            dFr_dVol_00 = dFr_dK_00*dK_dVol(I)            
         end if
         
         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = Fr_00
            s% solver_test_partials_var = i_var_T
            s% solver_test_partials_dval_dx = dFr_dT_out
            write(*,*) 'T_form_of_calc_Fr', s% solver_test_partials_var
         end if
      end subroutine T_form_of_calc_Fr
      
      
      subroutine T_form_of_Fr_eqn(s,i)
         type (star_info), pointer :: s
         integer, intent(in) :: i
         integer :: IL, k
         real(dp) :: residual, Fr_00, &
            dFr_00_dr_out, dFr_00_dr_in, dFr_00_dr_00, &
            dFr_00_dT_out, dFr_00_dT_00, dFr_dVol_00
         logical :: test_partials

         include 'formats'

         k = NZN+1-i
         IL = i_var_Fr + NV*(i-1)
         HD(1:LD_HD,IL) = 0d0        

         call T_form_of_calc_Fr(s, i, Fr_00, &
            dFr_00_dr_out, dFr_00_dr_in, dFr_00_dr_00, &
            dFr_00_dT_out, dFr_00_dT_00, dFr_dVol_00)
             
         residual = Fr_00 - s% Fr(k) 
         HR(IL) = -residual
                   
         HD(i_Fr_dFr_00,IL) = -1d0
         HD(i_Fr_dr_in,IL) = dFr_00_dr_in
         HD(i_Fr_dr_00,IL) = dFr_00_dr_00
         HD(i_Fr_dr_out,IL) = dFr_00_dr_out
         HD(i_Fr_dT_00,IL) = dFr_00_dT_00
         HD(i_Fr_dT_out,IL) = dFr_00_dT_out
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = residual
            s% solver_test_partials_var = i_var_T
            s% solver_test_partials_dval_dx = HD(i_Fr_dT_00,IL)
            write(*,*) 'T_form_of_Fr_eqn', s% solver_test_partials_var
         end if        
      end subroutine T_form_of_Fr_eqn


      end module rsp_step
      
      
            
