! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

      module solve_omega_mix

      use star_private_def
      use const_def

      implicit none

      private
      public :: do_solve_omega_mix

      contains

      integer function do_solve_omega_mix(s, dt_total)
         use star_utils, only: start_time, update_time, total_angular_momentum
         use mix_info, only: update_rotation_mixing_info
         use hydro_rotation, only: get_rotation_sigmas, set_omega, set_i_rot

         type (star_info), pointer :: s
         real(dp), intent(in) :: dt_total

         integer :: ierr, nz, i, j, k, max_iters_per_substep, &
            max_iters_total, total_num_iters, num_iters
         integer(8) :: time0, clock_rate
         integer :: steps_used, max_steps, min_steps
         real(qp) :: remaining_time, total_time, time, dt, &
            J_tot0, J_tot1, max_del, avg_del, &
            tol_correction_max, tol_correction_norm
         real(dp) :: total
         real(dp), pointer, dimension(:) :: am_sig_omega, am_sig_j
         real(qp), pointer, dimension(:) :: &
            du, d, dl, x, b, bp, vp, xp, dX, X_0, X_1, rhs, del
         logical :: okay, recalc_mixing_info_each_substep
         logical, parameter :: dbg = .false.

         include 'formats'

         do_solve_omega_mix = keep_going

         ierr = 0

         nz = s% nz
         total_time = dt_total
         time = 0
         max_steps = 20
         min_steps = 4
         max_iters_per_substep = 4
         max_iters_total = 40
         total_num_iters = 0
         tol_correction_max = 1d-4
         tol_correction_norm = 1d-7

         ! update omega for new i_rot and previous j_rot to conserve angular momentum
         call set_i_rot(s, .false.)
         call set_omega(s, 'solve_omega_mix')

         if (dt_total <= 0d0) return

         if (s% doing_timing) call start_time(s, time0, total)

         s% extra_jdot(1:nz) = 0
         s% extra_omegadot(1:nz) = 0

         if (s% use_other_torque) then
            call s% other_torque(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*, *) 'solve_omega_mix: other_torque returned ierr', ierr
               return
            end if
         end if

         if (associated(s% binary_other_torque)) then
            call s% binary_other_torque(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*, *) 'solve_omega_mix: binary_other_torque returned ierr', ierr
               return
            end if
         end if

         call do_alloc(ierr)
         if (ierr /= 0) then
            do_solve_omega_mix = terminate
            s% termination_code = t_solve_omega_mix
            s% result_reason = nonzero_ierr
            if (s% report_ierr) write(*,*) 'allocate failed in do_solve_omega_mix'
            return
         end if

         J_tot0 = total_angular_momentum(s)

         okay = .true.
         do k=1,nz
            if (is_bad_num(s% omega(k)) .or. abs(s% omega(k)) > 1d50) then
               if (s% stop_for_bad_nums) then
                  write(*,2) 's% omega(k)', k, s% omega(k)
                  stop 'solve omega'
               end if
               okay = .false.
               exit
            end if
         end do
         if (.not. okay) then
            write(*,2) 'model_number', s% model_number
            stop 'start solve omega: bad num omega'
         end if

         steps_used = 0
         recalc_mixing_info_each_substep = s% recalc_mixing_info_each_substep

      step_loop: do while &
               (total_time - time > 1d-10*total_time .and. &
                  steps_used < max_steps)

            if (steps_used > 0 .and. recalc_mixing_info_each_substep) then
               call update_rotation_mixing_info(s,ierr)
               if (ierr /= 0) then
                  do_solve_omega_mix = terminate
                  s% termination_code = t_solve_omega_mix
                  s% result_reason = nonzero_ierr
                  if (s% report_ierr) write(*,*) 'update_rotation_mixing_info failed in do_solve_omega_mix'
                  return
               end if
            end if

            if (steps_used == 0 .or. recalc_mixing_info_each_substep) then
               if (steps_used > 0) then
                  do k=1,nz
                     am_sig_omega(k) = s% am_sig_omega(k)
                     am_sig_j(k) = s% am_sig_j(k)
                  end do
               end if
               call get_rotation_sigmas(s, 1, nz, dt_total, ierr)
               if (ierr /= 0) then
                  do_solve_omega_mix = terminate
                  s% termination_code = t_solve_omega_mix
                  s% result_reason = nonzero_ierr
                  if (s% report_ierr) write(*,*) 'get_rotation_sigmas failed in do_solve_omega_mix'
                  return
               end if
               if (steps_used > 0) then
                  do k=1,nz
                     s% am_sig_omega(k) = 0.5d0*(s% am_sig_omega(k) + am_sig_omega(k))
                     s% am_sig_j(k) = 0.5d0*(s% am_sig_j(k) + am_sig_j(k))
                  end do
               end if
            end if

            steps_used = steps_used + 1

            dt = 0.5d0*min_mixing_timescale()
            remaining_time = total_time - time
            dt = max(dt, 1d-6*remaining_time)
            dt = min(dt, real(dt_total/min_steps,kind=qp))
            if (dt >= remaining_time) then
               dt = remaining_time
            else
               dt = min(dt, 0.5d0*remaining_time)
            end if
            if (steps_used >= max_steps) dt = remaining_time ! just go for it
            if (dbg) write(*,3) 'mix dt', &
                  s% model_number, steps_used, dt, dt/remaining_time

            ! X_0 is omega at start of substep
            ! X_1 is current candidate for omega at end of substep
            ! dX = X_1 - X_0
            do k=1,nz
               X_0(k) = s% omega(k)
               X_1(k) = X_0(k)
               dX(k) = 0d0
            end do

         solve_loop: do num_iters = 1, max_iters_per_substep

               if (total_num_iters >= max_iters_total) then
                  s% retry_message = 'solve omega mix failed to converge in allowed number of steps'
                  do_solve_omega_mix = retry
                  exit step_loop
               end if

               total_num_iters = total_num_iters+1

               if (s% use_other_torque_implicit) then
                  call s% other_torque_implicit(s% id, ierr)
                  if (ierr /= 0) then
                     s% retry_message = 'other_torque_implicit returned ierr'
                     do_solve_omega_mix = retry
                     exit step_loop
                  end if
               end if

               if (associated(s% binary_other_torque_implicit)) then
                  call s% binary_other_torque_implicit(s% id, ierr)
                  if (ierr /= 0) then
                     s% retry_message = 'binary_other_torque_implicit returned ierr'
                     do_solve_omega_mix = retry
                     exit step_loop
                  end if
               end if

               call create_matrix_and_rhs(dt)

               ! solve for del
               call solve_tridiag(dl, d, du, rhs(1:nz), del(1:nz), nz, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'matrix solve failed in solve mix'
                  do_solve_omega_mix = retry
                  exit step_loop
               end if

               ! apply the correction dX = dX + del
               ! X_1 = X_0 + dX
               ! X_0 is omega at start of substep
               ! X_1 is candidate for omega at end of substep
               do k=2,nz
                  dX(k) = dX(k) + del(k)
                  X_1(k) = X_0(k) + dX(k)
                  s% omega(k) = X_1(k)
               end do
               s% omega(1) = s% omega(2)

               ! if correction small enough, exit solve_loop
               max_del = maxval(abs(del(1:nz)))
               avg_del = sum(abs(del(1:nz)))/nz
               if (max_del <= tol_correction_max .and. avg_del <= tol_correction_norm) then
                  if (dbg) &
                     write(*,3) 'substep converged: iters max_del avg_del dt/total', &
                        steps_used, num_iters, max_del, avg_del, dt/total_time
                  exit solve_loop ! this substep is done
               end if

               if (num_iters == max_iters_per_substep) then
                  s% retry_message = 'num_iters == max_iters_per_substep in solve mix'
                  do_solve_omega_mix = retry
                  exit step_loop
               end if

            end do solve_loop

            time = time + dt

         end do step_loop

         !if (recalc_mixing_info_each_substep) &
         !   write(*,3) 'omega mix steps_used', steps_used, s% model_number

         if (dbg) write(*,2) 'omega mix steps_used', steps_used

         s% num_rotation_solver_steps = max(steps_used, s% num_rotation_solver_steps)

         if (do_solve_omega_mix == keep_going .and. total_time - time > 1d-10*total_time) then
            do_solve_omega_mix = retry
            s% retry_message = 'failed in mixing angular momentum'
         end if

         if (do_solve_omega_mix == keep_going) then

            okay = .true.
            do k=1,nz
               if (is_bad_num(s% omega(k)) .or. abs(s% omega(k)) > 1d50) then
                  write(*,2) 's% omega(k)', k, s% omega(k)
                  okay = .false.
                  exit
               end if
            end do
            if (.not. okay) then
               write(*,2) 'model_number', s% model_number
               stop 'end solve omega'
            end if

            do k=1,nz
               s% j_rot(k) = s% i_rot(k)*s% omega(k)
            end do

            if (.not. (s% use_other_torque .or. s% use_other_torque_implicit .or. &
                  associated(s% binary_other_torque))) then

               ! check conservation for cases with no extra torque
               J_tot1 = total_angular_momentum(s) ! what we have

               if (abs(J_tot0 - J_tot1) > s% angular_momentum_error_retry*abs(J_tot0)) then
                  s% retry_message = 'retry: failed to conserve angular momentum in mixing'
                  write(*,*) "angular momentum error larger than angular_momentum_error_retry", abs(J_tot0 - J_tot1)/abs(J_tot0)
                  do_solve_omega_mix = retry
               else if (abs(J_tot0 - J_tot1) > s% angular_momentum_error_warn*abs(J_tot0)) then
                  write(*,*) "angular momentum error larger than angular_momentum_error_warn", abs(J_tot0 - J_tot1)/abs(J_tot0)
               end if
               if (dbg) then
                  write(*,2) 'final J_tot1', s% model_number, J_tot1
                  write(*,2) '(J_tot1 - J_tot0)/J_tot0', &
                     steps_used, (J_tot1 - J_tot0)/J_tot0, J_tot0, J_tot1
               end if
            end if

         end if

         if (dbg) write(*,*)

         call dealloc

         if (s% doing_timing) &
            call update_time(s, time0, total, s% time_solve_omega_mix)


         contains


         subroutine do_alloc(ierr)
            use alloc, only: non_crit_get_quad_array
            integer, intent(out) :: ierr            
            call do_work_arrays(.true.,ierr)

            call non_crit_get_quad_array(s, du, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, d, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, dl, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, x, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, b, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, bp, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, vp, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, xp, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, dX, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, X_0, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, X_1, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, rhs, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, del, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return

         end subroutine do_alloc


         subroutine dealloc
            use alloc, only: non_crit_return_quad_array
            call do_work_arrays(.false.,ierr)
            
            call non_crit_return_quad_array(s, du, 'solve_omega_mix')
            call non_crit_return_quad_array(s, d, 'solve_omega_mix')
            call non_crit_return_quad_array(s, dl, 'solve_omega_mix')
            call non_crit_return_quad_array(s, x, 'solve_omega_mix')
            call non_crit_return_quad_array(s, b, 'solve_omega_mix')
            call non_crit_return_quad_array(s, bp, 'solve_omega_mix')
            call non_crit_return_quad_array(s, vp, 'solve_omega_mix')
            call non_crit_return_quad_array(s, xp, 'solve_omega_mix')

            call non_crit_return_quad_array(s, dX, 'solve_omega_mix')
            call non_crit_return_quad_array(s, X_0, 'solve_omega_mix')
            call non_crit_return_quad_array(s, X_1, 'solve_omega_mix')
            call non_crit_return_quad_array(s, rhs, 'solve_omega_mix')
            call non_crit_return_quad_array(s, del, 'solve_omega_mix')
            
         end subroutine dealloc
         
         
         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               am_sig_omega, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               am_sig_j, nz, nz_alloc_extra, 'solve_omega_mix', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays


         subroutine solve_tridiag(sub, diag, sup, rhs, x, n, ierr)
            implicit none
            !      sub - sub-diagonal
            !      diag - the main diagonal
            !      sup - sup-diagonal
            !      rhs - right hand side
            !      x - the answer
            !      n - number of equations

            integer, intent(in) :: n
            real(qp), dimension(:), intent(in) :: sup, diag, sub
            real(qp), dimension(:), intent(in) :: rhs
            real(qp), dimension(:), intent(out) :: x
            integer, intent(out) :: ierr

            real(qp) :: m
            integer i

            ierr = 0

            bp(1) = diag(1)
            vp(1) = rhs(1)

            do i = 2,n
               m = sub(i-1)/bp(i-1)
               bp(i) = diag(i) - m*sup(i-1)
               vp(i) = rhs(i) - m*vp(i-1)
            end do

            xp(n) = vp(n)/bp(n)
            x(n) = xp(n)
            do i = n-1, 1, -1
               xp(i) = (vp(i) - sup(i)*xp(i+1))/bp(i)
               x(i) = xp(i)
            end do

         end subroutine solve_tridiag


         real(dp) function min_mixing_timescale() result(dt)
            integer :: k
            real(dp) :: & ! use dp instead of qp to get same answer in ifort and gfortran
               omega, irot, irot_mid_00, am_sig_omega_00, c_omega_00, del00_omega, &
               omega_mid_00, am_sig_irot_00, c_irot_00, del00_irot, &
               dmbar, irot_mid_m1, am_sig_omega_m1, c_omega_m1, delm1_omega, &
               omega_mid_m1, am_sig_irot_m1, c_irot_m1, delm1_irot, &
               d2omega, d2irot, dt00

            include 'formats'

            dt = 1d99

            do k = 1, nz

               omega = s% omega(k)
               irot = s% i_rot(k)

               if (k < nz) then

                  irot_mid_00 = 0.5d0*(irot + s% i_rot(k+1))
                  am_sig_omega_00 = s% am_sig_omega(k) + s% am_sig_j(k)
                  c_omega_00 = am_sig_omega_00*irot_mid_00
                  del00_omega = omega - s% omega(k+1)

                  omega_mid_00 = 0.5d0*(omega + s% omega(k+1))
                  am_sig_irot_00 = s% am_sig_j(k)
                  c_irot_00 = am_sig_irot_00*omega_mid_00
                  del00_irot = irot - s% i_rot(k+1)

               else

                  c_omega_00 = 0
                  del00_omega = 0
                  c_irot_00 = 0
                  del00_irot = 0

               end if

               if (k > 1) then

                  if (k < nz) then
                     dmbar = 0.5d0*(s% dm(k-1) + s% dm(k))
                  else
                     dmbar = 0.5d0*s% dm(k-1) + s% dm(k)
                  end if

                  irot_mid_m1 = 0.5d0*(s% i_rot(k-1) + irot)
                  am_sig_omega_m1 = s% am_sig_omega(k-1) + s% am_sig_j(k-1)
                  c_omega_m1 = am_sig_omega_m1*irot_mid_m1
                  delm1_omega = s% omega(k-1) - omega

                  omega_mid_m1 = 0.5d0*(s% omega(k-1) + omega)
                  am_sig_irot_m1 = s% am_sig_j(k-1)
                  c_irot_m1 = am_sig_irot_m1*omega_mid_m1
                  delm1_irot = s% i_rot(k-1) - irot

               else

                  dmbar = 0.5d0*s% dm(k)
                  c_omega_m1 = 0
                  delm1_omega = 0
                  c_irot_m1 = 0
                  delm1_irot = 0

               end if

               if (k == 1) then
                  d2omega = -c_omega_00*del00_omega
                  d2irot = -c_irot_00*del00_irot
               else if (k == nz) then
                  d2omega = c_omega_m1*delm1_omega
                  d2irot = c_irot_m1*delm1_irot
               else
                  d2omega = c_omega_m1*delm1_omega - c_omega_00*del00_omega
                  d2irot = c_irot_m1*delm1_irot - c_irot_00*del00_irot
               end if

               dt00 = max(1d-12,abs(omega))*irot/ &
                  max(1d-50, abs((d2omega + d2irot)/dmbar + &
                                 s% extra_omegadot(k)*irot + s% extra_jdot(k)))
               if (dt00 < dt) dt = dt00

            end do

         end function min_mixing_timescale


         subroutine create_matrix_and_rhs(dt)
            ! basic equation from Heger, Langer, & Woosley, 2000, eqn 46.
            ! with source terms added.
            ! and term for j curvature as well as omega curvature
            real(qp), intent(in) :: dt
            integer :: k
            real(qp) :: &
               dmbar, f, &
               omega, omega_mid_00, omega_mid_m1, &
               irot, irot_mid_00, irot_mid_m1, &
               am_sig_omega_00, am_sig_omega_m1, c_omega_00, c_omega_m1, &
               am_sig_irot_00, am_sig_irot_m1, c_irot_00, c_irot_m1, &
               d_c_irot_00_domega_p1, d_c_irot_00_domega_00, &
               d_c_irot_m1_domega_00, d_c_irot_m1_domega_m1, &
               del00_omega, delm1_omega, &
               del00_irot, delm1_irot, &
               d2omega, d_d2omega_domega_p1, d_d2omega_domega_m1, d_d2omega_domega_00, &
               d2irot, d_d2irot_domega_p1, d_d2irot_domega_m1, d_d2irot_domega_00

            include 'formats'

            do k = 1, nz

               omega = s% omega(k)
               irot = s% i_rot(k)

               if (k < nz) then

                  irot_mid_00 = 0.5d0*(irot + s% i_rot(k+1))
                  am_sig_omega_00 = s% am_sig_omega(k) + s% am_sig_j(k)
                  c_omega_00 = am_sig_omega_00*irot_mid_00
                  del00_omega = omega - s% omega(k+1)

                  omega_mid_00 = 0.5d0*(omega + s% omega(k+1))
                  am_sig_irot_00 = s% am_sig_j(k)
                  c_irot_00 = am_sig_irot_00*omega_mid_00
                  d_c_irot_00_domega_p1 = 0.5d0*am_sig_irot_00
                  d_c_irot_00_domega_00 = 0.5d0*am_sig_irot_00
                  del00_irot = irot - s% i_rot(k+1)

               else

                  c_omega_00 = 0
                  del00_omega = 0
                  c_irot_00 = 0
                  d_c_irot_00_domega_p1 = 0
                  d_c_irot_00_domega_00 = 0
                  del00_irot = 0

               end if

               if (k > 1) then

                  if (k < nz) then
                     dmbar = 0.5d0*(s% dm(k-1) + s% dm(k))
                  else
                     dmbar = 0.5d0*s% dm(k-1) + s% dm(k)
                  end if

                  irot_mid_m1 = 0.5d0*(s% i_rot(k-1) + irot)
                  am_sig_omega_m1 = s% am_sig_omega(k-1) + s% am_sig_j(k-1)
                  c_omega_m1 = am_sig_omega_m1*irot_mid_m1
                  delm1_omega = s% omega(k-1) - omega

                  omega_mid_m1 = 0.5d0*(s% omega(k-1) + omega)
                  am_sig_irot_m1 = s% am_sig_j(k-1)
                  c_irot_m1 = am_sig_irot_m1*omega_mid_m1
                  d_c_irot_m1_domega_00 = 0.5d0*am_sig_irot_m1
                  d_c_irot_m1_domega_m1 = 0.5d0*am_sig_irot_m1
                  delm1_irot = s% i_rot(k-1) - irot

               else

                  dmbar = 0.5d0*s% dm(k)
                  c_omega_m1 = 0
                  delm1_omega = 0
                  c_irot_m1 = 0
                  d_c_irot_m1_domega_00 = 0
                  d_c_irot_m1_domega_m1 = 0
                  delm1_irot = 0

               end if

               if (k == 1) then
                  d2omega = -c_omega_00*del00_omega
                  d2irot = -c_irot_00*del00_irot
               else if (k == nz) then
                  d2omega = c_omega_m1*delm1_omega
                  d2irot = c_irot_m1*delm1_irot
               else
                  d2omega = c_omega_m1*delm1_omega - c_omega_00*del00_omega
                  d2irot = c_irot_m1*delm1_irot - c_irot_00*del00_irot
               end if
               d_d2omega_domega_00 = -(c_omega_m1 + c_omega_00)
               d_d2irot_domega_00 = &
                  d_c_irot_m1_domega_00*delm1_irot - d_c_irot_00_domega_00*del00_irot

               ! X_1 = X_0 + dX
               ! X_0 is omega at start of substep
               ! X_1 is candidate for omega at end of substep

               ! residual = dX - dt*(((d2omega+d2irot)/dmbar + extra_jdot)/irot + extra_omegadot)
               ! J = d(residual)/d(omega)
               ! del is linear estimate of change to dX to make residual = 0
               ! solve J*del = -residual == rhs

               rhs(k) = -dX(k) + &
                  dt*(((d2omega + d2irot)/dmbar + &
                        s% extra_jdot(k))/irot + s% extra_omegadot(k))

               f = dt/(dmbar*irot)
               d(k) = 1d0 - (d_d2omega_domega_00 + d_d2irot_domega_00)*f

               if (k < nz) then
                  d_d2omega_domega_p1 = c_omega_00
                  d_d2irot_domega_p1 = -d_c_irot_00_domega_p1*del00_irot
                  du(k) = -(d_d2omega_domega_p1 + d_d2irot_domega_p1)*f
               end if

               if (k > 1) then
                  d_d2omega_domega_m1 = c_omega_m1
                  d_d2irot_domega_m1 = d_c_irot_m1_domega_m1*delm1_irot
                  dl(k-1) = -(d_d2omega_domega_m1 + d_d2irot_domega_m1)*f
               end if

               if (s% use_other_torque_implicit) then
                  d(k) = d(k) - &
                     dt*(s% d_extra_jdot_domega_00(k)/irot + &
                           s% d_extra_omegadot_domega_00(k))
                  if (k < nz) &
                     du(k) = du(k) - &
                        dt*(s% d_extra_jdot_domega_p1(k)/irot + &
                              s% d_extra_omegadot_domega_p1(k))
                  if (k > 1) &
                     dl(k-1) = dl(k-1) - &
                        dt*(s% d_extra_jdot_domega_m1(k)/irot + &
                              s% d_extra_omegadot_domega_m1(k))
               end if

            end do

         end subroutine create_matrix_and_rhs


      end function do_solve_omega_mix


      end module solve_omega_mix


