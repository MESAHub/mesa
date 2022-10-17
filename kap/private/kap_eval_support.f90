! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module kap_eval_support
   
   use const_def, only : dp, one_sixth
   use math_lib
   
   implicit none

contains
   
   
   subroutine Locate_log(&
      rq, num_logs, log_min_in, log_max_in, ili_logs, logs, log_find, i, log0, log1, ierr)
      use kap_def
      use utils_lib, only : is_bad
      use num_lib, only : binary_search
      type (Kap_General_Info), pointer :: rq
      integer, intent(in) :: num_logs, ili_logs
      real(dp), intent(in) :: log_min_in, log_max_in
      real(dp), intent(in), pointer :: logs(:) ! (num_logs)
      real(dp), intent(inout) :: log_find ! can change log_find if clipping to table boundaries
      integer, intent(out) :: i ! index in logs s.t. logs(i) <= log_find < logs(i+1)
      ! one exception: if log_find == log_max then will get i = num_logs-1
      real(dp), intent(out) :: log0, log1
      integer, intent(out) :: ierr
      real(dp) :: dlog, log_min, log_max
      integer :: j
      ierr = 0
      log_min = max(log_min_in, logs(1))
      log_max = min(log_max_in, logs(num_logs))
      if (num_logs == 1) then
         i = 1
         log_find = log_min
         log0 = log_find
         log1 = log_find
         return
      end if
      if (log_find < log_min .or. log_find > log_max) then
         if (.not. clip_to_kap_table_boundaries) then
            ierr = -1
            return
         end if
         if (log_find < log_min) then
            i = 1
            log_find = log_min
         else
            i = num_logs - 1
            log_find = log_max
         end if
      else if (abs(log_find - log_max) < 1d-7) then
         i = num_logs - 1
         log_find = log_max
      else if (ili_logs == 1) then ! logs equally spaced
         dlog = (log_max - log_min) / (num_logs - 1)
         i = int((log_find - log_min) / dlog) + 1
         ! might not be exactly evenly spaced, so minor fixup if necessary
         if (logs(i) > log_find .and. i > 1) then
            i = i - 1
         else if (log_find >= logs(i + 1) .and. i + 1 < num_logs) then
            i = i + 1
         end if
      else
         i = binary_search(num_logs, logs, 0, log_find)
         if (i >= num_logs) then
            ierr = -1
            return
            !$OMP critical (kap_eval_crit1)
            write(*, *) 'i', i
            write(*, *) 'num_logs', num_logs
            call mesa_error(__FILE__, __LINE__, 'Locate_log')
            !$OMP end critical (kap_eval_crit1)
         end if
      end if
      
      if (i < 1 .or. i >= num_logs) then
         ierr = -1
         return
         write(*, *) 'i', i
         write(*, *) 'num_logs', num_logs
         call mesa_error(__FILE__, __LINE__, 'Locate_log')
      end if
      
      if (logs(i) > log_find .or. log_find > logs(i + 1)) then
         ierr = -1
         !$OMP critical (kap_eval_crit2)
         write(*, *) 'dlog', (log_max - log_min) / (num_logs - 1)
         write(*, *) 'log_max', log_max
         write(*, *) 'log_min', log_min
         write(*, *) 'log_max_in', log_max_in
         write(*, *) 'log_min_in', log_min_in
         write(*, *) 'logs(i)', logs(i)
         write(*, *) 'log_find', log_find
         write(*, *) 'logs(i+1)', logs(i + 1)
         write(*, *) 'i', i
         write(*, *) 'num_logs', num_logs
         write(*, *) 'ili_logs', ili_logs
         !$OMP end critical (kap_eval_crit2)
         return
         call mesa_error(__FILE__, __LINE__, 'error in Locate_log')
      end if
      log0 = logs(i)
      log1 = logs(i + 1)
      if (is_bad(log0) .or. is_bad(log1)) then
         ierr = -1
         return
         !$OMP critical (kap_eval_crit3)
         write(*, *) 'logs(i)', logs(i)
         write(*, *) 'log_find', log_find
         write(*, *) 'logs(i+1)', logs(i + 1)
         write(*, *) 'i', i
         write(*, *) 'num_logs', num_logs
         call mesa_error(__FILE__, __LINE__, 'error in Locate_log')
         !$OMP end critical (kap_eval_crit3)
      end if
   end subroutine Locate_log
   
   
   subroutine Locate_logT(&
      rq, num_logTs, logT_min, logT_max, ili_logTs, logTs, logT, iT, logT0, logT1, ierr)
      use kap_def
      type (Kap_General_Info), pointer :: rq
      integer, intent(in) :: num_logTs, ili_logTs
      real(dp), intent(in) :: logT_min, logT_max
      real(dp), intent(in), pointer :: logTs(:) ! (num_logTs)
      real(dp), intent(inout) :: logT ! can change logT if clipping to table boundaries
      integer, intent(out) :: iT ! index in logTs s.t. logTs(i) <= logT < logTs(i+1)
      real(dp), intent(out) :: logT0, logT1
      integer, intent(out) :: ierr
      call Locate_log(&
         rq, num_logTs, logT_min, logT_max, ili_logTs, logTs, logT, iT, logT0, logT1, ierr)
   end subroutine Locate_logT
   
   
   subroutine Locate_logR(&
      rq, num_logRs, logR_min, logR_max, ili_logRs, logRs, logR, iR, logR0, logR1, ierr)
      use kap_def
      type (Kap_General_Info), pointer :: rq
      integer, intent(in) :: num_logRs, ili_logRs
      real(dp), intent(in) :: logR_min, logR_max
      real(dp), intent(in), pointer :: logRs(:) ! (num_logRs)
      real(dp), intent(inout) :: logR ! can change logR if clipping to table boundaries
      integer, intent(out) :: iR ! index in logRs s.t. logRs(i) <= logR < logRs(i+1)
      real(dp), intent(out) :: logR0, logR1
      integer, intent(out) :: ierr
      call Locate_log(&
         rq, num_logRs, logR_min, logR_max, ili_logRs, logRs, logR, iR, logR0, logR1, ierr)
   end subroutine Locate_logR
   
   
   subroutine Do_Kap_Interpolations(&
      fin1, nx, ny, i, j, x0, xget, x1, y0, yget, y1, fval, df_dx, df_dy)
      ! derived from routines in the PSPLINE package written by Doug McCune
      
      real(dp), dimension(:), pointer :: fin1 ! the spline data array, dimensions (4, nx, ny)
      integer, intent(in) :: nx, ny, i, j           ! target cell in the spline data
      real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
      real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
      real(dp), intent(out) :: fval, df_dx, df_dy
      
      real(dp), parameter :: z36th = 1d0 / 36d0
      
      real(dp), pointer :: fin(:, :, :)
      
      real(dp) :: xp, xpi, xp2, xpi2, cx, cxi, hx2, cxd, cxdi, hx, hxi
      real(dp) :: yp, ypi, yp2, ypi2, cy, cyi, hy2, cyd, cydi, hy, hyi
      
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      fin(1:4, 1:nx, 1:ny) => fin1(1:4 * nx * ny)
      
      hx = x1 - x0
      hxi = 1d0 / hx
      hx2 = hx * hx
      
      xp = (xget - x0) * hxi
      xpi = 1d0 - xp
      xp2 = xp * xp
      xpi2 = xpi * xpi
      
      cx = xp * (xp2 - 1d0)
      cxi = xpi * (xpi2 - 1d0)
      cxd = 3d0 * xp2 - 1d0
      cxdi = -3d0 * xpi2 + 1d0
      
      hy = y1 - y0
      hyi = 1d0 / hy
      hy2 = hy * hy
      
      yp = (yget - y0) * hyi
      ypi = 1d0 - yp
      yp2 = yp * yp
      ypi2 = ypi * ypi
      
      cy = yp * (yp2 - 1d0)
      cyi = ypi * (ypi2 - 1d0)
      cyd = 3d0 * yp2 - 1d0
      cydi = -3d0 * ypi2 + 1d0
      
      ! bicubic spline interpolation
      fval = &
         xpi * (ypi * fin(1, i, j) + yp * fin(1, i, j + 1)) &
            + xp * (ypi * fin(1, i + 1, j) + yp * fin(1, i + 1, j + 1)) &
            + one_sixth * hx2 * (&
            cxi * (ypi * fin(2, i, j) + yp * fin(2, i, j + 1)) + &
               cx * (ypi * fin(2, i + 1, j) + yp * fin(2, i + 1, j + 1))) &
            + one_sixth * hy2 * (&
            xpi * (cyi * fin(3, i, j) + cy * fin(3, i, j + 1)) + &
               xp * (cyi * fin(3, i + 1, j) + cy * fin(3, i + 1, j + 1))) &
            + z36th * hx2 * hy2 * (&
            cxi * (cyi * fin(4, i, j) + cy * fin(4, i, j + 1)) + &
               cx * (cyi * fin(4, i + 1, j) + cy * fin(4, i + 1, j + 1)))
      
      if (.false.) then
         write(*, 3) 'fin(1,i,j)', i, j, fin(1, i, j)
         write(*, 3) 'fin(1,i,j+1)', i, j + 1, fin(1, i, j + 1)
         write(*, 3) 'fin(1,i+1,j)', i + 1, j, fin(1, i + 1, j)
         write(*, 3) 'fin(1,i+1,j+1)', i + 1, j + 1, fin(1, i + 1, j + 1)
         write(*, 3) 'fin(2,i,j)', i, j, fin(2, i, j)
         write(*, 3) 'fin(2,i,j+1)', i, j + 1, fin(2, i, j + 1)
         write(*, 3) 'fin(2,i+1,j)', i + 1, j, fin(2, i + 1, j)
         write(*, 3) 'fin(2,i+1,j+1)', i + 1, j + 1, fin(2, i + 1, j + 1)
         write(*, 3) 'fin(3,i,j)', i, j, fin(3, i, j)
         write(*, 3) 'fin(3,i,j+1)', i, j + 1, fin(3, i, j + 1)
         write(*, 3) 'fin(3,i+1,j)', i + 1, j, fin(3, i + 1, j)
         write(*, 3) 'fin(3,i+1,j+1)', i + 1, j + 1, fin(3, i + 1, j + 1)
         write(*, 3) 'fin(4,i,j)', i, j, fin(4, i, j)
         write(*, 3) 'fin(4,i,j+1)', i, j + 1, fin(4, i, j + 1)
         write(*, 3) 'fin(4,i+1,j)', i + 1, j, fin(4, i + 1, j)
         write(*, 3) 'fin(4,i+1,j+1)', i + 1, j + 1, fin(4, i + 1, j + 1)
         write(*, 1) 'ypi', ypi
         write(*, 1) 'yp', yp
         write(*, 1) 'xp', xp
         write(*, 1) 'yp', yp
         write(*, 1) 'hx2', hx2
         write(*, 1) 'cxi', cxi
         write(*, 1) 'cx', cx
         write(*, 1) 'hy2', hy2
         write(*, 1) 'cy', cy
         write(*, 1) 'cyi', cyi
         write(*, 1) 'one_sixth', one_sixth
         write(*, 1) 'z36th', z36th
         write(*, 1) 't1', xpi * (ypi * fin(1, i, j) + yp * fin(1, i, j + 1))
         write(*, 1) 't2', xp * (ypi * fin(1, i + 1, j) + yp * fin(1, i + 1, j + 1))
         write(*, 1) 't3', one_sixth * hx2 * (&
            cxi * (ypi * fin(2, i, j) + yp * fin(2, i, j + 1)) + &
               cx * (ypi * fin(2, i + 1, j) + yp * fin(2, i + 1, j + 1)))
         write(*, 1) 't4', one_sixth * hy2 * (&
            xpi * (cyi * fin(3, i, j) + cy * fin(3, i, j + 1)) + &
               xp * (cyi * fin(3, i + 1, j) + cy * fin(3, i + 1, j + 1)))
         write(*, 1) 't5', z36th * hx2 * hy2 * (&
            cxi * (cyi * fin(4, i, j) + cy * fin(4, i, j + 1)) + &
               cx * (cyi * fin(4, i + 1, j) + cy * fin(4, i + 1, j + 1)))
         write(*, 1) 'fval', fval
      end if
      
      ! derivatives of bicubic splines
      df_dx = &
         hxi * (&
            -(ypi * fin(1, i, j) + yp * fin(1, i, j + 1)) &
               + (ypi * fin(1, i + 1, j) + yp * fin(1, i + 1, j + 1))) &
            + one_sixth * hx * (&
            cxdi * (ypi * fin(2, i, j) + yp * fin(2, i, j + 1)) + &
               cxd * (ypi * fin(2, i + 1, j) + yp * fin(2, i + 1, j + 1))) &
            + one_sixth * hxi * hy2 * (&
            -(cyi * fin(3, i, j) + cy * fin(3, i, j + 1)) &
               + (cyi * fin(3, i + 1, j) + cy * fin(3, i + 1, j + 1))) &
            + z36th * hx * hy2 * (&
            cxdi * (cyi * fin(4, i, j) + cy * fin(4, i, j + 1)) + &
               cxd * (cyi * fin(4, i + 1, j) + cy * fin(4, i + 1, j + 1)))
      
      df_dy = &
         hyi * (&
            xpi * (-fin(1, i, j) + fin(1, i, j + 1)) + &
               xp * (-fin(1, i + 1, j) + fin(1, i + 1, j + 1))) &
            + one_sixth * hx2 * hyi * (&
            cxi * (-fin(2, i, j) + fin(2, i, j + 1)) + &
               cx * (-fin(2, i + 1, j) + fin(2, i + 1, j + 1))) &
            + one_sixth * hy * (&
            xpi * (cydi * fin(3, i, j) + cyd * fin(3, i, j + 1)) + &
               xp * (cydi * fin(3, i + 1, j) + cyd * fin(3, i + 1, j + 1))) &
            + z36th * hx2 * hy * (&
            cxi * (cydi * fin(4, i, j) + cyd * fin(4, i, j + 1)) + &
               cx * (cydi * fin(4, i + 1, j) + cyd * fin(4, i + 1, j + 1)))
   
   end subroutine Do_Kap_Interpolations


end module kap_eval_support
      
