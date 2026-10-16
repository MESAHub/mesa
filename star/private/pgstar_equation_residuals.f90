! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module pgstar_equation_residuals
   ! Composition equations are not included.
   use star_private_def
   use const_def, only: dp
   use pgstar_support
   use star_pgstar

   implicit none

   private
   public :: Max_eq_resid_plot, do_Max_eq_resid_plot, Max_eq_resid_record

   real(dp), parameter :: resid_floor = 1d-20, resid_ceil = 1d30
   integer, parameter :: default_capacity = 128

   ! per-star history; indexed by s% id so multiple stars (binaries) don't mix
   type resid_hist_t
      integer :: nvar = 0
      integer :: count = 0
      integer, allocatable :: model(:)
      real, allocatable :: vals(:, :)  ! (nvar, capacity)
      character(len=name_len), allocatable :: names(:)
   end type resid_hist_t

   type(resid_hist_t), target, save :: hist(max_star_handles)

contains

   subroutine Max_eq_resid_plot(id, device_id, ierr)
      integer, intent(in)  :: id, device_id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call do_Max_eq_resid_plot(s, id, &
         s% pg% Max_eq_resid_xleft, s% pg% Max_eq_resid_xright, &
         s% pg% Max_eq_resid_ybot, s% pg% Max_eq_resid_ytop, .false., &
         s% pg% Max_eq_resid_title, s% pg% Max_eq_resid_txt_scale, &
         s% pg% Max_eq_resid_max_width, ierr)

      call pgebuf()

   end subroutine Max_eq_resid_plot


   subroutine Max_eq_resid_record(s, max_width, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: max_width  ! <= 0: unlimited
      integer, intent(out) :: ierr

      type(resid_hist_t), pointer :: h
      integer :: nvar, nz, i_eq, k, c

      ierr = 0
      if (s% id < 1 .or. s% id > max_star_handles) then
         ierr = -1
         return
      end if
      h => hist(s% id)

      nvar = s% nvar_hydro
      nz = s% nz
      if (nvar <= 0 .or. nz <= 0 .or. .not. associated(s% equ)) return
      if (size(s% equ, 1) < nvar .or. size(s% equ, 2) < nz) return

      ! (re)initialize on first call or if the set of equations changed
      if (.not. allocated(h% model) .or. h% nvar /= nvar) then
         if (max_width > 0) then
            c = max_width
         else
            c = default_capacity
         end if
         call reset_hist(h, nvar, c, ierr)
         if (ierr /= 0) return
         do i_eq = 1, nvar
            h% names(i_eq) = s% nameofequ(i_eq)
         end do
      end if

      ! drop entries that are not older than the current model
      do while (h% count > 0)
         if (h% model(h% count) < s% model_number) exit
         h% count = h% count - 1
      end do

      if (max_width > 0 .and. h% count >= max_width) then
         ! roll: keep the newest max_width-1 entries
         k = h% count - max_width + 1
         h% model(1:h% count - k) = h% model(k + 1:h% count)
         h% vals(:, 1:h% count - k) = h% vals(:, k + 1:h% count)
         h% count = h% count - k
      end if

      if (h% count >= size(h% model)) then
         call grow_hist(h, 2 * size(h% model), ierr)
         if (ierr /= 0) return
      end if

      h% count = h% count + 1
      h% model(h% count) = s% model_number
      do i_eq = 1, nvar
         h% vals(i_eq, h% count) = real(min(resid_ceil, &
            max(resid_floor, maxval(abs(s% equ(i_eq, 1:nz))))))
      end do

   end subroutine Max_eq_resid_record


   subroutine do_Max_eq_resid_plot(s, id, &
         win_xleft, win_xright, win_ybot, win_ytop, subplot, title, txt_scale, &
         max_width, ierr)

      type (star_info), pointer :: s
      integer, intent(in) :: id, max_width
      logical, intent(in) :: subplot
      real, intent(in) :: win_xleft, win_xright, win_ybot, win_ytop, txt_scale
      character (len=*), intent(in) :: title
      integer, intent(out) :: ierr

      type(resid_hist_t), pointer :: h
      real, allocatable :: xvec(:), yvec(:)
      real :: xmin, xmax, ymin, ymax, dy_leg
      integer :: i_eq, n

      ierr = 0
      call Max_eq_resid_record(s, max_width, ierr)
      if (ierr /= 0) return

      h => hist(s% id)
      n = h% count
      if (n < 1) return

      allocate(xvec(n), yvec(n))
      xvec = real(h% model(1:n))

      xmin = xvec(1)
      xmax = xvec(n)
      if (xmin >= xmax) xmax = xmin + 1.0

      ymin = log10(minval(h% vals(:, 1:n)))
      ymax = log10(maxval(h% vals(:, 1:n)))
      if (ymax - ymin < 1.0) then
         ymin = ymin - 0.5
         ymax = ymax + 0.5
      end if
      ymin = ymin - 0.05 * (ymax - ymin)
      ymax = ymax + 0.05 * (ymax - ymin)

      call pgsave
      call pgsci(0)                       ! background
      call pgsvp(win_xleft, win_xright, win_ybot, win_ytop)
      call pgrect(0.0, 1.0, 0.0, 1.0)     ! pseudo erase
      call pgswin(xmin, xmax, ymin, ymax)
      call pgsch(txt_scale)
      call pgscf(1)
      call pgsci(1)
      call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      call pgmtxt('B', 2.5, 0.5, 0.5, 'Model Number')
      call pgmtxt('L', 3.0, 0.5, 0.5, 'log\d10\u max\dk\u |residual|')
      call pgmtxt('T', 1.0, 0.5, 0.5, title)

      do i_eq = 1, h% nvar
         yvec = log10(h% vals(i_eq, 1:n))
         call pgsci(mod(i_eq - 1, 13) + 2)
         if (n == 1) then
            call pgpt(1, xvec, yvec, 17)
         else
            call pgline(n, xvec, yvec)
         end if
      end do

      call pgsch(0.75 * txt_scale)
      dy_leg = min(0.08, 1.0 / real(h% nvar))
      do i_eq = 1, h% nvar
         call pgsci(mod(i_eq - 1, 13) + 2)
         call pgmtxt('RV', 0.5, 1.0 - (real(i_eq) - 0.5) * dy_leg, 0.0, &
            trim(h% names(i_eq)))
      end do

      call pgunsa
      deallocate(xvec, yvec)

   end subroutine do_Max_eq_resid_plot


   subroutine reset_hist(h, nvar, capacity, ierr)
      type(resid_hist_t), intent(inout) :: h
      integer, intent(in) :: nvar, capacity
      integer, intent(out) :: ierr

      if (allocated(h% model)) deallocate(h% model)
      if (allocated(h% vals)) deallocate(h% vals)
      if (allocated(h% names)) deallocate(h% names)
      allocate(h% model(capacity), h% vals(nvar, capacity), h% names(nvar), stat=ierr)
      if (ierr /= 0) then
         write(*, *) 'Max_eq_resid: allocate failed'
         return
      end if
      h% nvar = nvar
      h% count = 0
   end subroutine reset_hist


   subroutine grow_hist(h, new_capacity, ierr)
      type(resid_hist_t), intent(inout) :: h
      integer, intent(in) :: new_capacity
      integer, intent(out) :: ierr
      integer, allocatable :: tmp_model(:)
      real, allocatable :: tmp_vals(:, :)

      call move_alloc(h% model, tmp_model)
      call move_alloc(h% vals, tmp_vals)
      allocate(h% model(new_capacity), h% vals(h% nvar, new_capacity), stat=ierr)
      if (ierr /= 0) then
         write(*, *) 'Max_eq_resid: allocate failed'
         return
      end if
      h% model(1:h% count) = tmp_model(1:h% count)
      h% vals(:, 1:h% count) = tmp_vals(:, 1:h% count)
   end subroutine grow_hist

end module pgstar_equation_residuals

