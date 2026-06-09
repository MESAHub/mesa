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

   use star_private_def
   use const_def, only: dp
   use pgstar_support
   use star_pgstar

   implicit none

contains

   subroutine Max_eq_resid_plot(id, device_id, ierr)

      integer, intent(in)  :: id
      integer, intent(in)  :: device_id
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
      if (ierr /= 0) return

      call pgebuf()

   end subroutine Max_eq_resid_plot



   ! history of residuals for each structure equation
   subroutine do_Max_eq_resid_plot(s, id, &
         win_xleft, win_xright, win_ybot, win_ytop, subplot, title, txt_scale, &
         max_width, ierr)

      integer, intent(in) :: id, max_width
      logical, intent(in) :: subplot
      real, intent(in) :: &
         win_xleft, win_xright, win_ybot, win_ytop, txt_scale
      character (len=*), intent(in) :: title
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      integer :: i_eq, initial_width
      integer :: n
      integer, save :: resid_hist_nvar
      integer, save :: init_model
      real :: xmin, xmax
      real :: ymin, ymax

      ! use save, we set these once, and they are reused in later calls of this function
      real, allocatable, dimension(:), save :: xvec, yvec, resid_hist_model
      real, allocatable, dimension(:, :), save :: resid_hist_vals
      character(len=strlen), allocatable, dimension(:), save :: resid_equ_names

      ierr = 0

      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (.not. allocated(resid_hist_model)) then
         init_model = s% model_number - 1
      end if

      if (max_width < 0) then
         n = s% model_number - init_model
      else
         if (s% model_number < max_width) then
            n = s% model_number - init_model
         else
            n = max_width
         end if
      end if

      if (.not. allocated(resid_hist_model)) then
         resid_hist_nvar = s% nvar_hydro
         if (max_width < 0) then
            initial_width = 10
         else
            initial_width = max_width
         end if
         allocate(resid_hist_model(initial_width))
         allocate(resid_hist_vals(resid_hist_nvar, initial_width))
         allocate(resid_equ_names(resid_hist_nvar))
         allocate(xvec(initial_width), yvec(initial_width), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed for PGSTAR'
            return
         end if

         do i_eq = 1, resid_hist_nvar
            resid_equ_names(i_eq) = trim(s% nameofequ(i_eq))
         end do

      else if (s% model_number - init_model >= size(resid_hist_model)) then
         if (max_width < 0) then
            ! grow the arrays
            call realloc_resid_hist(2*size(resid_hist_model))
         end if
      end if

      if (max_width > 1 .and. s% model_number > max_width) then  ! displace values back one step
         resid_hist_model(1:n-1) = resid_hist_model(2:n)
         resid_hist_vals(:,1:n-1) = resid_hist_vals(:,2:n)
      end if

      resid_hist_model(n) = s% model_number

      ! get the max resid values for this step
      do i_eq = 1, resid_hist_nvar
         resid_hist_vals(i_eq, n) = &
            real(max(1d-40, maxval(abs(s% equ(i_eq, 1:s% nz)))))
      end do

      xvec(1:n) = real(resid_hist_model(1:n))

      xmin = xvec(1)
      xmax = xvec(n)

      if (xmin >= xmax) xmax = xmin + 1.0

      ymin = minval(resid_hist_vals(:, 1:n))
      ymax = maxval(resid_hist_vals(:, 1:n))
      ymin = max(ymin, 1e-20)
      ymax = max(ymax, 10.0*ymin)

      call pgsave
      call pgsci(0)                       ! color index 0 = background
      call pgsvp(win_xleft, win_xright, win_ybot, win_ytop)
      call pgrect(0.0, 1.0, 0.0, 1.0)   ! pseudo erase
      call pgswin(xmin, xmax, log10(ymin), log10(ymax))
      call pgsch(txt_scale)
      call pgscf(1)
      call pgsci(1)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pgmtxt('B',2.5,0.5,0.5,'Model Number')
      call pgmtxt('L',3.0,0.5,0.5,'log10(max |residual|)')
      call pgmtxt('T',1.0,0.5,0.5, title)

      ! draw lines
      do i_eq = 1, resid_hist_nvar
         yvec(1:n) = log10(resid_hist_vals(i_eq, 1:n))
         call pgsci(mod(i_eq-1,13) + 2)
         call pgline(n, xvec(1:n), yvec(1:n))
      end do

      ! legend
      call pgsch(0.75 * txt_scale)

      do i_eq = 1, resid_hist_nvar
         call pgsci(mod(i_eq-1,13)+2)
         call pgptxt( &
            xmax + 0.03, &
            log10(ymax) - (0.03 + 0.06*real(i_eq-1))*(log10(ymax)-log10(ymin)), &
            0.0, 0.0, &
            trim(resid_equ_names(i_eq)))
      end do

      call pgsci(1)
      call pgunsa

      contains

      subroutine realloc_resid_hist(new_size)
         integer, intent(in) :: new_size
         real, allocatable :: tmp_model(:)
         real, allocatable :: tmp_vals(:,:)

         call move_alloc(resid_hist_model, tmp_model)
         call move_alloc(resid_hist_vals,  tmp_vals)

         allocate(resid_hist_model(new_size))
         allocate(resid_hist_vals(resid_hist_nvar, new_size))

         deallocate(xvec, yvec)
         allocate(xvec(new_size), yvec(new_size))

         resid_hist_model(1:n) = tmp_model(1:n)
         resid_hist_vals(:,1:n) = tmp_vals(:,1:n)

         deallocate(tmp_model, tmp_vals)
      end subroutine

   end subroutine do_Max_eq_resid_plot

end module pgstar_equation_residuals
