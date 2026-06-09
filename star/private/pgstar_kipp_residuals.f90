! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module pgstar_kipp_residuals

   use const_def
   use star_private_def
   use star_pgstar
   use pgstar_colors
   use const_def, only: dp, Msun

   implicit none

contains

   ! Kippenhahn-style log10(max(|res|)) across all equations for a
   ! given space-time location
   subroutine Kipp_residuals_plot(id, device_id, ierr)
      integer, intent(in)  :: id, device_id
      integer, intent(out) :: ierr

      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call do_Kipp_residuals_Plot(s, id, &
         s% pg% Kipp_residuals_xleft, s% pg% Kipp_residuals_xright, &
         s% pg% Kipp_residuals_ybot, s% pg% Kipp_residuals_ytop, .false., &
         s% pg% Kipp_residuals_title, s% pg% Kipp_residuals_txt_scale, &
         s% pg% Kipp_residuals_max_width, s% pg% Kipp_residuals_yaxis_name, ierr)
      if (ierr /= 0) return

      call pgebuf()

   end subroutine Kipp_residuals_plot


   subroutine do_Kipp_residuals_plot(s, id, &
         win_xleft, win_xright, win_ybot, win_ytop, subplot, title, txt_scale, &
         max_width, yaxis_name, &
         ierr)
      integer, intent(in)  :: id, max_width
      real, intent(in) :: win_xleft, win_xright, win_ybot, win_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title, yaxis_name
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      integer :: k, neq, initial_width, n, x_size, y_size, nr_n_cells
      integer, save :: init_model

      real :: resid_hi, resid_lo

      real, allocatable, dimension(:, :), save :: nr_resid_buf, nr_ycoord_buf
      integer, allocatable, dimension(:), save :: nr_model_buf, nr_zone_buf

      real, allocatable :: tmp_resid_buf(:, :), tmp_ycoord_buf(:, :)
      integer, allocatable :: tmp_model_buf(:), tmp_zone_buf(:)

      neq = size(s% equ, 1) ! number of equations per cell

      if (.not. allocated(nr_resid_buf)) then
         if (max_width < 0) then
            initial_width = 10
         else
            initial_width = max_width
         end if
         ! make 5 times bigger than current nz
         allocate(nr_resid_buf(initial_width, int(1.5 * s% nz)))
         allocate(nr_ycoord_buf(initial_width, int(1.5 * s% nz)))
         allocate(nr_model_buf(initial_width))
         allocate(nr_zone_buf(initial_width))
         nr_resid_buf(:, :) = -99.0
         nr_ycoord_buf(:, :) = 0.0   ! zero-init so unused cells are safely masked

         init_model = s% model_number - 1
      else if (s% nz > 0.95 * size(nr_resid_buf, 2)) then
         write(*, *) "reallocating due to nz"
         ! check if remeshing has increased nz beyond 95% of the
         ! current array size
         x_size = size(nr_resid_buf, 1)
         call move_alloc(nr_resid_buf, tmp_resid_buf)
         call move_alloc(nr_ycoord_buf, tmp_ycoord_buf)

         allocate(nr_resid_buf(x_size, 5 * s% nz))
         allocate(nr_ycoord_buf(x_size, 5 * s% nz))
         nr_resid_buf(:, :) = -99.0
         nr_ycoord_buf(:, :) = 0.0   ! zero-init so unused cells are safely masked

         nr_resid_buf(:, 1:nr_n_cells) = tmp_resid_buf(:, 1:nr_n_cells)
         nr_ycoord_buf(:, 1:nr_n_cells)  = tmp_ycoord_buf(:, 1:nr_n_cells)

         deallocate(tmp_resid_buf, tmp_ycoord_buf)
      end if

      nr_n_cells = s% nz
      if (max_width < 0) then
         n = s% model_number - init_model
      else
         n = min(s% model_number - init_model, max_width)
      end if

      ! resize n_model dimension
      if (n > size(nr_resid_buf, 1)) then
         write(*, *) "reallocating due to model_number"
         x_size = size(nr_resid_buf, 1)
         y_size = size(nr_resid_buf, 2)

         call move_alloc(nr_resid_buf, tmp_resid_buf)
         call move_alloc(nr_ycoord_buf, tmp_ycoord_buf)
         call move_alloc(nr_zone_buf, tmp_zone_buf)
         call move_alloc(nr_model_buf, tmp_model_buf)

         allocate(nr_resid_buf(2 * x_size, y_size))
         allocate(nr_ycoord_buf(2 * x_size, y_size))
         allocate(nr_zone_buf(2 * x_size))
         allocate(nr_model_buf(2 * x_size))

         nr_resid_buf(:, :) = -99.0
         nr_ycoord_buf(:, :) = 0.0   ! zero-init so unused cells are safely masked

         nr_resid_buf(1:n, :) = tmp_resid_buf(1:n, :)
         nr_ycoord_buf(1:n, :) = tmp_ycoord_buf(1:n, :)
         nr_zone_buf(1:n) = tmp_zone_buf(1:n)
         nr_model_buf(1:n) = tmp_model_buf(1:n)

         deallocate(tmp_resid_buf, tmp_ycoord_buf, tmp_zone_buf, tmp_model_buf)
      end if

      if (max_width > 1 .and. s% model_number - init_model > max_width) then
         nr_model_buf(1:max_width-1) = nr_model_buf(2:max_width)
         nr_resid_buf(1:max_width-1, :) = nr_resid_buf(2:max_width, :)
         nr_ycoord_buf(1:max_width-1, :) = nr_ycoord_buf(2:max_width, :)
      end if

      nr_model_buf(n) = s% model_number
      nr_zone_buf(n) = s% nz

      ! fill arrays with new vals
      do k = 1, s% nz
         nr_resid_buf(n, k) = real(safe_log10(maxval(abs(s% equ(1:neq, k)))))
         select case (trim(yaxis_name))
         case ('mass')
            nr_ycoord_buf(n, k) = real(s% m(k)/Msun)
         case ('logR')
            ! cell-centered radius in log10(R/Rsun)
            if (k < s%nz) then
               nr_ycoord_buf(n, k) = real(safe_log10(0.5d0 * (s% r(k) + s% r(k+1)) / Rsun))
            else
               nr_ycoord_buf(n, k) = real(safe_log10(0.5d0 * s% r(k) / Rsun))
            end if
         case ('tau')
            ! optical depth
            nr_ycoord_buf(n, k) = real(safe_log10(s% tau(k)))
         case default
            ! cell-centered mass coordinate
            nr_ycoord_buf(n, k) = real((s% m(k) - 0.5d0 * s% dm(k))/Msun)
         end select
      end do

      ! so actual plotting
      call pgsave
      call pgsci(0)                       ! color index 0 = background
      call pgsvp(win_xleft, win_xright, win_ybot, win_ytop)
      call pgrect(0.0, 1.0, 0.0, 1.0)   ! pseudo erase
      call Kipp_residuals_render(ierr, id)
      call pgunsa

      contains

      subroutine Kipp_residuals_render(ierr, id)
         integer, intent(in)  :: id
         integer, intent(out) :: ierr
         real, parameter :: missing_val = TINY(1.0)
         real, allocatable :: img(:,:), coord_grid(:)
         real :: tr(6), fg, bg, xlo, xhi, dx
         real :: clo, chi, dcoord, coord_j, frac
         integer :: nx, ny, i, j, k, nz_i
         logical :: y_reverse

         ! ----------------------------------------- colorbar definition
         real :: bright, contra
         real, dimension(9) :: l, r, g, b

         l = (/ &
            0.00, 0.12, 0.25, 0.38, 0.50, &
            0.62, 0.75, 0.88, 1.00 /)

         r = (/ &
            0.00, 0.05, 0.15, 0.32, 0.55, &
            0.78, 0.93, 0.99, 1.00 /)

         g = (/ &
            0.00, 0.02, 0.04, 0.06, 0.10, &
            0.18, 0.32, 0.58, 0.95 /)

         b = (/ &
            0.00, 0.08, 0.18, 0.32, 0.36, &
            0.28, 0.16, 0.08, 0.80 /)

         contra = 1.0
         bright = 0.5

         call pgscir(16,255)
         call pgctab(l, r, g, b, 9, contra, bright)
         ! ----------------------------------------- end colorbar

         nx = n
         ny = nr_n_cells

         allocate(img(nx, ny))
         allocate(coord_grid(ny))

         chi = -1.0e30
         clo =  1.0e30
         do i = 1, nx
            nz_i = nr_zone_buf(i)
            do k = 1, nz_i
               if (nr_ycoord_buf(i,k) > chi) chi = nr_ycoord_buf(i,k)
               ! m(k) is positive; skip any zero-initialised slots
               if (nr_ycoord_buf(i,k) > 0.0 .and. nr_ycoord_buf(i,k) < clo) &
                  clo = nr_ycoord_buf(i,k)
            end do
         end do

         resid_hi = -1.0e30
         resid_lo =  1.0e30
         do i = 1, nx
            nz_i = nr_zone_buf(i)
            do k = 1, nz_i
               if (nr_resid_buf(i,k) > resid_hi) resid_hi = nr_resid_buf(i,k)
               if (nr_resid_buf(i,k) < resid_lo) resid_lo = nr_resid_buf(i,k)
            end do
         end do

         select case (yaxis_name)
         case ('mass')
            y_reverse = .false.
         case ('logR')
            y_reverse = .false.
         case ('tau')
            y_reverse = .true.
         case default
            y_reverse = .false.
         end select

         dcoord = (chi - clo) / real(ny - 1)

         do j = 1, ny
            if (.not. y_reverse) then
               coord_grid(j) = clo + real(j - 1) * dcoord
            else
               coord_grid(j) = chi - real(j - 1) * dcoord
            end if
         end do

         ! --- interpolate residuals onto the uniform grid ---
         img = missing_val
         do i = 1, nx
            nz_i = nr_zone_buf(i)
            do j = 1, ny
               coord_j = coord_grid(j)
               do k = 1, nz_i - 1
                  ! skip invalid cells
                  if (nr_resid_buf(i,k)   == missing_val) cycle
                  if (nr_resid_buf(i,k+1) == missing_val) cycle
                  ! works for either increasing or decreasing coordinates
                  if ((nr_ycoord_buf(i,k) - coord_j) * &
                     (nr_ycoord_buf(i,k+1) - coord_j) <= 0d0) then
                     if (abs(nr_ycoord_buf(i,k) - nr_ycoord_buf(i,k+1)) > 1d-99) then
                        frac = (coord_j - nr_ycoord_buf(i,k+1)) / &
                           (nr_ycoord_buf(i,k) - nr_ycoord_buf(i,k+1))
                     else
                        frac = 0.5d0
                     end if
                     img(i,j) = nr_resid_buf(i,k+1) + &
                        frac * (nr_resid_buf(i,k) - nr_resid_buf(i,k+1))
                     exit
                  end if
               end do
            end do
         end do
         ! --- x range centred on model numbers ---
         dx  = max(1.0, real(nr_model_buf(nx) - nr_model_buf(1)) / real(nx - 1))
         xlo = real(nr_model_buf(1))  - 0.5*dx
         xhi = real(nr_model_buf(nx)) + 0.5*dx

         ! --- transformation matrix for pgimag ---
         ! world_y(j) = tr(4) + tr(6)*j
         ! j=1  -> mhi  =>  tr(4) = mhi - tr(6) = mhi + dm
         ! j=ny -> mlo  =>  confirms tr(6) = -dm
         tr(1) = xlo
         tr(2) = (xhi - xlo) / real(nx)
         tr(3) = 0.0
         tr(4) = chi + dcoord
         tr(5) = 0.0
         tr(6) = -dcoord

         ! color scale
         fg = resid_hi
         bg = resid_lo
         if (fg <= -100) fg = maxval(img)
         if (bg <= -100) bg = minval(img, mask=(img > -19.9))
         if (bg >= fg)   bg = fg - 6.0

         ! --- main panel ---
         ! y-axis: mlo at bottom (centre), mhi at top (surface)
         call pgsvp(0.10, 0.82, 0.12, 0.92)
         if (.not. y_reverse) then
            call pgswin(xlo, xhi, clo - 0.5*dcoord, chi + 0.5*dcoord)
         else
            call pgswin(xlo, xhi, chi + 0.5*dcoord, clo - 0.5*dcoord)
         end if
         call pgimag(img, nx, ny, 1, nx, 1, ny, fg, bg, tr)

         call pgsci(clr_Foreground)
         call pgsch(txt_scale)
         call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

         select case (yaxis_name)
         case ('mass')
            y_reverse = .false.
            call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
               'Newton Raphson Solver Residuals of accepted step')
         case ('logR')
            y_reverse = .false.
            call pglab('Model Number', 'log(R/R\d\(2281)\u)', &
               'Newton Raphson Solver Residuals of accepted step')
         case ('tau')
            y_reverse = .true.
            call pglab('Model Number', 'log(tau)', &
               'Newton Raphson Solver Residuals of accepted step')
         case default
            y_reverse = .false.
            call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
               'Newton Raphson Solver Residuals of accepted step')
         end select

         ! --- color wedge ---
         call pgsvp(0.84, 0.90, 0.12, 0.92)
         call pgswin(0.0, 1.0, bg, fg)
         call pgwedg('RI', 0.0, 4.0, bg, fg, 'log10(max(|res|))')

         deallocate(img, coord_grid)

      end subroutine Kipp_residuals_render

   end subroutine do_Kipp_residuals_plot



end module pgstar_kipp_residuals
