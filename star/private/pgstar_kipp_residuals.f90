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

   use const_def, only: dp, Msun, Rsun
   use star_private_def
   use star_pgstar
   use pgstar_colors

   implicit none

   private
   public :: Kipp_residuals_plot, do_Kipp_residuals_plot

   real(dp), parameter :: resid_floor = 1d-20

   ! Grid points with no data (e.g, because of mass loss)
   ! are not drawn and show the background color.
   integer, parameter :: c1 = colormap_offset + 1
   integer, parameter :: c2 = colormap_offset + colormap_length

contains

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

      call do_Kipp_residuals_plot(s, id, &
         s% pg% Kipp_residuals_xleft, s% pg% Kipp_residuals_xright, &
         s% pg% Kipp_residuals_ybot, s% pg% Kipp_residuals_ytop, .false., &
         s% pg% Kipp_residuals_title, s% pg% Kipp_residuals_txt_scale, &
         s% pg% Kipp_residuals_max_width, s% pg% Kipp_residuals_yaxis_name, ierr)

      call pgebuf()

   end subroutine Kipp_residuals_plot


   subroutine do_Kipp_residuals_plot(s, id, &
         win_xleft, win_xright, win_ybot, win_ytop, subplot, title, txt_scale, &
         max_width, yaxis_name, ierr)
      type (star_info), pointer :: s
      integer, intent(in)  :: id, max_width
      real, intent(in) :: win_xleft, win_xright, win_ybot, win_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title, yaxis_name
      integer, intent(out) :: ierr

      integer :: k, neq, initial_width, n, x_size, y_size, x_new
      integer, save :: init_model

      real, allocatable, dimension(:, :), save :: nr_resid_buf, nr_ycoord_buf
      integer, allocatable, dimension(:), save :: nr_model_buf, nr_zone_buf

      real, allocatable :: tmp_resid_buf(:, :), tmp_ycoord_buf(:, :)
      integer, allocatable :: tmp_model_buf(:), tmp_zone_buf(:)

      ierr = 0
      neq = size(s% equ, 1)  ! number of equations per cell

      if (.not. allocated(nr_resid_buf)) then
         if (max_width < 0) then
            initial_width = 10
         else
            initial_width = max_width
         end if
         ! 1.5 times the current nz; grown below if remeshing increases nz
         allocate(nr_resid_buf(initial_width, int(1.5 * s% nz)))
         allocate(nr_ycoord_buf(initial_width, int(1.5 * s% nz)))
         allocate(nr_model_buf(initial_width))
         allocate(nr_zone_buf(initial_width))
         nr_resid_buf = 0.0
         nr_ycoord_buf = 0.0
         nr_model_buf = 0
         nr_zone_buf = 0   ! 0 zones = empty column, skipped when plotting

         init_model = s% model_number - 1
      else if (s% nz > 0.95 * size(nr_resid_buf, 2)) then
         ! remeshing increased nz beyond 95% of the current array size
         x_size = size(nr_resid_buf, 1)
         y_size = size(nr_resid_buf, 2)
         call move_alloc(nr_resid_buf, tmp_resid_buf)
         call move_alloc(nr_ycoord_buf, tmp_ycoord_buf)

         allocate(nr_resid_buf(x_size, 5 * s% nz))
         allocate(nr_ycoord_buf(x_size, 5 * s% nz))
         nr_resid_buf = 0.0
         nr_ycoord_buf = 0.0

         nr_resid_buf(:, 1:y_size) = tmp_resid_buf(:, 1:y_size)
         nr_ycoord_buf(:, 1:y_size) = tmp_ycoord_buf(:, 1:y_size)

         deallocate(tmp_resid_buf, tmp_ycoord_buf)
      end if

      if (max_width < 0) then
         n = s% model_number - init_model
      else
         n = min(s% model_number - init_model, max_width)
      end if
      if (n < 1) return

      ! grow dimension
      if (n > size(nr_resid_buf, 1)) then
         x_size = size(nr_resid_buf, 1)
         y_size = size(nr_resid_buf, 2)
         x_new = max(2 * x_size, n)

         call move_alloc(nr_resid_buf, tmp_resid_buf)
         call move_alloc(nr_ycoord_buf, tmp_ycoord_buf)
         call move_alloc(nr_zone_buf, tmp_zone_buf)
         call move_alloc(nr_model_buf, tmp_model_buf)

         allocate(nr_resid_buf(x_new, y_size))
         allocate(nr_ycoord_buf(x_new, y_size))
         allocate(nr_zone_buf(x_new))
         allocate(nr_model_buf(x_new))

         nr_resid_buf = 0.0
         nr_ycoord_buf = 0.0
         nr_zone_buf = 0
         nr_model_buf = 0

         nr_resid_buf(1:x_size, :) = tmp_resid_buf(1:x_size, :)
         nr_ycoord_buf(1:x_size, :) = tmp_ycoord_buf(1:x_size, :)
         nr_zone_buf(1:x_size) = tmp_zone_buf(1:x_size)
         nr_model_buf(1:x_size) = tmp_model_buf(1:x_size)

         deallocate(tmp_resid_buf, tmp_ycoord_buf, tmp_zone_buf, tmp_model_buf)
      end if

      if (max_width > 1 .and. s% model_number - init_model > max_width) then
         nr_model_buf(1:max_width-1) = nr_model_buf(2:max_width)
         nr_resid_buf(1:max_width-1, :) = nr_resid_buf(2:max_width, :)
         nr_ycoord_buf(1:max_width-1, :) = nr_ycoord_buf(2:max_width, :)
         nr_zone_buf(1:max_width-1) = nr_zone_buf(2:max_width)
      end if

      nr_model_buf(n) = s% model_number
      nr_zone_buf(n) = s% nz

      ! store the current model; all y coordinates are cell centered
      do k = 1, s% nz
         nr_resid_buf(n, k) = real(log10(max(resid_floor, &
            maxval(abs(s% equ(1:neq, k))))))
         select case (trim(yaxis_name))
         case ('logR')
            nr_ycoord_buf(n, k) = real(log10(max(1d-99, s% rmid(k)) / Rsun))
         case ('tau')
            nr_ycoord_buf(n, k) = real(log10(max(1d-99, s% tau(k))))
         case default  ! 'mass'
            nr_ycoord_buf(n, k) = real((s% m(k) - 0.5d0 * s% dm(k)) / Msun)
         end select
      end do

      call pgsave
      call pgsci(0)                       ! color index 0 = background
      call pgsvp(win_xleft, win_xright, win_ybot, win_ytop)
      call pgrect(0.0, 1.0, 0.0, 1.0)     ! pseudo erase
      call render(n, nr_model_buf, nr_zone_buf, nr_resid_buf, nr_ycoord_buf, &
         max(2, s% nz), yaxis_name, title, txt_scale)
      call pgunsa

   end subroutine do_Kipp_residuals_plot


   subroutine render(nx, model_buf, zone_buf, resid_buf, ycoord_buf, &
         ny, yaxis_name, title, txt_scale)
      integer, intent(in) :: nx, ny
      integer, intent(in) :: model_buf(:), zone_buf(:)
      real, intent(in) :: resid_buf(:, :), ycoord_buf(:, :)
      character(len=*), intent(in) :: yaxis_name, title
      real, intent(in) :: txt_scale

      real, allocatable :: col(:, :), coord_grid(:)
      real :: tr(6), fg, bg, missing, xlo, xhi, xl, xr
      real :: clo, chi, dcoord, resid_hi, resid_lo
      integer :: i, j, j1, j2, nz_i, ip, cir_lo_save, cir_hi_save
      logical :: y_reverse, have_data

      ! --- data ranges over all stored columns ---
      have_data = .false.
      chi = -huge(1.0); clo = huge(1.0)
      resid_hi = -huge(1.0); resid_lo = huge(1.0)
      do i = 1, nx
         nz_i = zone_buf(i)
         if (nz_i < 1) cycle
         have_data = .true.
         chi = max(chi, maxval(ycoord_buf(i, 1:nz_i)))
         clo = min(clo, minval(ycoord_buf(i, 1:nz_i)))
         resid_hi = max(resid_hi, maxval(resid_buf(i, 1:nz_i)))
         resid_lo = min(resid_lo, minval(resid_buf(i, 1:nz_i)))
      end do
      if (.not. have_data) return
      if (chi <= clo) chi = clo + 1.0

      y_reverse = (trim(yaxis_name) == 'tau')

      ! uniform y grid, ascending with j; axis reversal is done in pgswin
      allocate(coord_grid(ny), col(1, ny))
      dcoord = (chi - clo) / real(ny - 1)
      do j = 1, ny
         coord_grid(j) = clo + real(j - 1) * dcoord
      end do

      ! --- color scale ---
      fg = resid_hi
      bg = resid_lo
      if (fg - bg < 1.0) bg = fg - 1.0
      missing = bg - 1.0   ! marker only, never drawn

      call pgqcir(cir_lo_save, cir_hi_save)
      call pgscir(c1, c2)

      ! --- x range from stored model numbers ---
      xlo = huge(1.0); xhi = -huge(1.0)
      ip = 0
      do i = 1, nx
         if (zone_buf(i) < 1) cycle
         call column_edges(i, ip, nx, model_buf, zone_buf, xl, xr)
         xlo = min(xlo, xl)
         xhi = max(xhi, xr)
         ip = i
      end do

      if (.not. y_reverse) then
         call pgswin(xlo, xhi, clo - 0.5*dcoord, chi + 0.5*dcoord)
      else
         call pgswin(xlo, xhi, chi + 0.5*dcoord, clo - 0.5*dcoord)
      end if

      ! --- one pgimag call per stored model, placed at its model number ---
      tr(3) = 0.0
      tr(4) = clo - dcoord   ! j = 1 -> clo
      tr(5) = 0.0
      tr(6) = dcoord
      ip = 0
      do i = 1, nx
         nz_i = zone_buf(i)
         if (nz_i < 1) cycle
         call column_edges(i, ip, nx, model_buf, zone_buf, xl, xr)
         ip = i
         call interp_column(nz_i, ycoord_buf(i, 1:nz_i), resid_buf(i, 1:nz_i), &
            ny, coord_grid, dcoord, bg, fg, missing, col(1, :))
         tr(2) = xr - xl
         tr(1) = 0.5*(xl + xr) - tr(2)   ! array index 1 -> column center
         ! the column coordinate is monotone, so missing points can only be
         ! at the ends; draw just the contiguous range with data
         j1 = 1
         do while (j1 <= ny)
            if (col(1, j1) > missing) exit
            j1 = j1 + 1
         end do
         j2 = ny
         do while (j2 >= j1)
            if (col(1, j2) > missing) exit
            j2 = j2 - 1
         end do
         if (j1 > j2) cycle
         call pgimag(col, 1, ny, 1, 1, j1, j2, bg, fg, tr)
      end do

      ! --- axes and labels ---
      call pgsci(clr_Foreground)
      call pgsch(txt_scale)
      call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      call pgmtxt('B', 2.5, 0.5, 0.5, 'Model Number')
      select case (trim(yaxis_name))
      case ('logR')
         call pgmtxt('L', 3.0, 0.5, 0.5, 'log(R/R\d\(2281)\u)')
      case ('tau')
         call pgmtxt('L', 3.0, 0.5, 0.5, 'log(\gt)')
      case default
         call pgmtxt('L', 3.0, 0.5, 0.5, 'm/M\d\(2281)\u')
      end select
      call pgmtxt('T', 1.0, 0.5, 0.5, trim(title))

      ! --- color wedge ---
      call pgwedg('RI', 0.5, 3.0, bg, fg, 'log\d10\u max |res|')

      ! restore the index range; no color representations were modified
      call pgscir(cir_lo_save, cir_hi_save)

      deallocate(col, coord_grid)

   end subroutine render


   ! x extent of column i: halfway to the neighboring non-empty columns
   subroutine column_edges(i, iprev, nx, model_buf, zone_buf, xl, xr)
      integer, intent(in) :: i, iprev, nx
      integer, intent(in) :: model_buf(:), zone_buf(:)
      real, intent(out) :: xl, xr
      integer :: inext, j
      real :: xc

      xc = real(model_buf(i))
      inext = 0
      do j = i + 1, nx
         if (zone_buf(j) > 0) then
            inext = j
            exit
         end if
      end do

      if (iprev > 0) then
         xl = 0.5 * (real(model_buf(iprev)) + xc)
      else if (inext > 0) then
         xl = xc - 0.5 * (real(model_buf(inext)) - xc)
      else
         xl = xc - 0.5
      end if
      if (inext > 0) then
         xr = 0.5 * (real(model_buf(inext)) + xc)
      else if (iprev > 0) then
         xr = xc + 0.5 * (xc - real(model_buf(iprev)))
      else
         xr = xc + 0.5
      end if
   end subroutine column_edges

   subroutine interp_column(m, ycol, rcol, ny, coord_grid, dcoord, bg, fg, missing, out)
      integer, intent(in) :: m, ny
      real, intent(in) :: ycol(m), rcol(m), coord_grid(ny), dcoord, bg, fg, missing
      real, intent(out) :: out(ny)
      real :: yc(m), rv(m), y, v, dy
      integer :: j, kk

      if (ycol(1) <= ycol(m)) then
         yc = ycol; rv = rcol
      else
         yc = ycol(m:1:-1); rv = rcol(m:1:-1)
      end if

      kk = 1
      do j = 1, ny
         y = coord_grid(j)
         if (y < yc(1)) then
            ! allow up to one grid step beyond the outermost cell center
            if (y < yc(1) - dcoord) then
               out(j) = missing
               cycle
            end if
            v = rv(1)
         else if (y > yc(m)) then
            if (y > yc(m) + dcoord) then
               out(j) = missing
               cycle
            end if
            v = rv(m)
         else if (m == 1) then
            v = rv(1)
         else
            do while (kk < m - 1)
               if (yc(kk + 1) >= y) exit
               kk = kk + 1
            end do
            dy = yc(kk + 1) - yc(kk)
            if (dy > 0.0) then
               v = rv(kk) + (y - yc(kk)) / dy * (rv(kk + 1) - rv(kk))
            else
               v = 0.5 * (rv(kk) + rv(kk + 1))
            end if
         end if
         out(j) = min(fg, max(bg, v))
      end do
   end subroutine interp_column

end module pgstar_kipp_residuals
