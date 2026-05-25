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

   use star_pgstar
   use pgstar_colors
   use const_def, only: dp, Msun, Rsun

   implicit none

   integer, parameter :: NR_MAX_MODELS = 2000
   integer, allocatable, save :: nr_zone_buf(:)
   integer, allocatable, save :: nr_model_buf(:)
   integer,              save :: nr_n_stored = 0, nr_n_cells = 0
   real, allocatable, save :: nr_resid_buf(:,:), tmp_resid_buf(:,:)
   real, allocatable, save :: nr_ycoord_buf(:,:), tmp_mass_buf(:,:)
   real,              save :: nr_resid_min = -101.0, nr_resid_max = -101.0  ! auto-scale
   ! TODO: define pgstar variable to substitute x_integer_ctrl(1) used here
   integer, parameter :: NR_COORD_MASS   = 1
   integer, parameter :: NR_COORD_LOGR   = 2
   integer, parameter :: NR_COORD_TAU    = 3

contains

   ! -----------------------------------------------------------------------
   ! Top-level plot callback: called by do_pgstar_win_file via p%plot.
   ! Manages the history buffer then delegates all drawing to the render
   ! routine.  The PGPLOT device framing (pgslct/pgbbuf/pgeras/pgebuf) is
   ! the standard wrapper used across all pgstar_*.f90 plot callbacks.
   ! -----------------------------------------------------------------------
   subroutine Kipp_residuals_plot(id, device_id, ierr)
      integer, intent(in)  :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      integer  :: k, neq, col

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      neq = size(s%equ, 1) ! number of equations per cell

      ! --- allocate history buffer on first call ---
      if (.not. allocated(nr_resid_buf)) then
         ! make 5 times bigger than current nz
         allocate(nr_resid_buf(NR_MAX_MODELS, 5 * s%nz))
         allocate(nr_ycoord_buf(NR_MAX_MODELS, 5 * s%nz))
         allocate(nr_model_buf(NR_MAX_MODELS))
         allocate(nr_zone_buf(NR_MAX_MODELS))
         nr_resid_buf = -99.0
         nr_ycoord_buf  =  0.0   ! zero-init so unused cells are safely masked
         nr_n_cells   = s%nz
      else if (s%nz > 0.95 * size(nr_resid_buf, 2)) then
         ! check if remeshing has increased s%nz beyond 95% of the
         ! current array size, reallocate
         allocate(tmp_resid_buf(NR_MAX_MODELS, 5 * s%nz))
         allocate(tmp_mass_buf(NR_MAX_MODELS, 5  *s%nz))

         tmp_resid_buf = -99.0
         tmp_mass_buf  = 0.0

         ! Copy older compressed history into the newly expanded arrays
         tmp_resid_buf(:, 1:nr_n_cells) = nr_resid_buf(:, 1:nr_n_cells)
         tmp_mass_buf(:, 1:nr_n_cells)  = nr_ycoord_buf(:, 1:nr_n_cells)

         call move_alloc(tmp_resid_buf, nr_resid_buf)
         call move_alloc(tmp_mass_buf, nr_ycoord_buf)
         nr_n_cells = s%nz
      end if


      ! --- roll buffer if full ---
      if (nr_n_stored < NR_MAX_MODELS) then
         nr_n_stored = nr_n_stored + 1
         col = nr_n_stored
      else
         nr_model_buf(1:NR_MAX_MODELS-1)      = nr_model_buf(2:NR_MAX_MODELS)
         nr_resid_buf(1:NR_MAX_MODELS-1, :)   = nr_resid_buf(2:NR_MAX_MODELS, :)
         nr_ycoord_buf(1:NR_MAX_MODELS-1, :)    = nr_ycoord_buf(2:NR_MAX_MODELS, :)
         col = NR_MAX_MODELS
      end if

      nr_model_buf(col) = s%model_number
      nr_zone_buf(col) = s%nz

      do k = 1, s%nz
         nr_resid_buf(col, k) = &
            safe_log10(maxval(abs(s%equ(1:neq, k))))
         select case (s%x_integer_ctrl(1))
         case (NR_COORD_MASS)
            nr_ycoord_buf(col, k) = s%m(k)/Msun
         case (NR_COORD_LOGR)
            ! cell-centered radius in log10(R/Rsun)
            if (k < s%nz) then
               nr_ycoord_buf(col, k) = safe_log10(0.5d0 * (s%r(k) + s%r(k+1)) / Rsun)
            else
               nr_ycoord_buf(col, k) = safe_log10(0.5d0*s%r(k)/Rsun)
            end if
         case (NR_COORD_TAU)
            ! log10 optical depth
            nr_ycoord_buf(col, k) = safe_log10(s%tau(k))
         case default
            ! cell-centered mass coordinate
            nr_ycoord_buf(col, k) = (s%m(k) - 0.5d0*s%dm(k))/Msun
         end select
      end do


      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call Kipp_residuals_render(ierr, id)
      call pgebuf()

   end subroutine Kipp_residuals_plot

   ! -----------------------------------------------------------------------
   ! Pure rendering routine: builds the 2-D image from the history buffer
   ! and draws it using pgimag + a color wedge.  Not called directly by
   ! the pgstar harness; invoked only from Kipp_residuals_plot above.
   ! This subroutine is intentionally left unchanged from the original.
   ! -----------------------------------------------------------------------
   subroutine Kipp_residuals_render(ierr, id)
      integer, intent(in)  :: id
      integer, intent(out) :: ierr
      real, parameter :: missing_val = 1d-30
      real, allocatable :: img(:,:), coord_grid(:)
      real :: tr(6), fg, bg, xlo, xhi, dx
      real :: clo, chi, dcoord, coord_j, frac
      integer :: nx, ny, i, j, k, nz_i
      logical :: y_reverse
      type(star_info), pointer :: s

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

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (nr_n_stored < 2) return

      nx = nr_n_stored
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

      select case (s%x_integer_ctrl(1))
      case (NR_COORD_MASS)
         y_reverse = .false.
      case (NR_COORD_LOGR)
         y_reverse = .false.
      case (NR_COORD_TAU)
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
               if ( (nr_ycoord_buf(i,k)   - coord_j) * &
                  (nr_ycoord_buf(i,k+1) - coord_j) <= 0d0 ) then
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

      ! color scale: use inlist limits; negative sentinel means auto-scale
      fg = nr_resid_max
      bg = nr_resid_min
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
      call pgsch(1.2)
      call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      select case (s%x_integer_ctrl(1))
      case (NR_COORD_MASS)
         y_reverse = .false.
         call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
            'Newton Raphson Solver Residuals of accepted step')
      case (NR_COORD_LOGR)
         y_reverse = .false.
         call pglab('Model Number', 'log(R/R\d\(2281)\u)', &
            'Newton Raphson Solver Residuals of accepted step')
      case (NR_COORD_TAU)
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

      deallocate(img)
      deallocate(coord_grid)

   end subroutine Kipp_residuals_render

   ! -----------------------------------------------------------------------
   ! Registration routine: called once per step by the pgstar driver.
   ! Mirrors the do_*_plot pattern used throughout star/private/pgstar_*.f90:
   !   1. obtain the plot slot via s%pg%pgstar_win_file_ptr(i_Kipp_residuals)
   !   2. populate p%* from s%pg inlist controls (set once on first call,
   !      re-read every step so inlist changes take effect at runtime)
   !   3. hand off to do_pgstar_win_file which owns the window / file logic
   ! -----------------------------------------------------------------------
   subroutine do_Kipp_residuals_plot(id, ierr)
      integer, intent(in)  :: id
      integer, intent(out) :: ierr

      type(pgstar_win_file_data), pointer :: p
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      p => s%pg%pgstar_win_file_ptr(i_Kipp_residuals)

      p%plot => Kipp_residuals_plot
      p%id   =  i_Kipp_residuals
      p%name = 'Kipp_residuals'

      ! --- display settings read from s%pg so inlist edits take effect live ---
      p%win_flag          = s%pg%Kipp_residuals_win_flag
      p%win_width         = s%pg%Kipp_residuals_win_width
      p%win_aspect_ratio  = s%pg%Kipp_residuals_win_aspect_ratio

      p%file_flag         = s%pg%Kipp_residuals_file_flag
      p%file_dir          = s%pg%Kipp_residuals_file_dir
      p%file_prefix       = s%pg%Kipp_residuals_file_prefix
      p%file_interval     = s%pg%Kipp_residuals_file_interval
      p%file_width        = s%pg%Kipp_residuals_file_width
      p%file_aspect_ratio = s%pg%Kipp_residuals_file_aspect_ratio

      ! Sync color-scale limits into module-level vars used by the render routine.
      ! Negative sentinel (-101) means auto-scale (see Kipp_residuals_render).
      nr_resid_min = s%pg%Kipp_residuals_min
      nr_resid_max = s%pg%Kipp_residuals_max

      call do_pgstar_win_file(s, p, ierr)

   end subroutine do_Kipp_residuals_plot

end module pgstar_kipp_residuals


! ! ***********************************************************************
! !
! !   Copyright (C) 2010-2019  The MESA Team
! !
! !   This program is free software: you can redistribute it and/or modify
! !   it under the terms of the GNU Lesser General Public License
! !   as published by the Free Software Foundation,
! !   either version 3 of the License, or (at your option) any later version.
! !
! !   This program is distributed in the hope that it will be useful,
! !   but WITHOUT ANY WARRANTY; without even the implied warranty of
! !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! !   See the GNU Lesser General Public License for more details.
! !
! !   You should have received a copy of the GNU Lesser General Public License
! !   along with this program. If not, see <https://www.gnu.org/licenses/>.
! !
! ! ***********************************************************************

! module pgstar_kipp_residuals

!    use star_pgstar
!    use pgstar_colors
!    use const_def, only: dp, Msun

!    implicit none

! contains

!    ! Kippenhahn-style log10(max(|res|)) across all equations for a
!    ! given space-time location
!    subroutine Kipp_residuals_plot(id, device_id, ierr)
!       integer, intent(in)  :: id, device_id
!       integer, intent(out) :: ierr
!       type(star_info), pointer :: s
!       integer  :: k, neq, col

!       ierr = 0
!       call star_ptr(id, s, ierr)
!       if (ierr /= 0) return

!       neq = size(s%equ, 1) ! number of equations per cell

!       ! --- allocate history buffer on first call ---
!       if (.not. allocated(nr_resid_buf)) then
!          ! make 5 times bigger than current nz
!          allocate(nr_resid_buf(NR_MAX_MODELS, 5 * s%nz))
!          allocate(nr_ycoord_buf(NR_MAX_MODELS, 5 * s%nz))
!          allocate(nr_model_buf(NR_MAX_MODELS))
!          allocate(nr_zone_buf(NR_MAX_MODELS))
!          nr_resid_buf = -99.0
!          nr_ycoord_buf  =  0.0   ! zero-init so unused cells are safely masked
!          nr_n_cells   = s%nz
!       else if (s%nz > 0.95 * size(nr_resid_buf, 2)) then
!          ! check if remeshing has increased nz beyond 95% of the
!          ! current array size, reallocate
!          allocate(tmp_resid_buf(NR_MAX_MODELS, 5 * s%nz))
!          allocate(tmp_mass_buf(NR_MAX_MODELS, 5  *s%nz))

!          tmp_resid_buf = -99.0
!          tmp_mass_buf  = 0.0

!          ! Copy older compressed history into the newly expanded arrays
!          tmp_resid_buf(:, 1:nr_n_cells) = nr_resid_buf(:, 1:nr_n_cells)
!          tmp_mass_buf(:, 1:nr_n_cells)  = nr_ycoord_buf(:, 1:nr_n_cells)

!          call move_alloc(tmp_resid_buf, nr_resid_buf)
!          call move_alloc(tmp_mass_buf, nr_ycoord_buf)
!          nr_n_cells = s%nz
!       end if


!       ! --- roll buffer if full ---
!       if (nr_n_stored < NR_MAX_MODELS) then
!          nr_n_stored = nr_n_stored + 1
!          col = nr_n_stored
!       else
!          nr_model_buf(1:NR_MAX_MODELS-1)      = nr_model_buf(2:NR_MAX_MODELS)
!          nr_resid_buf(1:NR_MAX_MODELS-1, :)   = nr_resid_buf(2:NR_MAX_MODELS, :)
!          nr_ycoord_buf(1:NR_MAX_MODELS-1, :)    = nr_ycoord_buf(2:NR_MAX_MODELS, :)
!          col = NR_MAX_MODELS
!       end if

!       nr_model_buf(col) = s%model_number
!       nr_zone_buf(col) = s%nz
!       ! do k = 1, s%nz
!       !    nr_resid_buf(col, k) = safe_log10(maxval(abs(s%equ(1:neq, k))))       ! residuals are cell-centered
!       !    nr_ycoord_buf(col, k) = s%m(k) - 0.5 * s%dm(k) ! outer face mass coord minus half mass of the cell
!       ! end do

!       do k = 1, s%nz
!          nr_resid_buf(col, k) = &
!             safe_log10(maxval(abs(s%equ(1:neq, k))))
!          select case (s%x_integer_ctrl(1))
!          case (NR_COORD_MASS)
!             nr_ycoord_buf(col, k) = s%m(k)/Msun
!          case (NR_COORD_LOGR)
!             ! cell-centered radius in log10(R/Rsun)
!             if (k < s%nz) then
!                nr_ycoord_buf(col, k) = safe_log10(0.5d0 * (s%r(k) + s%r(k+1)) / Rsun)
!             else
!                nr_ycoord_buf(col, k) = safe_log10(0.5d0*s%r(k)/Rsun)
!             end if
!          case (NR_COORD_TAU)
!             ! optical depth
!             nr_ycoord_buf(col, k) = safe_log10(s%tau(k))
!          case default
!             ! cell-centered mass coordinate
!             nr_ycoord_buf(col, k) = (s%m(k) - 0.5d0*s%dm(k))/Msun
!          end select
!       end do


!       call pgslct(device_id)
!       call pgbbuf()
!       call pgeras()
!       call Kipp_residuals_render(ierr, id)
!       call pgebuf()

!    end subroutine Kipp_residuals_plot

!    subroutine Kipp_residuals_render(ierr, id)
!       integer, intent(in)  :: id
!       integer, intent(out) :: ierr
!       real, parameter :: missing_val = 1d-30
!       real, allocatable :: img(:,:), coord_grid(:)
!       real :: tr(6), fg, bg, xlo, xhi, dx
!       real :: clo, chi, dcoord, coord_j, frac
!       integer :: nx, ny, i, j, k, nz_i
!       logical :: y_reverse
!       type(star_info), pointer :: s

!       ! ----------------------------------------- colorbar definition
!       real :: bright, contra
!       real, dimension(9) :: l, r, g, b

!       l = (/ &
!          0.00, 0.12, 0.25, 0.38, 0.50, &
!          0.62, 0.75, 0.88, 1.00 /)

!       r = (/ &
!          0.00, 0.05, 0.15, 0.32, 0.55, &
!          0.78, 0.93, 0.99, 1.00 /)

!       g = (/ &
!          0.00, 0.02, 0.04, 0.06, 0.10, &
!          0.18, 0.32, 0.58, 0.95 /)

!       b = (/ &
!          0.00, 0.08, 0.18, 0.32, 0.36, &
!          0.28, 0.16, 0.08, 0.80 /)

!       contra = 1.0
!       bright = 0.5

!       call pgscir(16,255)
!       call pgctab(l, r, g, b, 9, contra, bright)
!       ! ----------------------------------------- end colorbar

!       ierr = 0
!       call star_ptr(id, s, ierr)
!       if (ierr /= 0) return

!       if (nr_n_stored < 2) return

!       nx = nr_n_stored
!       ny = nr_n_cells

!       allocate(img(nx, ny))
!       allocate(coord_grid(ny))

!       chi = -1.0e30
!       clo =  1.0e30
!       do i = 1, nx
!          nz_i = nr_zone_buf(i)
!          do k = 1, nz_i
!             if (nr_ycoord_buf(i,k) > chi) chi = nr_ycoord_buf(i,k)
!             ! m(k) is positive; skip any zero-initialised slots
!             if (nr_ycoord_buf(i,k) > 0.0 .and. nr_ycoord_buf(i,k) < clo) &
!                clo = nr_ycoord_buf(i,k)
!          end do
!       end do

!       select case (s%x_integer_ctrl(1))
!       case (NR_COORD_MASS)
!          y_reverse = .false.
!       case (NR_COORD_LOGR)
!          y_reverse = .false.
!       case (NR_COORD_TAU)
!          y_reverse = .true.
!       case default
!          y_reverse = .false.
!       end select

!       dcoord = (chi - clo) / real(ny - 1)

!       do j = 1, ny
!          if (.not. y_reverse) then
!             coord_grid(j) = clo + real(j - 1) * dcoord
!          else
!             coord_grid(j) = chi - real(j - 1) * dcoord
!          end if
!       end do

!       ! --- interpolate residuals onto the uniform grid ---
!       img = missing_val
!       do i = 1, nx
!          nz_i = nr_zone_buf(i)
!          do j = 1, ny
!             coord_j = coord_grid(j)
!             do k = 1, nz_i - 1
!                ! skip invalid cells
!                if (nr_resid_buf(i,k)   == missing_val) cycle
!                if (nr_resid_buf(i,k+1) == missing_val) cycle
!                ! works for either increasing or decreasing coordinates
!                if ( (nr_ycoord_buf(i,k)   - coord_j) * &
!                   (nr_ycoord_buf(i,k+1) - coord_j) <= 0d0 ) then
!                   if (abs(nr_ycoord_buf(i,k) - nr_ycoord_buf(i,k+1)) > 1d-99) then
!                      frac = (coord_j - nr_ycoord_buf(i,k+1)) / &
!                         (nr_ycoord_buf(i,k) - nr_ycoord_buf(i,k+1))
!                   else
!                      frac = 0.5d0
!                   end if
!                   img(i,j) = nr_resid_buf(i,k+1) + &
!                      frac * (nr_resid_buf(i,k) - nr_resid_buf(i,k+1))
!                   exit
!                end if
!             end do
!          end do
!       end do
!       ! --- x range centred on model numbers ---
!       dx  = max(1.0, real(nr_model_buf(nx) - nr_model_buf(1)) / real(nx - 1))
!       xlo = real(nr_model_buf(1))  - 0.5*dx
!       xhi = real(nr_model_buf(nx)) + 0.5*dx

!       ! --- transformation matrix for pgimag ---
!       ! world_y(j) = tr(4) + tr(6)*j
!       ! j=1  -> mhi  =>  tr(4) = mhi - tr(6) = mhi + dm
!       ! j=ny -> mlo  =>  confirms tr(6) = -dm
!       tr(1) = xlo
!       tr(2) = (xhi - xlo) / real(nx)
!       tr(3) = 0.0
!       tr(4) = chi + dcoord
!       tr(5) = 0.0
!       tr(6) = -dcoord

!       ! color scale
!       fg = nr_resid_max
!       bg = nr_resid_min
!       if (fg <= -100) fg = maxval(img)
!       if (bg <= -100) bg = minval(img, mask=(img > -19.9))
!       if (bg >= fg)   bg = fg - 6.0

!       ! --- main panel ---
!       ! y-axis: mlo at bottom (centre), mhi at top (surface)
!       call pgsvp(0.10, 0.82, 0.12, 0.92)
!       if (.not. y_reverse) then
!          call pgswin(xlo, xhi, clo - 0.5*dcoord, chi + 0.5*dcoord)
!       else
!          call pgswin(xlo, xhi, chi + 0.5*dcoord, clo - 0.5*dcoord)
!       end if
!       call pgimag(img, nx, ny, 1, nx, 1, ny, fg, bg, tr)

!       call pgsci(clr_Foreground)
!       call pgsch(1.2)
!       call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
!       select case (s%x_integer_ctrl(1))
!       case (NR_COORD_MASS)
!          y_reverse = .false.
!          call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
!             'Newton Raphson Solver Residuals of accepted step')
!       case (NR_COORD_LOGR)
!          y_reverse = .false.
!          call pglab('Model Number', 'log(R/R\d\(2281)\u)', &
!             'Newton Raphson Solver Residuals of accepted step')
!       case (NR_COORD_TAU)
!          y_reverse = .true.
!          call pglab('Model Number', 'log(tau)', &
!             'Newton Raphson Solver Residuals of accepted step')
!       case default
!          y_reverse = .false.
!          call pglab('Model Number', 'Mass (M\d\(2281)\u)', &
!             'Newton Raphson Solver Residuals of accepted step')
!       end select

!       ! --- color wedge ---
!       call pgsvp(0.84, 0.90, 0.12, 0.92)
!       call pgswin(0.0, 1.0, bg, fg)
!       call pgwedg('RI', 0.0, 4.0, bg, fg, 'log10(max(|res|))')

!       deallocate(img)
!       deallocate(coord_grid)

!    end subroutine Kipp_residuals_render

!    subroutine do_Kipp_residuals_plot(id, ierr)
!       integer, intent(in)  :: id
!       integer, intent(out) :: ierr

!       type(pgstar_win_file_data), pointer :: p
!       type(star_info), pointer :: s

!       ierr = 0
!       call star_ptr(id, s, ierr)
!       if (ierr /= 0) return

!       p => s%pg%pgstar_win_file_ptr(i_Other)
!       p%plot              => Kipp_residuals_plot
!       p%id                =  i_Other
!       p%name              = 'Newton Raphson Residuals - accepted step'
!       p%win_flag          = .true.
!       p%win_width         =  12.0
!       p%win_aspect_ratio  =  0.6
!       p%file_flag         = .true.
!       p%file_dir          = 'png'         ! folder where to save files
!       p%file_prefix       = 'nr_resid_'
!       p%file_interval     =  5
!       p%file_width        = -1.0
!       p%file_aspect_ratio = -1.0

!    end subroutine do_Kipp_residuals_plot

! end module pgstar_kipp_residuals
