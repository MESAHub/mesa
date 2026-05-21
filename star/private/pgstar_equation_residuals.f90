! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

   use star_pgstar


   implicit none

contains

   subroutine equ_resid_plot(id, device_id, ierr)

      integer, intent(in)  :: id
      integer, intent(in)  :: device_id
      integer, intent(out) :: ierr

      type(star_info), pointer :: s

      integer :: i_eq
      integer :: n, last_model_plotted

      real :: xmin, xmax
      real :: ymin, ymax

      real, dimension(max_resid_hist) :: xvec
      real, dimension(max_resid_hist) :: yvec

      ierr = 0

      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (.not. allocated(resid_hist_model)) then
         resid_hist_nvar = s%nvar_hydro
         allocate(resid_hist_model(max_resid_hist))
         allocate(resid_hist_vals(resid_hist_nvar, max_resid_hist))
         allocate(resid_equ_names(resid_hist_nvar))

         do i_eq = 1, resid_hist_nvar
            resid_equ_names(i_eq) = trim(s%nameofequ(i_eq))
         end do

         n_resid_hist = 0
      end if

      if (s%model_number /= last_model_plotted) then
         last_model_plotted = s%model_number
         if (n_resid_hist < max_resid_hist) then
            n_resid_hist = n_resid_hist + 1
         else
            resid_hist_model(1:max_resid_hist-1) = resid_hist_model(2:max_resid_hist)
            resid_hist_vals(:,1:max_resid_hist-1) = resid_hist_vals(:,2:max_resid_hist)
         end if

         resid_hist_model(n_resid_hist) = s%model_number

         do i_eq = 1, resid_hist_nvar
            resid_hist_vals(i_eq,n_resid_hist) = &
               real(max(1d-40, &
               maxval(abs(s%equ(i_eq,1:s%nz)))))
         end do

      end if

      n = n_resid_hist
      if (n < 2) return ! nothing to plot yet

      xvec(1:n) = real(resid_hist_model(1:n))

      xmin = xvec(1)
      xmax = xvec(n)

      if (xmin >= xmax) xmax = xmin + 1.0

      ymin = minval(resid_hist_vals(:,1:n))
      ymax = maxval(resid_hist_vals(:,1:n))
      ymin = max(ymin, 1e-20)
      ymax = max(ymax, 10.0*ymin)

      call pgslct(device_id)
      call pgsave
      call pgeras
      call pgsvp(0.12, 0.88, 0.12, 0.88)
      call pgswin(xmin, xmax, log10(ymin), log10(ymax))
      call pgscf(1)
      call pgsch(1.0)
      call pgsci(1)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
      call pgmtxt('B',2.5,0.5,0.5,'Model Number')
      call pgmtxt('L',3.0,0.5,0.5, 'log10(max |residual|)')
      call pgmtxt('T',1.0,0.5,0.5,'Maximum structural equation residuals')

      ! draw lines
      do i_eq = 1, resid_hist_nvar
         yvec(1:n) = log10(resid_hist_vals(i_eq,1:n))
         call pgsci(mod(i_eq-1,13)+2)
         call pgline(n, xvec, yvec)
      end do

      ! legend
      call pgsch(0.65)

      do i_eq = 1, resid_hist_nvar
         call pgsci(mod(i_eq-1,13)+2)
         call pgptxt( &
            xmin + (0.55 + 0.06*real(i_eq-1))*(xmax-xmin), &
            log10(ymax) - 0.03*(log10(ymax)-log10(ymin)), &
            0.0, 0.0, &
            trim(resid_equ_names(i_eq)))
      end do

      call pgsci(1)
      call pgunsa

   end subroutine equ_resid_plot



   ! history of residuals for each structure equation
   subroutine do_equ_resid_plot(id, ierr)
      integer, intent(in)  :: id
      integer, intent(out) :: ierr

      type(pgstar_win_file_data), pointer :: p
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      p => s%pg%pgstar_win_file_ptr(i_Other)
      p%plot              => equ_resid_plot
      p%id                =  i_Other
      p%name              = 'Max Residual per equation across mesh points'
      p%win_flag          = .true.
      p%win_width         =  12.0
      p%win_aspect_ratio  =  0.6
      p%file_flag         = .true.
      p%file_dir          = 'png'         ! folder where to save files
      p%file_prefix       = 'eq_resid_'
      p%file_interval     =  5
      p%file_width        = -1.0
      p%file_aspect_ratio = -1.0

   end subroutine do_equ_resid_plot


end module pgstar_equation_residuals
