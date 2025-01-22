! ***********************************************************************
!
!   Copyright (C) 2025 Matthias Fabry & The MESA Team
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

module binary_roche_deformation

   ! computes the stellar structure correction factors fp/ft and the specific
   ! moment of inertia i_rot assuming shellularity in the Roche potential of
   ! a binary star. Follows Fabry, Marchant and Sana, 2022, A&A 661, A123.

   use interp_2d_lib_db
   use auto_diff
   use star_def
   use binary_def
   use binary_lib, only : binary_eval_rlobe


   implicit none

   real(dp), parameter :: nudge = 1d-4
   real(dp), pointer :: xvals(:), yvals(:), yvals_gtr_than_1(:), fpfunc1d(:), ftfunc1d(:), &
         irotfunc1d(:), otherrfunc1d(:), afunc1d(:)
   logical :: inter_ok = .false., dbg = .true.
   integer :: num_xpts, num_ypts, num_ypts_gtr_than_1


contains

   subroutine build_roche_interpolators
      use const_def, only: mesa_data_dir
      real(dp) :: xtest, ytest, testval
      integer :: ierr
      character(len=strlen) :: upstairs

      include 'formats'

      if (.not. inter_ok) then
         upstairs = trim(mesa_data_dir) // 'roche_data/'  ! where fp/ft data lives
         if (dbg) then
            write(*, 1) 'starting interpolator setup'
         end if
         call setup_interpolator(trim(upstairs) // 'fp_data.txt', xvals, num_xpts, yvals, &
               num_ypts, fpfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'ft_data.txt', xvals, num_xpts, yvals, &
               num_ypts, ftfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'irot_data.txt', xvals, num_xpts, &
               yvals, num_ypts, irotfunc1d, ierr)
         if (dbg) then
            xtest = -0.5
            ytest = 1.35
            ! test fp interpolator
            write(*, 11) 'grid size', num_xpts, num_ypts
            write(*, 1) 'setup interpolators succesful,'

            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  fpfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'fp   test gave should be close to 0.6', testval
            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  ftfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'ft   test gave should be close to 0.8', testval
            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  irotfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'irot test gave should be close to 0.4', testval
            inter_ok = .true.
         end if
      end if

      contains
      ! interpolator data reading
      subroutine setup_interpolator(filename, xs, num_xs, ys, num_ys, func1d, ierr)
         integer, intent(out) :: ierr, num_xs, num_ys
         character(len = *) :: filename
         integer :: k, iounit
         real(dp), pointer, intent(out) :: xs(:), ys(:), func1d(:)
         real(dp), pointer :: func(:,:,:)

         include 'formats'

         if (dbg) then
            write(*, 1) 'loading ' // filename
         end if
         ! open data to interpolate
         open(newunit = iounit, file = trim(filename), status = 'old', action = 'read',&
               iostat = ierr)

         read(iounit, *, iostat = ierr) num_xs
         allocate(xs(num_xs))
         do k = 1, num_xs
            read(iounit, *, iostat = ierr) xs(k)
         end do

         read(iounit, *, iostat = ierr) num_ys
         allocate(ys(num_ys))
         do k = 1, num_ys
            read(iounit, *, iostat = ierr) ys(k)
         end do

         ! create a 1d array with all the data, point func to it
         allocate(func1d(4 * num_xs * num_ys))
         func(1:4, 1:num_xs, 1:num_ys) => func1d(1:4 * num_xs * num_ys)
         do k = 1, num_xs
            read(iounit, *, iostat = ierr) func(1, k, :)
         end do

         if (ierr /= 0) then
            close(iounit)
         end if
         ! create interpolator
         call interp_mkbipm_db(xs, num_xs, ys, num_ys, func1d, num_xs, ierr)
      end subroutine setup_interpolator

   end subroutine build_roche_interpolators


   real(dp) function eval_fp(lq, ar, ierr) result(fp)
      ! evaluates fp of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            fpfunc1d, num_xpts, fp, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval fp", ar, fp
      end if
   end function eval_fp

   real(dp) function eval_ft(lq, ar, ierr) result(ft)
      ! evaluates ft of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            ftfunc1d, num_xpts, ft, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval ft", ar, ft
      end if
   end function eval_ft

   real(dp) function eval_irot(lq, ar, ierr) result(irot)
      ! evaluates moment of inertia irot of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), in units of separation^2
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            irotfunc1d, num_xpts, irot, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval irot", ar, irot
      end if
   end function eval_irot

   ! deformation
   subroutine roche_fp_ft(id, nz, r, fp, ft, r_polar, r_equatorial, report_ierr, ierr)
      integer, intent(in) :: id, nz
      real(dp), intent(in) :: r(:) ! (nz)
      logical, intent(in) :: report_ierr
      real(dp), intent(inout) :: r_polar(:), r_equatorial(:) ! (nz)
      type (auto_diff_real_star_order1), intent(out) :: fp(:), ft(:)
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: j, this_star=0, other_star=0
      real(dp) :: r_roche, lq, m1, m2, a, ar

      include 'formats'

      if (.not. inter_ok) then
         write(*, *) "interpolators not setup, should happen in startup"
         stop
      end if

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return
      call assign_stars(id, this_star, other_star, ierr)
      if (ierr /= 0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      ! if masses or sep are not yet set at the binary level, use these
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
            pow((b% initial_period_in_days) * 86400, 2) / (4 * pi2), one_third)
      if (r_roche <= 0) r_roche = binary_eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)

      !$OMP PARALLEL DO PRIVATE(j, ar) SCHEDULE(dynamic,2)
      do j = 1, s% nz  ! for every cell, compute fp, ft from an interpolating table
         ar = r(j) / r_roche
         ! set values
         fp(j) = eval_fp(lq, ar, ierr)
         ft(j) = eval_ft(lq, ar, ierr)
         ! set log derivatives
         fp(j)% d1Array(i_lnR_00) = (eval_fp(lq, ar + nudge, ierr) - fp(j)% val) / nudge * ar
         ft(j)% d1Array(i_lnR_00) = (eval_ft(lq, ar + nudge, ierr) - ft(j)% val) / nudge * ar
!            if (fp(j) == 0d0 .or. fp(j) == 0d0) write(*, *) j, ar, fp(j), ft(j), ierr  ! debug
         ! fix these to the current radius, they're only used for some wind mass loss enhancement
         r_polar(j) = r(j)
         r_equatorial(j) = r(j)
      end do
      !$OMP END PARALLEL DO
   end subroutine roche_fp_ft

   subroutine roche_irot(id, k, r00, w_div_w_crit_roche, i_rot)
      integer, intent(in) :: id, k
      real(dp), intent(in) :: r00, w_div_w_crit_roche
      type (auto_diff_real_star_order1), intent(out) :: i_rot
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: ierr, j, this_star=0, other_star=0
      real(dp) :: r_roche, a, m1, m2, ar, lq

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr/=0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr/=0) return
      call assign_stars(id, this_star, other_star, ierr)
      if (ierr/=0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
         pow((b% initial_period_in_days) * 86400, 2) / (4 * pi2), one_third)
      if (r_roche <= 0) r_roche = binary_eval_rlobe(m1, m2, a)
!
      lq = log10(m2 / m1)
      i_rot = 0d0

      if (inter_ok) then
         ar = r00 / r_roche
         ! set value
         i_rot = eval_irot(lq, ar, ierr) * a * a
         ! set radius derivative
         i_rot% d1Array(i_lnR_00) = (eval_irot(lq, ar + nudge, ierr) * a * a - i_rot% val) / nudge * ar
!         write(*, *) r00, r_roche, i_rot, w_div_w_crit_roche  ! debug
      end if
   end subroutine roche_irot

end module binary_roche_deformation
