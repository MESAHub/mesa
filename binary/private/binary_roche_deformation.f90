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

   use ieee_arithmetic, only: ieee_is_nan
   use interp_2d_lib_db
   use auto_diff
   use star_def
   use binary_def
   use binary_utils, only : eval_rlobe


   implicit none

   real(dp), parameter :: nudge = 1d-4
   real(dp), pointer :: xvals(:), yvals(:), yvals_gtr_than_1(:), fpfunc1d(:), ftfunc1d(:), &
         irotfunc1d(:), otherrfunc1d(:), afunc1d(:)
   logical :: inter_ok = .false., dbg = .false.
   integer :: num_xpts, num_ypts, num_ypts_gtr_than_1

contains

   subroutine build_roche_interpolators
      use const_def, only: mesa_data_dir
      real(dp) :: xtest, ytest, testval
      integer :: ierr
      character(len=strlen) :: upstairs
      include 'formats'
      if (.not. inter_ok) then
         upstairs = trim(mesa_data_dir) // '/roche_data/'  ! where fp/ft data lives
         if (dbg) write(*, 1) 'starting interpolator setup'
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
         end if
         inter_ok = .true.
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
            write(*, '(a100)') 'loading ' // filename
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
      ! evaluates fp of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            fpfunc1d, num_xpts, fp, ierr)
      if (ierr /= 0) write(*, 1) "error in eval fp", ar, fp
   end function eval_fp

   real(dp) function eval_ft(lq, ar, ierr) result(ft)
      ! evaluates ft of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            ftfunc1d, num_xpts, ft, ierr)
      if (ierr /= 0) write(*, 1) "error in eval ft", ar, ft

   end function eval_ft

   real(dp) function eval_irot(lq, ar, ierr) result(irot)
      ! evaluates moment of inertia irot divided by the spherical equivalent moi (2/3 r_\psi^2)
      ! of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            irotfunc1d, num_xpts, irot, ierr)
      if (ierr /= 0) write(*, 1) "error in eval irot", ar, irot

   end function eval_irot

   subroutine roche_fp_ft(id, r, fp, ft, r_polar, r_equatorial, report_ierr, ierr)
      integer, intent(in) :: id
      real(dp), intent(in) :: r
      logical, intent(in) :: report_ierr
      real(dp), intent(inout) :: r_polar, r_equatorial
      type (auto_diff_real_star_order1), intent(out) :: fp, ft
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: j, this_star=0, other_star=0
      real(dp) :: r_roche, lq, m1, m2, a, ar

      include 'formats'

      if (.not. inter_ok) then
         write(*, 1) "interpolators not setup, should happen in startup"
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
      if (r_roche <= 0) r_roche = eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)

      ar = r / r_roche
      ! set values
      fp = eval_fp(lq, ar, ierr)
      ft = eval_ft(lq, ar, ierr)
      ! set log derivatives
      fp% d1Array(i_lnR_00) = (eval_fp(lq, ar + nudge, ierr) - fp% val) / nudge * ar
      ft% d1Array(i_lnR_00) = (eval_ft(lq, ar + nudge, ierr) - ft% val) / nudge * ar
!            if (fp(j) == 0d0 .or. fp(j) == 0d0) write(*, *) j, ar, fp(j), ft(j), ierr  ! debug
      ! fix these to the current radius, they're only used for some wind mass loss enhancement
      r_polar = r
      r_equatorial = r
   end subroutine roche_fp_ft

   subroutine roche_irot(id, r00, i_rot)
      integer, intent(in) :: id
      real(dp), intent(in) :: r00
      type (auto_diff_real_star_order1), intent(out) :: i_rot
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: ierr, j, this_star=0, other_star=0
      real(dp) :: r_roche, a, m1, m2, ar, lq, eval, d_eval_d_ar

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
      if (m1 <= 0d0) m1 = b% m1 * Msun
      if (m2 <= 0d0) m2 = b% m2 * Msun
      if (a <= 0d0) a = pow(standard_cgrav * (m1 + m2) * &
         pow((b% initial_period_in_days) * secday, 2) / (4d0 * pi2), one_third)
      if (r_roche <= 0d0) r_roche = eval_rlobe(m1, m2, a)
!
      lq = log10(m2 / m1)
      i_rot = 0d0

      if (inter_ok) then
         ar = r00 / r_roche
         ! set value
         eval = eval_irot(lq, ar, ierr)
         d_eval_d_ar = (eval_irot(lq, ar + nudge, ierr) - eval) / nudge
         ! scale value with spherical moi (see eval_irot)
         i_rot = eval * (two_thirds * r00 * r00)
         ! set radius derivative
         i_rot% d1Array(i_lnR_00) = two_thirds * r00 * r00 * (ar * d_eval_d_ar + 2 * eval)
!         write(*, *) r00, r_roche, i_rot, w_div_w_crit_roche  ! debug
      end if
   end subroutine roche_irot

   subroutine assign_stars(id, this_star, other_star, ierr)
      ! determine which star is which in the binary
      integer, intent(in) :: id
      integer, intent(out) :: this_star, other_star, ierr
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      call binary_ptr(s% binary_id, b, ierr)

      if (ierr /= 0) return

      if (b% s1% id == s% id) then
         this_star = 1
         other_star = 2
      else if (b% s2% id == s% id) then
         this_star = 2
         other_star = 1
      else
         ierr = 1
      end if

   end subroutine assign_stars

   ! determines by what fraction the tidal deformation corrections should be used
   ! eg fp = f_switch * fp_tidal + (1-f_switch) * fp_single
   ! it uses the synchronicity parameter to estimate this. When the shell is quite synchronous,
   ! you probably want to use the tidal corrections, so f_switch -> 1, when not synchronous, the single
   ! star rotation deformation is likely more accurate, so f_switch -> 0 in that case.
   subroutine synchronicity(id, k, omega_in, f_switch, ierr)
      integer, intent(in) :: id, k
      real(dp), intent(in) :: omega_in
      integer, intent(out) :: ierr
      real(dp), intent(out) :: f_switch

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      real(dp) :: omega, omega_sync, f_sync, p

      include 'formats'

      call star_ptr(id, s, ierr)
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return

      if (b% use_other_tidal_deformation_switch_function) then
         call b% other_tidal_deformation_switch_function(id, k, omega_in, f_switch, ierr)
         return
      end if

      p = b% period
      if (p <= 0d0) p = b% initial_period_in_days * secday
      omega_sync = 2*pi/p
      if (ieee_is_nan(omega_in)) then
         omega = 0d0
      else if (omega_in < 0d0) then
         omega = abs(omega_in)
      else
         omega = omega_in
      end if
      f_sync = omega / omega_sync
      f_sync = min(f_sync, 1d0 / f_sync)  ! we could be super or sub synchronous

      f_switch = 1d0 / (1d0 + exp(-b% f_sync_switch_width * (f_sync - b% f_sync_switch_from_rot_defor)))
      ! apply limiting values if f_switch is far up (or down) the sigmoid
      if (f_switch < b% f_sync_switch_lim) then
         f_switch = 0d0
      else if (1d0 - f_switch < b% f_sync_switch_lim) then
         f_switch = 1d0
      end if

      if (ieee_is_nan(f_switch) .and. k == 1) then
         write(*, 1) "error in synchronicity", f_switch, f_sync, omega, omega_sync, p
         ierr = 1
      end if

   end subroutine synchronicity
end module binary_roche_deformation
