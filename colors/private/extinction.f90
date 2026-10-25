! ***********************************************************************
!
!   Copyright (C) 2026  Niall Miller & The MESA Team
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

! interstellar extinction, applied to the SED before it is integrated through the filters:
!
!    F_obs(lambda) = F(lambda) * 10**(-0.4 * A_V * k(lambda)),    k = A(lambda)/A(V)
!
! A_V sets the amount of dust and the extinction law sets the shape k, which for the
! Milky Way laws depends on R_V = A(V)/E(B-V).
!
! each law is a pure subroutine  law(x, R_V, k)  that fills k(:) at the wavenumbers
! x(:) = 1/lambda [micron^-1], the variable the published fits are written in.
! to add a law, write one below, give it a case in extinction_curve, and add its name
! (with the R_V range it was calibrated over, if it uses R_V) to check_extinction_controls.
!
! outside a law's fitted range: longward, the law's own functional form is continued
! (each falls to zero as lambda -> infinity); shortward, k is held at its value at the
! law's blue limit, since the far-UV fits diverge beyond it.
!
! coefficients follow the dust_extinction package (Gordon 2024, JOSS, 9, 7023), and for
! FM07 the extinction package (Barbary 2016), which serve as reference implementations.

module extinction

   use const_def, only: dp
   use math_lib, only: pow, pow2, pow3, exp, exp10
   use utils_lib, only: StrUpCase

   implicit none

   private
   public :: apply_extinction, extinction_curve, check_extinction_controls

contains

   ! reddens fluxes in place; wavelengths in Angstrom
   subroutine apply_extinction(law, A_V, R_V, wavelengths, fluxes, ierr)
      character(len=*), intent(in) :: law
      real(dp), intent(in) :: A_V, R_V
      real(dp), intent(in) :: wavelengths(:)
      real(dp), intent(inout) :: fluxes(:)
      integer, intent(out) :: ierr

      real(dp), allocatable :: k(:)

      ! checked on every call, since the controls can change after setup
      call check_extinction_controls(law, A_V, R_V, ierr)
      if (ierr /= 0) return
      if (StrUpCase(law) == 'NONE' .or. A_V == 0.0_dp) return

      allocate (k(size(wavelengths)))
      call extinction_curve(law, R_V, wavelengths, k, ierr)
      if (ierr /= 0) return

      k = exp10(-0.4_dp*A_V*k)
      fluxes = fluxes*k
   end subroutine apply_extinction

   ! k = A(lambda)/A(V) for the named law; wavelengths in Angstrom
   subroutine extinction_curve(law, R_V, wavelengths, k, ierr)
      character(len=*), intent(in) :: law
      real(dp), intent(in) :: R_V
      real(dp), intent(in) :: wavelengths(:)
      real(dp), intent(out) :: k(:)
      integer, intent(out) :: ierr

      real(dp), allocatable :: x(:)

      ierr = 0
      allocate (x(size(wavelengths)))

      ! a non-positive wavelength (never present in a valid SED) is sent to the far-UV,
      ! where every law is held constant
      x = 1.0d4
      where (wavelengths > 0.0_dp) x = 1.0d4/wavelengths

      select case (StrUpCase(law))
      case ('NONE')
         k = 0.0_dp
      case ('CCM89')
         call ccm89(x, R_V, k)
      case ('O94')
         call o94(x, R_V, k)
      case ('F99')
         call f99(x, R_V, k)
      case ('FM07')
         call fm07(x, k)
      case ('G23')
         call g23(x, R_V, k)
      case ('G03_SMCBAR')
         call g03_smcbar(x, k)
      case ('G03_LMCAVG')
         call g03_lmcavg(x, k)
      case default
         write (*, '(a)') 'colors error: unknown extinction_law '//trim(law)
         k = 0.0_dp
         ierr = -1
      end select
   end subroutine extinction_curve

   ! checks the extinction controls: a known law, A_V >= 0, and R_V inside the range the
   ! law was calibrated over (outside it the fits can go negative and brighten the star).
   ! with report, also says when a control has no effect
   subroutine check_extinction_controls(law, A_V, R_V, ierr, report)
      character(len=*), intent(in) :: law
      real(dp), intent(in) :: A_V, R_V
      integer, intent(out) :: ierr
      logical, intent(in), optional :: report

      real(dp) :: rv_lo, rv_hi
      logical :: uses_rv

      ierr = 0
      rv_lo = 0.0_dp
      rv_hi = 0.0_dp

      ! every law, with the R_V range it was calibrated over; fixed-shape laws ignore R_V
      uses_rv = .true.
      select case (StrUpCase(law))
      case ('NONE', 'FM07', 'G03_SMCBAR', 'G03_LMCAVG')
         uses_rv = .false.
      case ('CCM89', 'O94', 'F99')
         rv_lo = 2.0_dp
         rv_hi = 6.0_dp
      case ('G23')
         rv_lo = 2.3_dp
         rv_hi = 5.6_dp
      case default
         write (*, '(a)') 'colors error: unknown extinction_law '//trim(law)// &
            '; options are none, CCM89, O94, F99, FM07, G23, G03_SMCBar, G03_LMCAvg'
         ierr = -1
         return
      end select

      if (A_V < 0.0_dp) then
         write (*, '(a)') 'colors error: A_V must be >= 0, got '//num_str(A_V)
         ierr = -1
         return
      end if
      if (uses_rv .and. (R_V < rv_lo .or. R_V > rv_hi)) then
         write (*, '(a)') 'colors error: '//trim(law)//' is calibrated for '//num_str(rv_lo)// &
            ' <= R_V <= '//num_str(rv_hi)//', got R_V = '//num_str(R_V)
         ierr = -1
         return
      end if

      if (.not. present(report)) return
      if (.not. report) return
      if (StrUpCase(law) == 'NONE') then
         if (A_V /= 0.0_dp) write (*, '(a)') 'colors: A_V has no effect because extinction_law = ''none'''
      else if (A_V == 0.0_dp) then
         write (*, '(a)') 'colors: extinction_law = '//trim(law)//' has no effect with A_V = 0'
      else if (.not. uses_rv .and. R_V /= 3.1_dp) then
         write (*, '(a)') 'colors: '//trim(law)//' has a fixed shape, so the R_V control is ignored'
      end if
   end subroutine check_extinction_controls

   ! v as short text for messages, e.g. -0.500
   function num_str(v) result(s)
      real(dp), intent(in) :: v
      character(len=:), allocatable :: s

      character(len=32) :: buf

      write (buf, '(f12.3)') v
      s = trim(adjustl(buf))
   end function num_str

   ! ----------------------------------------------------------------------
   ! extinction laws
   ! ----------------------------------------------------------------------

   ! Cardelli, Clayton & Mathis (1989, ApJ, 345, 245); fitted for 0.3 <= x <= 10
   pure subroutine ccm89(x, R_V, k)
      real(dp), intent(in) :: x(:), R_V
      real(dp), intent(out) :: k(:)

      ! optical a(y) and b(y), y = x - 1.82, in ascending powers of y
      real(dp), parameter :: opt_a(8) = [1.0_dp, 0.17699_dp, -0.50447_dp, -0.02427_dp, &
                                         0.72085_dp, 0.01979_dp, -0.77530_dp, 0.32999_dp]
      real(dp), parameter :: opt_b(8) = [0.0_dp, 1.41338_dp, 2.28305_dp, 1.07233_dp, &
                                         -5.38434_dp, -0.62251_dp, 5.30260_dp, -2.09002_dp]

      call cardelli_form(x, R_V, opt_a, opt_b, k)
   end subroutine ccm89

   ! O'Donnell (1994, ApJ, 422, 158): CCM89 with a refitted optical section
   pure subroutine o94(x, R_V, k)
      real(dp), intent(in) :: x(:), R_V
      real(dp), intent(out) :: k(:)

      real(dp), parameter :: opt_a(9) = [1.0_dp, 0.104_dp, -0.609_dp, 0.701_dp, 1.137_dp, &
                                         -1.718_dp, -0.827_dp, 1.647_dp, -0.505_dp]
      real(dp), parameter :: opt_b(9) = [0.0_dp, 1.952_dp, 2.908_dp, -3.989_dp, -7.985_dp, &
                                         11.102_dp, 5.491_dp, -10.805_dp, 3.347_dp]

      call cardelli_form(x, R_V, opt_a, opt_b, k)
   end subroutine o94

   ! Fitzpatrick (1999, PASP, 111, 63); fitted for 0.3 <= x <= 10. the last optical
   ! anchor uses +1.208, which reproduces F99 table 3 (table 4 prints -1.208)
   pure subroutine f99(x, R_V, k)
      real(dp), intent(in) :: x(:), R_V
      real(dp), intent(out) :: k(:)

      real(dp), parameter :: knot_lambda(6) = &
         [26500.0_dp, 12200.0_dp, 6000.0_dp, 5470.0_dp, 4670.0_dp, 4110.0_dp]
      real(dp) :: knot_k(6), c2

      ! A(lambda)/E(B-V) at the anchors
      knot_k(1:2) = [0.265_dp, 0.829_dp]*R_V/3.1_dp
      knot_k(3) = -0.426_dp + 1.0044_dp*R_V
      knot_k(4) = -0.050_dp + 1.0016_dp*R_V
      knot_k(5) = 0.701_dp + 1.0016_dp*R_V
      knot_k(6) = 1.208_dp + 1.0032_dp*R_V - 0.00033_dp*pow2(R_V)

      c2 = -0.824_dp + 4.717_dp/R_V
      call f99_form(x, R_V, [2.030_dp - 3.007_dp*c2, c2, 3.23_dp, 0.41_dp, 4.596_dp, 0.99_dp], &
                    1.0d4/knot_lambda, knot_k/R_V, k)
   end subroutine f99

   ! Fitzpatrick & Massa (2007, ApJ, 663, 320) average Milky Way curve, fixed at R_V = 3.1
   ! (the R_V control is not used), as implemented in the extinction package
   ! (Barbary 2016); fitted from 910 A to 6 micron. FM90 in the UV, with a quadratic far-UV
   ! rise above x = c5, and a natural cubic spline in E(x-V)/E(B-V) longward of 2700 A
   pure subroutine fm07(x, k)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: k(:)

      real(dp), parameter :: R_V = 3.1_dp, x_uv = 1.0d4/2700.0_dp
      real(dp), parameter :: sx(10) = [0.0_dp, 0.25_dp, 0.50_dp, 0.75_dp, 1.0_dp, 1.0d4/5530.0_dp, &
                                       1.0d4/4000.0_dp, 1.0d4/3300.0_dp, 1.0d4/2700.0_dp, 1.0d4/2600.0_dp]
      real(dp) :: sk(10), m(10), xi
      integer :: i

      sk(1:5) = (0.63_dp*R_V - 0.83_dp)*pow(sx(1:5), 1.84_dp) - R_V
      sk(6:8) = [0.0_dp, 1.322_dp, 2.055_dp]
      sk(9:10) = [fm07_uv(sx(9)), fm07_uv(sx(10))]
      call natural_spline(sx, sk, m)

      do i = 1, size(x)
         xi = min(x(i), 1.0d4/910.0_dp)
         if (xi >= x_uv) then
            k(i) = 1.0_dp + fm07_uv(xi)/R_V
         else
            k(i) = 1.0_dp + spline_value(sx, sk, m, xi)/R_V
         end if
      end do

   contains

      pure real(dp) function fm07_uv(xu)
         real(dp), intent(in) :: xu

         fm07_uv = -0.175_dp + 0.807_dp*xu + 2.991_dp*pow2(xu)/(pow2(pow2(xu) - pow2(4.592_dp)) + pow2(0.922_dp*xu))
         if (xu > 6.097_dp) fm07_uv = fm07_uv + 0.319_dp*pow2(xu - 6.097_dp)
      end function fm07_uv

   end subroutine fm07

   ! Gordon et al. (2023, ApJ, 950, 86); fitted for 0.0912 <= lambda <= 32 micron,
   ! including the 10 and 20 micron silicate features. the 1/0.9854 renormalisation
   ! corrects the F19 photometry to the HST/STIS spectrophotometry
   pure subroutine g23(x, R_V, k)
      real(dp), intent(in) :: x(:), R_V
      real(dp), intent(out) :: k(:)

      real(dp) :: xi, lam, w, a, b, a2, b2
      integer :: i

      do i = 1, size(x)
         xi = min(x(i), 1.0_dp/0.0912_dp)
         lam = 1.0_dp/xi
         if (lam > 1.1_dp) then
            call g23_ir(xi, a, b)
         else if (lam >= 0.9_dp) then  ! optical blends into IR over 0.9-1.1 micron
            call g23_optical(xi, a, b)
            call g23_ir(xi, a2, b2)
            w = smoothstep(lam, 0.9_dp, 1.1_dp)
            a = (1.0_dp - w)*a + w*a2
            b = (1.0_dp - w)*b + w*b2
         else if (lam > 0.33_dp) then
            call g23_optical(xi, a, b)
         else if (lam >= 0.3_dp) then  ! UV blends into optical over 0.30-0.33 micron
            call g23_uv(xi, a, b)
            call g23_optical(xi, a2, b2)
            w = smoothstep(lam, 0.3_dp, 0.33_dp)
            a = (1.0_dp - w)*a + w*a2
            b = (1.0_dp - w)*b + w*b2
         else
            call g23_uv(xi, a, b)
         end if
         k(i) = (a + b*(1.0_dp/R_V - 1.0_dp/3.1_dp))/0.9854_dp
      end do
   end subroutine g23

   ! Gordon et al. (2003, ApJ, 594, 279) SMC bar average, R_V = 2.74 (the R_V control is
   ! not used); fitted for 0.3 <= x <= 10
   pure subroutine g03_smcbar(x, k)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: k(:)

      call f99_form(x, 2.74_dp, [-4.959_dp, 2.264_dp, 0.389_dp, 0.461_dp, 4.6_dp, 1.0_dp], &
                    1.0_dp/[2.198_dp, 1.65_dp, 1.25_dp, 0.81_dp, 0.65_dp, 0.55_dp, 0.44_dp, 0.37_dp], &
                    [0.110_dp, 0.169_dp, 0.250_dp, 0.567_dp, 0.801_dp, 1.000_dp, 1.374_dp, 1.672_dp], k)
   end subroutine g03_smcbar

   ! Gordon et al. (2003, ApJ, 594, 279) LMC average, R_V = 3.41 (the R_V control is not
   ! used); fitted for 0.3 <= x <= 10
   pure subroutine g03_lmcavg(x, k)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: k(:)

      call f99_form(x, 3.41_dp, [-0.890_dp, 0.998_dp, 2.719_dp, 0.400_dp, 4.579_dp, 0.934_dp], &
                    1.0_dp/[2.198_dp, 1.65_dp, 1.25_dp, 0.55_dp, 0.44_dp, 0.37_dp], &
                    [0.100_dp, 0.186_dp, 0.257_dp, 1.000_dp, 1.293_dp, 1.518_dp], k)
   end subroutine g03_lmcavg

   ! ----------------------------------------------------------------------
   ! building blocks shared by the laws
   ! ----------------------------------------------------------------------

   ! k = a(x) + b(x)/R_V in the CCM89 form, with the given optical polynomials
   pure subroutine cardelli_form(x, R_V, opt_a, opt_b, k)
      real(dp), intent(in) :: x(:), R_V, opt_a(:), opt_b(:)
      real(dp), intent(out) :: k(:)

      real(dp) :: xi, y, a, b
      integer :: i

      do i = 1, size(x)
         xi = min(x(i), 10.0_dp)
         if (xi < 1.1_dp) then  ! infrared power law
            a = 0.574_dp*pow(xi, 1.61_dp)
            b = -0.527_dp*pow(xi, 1.61_dp)
         else if (xi < 3.3_dp) then  ! optical
            y = xi - 1.82_dp
            a = polynomial(opt_a, y)
            b = polynomial(opt_b, y)
         else if (xi <= 8.0_dp) then  ! near-UV, with extra curvature above x = 5.9
            a = 1.752_dp - 0.316_dp*xi - 0.104_dp/(pow2(xi - 4.67_dp) + 0.341_dp)
            b = -3.090_dp + 1.825_dp*xi + 1.206_dp/(pow2(xi - 4.62_dp) + 0.263_dp)
            if (xi >= 5.9_dp) then
               y = xi - 5.9_dp
               a = a - 0.04473_dp*pow2(y) - 0.009779_dp*pow3(y)
               b = b + 0.2130_dp*pow2(y) + 0.1207_dp*pow3(y)
            end if
         else  ! far-UV
            y = xi - 8.0_dp
            a = -1.073_dp - 0.628_dp*y + 0.137_dp*pow2(y) - 0.070_dp*pow3(y)
            b = 13.670_dp + 4.257_dp*y - 0.420_dp*pow2(y) + 0.374_dp*pow3(y)
         end if
         k(i) = a + b/R_V
      end do
   end subroutine cardelli_form

   ! the F99 construction: the FM90 ultraviolet curve shortward of 2700 A and, longward,
   ! a natural cubic spline through k = 0 at x = 0, the optical/IR anchors
   ! (knot_x, knot_k), and the FM90 curve at 2700 and 2600 A. fm90_c holds the FM90
   ! parameters; FM90 gives E(x-V)/E(B-V), so k = 1 + FM90/R_V
   pure subroutine f99_form(x, R_V, fm90_c, knot_x, knot_k, k)
      real(dp), intent(in) :: x(:), R_V, fm90_c(6), knot_x(:), knot_k(:)
      real(dp), intent(out) :: k(:)

      real(dp), parameter :: x_uv = 1.0d4/2700.0_dp, x_uv2 = 1.0d4/2600.0_dp
      real(dp) :: sx(size(knot_x) + 3), sk(size(knot_x) + 3), m(size(knot_x) + 3), xi
      integer :: i

      sx = [0.0_dp, knot_x, x_uv, x_uv2]
      sk = [0.0_dp, knot_k, 1.0_dp + fm90(x_uv, fm90_c)/R_V, 1.0_dp + fm90(x_uv2, fm90_c)/R_V]
      call natural_spline(sx, sk, m)

      do i = 1, size(x)
         xi = min(x(i), 10.0_dp)
         if (xi >= x_uv) then
            k(i) = 1.0_dp + fm90(xi, fm90_c)/R_V
         else
            k(i) = spline_value(sx, sk, m, xi)
         end if
      end do
   end subroutine f99_form

   ! Fitzpatrick & Massa (1990, ApJS, 72, 163) ultraviolet curve E(x-V)/E(B-V);
   ! c = [c1, c2, c3, c4, x0, gamma]: linear term, 2175 A bump, far-UV rise
   pure real(dp) function fm90(x, c)
      real(dp), intent(in) :: x, c(6)

      fm90 = c(1) + c(2)*x + c(3)*pow2(x)/(pow2(pow2(x) - pow2(c(5))) + pow2(x*c(6)))
      if (x >= 5.9_dp) fm90 = fm90 + c(4)*(0.5392_dp*pow2(x - 5.9_dp) + 0.05644_dp*pow3(x - 5.9_dp))
   end function fm90

   ! G23 near/mid-IR: a is a power law with a smooth break near 4.3 micron plus the two
   ! silicate features, b a single power law
   pure subroutine g23_ir(x, a, b)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: a, b

      real(dp), parameter :: scale = 0.38526_dp, alpha1 = 1.68467_dp, alpha2 = 0.78791_dp
      real(dp), parameter :: break_lam = 4.30578_dp, break_width = 4.78338_dp
      real(dp) :: lam, w

      lam = 1.0_dp/x
      w = smoothstep(lam, break_lam - 0.5_dp*break_width, break_lam + 0.5_dp*break_width)
      a = scale*((1.0_dp - w)*pow(lam, -alpha1) + w*pow(break_lam, alpha2 - alpha1)*pow(lam, -alpha2)) &
          + modified_drude(lam, 0.06652_dp, 9.8434_dp, 2.21205_dp, -0.24703_dp) &
          + modified_drude(lam, 0.0267_dp, 19.58294_dp, 17.0_dp, -0.27_dp)
      b = -1.01251_dp*pow(x, 1.06099_dp)
   end subroutine g23_ir

   ! G23 optical: quartic in x plus three intermediate-scale (Drude) features
   pure subroutine g23_optical(x, a, b)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: a, b

      a = polynomial([-0.35848_dp, 0.7122_dp, 0.08746_dp, -0.05403_dp, 0.00674_dp], x) &
          + drude(x, 0.03893_dp, 2.288_dp, 0.243_dp) &
          + drude(x, 0.02965_dp, 2.054_dp, 0.179_dp) &
          + drude(x, 0.01747_dp, 1.587_dp, 0.243_dp)
      b = polynomial([0.12354_dp, -2.68335_dp, 2.01901_dp, -0.39299_dp, 0.03355_dp], x) &
          + drude(x, 0.18453_dp, 2.288_dp, 0.243_dp) &
          + drude(x, 0.19728_dp, 2.054_dp, 0.179_dp) &
          + drude(x, 0.1713_dp, 1.587_dp, 0.243_dp)
   end subroutine g23_optical

   ! G23 ultraviolet: FM90 curves for a and b
   pure subroutine g23_uv(x, a, b)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: a, b

      a = fm90(x, [0.81297_dp, 0.2775_dp, 1.06295_dp, 0.11303_dp, 4.60_dp, 0.99_dp])
      b = fm90(x, [-2.97868_dp, 1.89808_dp, 3.10334_dp, 0.65484_dp, 4.60_dp, 0.99_dp])
   end subroutine g23_uv

   ! Drude profile with central value amp at x0 and full width fwhm
   pure real(dp) function drude(x, amp, x0, fwhm)
      real(dp), intent(in) :: x, amp, x0, fwhm

      drude = amp*pow2(fwhm/x0)/(pow2(x/x0 - x0/x) + pow2(fwhm/x0))
   end function drude

   ! Drude profile whose width changes across the feature, with asymmetry asym
   pure real(dp) function modified_drude(x, amp, x0, fwhm, asym)
      real(dp), intent(in) :: x, amp, x0, fwhm, asym

      modified_drude = drude(x, amp, x0, 2.0_dp*fwhm/(1.0_dp + exp(asym*(x - x0))))
   end function modified_drude

   ! cubic step from 0 at v = lo to 1 at v = hi
   pure real(dp) function smoothstep(v, lo, hi)
      real(dp), intent(in) :: v, lo, hi

      real(dp) :: t

      t = min(max((v - lo)/(hi - lo), 0.0_dp), 1.0_dp)
      smoothstep = pow2(t)*(3.0_dp - 2.0_dp*t)
   end function smoothstep

   ! sum of c(i)*y**(i-1)
   pure real(dp) function polynomial(c, y)
      real(dp), intent(in) :: c(:), y

      integer :: i

      polynomial = c(size(c))
      do i = size(c) - 1, 1, -1
         polynomial = polynomial*y + c(i)
      end do
   end function polynomial

   ! second derivatives m of the natural cubic spline through (xs, ys), xs ascending
   pure subroutine natural_spline(xs, ys, m)
      real(dp), intent(in) :: xs(:), ys(:)
      real(dp), intent(out) :: m(:)

      real(dp) :: c(size(xs)), r(size(xs)), h_lo, h_hi, w
      integer :: i, n

      n = size(xs)
      m = 0.0_dp
      c = 0.0_dp
      r = 0.0_dp

      ! tridiagonal system for m(2:n-1), with m(1) = m(n) = 0
      do i = 2, n - 1
         h_lo = xs(i) - xs(i - 1)
         h_hi = xs(i + 1) - xs(i)
         w = 2.0_dp*(h_lo + h_hi) - h_lo*c(i - 1)
         c(i) = h_hi/w
         r(i) = (6.0_dp*((ys(i + 1) - ys(i))/h_hi - (ys(i) - ys(i - 1))/h_lo) - h_lo*r(i - 1))/w
      end do
      do i = n - 1, 2, -1
         m(i) = r(i) - c(i)*m(i + 1)
      end do
   end subroutine natural_spline

   ! value at x of the cubic spline with knots (xs, ys) and second derivatives m
   pure real(dp) function spline_value(xs, ys, m, x)
      real(dp), intent(in) :: xs(:), ys(:), m(:), x

      real(dp) :: h, t_lo, t_hi
      integer :: j

      j = 1
      do while (j < size(xs) - 1 .and. x > xs(j + 1))
         j = j + 1
      end do

      h = xs(j + 1) - xs(j)
      t_lo = xs(j + 1) - x
      t_hi = x - xs(j)
      spline_value = (m(j)*pow3(t_lo) + m(j + 1)*pow3(t_hi))/(6.0_dp*h) &
                     + (ys(j)/h - m(j)*h/6.0_dp)*t_lo + (ys(j + 1)/h - m(j + 1)*h/6.0_dp)*t_hi
   end function spline_value

end module extinction
