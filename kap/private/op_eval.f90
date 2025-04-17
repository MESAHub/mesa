! ***********************************************************************
!
!   Copyright (C) 2013-2019  Haili Hu & The MESA Team
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

! FORTRAN 90 module for calculation of radiative accelerations,
! based on the Opacity Project (OP) code "OPserver".
! See CHANGES_HU for changes made to the original code.

! Haili Hu 2010

module op_eval

   use op_load
   use const_def, only: dp
   use math_lib
   use kap_def, only: kap_test_partials, kap_test_partials_val, kap_test_partials_dval_dx

   implicit none

   private
   public :: eval_op_radacc
   public :: eval_op_ev
   public :: eval_alt_op

   logical, parameter :: dbg = .false.

contains

! HH: Based on "op_ax.f"
! Input:   kk = number of elements to calculate g_rad for
!          iz1(kk) = charge of element to calculate g_rad for
!          nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fac(nel) = scale factors for element opacity
!          flux = local radiative flux (Lrad/4*pi*r^2)
!          fltp = log T
!          flrhop = log rho
!          screening   if true, use screening corrections
! Output: g1 = log kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         gp1(kk) = d(log kappa)/d(log xi)
!         grl1(kk) = log grad
!         fx1(kk) = d(log grad)/d(log T)
!         fy1(kk) = d(log grad)/d(log rho)
!         grlp1(kk) = d(log grad)/d(log xi)
!         meanZ(nel) = average ionic charge of elements
!         zetx1(nel) = d(meanZ)/d(log T)
!         zety1(nel) = d(meanZ)/d(log rho)
!         ierr = 0 for correct use
   subroutine eval_op_radacc(kk, izk, nel, izzp, fap, fac, flux, fltp, flrhop, screening, g1, grl1, umesh, semesh, ff, ta, rs, ierr)
      use op_radacc
      use op_load, only: msh
      use op_common
      integer, intent(in) :: kk, nel
      integer, intent(in) :: izk(kk), izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: flux, fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(out) :: g1
      real(dp), intent(inout) :: grl1(kk)
      real, pointer :: umesh(:), semesh(:), ff(:, :, :, :), ta(:, :, :, :), rs(:, :, :)
      ! umesh(nptot)
      ! semesh(nptot)
      ! ff(nptot, ipe, 4, 4)
      ! ta(nptot, nrad, 4, 4),
      ! rs(nptot, 4, 4)
      integer, intent(out) :: ierr
      ! local variables
      integer :: n, i, k2, i3, ntot, jhmin, jhmax
      integer :: ih(4), jh(4), ilab(4), kzz(nrad), nkz(ipe), izz(ipe), iz1(nrad)
      real :: const, gx, gy, flt, flrho, flmu, dscat, dv, xi, flne, epa, eta, ux, uy, g
      real :: uf(0:100), rion(28, 4, 4), rossl(4, 4), flr(4, 4), rr(28, ipe, 4, 4), &
              fa(ipe), gaml(4, 4, nrad), f(nrad), am1(nrad), fmu1(nrad)

      include 'formats'

      ierr = 0
      if (nel <= 0 .or. nel > ipe) then
         write (6, *) 'OP - NUMBER OF ELEMENTS OUT OF RANGE:', nel
         ierr = 1
         return
      end if

      if (dbg) write (*, *) 'start eval_op_radacc'

      ! Get i3 for mesh type q='m'
      i3 = 2

      ! HH: k2 loops over elements for which to calculate grad.
      do k2 = 1, kk
         do n = 1, ipe
            if (izk(k2) == kz(n)) then
               iz1(k2) = izk(k2)
               exit
            end if
            if (n == ipe) then
               write (6, *) 'OP - SELECTED ELEMENT CANNOT BE TREATED: Z = ', izk(k2)
               write (6, *) 'OP supports C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni'
               ierr = 5
               return
            end if
         end do
      end do

      outer: do i = 1, nel
         inner: do n = 1, ipe
            if (izzp(i) == kz(n)) then
               izz(i) = izzp(i)
               fa(i) = fap(i)
               if (fa(i) < 0.0) then
                  write (6, *) 'OP - NEGATIVE FRACTIONAL ABUNDANCE:', fa(i)
                  ierr = 7
                  return
               end if
               cycle outer
            end if
         end do inner
         write (6, *) 'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ', izzp(i)
         ierr = 8
         return
      end do outer

      ! Calculate mean atomic weight (flmu) and
      ! array kzz indicating elements for which to calculate g_rad
      if (dbg) write (*, *) 'call abund'
      call abund(nel, izz, kk, iz1, fa, kzz, flmu, am1, fmu1, nkz)

      ! Other initializations
      !  dv = interval in frequency variable v
      !  ntot=number of frequency points
      !  umesh, values of u=(h*nu/k*T) on mesh points
      !  semesh, values of 1-exp(-u) on mesh points
      !  uf, dscat used in scattering correction
      if (dbg) write (*, *) 'call msh'
      call msh(dv, ntot, umesh, semesh, uf, dscat)  !output variables

      ! Start loop on temperature-density points
      ! flt=log10(T, K)
      ! flrho=log10(rho, cgs)

      flt = fltp
      flrho = flrhop
      ! Get temperature indices
      !  Let ite(i) be temperature index used in mono files
      !  Put ite(i)=2*ih(i)
      !  Use ih(i), i=1 to 4
      !  xi=interpolation variable
      !  log10(T)=flt=0.025*(ite(1)+xi+3)
      !  ilab(i) is temperature label
      if (dbg) write (*, *) 'call xindex'
      call xindex(flt, ilab, xi, ih, i3, ierr)
      if (ierr /= 0) then
         write (*, *) "xindex errored in radacc"
         return
      end if

      ! Get density indices
      !  Let jne(j) be density index used  in mono files
      !  Put jne(j)=2*jh(j)
      !  Use jh(j), j=1 to 4
      !  Get extreme range for jh
      if (dbg) write (*, *) 'call jrange'
      call jrange(ih, jhmin, jhmax, i3)

      ! Get electron density flne=log10(Ne) for specified mass density flrho
      !  Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)]
      !    epa=electrons per atom
      if (dbg) write (*, *) 'call findne'
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt, xi, flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr /= 0) then
         write (*, *) "findne encountered error in radacc"
         return
      end if

      ! Get density indices jh(j), j=1 to 4,
      !  Interpolation variable eta
      !  log10(Ne)=flne=0.25*(jne(1)+eta+3)
      if (dbg) write (*, *) 'call yindex'
      call yindex(jhmin, jhmax, flne, jh, i3, eta)

      ! Get ux=0.025*[d log10(rho)]/[d log10(T)]
      if (dbg) write (*, *) 'call findux'
      call findux(flr, xi, eta, ux)

      ! rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
      ! Get new mono opacities, ff(n,k,i,j)
      if (dbg) write (*, *) 'call rd'
      call rd(i3, kk, kzz, nel, nkz, izz, ilab, jh, ntot, umesh, semesh, ff, rr, ta, fac)

      ! Get rs = weighted sum of monochromatic opacity cross sections
      if (dbg) write (*, *) 'call mix'
      call mix(kk, kzz, ntot, nel, fa, ff, rr, rs, rion)

      ! Screening corrections
      if (screening) then
         ! Get Boercker scattering correction
         if (dbg) write (*, *) 'call scatt'
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr /= 0) return
         ! Get correction for Debye screening
         if (dbg) write (*, *) 'call screen1'
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      end if

      ! Get rossl, array of log10(Rosseland mean in cgs)
      if (dbg) write (*, *) 'call ross'
      call ross(kk, flmu, fmu1, dv, ntot, rs, rossl, gaml, ta)

      ! Interpolate to required flt, flrho
      ! g=log10(ross, cgs)
      if (dbg) write (*, *) 'call interp'
      call interp(nel, kk, rossl, gaml, xi, eta, g, i3, f)
      if (dbg) write (*, *) 'done interp'

      ! Write grad in terms of local radiative flux instead of (Teff, r/R*):
      const = 13.30295d0 + log10(flux)  ! = -log10(c) - log10(amu) + log10(flux)
      do k2 = 1, kk
         grl1(k2) = const + flmu - log10(dble(am1(k2))) + f(k2) + g  ! log g_rad
      end do

      g1 = g  ! log kappa

      if (dbg) write (*, *) 'done eval_op_radacc'

   end subroutine eval_op_radacc

! ***********************************************************************
! HH: Based on "op_mx.f", opacity calculations to be used for stellar evolution calculations
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fac(nel) = scale factors for element opacity
!          fltp = log (temperature)
!          flrhop = log (mass density)
!          screening   if true, use screening corrections
! Output: g1 = log kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
   subroutine eval_op_ev(nel, izzp, fap, fac, fltp, flrhop, screening, g1, gx1, gy1, umesh, semesh, ff, rs, ierr)
      use op_ev
      use op_load, only: msh
      use op_common
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(inout) :: g1, gx1, gy1
      real, pointer :: umesh(:), semesh(:), ff(:, :, :, :), rs(:, :, :)
      ! umesh(nptot)
      ! semesh(nptot)
      ! ff(nptot, ipe, 4, 4)
      ! rs(nptot, 4, 4)
      ! s(nptot, nrad, 4, 4)
      integer, intent(out) :: ierr
      ! local variables
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(4), jh(4), ilab(4), izz(ipe), nkz(ipe)
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g, eta, epa, xi, ux, uy
      real :: uf(0:100), rion(28, 1:4, 1:4), rossl(4, 4), flr(4, 4), fa(ipe), rr(28, ipe, 4, 4), fmu1(nrad)

      ! Initializations
      ierr = 0
      if (nel <= 0 .or. nel > ipe) then
         write (6, *) 'OP - NUMBER OF ELEMENTS OUT OF RANGE:', nel
         ierr = 1
         return
      end if

      ! Get i3 for mesh type q='m'
      i3 = 2

      outer: do i = 1, nel
         inner: do n = 1, ipe
            if (izzp(i) == kz(n)) then
               izz(i) = izzp(i)
               fa(i) = fap(i)
               if (fa(i) < 0.d0) then
                  write (6, *) 'OP - NEGATIVE FRACTIONAL ABUNDANCE:', fa(i)
                  ierr = 7
                  return
               end if
               cycle outer
            end if
         end do inner
         write (6, *) 'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ', izzp(i)
         ierr = 8
         return
      end do outer

      ! Calculate mean atomic weight (flmu)
      call abund(nel, izz, fa, flmu, fmu1, nkz)

!  Other initializations
!       dv = interval in frequency variable v
!       ntot=number of frequency points
!       umesh, values of u=(h*nu/k*T) on mesh points
!       semesh, values of 1-exp(-u) on mesh points
!       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, semesh, uf, dscat)

!  Start loop on temperature-density points
!  flt=log10(T, K)
!  flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
!     Get temperature indices
!       Let ite(i) be temperature index used in mono files
!       Put ite(i)=2*ih(i)
!       Use ih(i), i=1 to 4
!       xi=interpolation variable
!       log10(T)=flt=0.025*(ite(1)+xi+3)
!       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3, ierr)
      if (ierr /= 0) then
         write (*, *) "eval_op_ev failed in xindex"
         return
      end if

!     Get density indices
!       Let jne(j) be density index used  in mono files
!       Put jne(j)=2*jh(j)
!       Use jh(j), j=1 to 4
!       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)

!     Get electron density flne=log10(Ne) for specified mass density flrho
!       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)]
!              epa=electrons per atom
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt, xi, flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr /= 0) then
         write (*, *) "eval_op_ev failed in findne"
         return
      end if

!     Get density indices jh(j), j=1 to 4,
!       Interpolation variable eta
!       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)

!     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)

!    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
!     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, nkz, izz, ilab, jh, ntot, ff, rr, i3, umesh, fac)

!     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion)

      if (screening) then
!        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr /= 0) then
            write (*, *) "scattering correction failed in op_eval_ev"
            return
         end if
!        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      end if

!     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, fmu1, dv, ntot, rs, rossl)

!     Interpolate to required flt, flrho
!     g=log10(ross, cgs)
      call interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)

      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho

   end subroutine eval_op_ev

! ***********************************************************************

!     HH: Based on "op_mx.f", opacity calculations to be used for non-adiabatic pulsation calculations
! Special care is taken to ensure smoothness of opacity derivatives
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fac(nel) = scale factors for element opacity
!          fltp = log (temperature)
!          flrhop = log (mass density)
!          screening   if true, use screening corrections
! Output: g1 = log kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
   subroutine eval_alt_op(nel, izzp, fap, fac, fltp, flrhop, screening, g1, gx1, gy1, umesh, semesh, ff, rs, ierr)
      use op_osc
      use op_load, only: msh
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(out) :: g1, gx1, gy1
      ! real(dp), intent(out) :: meanZ(nel)
      real, pointer :: umesh(:), semesh(:), ff(:, :, :, :), rs(:, :, :)
      ! umesh(nptot)
      ! semesh(nptot)
      ! ff(nptot, ipe, 0:5, 0:5)
      ! rs(nptot, 0:5, 0:5)
      integer, intent(out) :: ierr
      ! local variables
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(0:5), jh(0:5), ilab(0:5), izz(ipe), nkz(ipe)
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g, eta, epa, xi, ux, uy
      real :: uf(0:100), rion(28, 0:5, 0:5), rossl(0:5, 0:5), flr(4, 4), fa(ipe), rr(28, ipe, 0:5, 0:5)

      ! Initializations
      ierr = 0
      if (nel <= 0 .or. nel > ipe) then
         write (6, *) 'OP - NUMBER OF ELEMENTS OUT OF RANGE:', nel
         ierr = 1
         return
      end if

      ! Get i3 for mesh type q='m'
      i3 = 2

      outer: do i = 1, nel
         inner: do n = 1, ipe
            if (izzp(i) == kz(n)) then
               izz(i) = izzp(i)
               fa(i) = fap(i)
               if (fa(i) < 0.0) then
                  write (6, *) 'OP - NEGATIVE FRACTIONAL ABUNDANCE:', fa(i)
                  ierr = 7
                  return
               end if
               cycle outer
            end if
         end do inner
         write (6, *) 'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ', izzp(i)
         ierr = 8
         return
      end do outer

      ! Calculate mean atomic weight (flmu)
      call abund(nel, izz, fa, flmu, nkz)

!  Other initializations
!       dv = interval in frequency variable v
!       ntot=number of frequency points
!       umesh, values of u=(h*nu/k*T) on mesh points
!       semesh, values of 1-exp(-u) on mesh points
!       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, semesh, uf, dscat)

      ! Start loop on temperature-density points
      ! flt=log10(T, K)
      ! flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
!     Get temperature indices
!       Let ite(i) be temperature index used in mono files
!       Put ite(i)=2*ih(i)
!       Use ih(i), i=1 to 4
!       xi=interpolation variable
!       log10(T)=flt=0.025*(ite(1)+xi+3)
!       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3, ierr)
      if (ierr /= 0) return

!     Get density indices
!       Let jne(j) be density index used  in mono files
!       Put jne(j)=2*jh(j)
!       Use jh(j), j=1 to 4
!       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)

!     Get electron density flne=log10(Ne) for specified mass density flrho
!       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)]
!              epa=electrons per atom
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt, xi, flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr /= 0) return

!     Get density indices jh(j), j=1 to 4,
!       Interpolation variable eta
!       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)

!     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)

!     rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
!     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, nkz, izz, ilab, jh, ntot, ff, rr, i3, umesh, fac)

!     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion)

      if (screening) then
!        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr /= 0) return
!        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      end if

!     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, dv, ntot, rs, rossl)

!     Interpolate to required flt, flrho
!     g=log10(ross, cgs)
      call interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)

      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho

   end subroutine eval_alt_op

end module op_eval
