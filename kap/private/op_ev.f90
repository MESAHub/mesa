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

module op_ev

   use math_lib
   use op_def
   use utils_lib
   use const_def, only: dp
   use kap_def, only: kap_test_partials, kap_test_partials_val, kap_test_partials_dval_dx

   implicit none

   private
   public :: abund
   public :: rd
   public :: ross
   public :: mix
   public :: interp

contains

   subroutine abund(nel, izz, fa, flmu, fmu1, nkz)
      integer, intent(in) :: nel, izz(ipe)
      real, intent(in) :: fa(ipe)
      real, intent(out) :: flmu, fmu1(nrad)
      integer, intent(out) :: nkz(ipe)
      integer :: k, k1, k2, m
      real :: amamu(ipe), fmu, a1, c1, fmu0
      logical :: found

      ! Get k1, get amamu(k)
      do k = 1, nel
         found = .false.
         do m = 1, ipe
            if (izz(k) == kz(m)) then
               amamu(k) = amass(m)
               nkz(k) = m
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            write (*, *) ' k=', k, ', izz(k)=', izz(k)
            call mesa_error(__FILE__, __LINE__, ' kz(m) not found')
         end if
      end do

      ! Mean atomic weight = fmu
      fmu = 0.0d0
      do k = 1, nel
         fmu = fmu + fa(k)*amamu(k)
      end do

      do k2 = 1, nel
         a1 = fa(k2)
         c1 = 1.0d0/(1.0d0 - a1)
         fmu0 = c1*(fmu - fa(k2)*amamu(k2))
         fmu1(k2) = a1*(amamu(k2) - fmu0)*1.660531d-24  ! dmu/dlog xi
      end do

      fmu = fmu*1.660531d-24  ! Convert to cgs
      flmu = log10(dble(fmu))

   end subroutine abund
   ! **********************************************************************

   subroutine rd(nel, nkz, izz, ilab, jh, n_tot, ff, rr, i3, umesh, fac)
      integer, intent(in) :: nel, nkz(ipe), izz(ipe), ilab(4), jh(4), n_tot, i3
      real, intent(in) :: umesh(:)  ! (nptot)
      real(dp), intent(in) :: fac(nel)
      real, intent(out) :: ff(:, :, :, :)  ! (nptot, ipe, 4, 4)
      real, intent(out) :: rr(28, ipe, 4, 4)
      integer :: i, j, k, l, m, n, itt, jnn, izp, ne1, ne2, ne, ib, ia
      real :: fion(-1:28), yb, ya, d
      ! declare variables in common block (instead of by default: real (a-h, o-z), integer (i-n))
      !       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int, ne1p, ne2p, np, kp1, kp2, kp3, npp, mx, nx
      !       real :: umin, umax, epatom, oplnck, fionp, yy1, yy2, yx
      !       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
      !      + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1p(17,91,25),
      !      + ne2p(17,91,25),fionp(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
      !      + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
      !      + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
      !       save /atomdata/

      !  i=temperature index
      !  j=density index
      !  k=frequency index
      !  n=element index
      !  Get:
      !    mono opacity cross-section ff(k,n,i,j)
      !    modified cross-section for selected element, ta(k,i,j)

      rr = 0.0d0

      !  Start loop on i (temperature index)
      do i = 1, 4
         itt = (ilab(i) - ite1)/2 + 1
         do j = 1, 4
            jnn = (jh(j)*i3 - jn1(itt))/2 + 1
            ! Read mono opacities
            do n = 1, nel
               izp = izz(n)
               ne1 = ne1p(nkz(n), itt, jnn)
               ne2 = ne2p(nkz(n), itt, jnn)
               do ne = ne1, ne2
                  fion(ne) = fionp(ne, nkz(n), itt, jnn)
               end do
               do ne = ne1, min(ne2, izp - 2)
                  rr(izp - 1 - ne, n, i, j) = fion(ne)
               end do
               ! call zetbarp(izp, ne1, ne2, fion, zet, i3)
               ! zetal(n, i, j) = zet
               do k = 1, n_tot
                  ff(k, n, i, j) = yy2(k + kp2(nkz(n), itt, jnn))
               end do

               if (fac(n) /= 1d0) then
                  do k = 1, size(ff, dim=1)
                     ff(k, n, i, j) = fac(n)*ff(k, n, i, j)
                  end do
               end if
            end do
         end do
      end do

   end subroutine rd

   subroutine ross(flmu, fmu1, dv, ntot, rs, rossl)
      integer, intent(in) :: ntot
      real, intent(in) :: flmu, dv, fmu1(nrad)
      real, intent(in) :: rs(:, :, :)  ! (nptot, 4, 4)
      real, intent(out) :: rossl(4, 4)
      integer :: i, j, n, k2
      real(dp) :: drs, dd, oross
      real :: fmu, tt, ss, dd2, drsp(nrad), exp10_flmu
      real, parameter :: log10_bohr_radius_sqr = -16.55280d0

      ! oross=cross-section in a.u.
      ! rossl=log10(ROSS in cgs)
      exp10_flmu = real(exp10(dble(flmu)))
      do i = 1, 4
         do j = 1, 4
            drs = 0.d0
            do n = 1, ntot
               dd = 1.d0/rs(n, i, j)
               dd2 = dd*dd
               drs = drs + dd
            end do
            oross = 1.d0/(drs*dv)
            rossl(i, j) = log10(dble(oross)) + log10_bohr_radius_sqr - flmu  !log10(fmu)
         end do
      end do

   end subroutine ross

   subroutine mix(ntot, nel, fa, ff, rs, rr, rion)
      integer, intent(in) :: ntot, nel
      real, intent(in) :: ff(:, :, :, :)  ! (nptot, ipe, 4, 4)
      real, intent(in) :: fa(ipe), rr(28, 17, 4, 4)
      real, intent(out) :: rs(:, :, :)  ! (nptot, 4, 4)
      real, intent(out) :: rion(28, 4, 4)
      integer :: i, j, k, n, m
      real :: rs_temp, rion_temp

      do j = 1, 4
         do i = 1, 4
            do n = 1, ntot
               rs_temp = ff(n, 1, i, j)*fa(1)
               do k = 2, nel
                  rs_temp = rs_temp + ff(n, k, i, j)*fa(k)
               end do
               rs(n, i, j) = rs_temp
               !rs(n, i, j) = dot_product(ff(n,1:nel,i,j),fa(1:nel))
            end do
            do m = 1, 28
               rion_temp = rr(m, 1, i, j)*fa(1)
               do k = 2, nel
                  rion_temp = rion_temp + rr(m, k, i, j)*fa(k)
               end do
               rion(m, i, j) = rion_temp
               !rion(m,i,j) = dot_product(rr(m,1:nel,i,j),fa(1:nel))
            end do
         end do
      end do

   end subroutine mix

   subroutine interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)
      use op_common, only: fint, fintp
      integer, intent(in) :: nel, i3
      real, intent(in) :: ux, uy, xi, eta, rossl(4, 4)
      real, intent(out) :: gx, gy, g
      integer :: i, j, l, k2
      real :: V(4), U(4), vyi(4)

      do i = 1, 4
         do j = 1, 4
            U(j) = rossl(i, j)
         end do
         V(i) = fint(U, eta)
         vyi(i) = fintp(U, eta)
      end do
      g = fint(V, xi)
      gy = fint(vyi, xi)
      gx = fintp(V, xi)

      gy = gy/uy
      gx = (80.0d0/real(i3))*(gx - gy*ux)

   end subroutine interp

end module op_ev
