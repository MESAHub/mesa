! ***********************************************************************
!
!   Copyright (C) 2013-2019  Haili Hu & The MESA Team
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

      module op_radacc

      use math_lib
      use op_def
      use const_def
    !  use op_load, only: nptot, ipe

      logical, parameter :: dbg = .false.



      contains
c***********************************************************************1
      subroutine abund(nel, izz, kk, iz1, fa, kzz, flmu, am1, fmu1, nkz)
      implicit none
      integer, intent(in) :: nel, izz(ipe), kk, iz1(nrad)
      real, intent(in) :: fa(ipe)
      integer,  intent(out) :: kzz(nrad), nkz(ipe)
      real, intent(out) :: flmu, am1(nrad), fmu1(nrad)
c local variables
      integer :: k, k1, k2, m
      real :: fmu, a1, c1, fmu0, amamu(ipe)
c
c Get k1,get amamu(k)
      do k = 1, nel
         do m = 1, ipe
            if(izz(k).eq.kz(m))then
               amamu(k) = amass(m)
               nkz(k) = m
            goto 1
            endif
         end do
         print*,' k=',k,', izz(k)=',izz(k)
         print*,' kz(m) not found'
         stop
    1    continue
      end do
c
      k1=-1
      do k2 = 1, kk
         do k = 1, nel
            if(izz(k).eq.iz1(k2))then
               k1 = k
            endif
         end do
         kzz(k2) = k1
         am1(k2) = amamu(k1)
      end do
c
c  Mean atomic weight = fmu
      fmu = 0d0
      do k = 1, nel
         fmu = fmu + fa(k)*amamu(k)
      end do
c
      do k2 = 1, kk
         a1 = fa(kzz(k2))
         c1 = 1d0/(1d0-a1)
         fmu0 = c1*(fmu - fa(kzz(k2))*amamu(kzz(k2)))
         fmu1(k2) = a1*(am1(k2)-fmu0)*1.660531d-24 !dmu/dlog xi
      end do
c
      fmu = fmu*1.660531d-24 ! Convert to cgs
      flmu = log10(dble(fmu))
c
      return
      end subroutine abund

c**********************************************************************
      subroutine rd(i3, kk, kzz, nel, nkz, izz, ilab, jh, ntot, umesh,
     : semesh, ff, rr, ta, fac)
      implicit none
      integer, intent(in) :: i3, kk, kzz(nrad), nel, nkz(ipe),
     : izz(ipe), ilab(4), jh(4), ntot
      real, intent(in) :: umesh(:), semesh(:) ! (nptot)
      real(dp), intent(in) :: fac(nel)
      real, intent(out) :: ff(:,:,:,:) ! (nptot, ipe, 4, 4)
      real, intent(out) :: ta(:,:,:,:) ! (nptot, nrad, 4, 4)
      real, intent(out) :: rr(28, ipe, 4, 4)
c local variables
      integer :: i, j, k, k2, l, n, m, itt, jnn, izp, ne1, ne2, ne, ib, ia
      real :: ya, yb, d, se, u, fion(-1:28) !, ff_temp(nptot), ta_temp(nptot)

c declare variables in common block, by default: real (a-h, o-z), integer (i-n)
!       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntotp, nc, nf, int,
!      : ne1p, ne2p, np, kp1, kp2, kp3, npp, mx, nx
!       real :: umin, umax, epatom, oplnck, fionp, yy1, yy2, yx
!       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntotp,
!      + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1p(17,91,25),
!      + ne2p(17,91,25),fionp(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
!      + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
!      + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
!       save /atomdata/

      include 'formats'

c
c  i=temperature inex
c  j=density index
c  k=frequency index
c  n=element index
c  Get:
c    mono opacity cross-section ff(k,n,i,j)
c    modified cross-section for selected element, ta(k,i,j)
c    zet(i,j) for diffusion coefficient
c
c     Initialisations
      if (dbg) then
         write(*,2) 'size(ff,dim=1)', size(ff,dim=1)
         write(*,2) 'size(ff,dim=2)', size(ff,dim=2)
         write(*,2) 'size(ff,dim=3)', size(ff,dim=3)
         write(*,2) 'size(ff,dim=4)', size(ff,dim=4)
         write(*,2) 'size(ta,dim=1)', size(ta,dim=1)
         write(*,2) 'size(ta,dim=2)', size(ta,dim=2)
         write(*,2) 'size(ta,dim=3)', size(ta,dim=3)
         write(*,2) 'size(ta,dim=4)', size(ta,dim=4)
         write(*,2) 'size(rr,dim=1)', size(rr,dim=1)
         write(*,2) 'size(rr,dim=2)', size(rr,dim=2)
         write(*,2) 'size(rr,dim=3)', size(rr,dim=3)
         write(*,2) 'size(rr,dim=4)', size(rr,dim=4)
         write(*,2) 'ipe', ipe
         write(*,*) 'rd: ff = 0'
      end if

      ff = 0d0
      if (dbg) write(*,*) 'rd: rr = 0'
      rr = 0d0
      if (dbg) write(*,*) 'rd: ta = 0'
      ta = 0d0
c
c  Start loop on i (temperature index)
      do i = 1, 4
         if (dbg) write(*,2) 'rd: i', i
         itt = (ilab(i) - ite1)/2 + 1
c        Read mono opacities
         do j = 1, 4
            if (dbg) write(*,2) 'rd: j', j
            jnn = (jh(j)*i3 - jn1(itt))/2 + 1
            do n = 1, nel
               if (dbg) write(*,2) 'rd: n', n
               izp = izz(n)
               ne1 = ne1p(nkz(n), itt, jnn)
               ne2 = ne2p(nkz(n), itt, jnn)
               do ne = ne1, ne2
                  fion(ne) = fionp(ne, nkz(n), itt, jnn)
                  if (ne .le. min(ne2, izp-2)) rr(izp-1-ne, n, i, j) = fion(ne)
               end do

               do k = 1, ntot
                  ff(k, n, i, j) = yy2(k+kp2(nkz(n), itt, jnn))
               end do

               if (fac(n) /= 1d0) then
                  do k = 1, size(ff,dim=1)
                     ff(k,n,i,j) = fac(n)*ff(k,n,i,j)
                  end do
               end if

               do k2 = 1, kk
                  if (dbg) write(*,2) 'rd: k2', k2
                  if (kzz(k2)== n) then                                    
                     ib = 1
                     yb = yx(1+kp3(nkz(kzz(k2)), itt, jnn))
                     u = umesh(ib)
                     se = semesh(ib)
                     ta(1, k2, i, j) = se*ff(1, n, i, j) - yb
                     do m = 2, npp(nkz(kzz(k2)), itt, jnn)
                        ia = ib
                        ya = yb
                        ib = nx(m+kp3(nkz(kzz(k2)), itt, jnn))
                        yb = yx(m+kp3(nkz(kzz(k2)), itt, jnn))
                        d = (yb-ya)/float(ib-ia)
                        do l = ia+1, ib-1
                           u = umesh(l)
                           se = semesh(l)
                           ta(l, k2, i, j) = se*ff(l, n, i, j) -(ya + (l-ia)*d)
                        end do
                        u = umesh(ib)
                        se = semesh(ib)
                        ta(ib, k2, i, j) = se*ff(ib, n, i, j) - yb                        
                     end do
                     goto 101  ! get out of k2-loop
                  endif
               end do !k2
 101           continue
            end do !n
         end do !j
      end do !i
c
      return
c
      end subroutine rd


c***********************************************************************
      subroutine mix(kk, kzz, ntot, nel, fa, ff, rr, rs, rion)
      implicit none
      integer, intent(in) :: kk, kzz(nrad), ntot, nel
      real, intent(in) :: fa(ipe), rr(28, 17, 4, 4)
      real, intent(in) :: ff(:,:,:,:) ! (nptot, ipe, 4, 4)
      real, intent(out) :: rs(:,:,:) ! (nptot, 4, 4)
      real, intent(out) :: rion(28, 4, 4)
c local variables
      integer :: i, j, k, n, m, k2

c
      do i = 1, 4
         do j = 1, 4
            do n = 1, ntot
               rs(n, i, j) = dot_product(ff(n,1:nel,i,j),fa(1:nel))
            end do

            do m = 1, 28
               rion(m,i,j) = dot_product(rr(m,1:nel,i,j),fa(1:nel))
            end do
         end do
      end do
c
      return
      end subroutine mix
c***********************************************************************
      subroutine ross(kk, flmu, fmu1, dv, ntot, rs, rossl, gaml, ta)
      implicit none
      integer, intent(in) :: kk, ntot
      real, intent(in) :: rs(:,:,:) ! (nptot, 4, 4)
      real, intent(in) :: ta(:,:,:,:) ! (nptot, nrad, 4, 4)
      real, intent(in) :: flmu, dv, fmu1(nrad)
      real, intent(out) :: rossl(4, 4), gaml(4, 4, nrad)

c local variables
      integer :: k2, i, j, n
      real(dp) :: drs, dd,  dgm(nrad)
      real :: fmu, oross, tt, exp10_flmu
c
c  oross=cross-section in a.u.
c  rossl=log10(ROSS in cgs)
         exp10_flmu = real(exp10(dble(flmu)))
         do i = 1, 4
            do j = 1, 4
               drs = 0d0
               dgm(:) = 0d0
               do n = 1, ntot   !10000
                  dd = 1d0/rs(n, i, j)
                  drs = drs + dd
                  do k2 = 1, kk
                     tt = ta(n, k2, i, j)
                     dgm(k2) = dgm(k2) + tt*dd
                  end do
               end do
               oross = 1d0/(drs*dv)
               rossl(i, j) = log10(dble(oross)) - 16.55280d0 - flmu
               do k2 = 1, kk
                  if(dgm(k2).gt.0) then
                     dgm(k2) = dgm(k2)*dv
                     gaml(i, j, k2) = log10(dgm(k2))
                  else
                     gaml(i, j, k2) = -30d0
                  endif
               end do
            end do !j
         end do !i
c
      return
      end subroutine ross
c***********************************************************************
      subroutine interp(nel, kk, rossl, gaml, xi, eta, g, i3, f)
      use op_common, only: fint, fintp
      implicit none
      integer, intent(in) :: nel, kk, i3
      real, intent(in) :: gaml(4, 4, nrad), rossl(4, 4), eta
      real, intent(out) :: f(nrad), g
      integer :: l, i, j, k2
      real ::  v(4), u(4), vyi(4), xi

c     interpolation of g (=rosseland mean opacity)
      do i = 1, 4
         do j = 1, 4
            u(j) = rossl(i, j)
         end do
         v(i) = fint(u, eta)
         vyi(i) = fintp(u, eta)
      end do
      g = fint(v, xi)

      do k2 = 1, kk
c     Interpolation of gam (=dimensionless parameter in g_rad)
c     f=log10(gamma)
         do i = 1, 4
            do j = 1, 4
! HH: This gives irregularities, perhaps it is preferable to assign nonnegative values of 
! neighbouring interpolation points (in subroutine "ross") ?
               u(j) = gaml(i, j, k2)
            end do
            v(i) = fint(u, eta)
            vyi(i) = fintp(u, eta)
         end do

         f(k2) = fint(v, xi)

    1 continue
      end do !k2
c
      return
      end subroutine interp

c***********************************************************************
      end module op_radacc
