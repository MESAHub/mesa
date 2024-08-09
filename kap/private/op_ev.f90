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
 
module op_ev
   use math_lib
   use op_def
   use const_def
   use kap_def, only: kap_test_partials, kap_test_partials_val, kap_test_partials_dval_dx

   implicit none

contains

  subroutine abund(nel, izz, fa, flmu, fmu1, nkz)
     integer, intent(in) :: nel, izz(ipe)
     real, intent(in) :: fa(ipe)
     real, intent(out) :: flmu, fmu1(nrad)
     integer, intent(out) :: nkz(ipe)

      integer :: k, k1, k2, m
      real :: amamu(ipe), fmu, a1, c1, fmu0

      ! Get k1,get amamu(k)
      do k = 1, nel              
         do m = 1, ipe
            if(izz(k).eq.kz(m))then
               amamu(k) = amass(m)
               nkz(k) = m
               goto 1
            end if
         end do  
         write(*,*) ' k=',k,', izz(k)=',izz(k)
         write(*,*) ' kz(m) not found'
         stop
    1    continue           
      end do
   
      ! Mean atomic weight = fmu
      fmu = 0d0
      do k = 1, nel
         fmu = fmu + fa(k)*amamu(k)
      end do

      do k2 = 1, nel
         a1 = fa(k2)
         c1 = 1d0/(1d0-a1)
         fmu0 = c1*(fmu - fa(k2)*amamu(k2))
         fmu1(k2) = a1*(amamu(k2)-fmu0)*1.660531d-24 !dmu/dlog xi
      end do   

      fmu = fmu*1.660531d-24 ! Convert to cgs
      flmu = log10(dble(fmu))

      return
   end subroutine abund


   subroutine rd(nel, nkz, izz, ilab, jh, n_tot, ff, rr, i3, umesh, fac)
      integer, intent(in) :: nel, nkz(ipe), izz(ipe), ilab(4), jh(4), n_tot, i3
      real, intent(in) :: umesh(:) ! (nptot)
      real(dp), intent(in) :: fac(nel)
      real, intent(out) :: ff(:,:,:,:) ! (nptot, ipe, 4, 4)
      real, intent(out) :: rr(28, ipe, 4, 4)
      ! local variables      
      integer :: i, j, k, l, m, n, itt, jnn, izp, ne1, ne2, ne, ib, ia
      real :: fion(-1:28), yb, ya, d
      ! declare variables in common block (instead of by default: real (a-h, o-z), integer (i-n))   
      !       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int,
      !      : ne1p, ne2p, np, kp1, kp2, kp3, npp, mx, nx   
      !       real :: umin, umax, epatom, oplnck, fionp, yy1, yy2, yx    
      !       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
      !      + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1p(17,91,25),
      !      + ne2p(17,91,25),fionp(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
      !      + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
      !      + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)  
      !       save /atomdata/

      !  i=temperature inex
      !  j=density index
      !  k=frequency index
      !  n=element index
      !  Get:
      !    mono opacity cross-section ff(k,n,i,j)
      !    modified cross-section for selected element, ta(k,i,j)
      !
      !     Initialisations    
      rr = 0d0

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
               do ne = ne1, min(ne2, izp-2)
                  rr(izp-1-ne, n, i, j) = fion(ne)
               end do                                       
               ! call zetbarp(izp, ne1, ne2, fion, zet, i3)
               ! zetal(n, i, j) = zet            
               do k = 1, n_tot
                  ff(k, n, i, j) = yy2(k+kp2(nkz(n), itt, jnn))
               end do  
              
              if (fac(n) /= 1d0) then
                 do k = 1, size(ff,dim=1)
                    ff(k,n,i,j) = fac(n)*ff(k,n,i,j)
                 end do
              end if
            end do !n
         end do !j
      end do !i

      return

   end subroutine rd


   subroutine ross(flmu, fmu1, dv, ntot, rs, rossl) 
      integer, intent(in) :: ntot
      real, intent(in) :: flmu, dv, fmu1(nrad)
      real, intent(in) :: rs(:,:,:) ! (nptot, 4, 4)
      real, intent(out) :: rossl(4, 4)
      integer :: i, j, n, k2
      real(dp) :: drs, dd, oross
      real :: fmu, tt, ss, dd2, drsp(nrad), exp10_flmu
      real :: log10_bohr_radius_sqr = -16.55280d0

      !  oross=cross-section in a.u.
      !  rossl=log10(ROSS in cgs)
      exp10_flmu = real(exp10(dble(flmu)))
      do i = 1, 4
         do j = 1, 4
            drs = 0.d0
            do n = 1, ntot
               dd = 1d0/rs(n, i, j)
               dd2 = dd*dd   
               drs = drs + dd
            end do            
            oross = 1d0/(drs*dv)
            rossl(i, j) = log10(dble(oross)) + log10_bohr_radius_sqr - flmu !log10(fmu) 
         end do !j
      end do !i

      return
   end subroutine ross


   subroutine mix(ntot, nel, fa, ff, rs, rr, rion)
      integer, intent(in) :: ntot, nel
      real, intent(in) :: ff(:,:,:,:) ! (nptot, ipe, 4, 4)
      real, intent(in) :: fa(ipe), rr(28, 17, 4, 4)
      real, intent(out) :: rs(:,:,:) ! (nptot, 4, 4)
      real, intent(out) :: rion(28, 4, 4)

      ! local variables      
      integer :: i, j, k, n, m
      real :: rs_temp, rion_temp

      do j = 1, 4
         do i = 1, 4
            do n = 1, ntot
               rs_temp = ff(n,1,i,j)*fa(1)
               do k = 2, nel
                  rs_temp = rs_temp + ff(n,k,i,j)*fa(k)
               end do
               rs(n, i, j) = rs_temp  
               !rs(n, i, j) = dot_product(ff(n,1:nel,i,j),fa(1:nel))
            end do            
            do m = 1, 28
               rion_temp = rr(m, 1, i, j)*fa(1)
               do k = 2, nel
                  rion_temp = rion_temp + rr(m,k,i,j)*fa(k)
               end do
               rion(m,i,j) = rion_temp
               !rion(m,i,j) = dot_product(rr(m,1:nel,i,j),fa(1:nel))
            end do   
         end do
      end do

      return
   end subroutine mix


   subroutine interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)
      use op_common, only: fint, fintp

      integer, intent(in) :: nel, i3
      real, intent(in) :: ux, uy, xi, eta, rossl(4, 4)
      real, intent(out) :: gx, gy, g
      ! local variables
      integer :: i, j, l, k2
      real :: V(4), U(4),  vyi(4)  

      do i = 1, 4
         do j = 1, 4
            U(J) = rossl(i,j)
         end do
         V(I) = fint(U,eta)
         vyi(i) = fintp(u,eta)
      end do
      g = fint(V, xi)
      gy = fint(vyi, xi)
      gx = fintp(v, xi)

      gy=gy/uy
      gx=(80d0/real(I3))*(gx-gy*ux)
      
      return
   end subroutine interp

end module op_ev
