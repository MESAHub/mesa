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
 
      module op_osc
      use math_lib
      use op_def
      use const_def
      use kap_def, only: kap_test_partials, kap_test_partials_val, kap_test_partials_dval_dx

      contains
c***********************************************************************
      subroutine abund(nel, izz, fa, flmu, nkz)
      implicit none
      integer, intent(in) :: nel, izz(ipe)
      real, intent(in) :: fa(ipe)
      real, intent(out) :: flmu
      integer, intent(out) :: nkz(ipe)
c local variables     
      integer :: k, k1, k2, m
      real :: amamu(ipe), fmu
c
c Get k1,get amamu(k)
      do k = 1, nel              
         do m = 1, ipe
            if(izz(k).eq.kz(m))then
               amamu(k) = amass(m)
               nkz(k) = m
            goto 1
            endif
         enddo  
         print*,' k=',k,', izz(k)=',izz(k)
         print*,' kz(m) not found'
         stop
    1    continue           
      enddo      
c   
c  Mean atomic weight = fmu
      fmu = 0.
      do k = 1, nel
         fmu = fmu + fa(k)*amamu(k)
      enddo
c
      fmu = fmu*1.660531e-24 ! Convert to cgs
      flmu = log10(dble(fmu))
c  
      return
      end subroutine abund
c**********************************************************************
      subroutine xindex(flt, ilab, xi, ih, i3, ierr)
      implicit none
      integer, intent(in) :: i3
      real, intent(in) :: flt
      integer, intent(out) :: ih(0:5), ilab(0:5)
      real, intent(out) :: xi
      integer, intent(out) :: ierr
      integer :: i, ih2
      real :: x
c
      ierr = 0
      if(flt.lt.3.5) then
        ierr = 102
        return
      elseif(flt.gt.8.) then
        ierr = 102
        return
      endif
c      
      x = 40.*flt/real(i3)
      ih2 = x
      ih2 = max(ih2, 140/i3+2)
      ih2 = min(ih2, 320/i3-3)
      do i = 0, 5
         ih(i) = ih2 + i - 2
         ilab(i) = i3*ih(i)
      enddo         
      xi = 2.*(x-ih2) - 1
c
      return
      end subroutine xindex
c**********************************************************************
      subroutine jrange(ih, jhmin, jhmax, i3) 
      implicit none
      integer, intent(in) :: ih(0:5), i3
      integer, intent(out) :: jhmin, jhmax  
      integer :: i
c      
      jhmin = 0
      jhmax = 1000
      do i = 0, 5
         jhmin = max(jhmin, js(ih(i)*i3)/i3)
         jhmax = min(jhmax, je(ih(i)*i3)/i3)
      enddo
c      
      return
      end subroutine jrange
c**********************************************************************
      subroutine findne(ilab, fa, nel, nkz, jhmin, jhmax, ih,
     + flrho, flt, xi, flne, flmu, flr, epa, uy, i3, ierr)
      use op_load, only : solve
      implicit none
      integer, intent(in) :: ilab(0:5), nel, nkz(ipe), jhmin, 
     >      ih(0:5), i3
      integer, intent(inout) :: jhmax
      integer, intent(out) :: ierr
      real, intent(in) :: fa(ipe), flt, xi, flmu
      real,intent(out) :: flne, uy, epa
      real, intent(inout) :: flrho
c local variables
      integer :: i, j, n, jh, jm, itt, jne, jnn
      real :: flrmin, flrmax, flr(4,4), uyi(4), efa(0:5, 7:118), 
     : flrh(0:5, 7:118), u(4), flnei(4), y, zeta, efa_temp
c declare variables in common block, by default: real (a-h, o-z), integer (i-n)      
!       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int,
!      : ne1, ne2, np, kp1, kp2, kp3, npp, mx, nx   
!       real :: umin, umax, epatom, oplnck, fion, yy1, yy2, yx
!       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
!      + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25),
!      + ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
!      + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
!      + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
!       save /atomdata/
c
c  efa(i,jh)=sum_n epa(i,jh,n)*fa(n)
c  flrh(i,jh)=log10(rho(i,jh))
c
c  Get efa    
      do i = 0, 5
         itt = (ilab(i)-ite1)/2 + 1
         do jne = jn1(itt), jn2(itt), i3
            jnn = (jne-jn1(itt))/2 + 1
            jh = jne/i3
            efa_temp = 0.
            do n = 1, nel               
               efa_temp = efa_temp + epatom(nkz(n), itt, jnn)*fa(n)
            enddo !n    
            efa(i, jh) = efa_temp             
         enddo !jne
      enddo !i 
c        
c Get range for efa.gt.0
      do i = 0, 5
         do jh = jhmin, jhmax
            if(efa(i, jh) .le. 0.)then
               jm = jh - 1
               goto 3
            endif
         enddo
         goto 4
    3    jhmax = MIN(jhmax, jm)
    4    continue
      enddo
c    
c Get flrh
      do jh = jhmin,jhmax
        do i = 0,5
         flrh(i, jh) = flmu + 0.25*i3*jh - log10(dble(efa(i,jh)))
        enddo
      enddo  
c
c  Find flrmin and flrmax
      flrmin = -1000
      flrmax = 1000
      do i = 0, 5
         flrmin = max(flrmin, flrh(i,jhmin))
         flrmax = min(flrmax, flrh(i,jhmax))
      enddo
c      
c  Check range of flrho
      if(flrho .lt. flrmin .or. flrho .gt. flrmax)then
         ierr = 101
         return
      endif
c
c  Interpolations in j for flne
      do jh = jhmin, jhmax
         if(flrh(2,jh) .gt. flrho)then
            jm = jh - 1
            goto 5
         endif
      enddo
      print*,' Interpolations in j for flne'
      print*,' Not found, i=',i
      stop
    5 jm=max(jm,jhmin+1)
      jm=min(jm,jhmax-2)
c      
      do i = 1, 4
         do j = 1, 4
            u(j) = flrh(i, jm+j-2)
            flr(i,j) = flrh(i, jm+j-2)
         enddo
         call solve(u, flrho, zeta, uyi(i), ierr)
         if (ierr /= 0) return
         y = jm + 0.5*(zeta+1)
         flnei(i) = .25*i3*y
      enddo  
c
c  Interpolations in i
      flne = fint(flnei, xi)
      uy = fint(uyi, xi)
c Get epa      
      epa = exp10(dble(flne + flmu - flrho))  
c      
      return
c   
  601 format(' For flt=',1p,e11.3,', flrho=',e11.3,' is out of range'/
     +  ' Allowed range for flrho is ',e11.3,' to ',e11.3) 
      end subroutine findne
c***********************************************************************
      subroutine yindex(jhmin, jhmax, flne, jh, i3, eta)
      implicit none
      integer, intent(in) :: jhmin, jhmax, i3
      real, intent(in) :: flne
      integer, intent(out) :: jh(0:5)
      real, intent(out) :: eta
c local variables      
      integer :: j, k
      real :: y
c
      y = 4.*flne/real(i3)
      j = y
      j = max(j,jhmin+2)
      j = min(j,jhmax-3)
      do k = 0, 5
         jh(k) = j + k - 2
      enddo
      eta = 2.*(y-j)-1
c      
      return
      end subroutine yindex
c***********************************************************************
      subroutine findux(flr, xi, eta, ux)
      implicit none
      real, intent(in) :: flr(4, 4), xi, eta
      real, intent(out) :: ux
c local variables      
      integer :: i, j      
      real :: uxj(4), u(4)
c
      do j = 1, 4
         do i = 1, 4
            u(i) = flr(i, j)
         enddo
         uxj(j) = fintp(u, xi)
      enddo
      ux = fint(uxj, eta)
c
      return
      end subroutine findux   
c**********************************************************************
      subroutine rd(nel, nkz, izz, ilab, jh, n_tot, ff, rr, i3, umesh, fac)
      implicit none
      integer, intent(in) :: nel, nkz(ipe), izz(ipe), ilab(0:5), 
     >      jh(0:5), n_tot, i3
      real, intent(in) :: umesh(nptot)
      real(dp), intent(in) :: fac(nel)
      real, intent(out) :: ff(:,:,0:,0:)  ! (nptot, ipe, 6, 6)
      real, intent(out) :: rr(28, ipe, 0:5, 0:5)
c local variables      
      integer :: i, j, k, l, m, n, itt, jnn, izp, ne1, ne2, ne, ib, ia
      real :: fion(-1:28), yb, ya, d
c declare variables in common block (instead of by default: real (a-h, o-z), integer (i-n))   
!       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int,
!      : ne1p, ne2p, np, kp1, kp2, kp3, npp, mx, nx   
!       real :: umin, umax, epatom, oplnck, fionp, yy1, yy2, yx    
!       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
!      + nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1p(17,91,25),
!      + ne2p(17,91,25),fionp(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
!      + kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
!      + yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)  
!       save /atomdata/
c
c  i=temperature inex
c  j=density index
c  k=frequency index
c  n=element index
c  Get:
c    mono opacity cross-section ff(k,n,i,j)
c    modified cross-section for selected element, ta(k,i,j)
c
c     Initialisations    
      rr=0.
      ff=0.
c     
c  Start loop on i (temperature index)
      do i = 0, 5
         itt = (ilab(i) - ite1)/2 + 1        
         do j = 0, 5
            jnn = (jh(j)*i3 - jn1(itt))/2 + 1                  
c        Read mono opacities  
            do n = 1, nel
               izp = izz(n)
               ne1 = ne1p(nkz(n), itt, jnn)
               ne2 = ne2p(nkz(n), itt, jnn)
               do ne = ne1, ne2
                  fion(ne) = fionp(ne, nkz(n), itt, jnn)
               enddo    
               do ne = ne1, min(ne2, izp-2)
                  rr(izp-1-ne, n, i, j) = fion(ne)
               enddo
               
               do k = 1, n_tot
                  ff(k, n, i, j) = yy2(k+kp2(nkz(n), itt, jnn))
               enddo
                 
               if (fac(n) /= 1d0) then
                  do k = 1, size(ff,dim=1)
                     ff(k,n,i,j) = fac(n)*ff(k,n,i,j)
                  end do
               end if
            enddo !n
         enddo !j
      enddo !i
c
      return
c
      end subroutine rd
c***********************************************************************
      subroutine ross(flmu, dv, ntot,rs, rossl) 
      implicit none
      integer, intent(in) :: ntot
      real, intent(in) :: flmu, dv, rs(nptot, 0:5, 0:5)
      real, intent(out) :: rossl(0:5, 0:5)
      integer :: i, j, n
      real(dp) :: drs, dd, oross
      real :: fmu, tt
c
c  oross=cross-section in a.u.
c  rossl=log10(ROSS in cgs)
      do i = 0, 5
         do j = 0, 5
            drs = 0.d0
            do n = 1, ntot
               dd = 1.d0/rs(n, i, j)
               drs = drs + dd
            enddo
            oross = 1.d0/(drs*dv)
            rossl(i, j) = log10(oross) - 16.55280d0 - flmu !log10(fmu)                 
         enddo !j
      enddo !i
c
      return
      end subroutine ross
c***********************************************************************
      subroutine mix(ntot, nel, fa, ff, rs, rr, rion)
      implicit none
      integer, intent(in) :: ntot, nel
      real, intent(in) :: ff(nptot, ipe, 0:5, 0:5), fa(ipe), rr(28, 17, 0:5, 0:5)
      real, intent(out) :: rs(nptot, 0:5, 0:5), rion(28, 0:5, 0:5)
c local variables      
      integer :: i, j, k, n, m
      real :: rs_temp, rion_temp
c
      do i = 0, 5
         do j = 0, 5
            do n = 1, ntot
               !rs_temp = ff(n,1,i,j)*fa(1)
               !do k = 2, nel
               !   rs_temp = rs_temp + ff(n,k,i,j)*fa(k)
               !enddo
               !rs(n,i,j) = rs_temp 
               rs(n, i, j) = dot_product(ff(n,1:nel,i,j),fa(1:nel))
            enddo            
            do m = 1, 28
               !rion_temp = rr(m, 1, i, j)*fa(1)
               !do k = 2, nel
               !   rion_temp = rion_temp + rr(m,k,i,j)*fa(k)
               !enddo
               !rion(m,i,j) = rion_temp
               rion(m,i,j) = dot_product(rr(m,1:nel,i,j),fa(1:nel))
            enddo    
         enddo
      enddo
c
      return
      end subroutine mix
C***********************************************************************    
      subroutine interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)
      implicit none
      integer, intent(in) :: nel, i3
      real, intent(in) :: ux, uy, xi, eta
      real, intent(out) :: gx, gy, g
c local variables
      integer :: i, j, l
      real :: V(4), U(4),  vyi(4)
      real :: x3(3), fx3!, fxxy(0:5, 0:5), fyxy(0:5, 0:5)
c pointers and targets      
      real, target :: fx(0:5, 0:5), fy(0:5, 0:5), fxy(0:5, 0:5), fyx(0:5, 0:5), 
     : fxx(0:5, 0:5), fyy(0:5, 0:5), rossl(0:5, 0:5)
      real, pointer :: f3(:), fin(:, :), finx(:, :), finy(:, :)   
c        
c     interpolation of g (=rosseland mean opacity)
c Use refined techniques of bi-cubic spline interpolation (Seaton 1993):
!      call deriv(rossl, fx, fy, fxy)
!      call interp2(rossl, fx, fy, fxy, xi, eta, g, gx, gy)
!      gy = 0.5*gy/uy
!      gx = (80./real(i3))*(0.5*gx-gy*ux)
c Alternatively, use interpolation techniques by M.-A. Dupret to ensure smooothness required
c for pulsation studies:      
      do i = 0, 3
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2 
         do j = 1, 4
            f3 => rossl(i:i+2, j)
            call deriv3(f3, x3, fx3)
            fx(i+1, j) = fx3
         enddo
      enddo   
      do i = 0, 3
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2       
         do j = 1, 4        
            f3 => rossl(j, i:i+2)
            call deriv3(f3, x3, fx3)
            fy(j, i+1) = fx3
         enddo
      enddo   
      do i = 1, 2
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2 
         do j = 2, 3
            f3 => fx(i:i+2, j)
            call deriv3(f3, x3, fx3)
            fxx(i+1, j) = fx3
         enddo
      enddo 
      do i = 1, 2
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2 
         do j = 2, 3
            f3 => fx(j, i:i+2)
            call deriv3(f3, x3, fx3)
            fxy(j, i+1) = fx3
         enddo
      enddo      
       do i = 1, 2
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2        
         do j = 2, 3
            f3 => fy(i:i+2, j)
            call deriv3(f3, x3, fx3)
            fyx(i+1, j) = fx3
         enddo
      enddo 
      do i = 1, 2
         x3(1) = i
         x3(2) = i+1
         x3(3) = i+2 
         do j = 2, 3
            f3 => fy(j, i:i+2)
            call deriv3(f3, x3, fx3)
            fyy(j, i+1) = fx3
         enddo
      enddo   
!      call deriv(rossl, fx, fy, fxy) 
!      call deriv(fx, fxx, fxy, fxxy) 
!      call deriv(fy, fyx, fyy, fyxy)             
      fin => rossl(2:3, 2:3)
      finx => fx(2:3, 2:3)
      finy => fy(2:3, 2:3)
      call interp3(fin, finx, finy, xi, eta, g)            
      fin => fx(2:3, 2:3)
      finx => fxx(2:3, 2:3)
      finy => fxy(2:3, 2:3)      
      call interp3(fin, finx, finy, xi, eta, gx)      
      fin => fy(2:3, 2:3)
      finx => fyx(2:3, 2:3)
      finy => fyy(2:3, 2:3)      
      call interp3(fin, finx, finy, xi, eta, gy)             
      gy = 0.5*gy/uy
      gx = (80./real(i3))*(0.5*gx-gy*ux)
c
      RETURN
      end subroutine interp  
C**************************************
      function fint(u,r)
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1     1     3
c  then a cubic fit is:
      P(R)=( 
     +  27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(
     +  27*(u(3)-u(2))-(u(4)-u(1))   +R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
c
        fint=p(r)
c
      return
      end function fint
c***********************************************************************
      function fintp(u,r)
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1     1     3
c  then a cubic fit to the derivative is:
      PP(R)=( 
     +  27*(u(3)-u(2))-(u(4)-u(1))   +2.*R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +3.*R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.
c
        fintp=pp(r)
c
      return
      end function fintp
C
c***********************************************************************
      SUBROUTINE DERIV(f, fx, fy, fxy)
C
      real, intent(in) :: f(0:5, 0:5)
      real, intent(out) :: fx(0:5, 0:5), fy(0:5, 0:5), fxy(0:5, 0:5) 
      real :: C(6)
C
C  GET FX
      DO 70 J = 0, 5
         L=0
         DO 50 I = 0, 5
            L=L+1
            C(L)=F(I,J)
   50    CONTINUE
         CALL GET(C,L)
         L=0
         DO 60 I = 0, 5
            L = L + 1
            FX(I, J) = C(L)
   60    CONTINUE
   70 CONTINUE
C
C  GET FY
      DO 100 I = 0, 5
         L=0
         DO 80 J = 0, 5
            L = L + 1
            C(L) = F(I, J)
   80    CONTINUE
         CALL GET(C,L)
         L=0
         DO 90 J = 0, 5
            L = L + 1
            FY(I,J) = C(L)
   90    CONTINUE
  100 CONTINUE
C
C  GET FXY
      DO 130 I = 0, 5
         L = 0
         DO 110 J = 0, 5
            L = L + 1
            C(L) = FX(I, J)
  110    CONTINUE
         CALL GET(C,L)
         L=0
         DO 120 J = 0, 5
            L = L + 1
            FXY(I,J) = C(L)
  120    CONTINUE
  130 CONTINUE
c
      RETURN
C
      END SUBROUTINE DERIV
CC******************************************************************
C
      SUBROUTINE GET(F,N)
C
C  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS
C  OF UNITY.
C  RETURNS DERIVATIVES OF ORIGINAL F IN LOCATION F
C
c  REVISED 5.5.95
C
      PARAMETER (IPI=6)
      DIMENSION F(IPI),D(IPI),T(IPI)
C
      IF(N.LE.0)THEN
         WRITE(6,*)' Error in SUBROUTINE GET: N=',N
         STOP
      ELSEIF(N.EQ.1)THEN
         F(1)=0.
         RETURN
      ELSEIF(N.EQ.2)THEN
         F(1)=F(2)-F(1)
         F(2)=F(1)
         RETURN
      ELSEIF(N.EQ.3)THEN
         FP1=.5*(-3.*F(1)+4.*F(2)-F(3))
         FPN=.5*(F(1)-4.*F(2)+3.*F(3))
      ELSE
         FP1=(-11.*F(1)+18.*F(2)-9.*F(3)+2.*F(4))/6.
         FPN=(11.*F(N)-18.*F(N-1)+9.*F(N-2)-2.*F(N-3))/6.
      ENDIF
C
      D(1)=-.5
      T(1)=.5*(-F(1)+F(2)-FP1)
C
      DO 10 J=2,N-1
         D(J)=-1./(4.+D(J-1))
         T(J)=-D(J)*(F(J-1)-2.*F(J)+F(J+1)-T(J-1))
   10 CONTINUE
C
      D(N)=(FPN+F(N-1)-F(N)-T(N-1))/(2.+D(N-1))
C
      DO 20 J=N-1,1,-1
         D(J)=D(J)*D(J+1)+T(J)
   20 CONTINUE
C
      DO 30 J=2,N-1
         F(J)=-F(J)+F(J+1)-2.*D(J)-D(J+1)
   30 CONTINUE
      F(1)=FP1
      F(N)=FPN
C
      RETURN
      END SUBROUTINE GET
C
c***********************************************************************
      subroutine INTERP2(f, fx, fy, fxy, xi, eta, g, gx, gy)
      real, intent(in) :: eta, xi, f(0:5, 0:5), fx(0:5, 0:5), fy(0:5, 0:5), fxy(0:5, 0:5)
      real, intent(out) :: g, gx, gy
      integer :: i , j
      real :: x, y, B(16)
C
C  FUNCTION DEFINITIONS FOR CUBIC EXPANSION
C
      FF(S,T)=    B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))
     +   +S*(     B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(     B( 9)+T*(B(10)+T*(B(11)+T*B(12)))
     +   +S*(     B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))
C
      FFX(S,T)=   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(  2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))
     +   +S*(  3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))
C
      FFY(S,T)=   B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))
     +   +T*(  2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))
     +   +T*(  3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))
C
C
      Y = (eta + 5.)/2.
      X = (xi + 5.)/2.    
!      i = floor(x)
!      j = floor(y)
      I = X + 1.E-5
      IF(ABS(X-I).LE.1.E-5) X = I
      J = Y + 1.E-5
      IF(ABS(Y-J).LE.1.E-5) Y = J                
C
C     INTERPOLATE
C
C  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
      B(1)=F(I,J)
      B(2)=FY(I,J)
      B(3)=3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1)
      B(4)=2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)
C
      B(5)=FX(I,J)
      B(6)=FXY(I,J)
      B(7)=3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1)
      B(8)=2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)
C
      B(9)=3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J)
      B(10)=3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J)
      B(11)=9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     + +4*FXY(I,J)
     + +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     + +2*(FXY(I,J+1)+FXY(I+1,J))
     + +FXY(I+1,J+1)
      B(12)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FX(I,J)+FX(I,J+1))
     + +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     + +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     + -FXY(I+1,J)-FXY(I+1,J+1)
C
      B(13)=2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14)=2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FY(I,J)+FY(I+1,J))
     + +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     + +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     + -FXY(I+1,J+1)-FXY(I,J+1)
      B(16)=4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
     +    +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))
     + +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)
C
C  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
C      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)
      U = X - I
      V = Y - J
      G = FF(U, V)
      gx = FFX(U, V)
      gy = FFY(U, V)
!      DGDT=(1./CT)*FFX(U,V)-(3./CN)*FFY(U,V)
!      DGDRHO=(1./CN)*FFY(U,V)
      RETURN
C
      END SUBROUTINE INTERP2
c*************************************************************
c This subroutine estimates the partial derivative of a function 
c as described in the PhD thesis by M.-A. Dupret.
      subroutine deriv3(f, x, fx)
      real, intent(in) :: f(3), x(3)
      real, intent(out) :: fx
      real :: a, a1, a2, b, x1, x2
c            
      x1 = (f(3)- f(2))/(x(3)-x(2))
      x2 = (f(2)-f(1))/(x(2)-x(1))
      b = abs(x1/x2)
      a = (2.*x(2) - x(1) - x(3))/((x(2)-x(1))*(x(2)-x(3)))
      a1 = (x(2)-x(3))/((x(1)-x(3))*(x(1)-x(2)))
      a2 = (x(2)-x(1))/((x(3)-x(2))*(x(3)-x(1)))
      if (b .ge. 0.2 .and. b .le. 5.) then            
         fx = a*f(2) + a1*f(1) + a2*f(3)
      else
         if (abs(x1)<abs(x2)) fx = sign(1.0,x1)*min(abs(x1), abs(x2))
         if (abs(x2)<abs(x1)) fx = sign(1.0,x2)*min(abs(x1), abs(x2))         
      endif
c     
      end subroutine deriv3            
c***********************************************************************      
      subroutine interp3(f, fx, fy, xi, eta, g)
      real, intent(in) :: eta, xi, f(2:3, 2:3), fx(2:3, 2:3), fy(2:3, 2:3)
      real, intent(out) :: g
      real :: x, y
      
      x = (xi + 5.)/2.        
      y = (eta + 5.)/2.
      
      P11 = P1(x,2.,3.)*P1(y,2.,3.)
      P21 = P2(x,2.,3.)*P1(y,2.,3.)
      P12 = P1(x,2.,3.)*P2(y,2.,3.)
      P22 = P2(x,2.,3.)*P2(y,2.,3.)
      Px11 = Px1(x,2.,3.)*P1(y,2.,3.)
      Px21 = Px2(x,2.,3.)*P1(y,2.,3.)
      Px12 = Px1(x,2.,3.)*P2(y,2.,3.)
      Px22 = Px2(x,2.,3.)*P2(y,2.,3.)
      Py11 = P1(x,2.,3.)*Px1(y,2.,3.)
      Py21 = P2(x,2.,3.)*Px1(y,2.,3.)
      Py12 = P1(x,2.,3.)*Px2(y,2.,3.)
      Py22 = P2(x,2.,3.)*Px2(y,2.,3.)
      
      g = f(2,2)*P11 + f(3,2)*P21 + f(2,3)*P12 + f(3,3)*P22 
     : + fx(2,2)*Px11 + fx(3,2)*Px21 + fx(2,3)*Px12 + fx(3,3)*Px22
     : + fy(2,2)*Py11 + fy(3,2)*Py21 + fy(2,3)*Py12 + fy(3,3)*Py22                
            
      end subroutine interp3
c***********************************************************************      
      function P1(x, xi, xj)
c
      P1 = (x - xj)*(x - xj) * (2*x - 3*xi + xj)/((xj - xi)*(xj - xi)*(xj - xi))
c      
      end function P1      
c***********************************************************************
      function P2(x, xi, xj)
      real, intent(in) :: x, xi, xj
c
      P2 = -(x - xi)*(x - xi) * (2*x - 3*xj + xi)/((xj - xi)*(xj - xi)*(xj - xi))
c            
      end function P2
c***********************************************************************      
      function Px1(x, xi, xj)
c
      Px1 = (x - xi) * (x - xj)*(x - xj)/((xj - xi)*(xj - xi))      
c
      end function Px1
c***********************************************************************
      function Px2(x, xi, xj)
c
      Px2 = (x - xi)*(x - xi) * (x - xj)/((xj - xi)*(xj - xi))    
c      
      end function Px2      
c***********************************************************************
      subroutine scatt(ih, jh, rion, uf, f, umesh, semesh, dscat, ntot, epa, ierr)
      use op_load, only: BRCKR
      integer, intent(inout) :: ierr
      dimension rion(28, 0:5, 0:5),uf(0:100),f(nptot,0:5,0:5),umesh(nptot),
     +  semesh(nptot),fscat(0:100),p(nptot),rr(28),ih(0:5),jh(0:5)
        integer i,j,k,n
c HH: always use meshtype q='m'
      ite3=2
      umin=umesh(1)
      CSCAT=EPA*2.37567E-8
      ierr = 0
c
      do i = 0, 5
         ft = exp10(ITE3*dble(ih(i))/40d0)
         do j = 0, 5
            fne = exp10(ITE3*dble(jh(j))/4d0)
            do k = 1, ntot
               p(k) = f(k, i, j)
            enddo
            do m = 1, 28
              rr(m) = rion(m, i, j)
            enddo      
            CALL BRCKR(FT,FNE,RR,28,UF,100,FSCAT,ierr)
            if (ierr /= 0) return
            do n = 0, 100
              u = uf(n)
              fscat(n) = cscat*(fscat(n)-1)
            enddo
            do n = 2, ntot-1
              u = umesh(n)
              se = semesh(n)
              m=(u-umin)/dscat
              ua=umin+dscat*m
              ub=ua+dscat
              p(n)=p(n)+((ub-u)*fscat(m)+(u-ua)*fscat(m+1))/(dscat*se)
            enddo
            u=umesh(ntot)
            se=semesh(ntot)
            p(ntot)=p(ntot)+fscat(100)/se
            p(1)=p(1)+fscat(1)/(1.-.5*umin)
            do k=1,ntot
              f(k,i,j)=p(k)
            enddo  
         enddo
      enddo  
C
      return
      end subroutine scatt
c***********************************************************************      
      subroutine screen1(ih,jh,rion,umesh,ntot,epa,f)
      use op_load, only: screen2
      dimension uf(0:100), umesh(nptot),
     +  fscat(0:100), ih(0:5), jh(0:5)    
      real, target :: f(nptot, 0:5, 0:5), rion(28, 0:5, 0:5)   
      integer i, j, k, m
      real, pointer :: p(:), rr(:)
c
      ite3=2
      umin=umesh(1)
      umax=umesh(ntot)
c
      do i = 0, 5
         ft=exp10(ITE3*dble(ih(i))/40d0)
         do j = 0, 5
            fne=exp10(ITE3*dble(jh(j))/4d0)
!            do k=1,ntot
              p => f(1:ntot,i,j)
!            enddo
!            do m=1,28
              rr => rion(1:28,i,j)
!            enddo      
            call screen2(ft,fne,rr,epa,ntot,umin,umax,umesh,p)
!            do k=1,ntot
!              f(k,i,j)=p(k)
!            enddo            
        enddo
      enddo  
C
      return
      end subroutine screen1
c***********************************************************************
      end module op_osc
