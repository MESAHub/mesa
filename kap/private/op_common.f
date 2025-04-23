      module op_common

      use math_lib
      use op_def

      contains
! **********************************************************************
      subroutine xindex(flt, ilab, xi, ih, i3, ierr)
      implicit none
      integer, intent(in) :: i3
      real, intent(in) :: flt
      integer, intent(out) :: ih(4), ilab(4)
      real, intent(out) :: xi
      integer, intent(out) :: ierr
      integer :: i, ih2
      real :: x

      ierr = 0
      if(flt<3.5) then
        ierr = 102
        return
      else if(flt>8.) then
        ierr = 102
        return
      end if

      x = 40.*flt/real(i3)
      ih2 = x
      ih2 = max(ih2, 140/i3+2)
      ih2 = min(ih2, 320/i3-3)
      do i = 1, 4
         ih(i) = ih2 + i - 2
         ilab(i) = i3*ih(i)
      end do
      xi = 2.*(x-ih2) - 1

      end subroutine xindex

! *********************************************************************
      subroutine jrange(ih, jhmin, jhmax, i3)
      implicit none
      integer, intent(in) :: ih(4), i3
      integer, intent(out) :: jhmin, jhmax
      integer :: i

      jhmin = 0
      jhmax = 1000
      do i = 1, 4
         jhmin = max(jhmin, js(ih(i)*i3)/i3)
         jhmax = min(jhmax, je(ih(i)*i3)/i3)
      end do

      end subroutine jrange
! *********************************************************************

      subroutine findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt, xi, flne, flmu, flr, epa, uy, i3, ierr)
      use op_load, only : solve
      implicit none
      integer, intent(in) :: ilab(4), nel, nkz(ipe), jhmin, ih(4), i3
      integer, intent(inout) :: jhmax
      integer, intent(out) :: ierr
      real, intent(in) :: fa(ipe), flt, xi, flmu
      real,intent(out) :: flne, uy, epa
      real, intent(inout) :: flrho
      integer :: i, j, n, jh, jm, itt, jne, jnn
      real :: flrmin, flrmax, flr(4,4), uyi(4), efa(4, 7:118), flrh(4, 7:118), u(4), flnei(4), y, zeta, efa_temp
! declare variables in common block, by default: real (a-h, o-z), integer (i-n)
!       integer :: ite1, ite2, ite3, jn1, jn2, jne3, ntot, nc, nf, int, ne1, ne2, np, kp1, kp2, kp3, npp, mx, nx
!       real :: umin, umax, epatom, oplnck, fion, yy1, yy2, yx
!       common /atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot, &
!        nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25), &
!        ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25), &
!        kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000), &
!        yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
!       save /atomdata/

!  efa(i,jh)=sum_n epa(i,jh,n)*fa(n)
!  flrh(i,jh)=log10(rho(i,jh))

!  Get efa
      do i = 1, 4
         itt = (ilab(i)-ite1)/2 + 1
         do jne = jn1(itt), jn2(itt), i3
            jnn = (jne-jn1(itt))/2 + 1
            jh = jne/i3
            efa_temp = 0.
            do n = 1, nel
               efa_temp = efa_temp + epatom(nkz(n), itt, jnn)*fa(n)
            end do !n
            efa(i, jh) = efa_temp
         end do !jne
      end do !i

! Get range for efa > 0
      do i = 1, 4
         do jh = jhmin, jhmax
            if(efa(i, jh) <= 0.)then
               jm = jh - 1
               jhmax = min(jhmax, jm)
               cycle
            end if
         end do
      end do

! Get flrh
      do jh=jhmin,jhmax
         do i=1,4
            flrh(i, jh) = flmu + 0.25*i3*jh - log10(dble(efa(i, jh)))
         end do
      end do

!  Find flrmin and flrmax
      flrmin = -1000
      flrmax = 1000
      do i = 1, 4
         flrmin = max(flrmin, flrh(i,jhmin))
         flrmax = min(flrmax, flrh(i,jhmax))
      end do

!  Check range of flrho
      if(flrho < flrmin .or. flrho > flrmax)then
         write(*,*) "findne failed because density is out of range for logT, logRho", flt, flrho
         write(*,*) "Allowed range for logRho is",flrmin," to ",flrmax
         ierr = 101
         return
      end if

!  Interpolations in j for flne
      do jh = jhmin, jhmax
         if (flrh(2,jh) > flrho) then
            jm = jh - 1
            cycle
         end if
         write(*,*) ' Interpolations in j for flne'
         write(*,*) ' Not found, i=',i
         stop
      end do
      jm=max(jm,jhmin+1)
      jm=min(jm,jhmax-2)

      do i = 1, 4
         do j = 1, 4
            u(j) = flrh(i, jm+j-2)
            flr(i,j) = flrh(i, jm+j-2)
         end do
         call solve(u, flrho, zeta, uyi(i), ierr)
         if (ierr /= 0) return
         y = jm + 0.5*(zeta+1)
         flnei(i) = .25*i3*y
      end do

!  Interpolations in i
      flne = fint(flnei, xi)
      uy = fint(uyi, xi)
! Get epa
      epa = exp10(dble(flne + flmu - flrho))

      return

!  601 format(' For flt=',1p,e11.3,', flrho=',e11.3,' is out of range. Allowed range for flrho is ',e11.3,' to ',e11.3)
      end subroutine findne

! **********************************************************************
      subroutine yindex(jhmin, jhmax, flne, jh, i3, eta)
      implicit none
      integer, intent(in) :: jhmin, jhmax, i3
      real, intent(in) :: flne
      integer, intent(out) :: jh(4)
      real, intent(out) :: eta
      integer :: j, k
      real :: y

      y = 4.*flne/real(i3)
      j = y
      j = max(j,jhmin+2)
      j = min(j,jhmax-3)
      do k = 1, 4
         jh(k) = j + k - 2
      end do
      eta = 2.*(y-j)-1

      end subroutine yindex

! **********************************************************************
      subroutine findux(flr, xi, eta, ux)
      implicit none
      real, intent(in) :: flr(4, 4), xi, eta
      real, intent(out) :: ux
      integer :: i, j
      real :: uxj(4), u(4)

      do j = 1, 4
         do i = 1, 4
            u(i) = flr(i, j)
         end do
         uxj(j) = fintp(u, xi)
      end do
      ux = fint(uxj, eta)

      end subroutine findux

! *************************************
      function fint(u,r)
      dimension u(4)

!  If  P(R) =   u(1)  u(2)  u(3)  u(4)
!  for   R  =    -3    -1     1     3
!  then a cubic fit is:
      P(R)=(27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(27*(u(3)-u(2))-(u(4)-u(1))
     &      + R*(-3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(-3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.

        fint=p(r)

      end function fint
! **********************************************************************
      function fintp(u,r)
      dimension u(4)

!  If  P(R) =   u(1)  u(2)  u(3)  u(4)
!  for   R  =    -3    -1     1     3
!  then a cubic fit to the derivative is:
      PP(R)=(27*(u(3)-u(2))-(u(4)-u(1))  +2.*R*(-3*(u(2)+u(3))+3*(u(4)+u(1)) +3.*R*(-3*(u(3)-u(2))+(u(4)-u(1)) )))/48.

        fintp=pp(r)

      end function fintp

! **********************************************************************
      subroutine scatt(ih, jh, rion, uf, f, umesh, semesh, dscat, ntot, epa, ierr)
      use op_load, only: BRCKR
      integer, intent(inout) :: ierr
      real :: umesh(:), semesh(:) ! (nptot)
      real :: f(:,:,:) ! (nptot,4,4)
      dimension rion(28, 4, 4),uf(0:100),fscat(0:100),p(nptot),rr(28),ih(4),jh(4)
        integer :: i,j,k,n
! HH: always use meshtype q='m'
      ite3=2
      umin=umesh(1)
      CSCAT=EPA*2.37567E-8
      ierr = 0

      do i = 1, 4
         ft = real(exp10(ITE3*dble(ih(i))/40d0))
         do j = 1, 4
            fne = real(exp10(ITE3*dble(jh(j))/4d0))
            do k = 1, ntot
               p(k) = f(k, i, j)
            end do
            do m = 1, 28
              rr(m) = rion(m, i, j)
            end do
            call BRCKR(FT,FNE,RR,28,UF,100,FSCAT,ierr)
            if (ierr /= 0) return
            do n = 0, 100
              u = uf(n)
              fscat(n) = cscat*(fscat(n)-1)
            end do
            do n = 2, ntot-1
              u=umesh(n)
              se=semesh(n)
              m=(u-umin)/dscat
              ua=umin+dscat*m
              ub=ua+dscat
              p(n)=p(n)+((ub-u)*fscat(m)+(u-ua)*fscat(m+1))/(dscat*se)
            end do
            u=umesh(ntot)
            se=semesh(ntot)
            p(ntot)=p(ntot)+fscat(100)/se
            p(1)=p(1)+fscat(1)/(1.-.5*umin)
            do k=1,ntot
              f(k,i,j)=p(k)
            end do
         end do
      end do

      end subroutine scatt
! **********************************************************************
      subroutine screen1(ih,jh,rion,umesh,ntot,epa,f)
      use op_load, only: screen2
      real :: umesh(:) ! (nptot)
      real, pointer :: f(:,:,:) ! (nptot,4,4)
      dimension uf(0:100),fscat(0:100), ih(4), jh(4)
      real, target :: rion(28, 4, 4)
      integer :: i, j, k, m
      real, pointer :: p(:), rr(:)

      ite3=2
      umin=umesh(1)
      umax=umesh(ntot)

      do i = 1, 4
         ft = real(exp10(ITE3*dble(ih(i))/40d0))
         do j = 1, 4
            fne = real(exp10(ITE3*dble(jh(j))/4d0))
!            do k=1,ntot
              p => f(1:ntot,i,j)
!            end do
!            do m=1,28
              rr => rion(1:28,i,j)
!            end do
            call screen2(ft,fne,rr,epa,ntot,umin,umax,umesh,p)
!            do k=1,ntot
!              f(k,i,j)=p(k)
!            end do
        end do
      end do

      end subroutine screen1

      end module op_common
