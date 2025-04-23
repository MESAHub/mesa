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
!
! Haili Hu 2010

      module op_load
      use math_lib
      use op_def
      logical :: have_loaded_op = .false.

      contains
! *****************************************************************
      subroutine op_dload(path, cache_filename, ierr)
      implicit none
      character (len=*), intent(in) :: path, cache_filename
      integer, intent(out) :: ierr

      real :: am,amm,delp,dpack
      integer :: ios,it,ite11,ite22,ite33,itt,itte1,itte2,itte3,izz,jne,ite
      integer :: jne1,jne22,jne33,jnn,k,n,ncount2,ncount3,ja,jn
      integer :: jne11,jne2,nccc,ne,nfff,ntott,nn
      real :: orss,um,ux,umaxx,uminn,u
      real,dimension(nptot):: umesh, semesh
      integer :: ntotv
      real :: dv,dv1

      integer :: cache_version

       common /mesh/ ntotv,dv,dv1,umesh,semesh
!    common /atomdata/
!       common/atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot, &
!         nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25), &
!         ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25), &
!         kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000), &
!         yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)

      integer,dimension(ipe) :: ifl,iflp
      character :: num(0:9)*1,zlab(ipe)*3,tlab*6,zlabp(ipe)*3
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/

      integer :: kz(17)
      data kz/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28/
      save  /mesh/    !HH: put common block in static memory

      integer :: nx_temp(nptot)
      real :: y_temp(nptot)
      integer :: nx_index, left_n, right_n, n_index
      real :: left_val, right_val, cross_section, slope

      if (allocated(yy2) .eqv. .false.) then
         ! yy2 actually needs 29,563 x 10,000 length
         ALLOCATE(yy2(30000*10000),nx(19305000),yx(19305000),stat=ierr)
         if (ierr/=0) return
         yy2=0.0
         nx=0.0
         yx=0.0
   !     write(*,*) "ierr",ierr
      end if


      ierr=0
      if (have_loaded_op) return

!$omp critical (critial_do_op_dload)

      if (have_loaded_op) GOTO 1001

      !path = '../OP4STARS_1.3'
      !call getenv("oppath", path)
      !if (len(trim(path)) == 0) then
      !   write(6,*) 'Define environmental variable oppath (directory of OP data)'
      !   stop
      !end if

      ios = 0
      open(1,file=trim(cache_filename),action='read',status='old',iostat=ios,form='unformatted')
      if (ios == 0) then
         write(*,*) 'reading OP cache file ' // trim(cache_filename)
         read(1,iostat=ios) cache_version,ntotv,dv,dv1,umesh,ite1,ite2,ite3,jn1,jn2,jne3,
     &                      umin,umax,ntotp,nc,nf,int,epatom,oplnck,ne1p,ne2p,fionp,np,kp1,kp2,kp3,npp,yy2,nx,yx
         write(*,*) 'done reading OP cache file'
         close(1)
         if (cache_version /= op_cache_version) then
            write(*,*) 'wrong version of OP cache'
            write(*,*) 'cache file path is set by op_mono_data_cache_filename'
            write(*,*) 'perhaps cache is shared between different MESA versions'
            write(*,*) 'please delete cache file and try again'
            stop 'op_load'
         end if
         if (ios == 0) then
            have_loaded_op = .true.
            call IMESH(UMESH,NTOTV)
            GOTO 1001
         end if
         write(*,*) 'failed in reading cache file'
         write(*,*) 'please delete cache file and try again'
         stop 'op_load'
      end if

      ite1=140
      ite2=320
      ite3=2
      do n=1,ipe
        ifl(n)=50+n
        zlab(n)='m'//num(kz(n)/10)//num(kz(n)-10*(kz(n)/10))
        iflp(n)=70+n
        zlabp(n)='a'//num(kz(n)/10)//num(kz(n)-10*(kz(n)/10))
      end do

      write(*,*) 'loading OP mono data...'

!  READ INDEX FILES
!     FIRST FILE
      NN=1
!     write(*,*) ' Opening '//'./'//zlab(1)//'.index'
      OPEN(1,FILE=trim(path)//'/'//ZLAB(1)//'.index',STATUS='OLD',iostat=ios)
      if (ios /= 0) then
         write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(1) // '.index'
         ierr = -1
         GOTO 1001
      end if
      read(1,*)IZZ,AMM
      read(1,*)ITTE1,ITTE2,ITTE3
      read(1,*)UMIN,UMAX
      read(1,*)NC,NF
      read(1,*)DPACK
      CLOSE(1)
      if (IZZ /= KZ(1)) then
         write(6,6001)zlab(1),izz,nn,kz(1)
         ierr=1
         GOTO 1001
      end if
      NTOTP=NF
      if (NTOTP > nptot) then
         write(6,6002)ntotp,nptot
         ierr=2
         GOTO 1001
      end if
      INT(1)=1
      if (ITTE3 /= ITE3) then
         write(6,6077)ite3,itte3,nn
         ierr=3
         GOTO 1001
      end if
!      ITE1=MAX(ITE1,ITTE1)
!      ITE2=MIN(ITE2,ITTE2)
!
!  READ MESH FILES
      OPEN(1,FILE=trim(path)//'/'//ZLAB(1)//'.mesh',status='old',form='unformatted',iostat=ios)
      if (ios /= 0) then
         write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(1) // '.mesh'
         ierr = -1
         GOTO 1001
      end if
      read(1)DV,NTOTV,(UMESH(N),N=1,NTOTV)
      umin=umesh(1)
      umax=umesh(ntotv)
      DV1=DV
      CLOSE(1)
!
!  GET MESH FOR SCREEN
      call IMESH(UMESH,NTOTV)
!
!     SUBSEQUENT FILES
      do N=2,ipe
         NN=N
         OPEN(1,FILE=trim(path)//'/'//ZLAB(N)//'.index',STATUS='OLD')
         read(1,*)IZZ,AMM
         read(1,*)ITE11,ITE22,ITE33
         read(1,*)UMINN,UMAXX
         read(1,*)NC,NF
         read(1,*)DPACK
         CLOSE(1)
         if (ITE33 /= ITE3) then
            write(6,6077)ite3,ite33,nn
            ierr=4
            GOTO 1001
         end if
!         ITE1=MAX(ITE1,ITE11)
!         ITE2=MIN(ITE2,ITE22)
         if (IZZ /= KZ(N)) then
            write(6,6001)zlab(n),izz,nn,kz(nn)
            ierr=5
            GOTO 1001
         end if
         NTOTT=NF
         if (NTOTT > NTOTP) then
            write(6,6006)nn,ntott,ntotp
            ierr=6
            GOTO 1001
         end if
! !!        if (UMIN /= UMINN.OR.UMAX /= UMAXX) GOTO 1003  !!
         INT(N)=NTOTP/NTOTT
         if (INT(N)*NTOTT /= NTOTP) then
            WRITE(6,6009)NN,NTOTT,NTOTP
            ierr=7
            GOTO 1001
         end if
         if (INT(N) /= 1)WRITE(6,6007)N,INT(N)
!
!        READ MESH FILES
!
         OPEN(1,FILE=trim(path)//'/'//ZLAB(N)//'.mesh', status='old',form='unformatted',iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(N) // '.mesh'
            ierr = -1
            GOTO 1001
         end if
         read(1)DV
       if (DV /= DV1) then
!          write(*,*) ' OP: N=',N,', DV=',DV,' NOT EQUAL TO DV1=',DV1
            ierr=8
          GOTO 1001
       end if
         CLOSE(1)
      end do
!
!  START TEMPERATURE LOOP
!
      ncount2=0
      ncount3=0
      do it=ite1,ite2,ite3
!
!        OPEN FILES
!
         TLAB='.'//NUM(IT/100)//NUM(IT/10-10*(IT/100))//NUM(IT-10*(IT/10))
         do n=1,ipe
!            if (SKIP(N))GOTO 70
            NN=N
            OPEN(IFL(N),FILE=trim(path)//'/'//ZLAB(N)//TLAB,FORM='UNFORMATTED',STATUS='OLD',iostat=ios)
            if (ios /= 0) then
               write(*,*) 'failed to open ' // trim(path)//'/'//ZLAB(N)//TLAB
               ierr = -1
               GOTO 1001
            end if
            if (n > 2) then
            OPEN(IFLP(N),FILE=trim(path)//'/'//ZLABP(N)//TLAB,FORM='UNFORMATTED',STATUS='OLD',iostat=ios)
            if (ios /= 0) then
               write(*,*) 'failed to open ' // trim(path)//'/'//ZLABP(N)//TLAB
               ierr = -1
               GOTO 1001
            end if
            end if
         end do
!        READ HEADINGS
         NN=1
         read(IFL(1))IZZ,ITE,AM,UM,UX,NCCC,NFFF,DelP,JNE1,JNE2,JNE3
         do n=2,ipe
!            if (SKIP(N))GOTO 80
            NN=N
            read(IFL(N))IZZ,ITE,AM,UM,UX,NC,NF,DelP,JNE11,JNE22,JNE33
            if (n > 2) read(iflp(n))
            if (JNE33 /= JNE3) then
              write(6,6099)jne3,jne33,nn
              ierr=9
              GOTO 1001
            end if
            JNE1=MAX(JNE1,JNE11)
            JNE2=MIN(JNE2,JNE22)
         end do
         itt=(it-ite1)/2+1
         jn1(itt)=jne1
         jn2(itt)=jne2
!
!         WRITE(98,9802)ITE,JNE1,JNE2,JNE3
!
!        START DENSITY LOOP
!
         do n=1,ipe
           do jn=jne1,jne2,jne3
             jnn=(jn-jne1)/2+1
!
!           START LOOP ON ELEMENTS
!
   95        read(IFL(N))JNE,EPATOM(n,itt,jnn),OPLNCK(n,itt,jnn),ORSS,
     &         NE1P(n,itt,jnn),NE2P(n,itt,jnn),
     &         (FIONP(NE,n,itt,jnn),NE=NE1P(n,itt,jnn),NE2P(n,itt,jnn))
             read(ifl(n))np(n,itt,jnn)
             if (np(n,itt,jnn) > 0) then
                read(ifl(n))(nx_temp(k),y_temp(k),k=1,np(n,itt,jnn))
                do nx_index = 2, np(n,itt,jnn)
                   left_val = y_temp(nx_index-1)
                   right_val = y_temp(nx_index)
                   left_n = nx_temp(nx_index-1)
                   right_n = nx_temp(nx_index)
                   slope = (right_val - left_val)/float(right_n - left_n)

                   do n_index = left_n, right_n
                      cross_section = left_val + (n_index-left_n)*slope
                      yy2(ncount2 + n_index) = cross_section
                   end do
                   yy2(ncount2 + left_n) = left_val
                   yy2(ncount2 + right_n) = right_val
                end do
                kp2(n, itt, jnn) = ncount2
                ncount2 = ncount2 + ntotp
             else
               read(ifl(n))(yy2(k+ncount2),k=1,ntotp)
                 kp2(n,itt,jnn)=ncount2
                 ncount2=ncount2+ntotp
             end if
               if (n > 2) then
                 read(iflp(n))ja,npp(n,itt,jnn)
                 if (npp(n,itt,jnn) > 0) then
           read(iflp(n))(nx(k+ncount3),yx(k+ncount3),k=1,npp(n,itt,jnn))
                   kp3(n,itt,jnn)=ncount3
                   ncount3=ncount3+npp(n,itt,jnn)
                 end if
               end if
            end do
          end do

!     write(6,610)it
!     write(6,*)'ncount1 = ',ncount1
!     write(6,*)'ncount2 = ',ncount2
!     write(6,*)'ncount3 = ',ncount3
!
!        CLOSE FILES

         do N=1,ipe
            close(IFL(N))
            close(iflp(n))
         end do

      end do

      write(*,*) 'done loading OP mono data'
      have_loaded_op = .true.

      !write(*,*)'ncount1 = ',ncount1
      !write(6,*)'ncount2 = ',ncount2
      !write(6,*)'ncount3 = ',ncount3
      ios = 0
      open(1, file=trim(cache_filename), iostat=ios, action='write', form='unformatted')
      if (ios == 0) then
         write(*,*) 'write ' // trim(cache_filename)
         write(1) op_cache_version, ntotv,dv,dv1,umesh,
     &     ite1,ite2,ite3,jn1,jn2,jne3,umin,umax,ntotp,nc,nf,int,epatom,oplnck, ne1p,
     &     ne2p,fionp,np,kp1,kp2,kp3,npp,yy2,nx,yx
         close(1)
      end if

1001  continue

!     pre-calculate semesh
      do n = 1, nptot
         u = umesh(n)
         semesh(n) = 1.d0 - exp(dble(-u))
      end do

!$omp end critical (critial_do_op_dload)


      return
610   format(10x,'Done IT= ',i3)
1004  WRITE(6,6004)ZLAB(NN),TLAB
      STOP
6001  FORMAT(//5X,'*** OP: FILE ',A3,' GIVES IZZ=',I3,'NOT EQUAL TO IZ(',I2,')=',I2,' ***')
6002  FORMAT(//5X,'*** OP: NTOT=',I7,' GREATER THAN nptot=',I7,' ***')
6003  FORMAT(//5X,'*** OP: DISCREPANCY BETWEEN DATA ON FILES ',A3,' AND ',A3,' ***')
6004  FORMAT(//5X,'*** OP: ERROR OPENING FILE ',A3,A6,'  ***')
6006  FORMAT(//5X,'OP: N=',I2,', NTOTT=',I7,', GREATER THAN NTOT=',I7)
6007  FORMAT(/5X,'OP: N=',I2,', INT(N)=',I4)
6009  FORMAT(' OP: N=',I5,', NTOTT=',I10,', NTOT=',I10/'   NTOT NOT MULTIPLE OF NTOTT')
! 6012  FORMAT(/10X,'ERROR, SEE WRITE(6,6012)'/10X,'IT=',I3,', JN=',I3,', N=',I3,', JNE=',I3/)
6077  FORMAT(//5X,'OP: DISCREPANCY IN ITE3'/10X,I5,' READ FROM UNIT 5'/10X,I5,' FROM INDEX FILE ELEMENT',I5)
6099  FORMAT(//5X,'OP: DISCREPANCY IN JNE3'/10X,I5,' READ FOR N=1'/10X,I5,' READ FOR N=',I5)

! 8000  FORMAT(5X,I5,F10.4/5X,3I5/2E10.2/2I10/10X,E10.2)
1010  write(*,*) ' OP: ERROR OPENING FILE '//'./'//ZLAB(1)//'.index'
      stop
1011  write(*,*) ' OP: ERROR OPENING FILE '//'./'//ZLAB(1)//'.mesh'
      stop
      end subroutine op_dload

! **********************************************************************
        subroutine IMESH(UMESH,NTOT)

      DIMENSION UMESH(nptot)
      COMMON/CIMESH/U(100),AA(nptot),BB(nptot),IN(nptot),ITOT,NN
      save /cimesh/

      UMIN=UMESH(1)
      UMAX=UMESH(NTOT)

      II=100
      A=(II*UMIN-UMAX)/REAL(II-1)
      B=(UMAX-UMIN)/REAL(II-1)
      do I=1,II
        U(I)=A+B*I
      end do

      ib=2
      ub=u(ib)
      ua=u(ib-1)
      d=ub-ua
      ibb=0
      do n=2,ntot
        if (umesh(n) > ub) then
          ua=ub
          ib=ib+1
          ub=u(ib)
          d=ub-ua
          if (umesh(n) > ub) then
            nn=n-1
            ibb=ib-1
            GOTO 1
          end if
        end if
        in(n)=ib
        aa(n)=(ub-umesh(n))/d
        bb(n)=(umesh(n)-ua)/d
      end do

    1      ib=ibb
      do n=nn+1,ntot
        ib=ib+1
        in(n)=ib
        u(ib)=umesh(n)
      end do
      itot=ib

      end subroutine IMESH


      subroutine msh(dv, ntot, umesh, semesh, uf, dscat)
      implicit none
      integer, intent(out) :: ntot
      real, intent(out) :: dv, uf(0:100), dscat
      real, intent(out) :: umesh(:), semesh(:) ! (nptot)
      integer :: i, ntotv
      real :: dvp, dv1, umin, umax, umeshp(nptot), semeshp(nptot)
      common /mesh/ ntotv, dvp, dv1, umeshp, semeshp
      save /mesh/

      ntot = ntotv
      dv = dvp
      do i=1,ntot
         umesh(i) = umeshp(i)
      end do
      do i=1,ntot
         semesh(i) = semeshp(i)
      end do

      umin = umesh(1)
      umax = umesh(ntot)
      dscat = (umax - umin)*0.01
      do i = 0, 100
         uf(i) = umin + i*dscat
      end do

      end subroutine msh


      subroutine solve(u,v,z,uz,ierr)
      integer, intent(inout) :: ierr
      dimension u(4)

!  If  P(R) =   u(1)  u(2)  u(3)  u(4)
!  for   R  =    -3    -1    1     3
!  then a cubic fit is:
      P(R)=(
     &  27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(
     &  27*(u(3)-u(2))-(u(4)-u(1))   +R*(
     &  -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(
     &  -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
!  First derivative is:
      PP(R)=(
     &  27*(u(3)-u(2))-(u(4)-u(1))+ 2*R*(
     &  -3*(u(2)+u(3))+3*(u(4)+u(1)) +3*R*(
     &  -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.

!      ierr = 0
!  Find value of z giving P(z)=v
!  First estimate
      z=(2.*v-u(3)-u(2))/(u(3)-u(2))
!  Newton-Raphson iterations
      do k=1,10
         uz=pp(z)
         d=(v-p(z))/uz
         z=z+d
         if (abs(d) < 1.e-4) return
      end do

!      write(*,*) ' Not converged after 10 iterations in SOLVE'
!      write(*,*) ' v=',v
!      do N=1,4
!         write(*,*) ' N, U(N)=',N,U(N)
!      end do
      ierr = 10
      return
!      stop

      end subroutine solve
! **********************************************************************

      subroutine BRCKR(T,FNE,RION,NION,U,NFREQ,SF, ierr)
      integer, intent(inout) :: ierr
!
!  CODE FOR COLLECTIVE EFFECTS ON THOMSON SCATTERING.
!  METHOD OF D.B. BOERCKER, AP. J., 316, L98, 1987.
!
!  INPUT:-
!     T=TEMPERATURTE IN K
!     FNE=ELECTRON DENSITY IN CM**(-3)
!     ARRAY RION (DIMENSIONED FOR 30 IONS).
!        RION(IZ) IS NUMBER OF IONS WITH NET CHARGE IZ.
!        NORMALISATION OF RION IS OF NO CONSEQUENCE.
!     NION=NUMBER OF IONS INCLUDED.
!     ARRAY U (DIMENSIONED FOR 1000). VALUES OF (H*NU/K*T).
!     NFREQ=NUMBER OF FREQUENCY POINTS.
!
!  OUTPUT:-
!     ARRAY SF, GIVING FACTORS BY WHICH THOMSON CROSS SECTION
!     SHOULD BE MULTIPLIED TO ALLOW FOR COLLECTIVE EFFECTS.
!
!  MODIFFICATIONS:-
!     (1) REPLACE (1.-Y) BY EXP(-Y) TO AVOID NEGATIVE FACTORS FOR
!         HIGHLY-DEGENERATE CASES.
!     (2) INCLUDE RELATIVISTIC CORRECTION.
!
      parameter (IPZ=28,IPNC=100)
      DIMENSION RION(IPZ),U(0:IPNC),SF(0:IPNC)

      AUNE=1.48185E-25*FNE
      AUT=3.16668E-6*T
      C1=-1.0650E-4*AUT
      C2=+1.4746E-8*AUT**2
      C3=-2.0084E-12*AUT*AUT*AUT
      V=7.8748*AUNE/(AUT*sqrt(AUT))
      call FDETA(V,ETA, ierr)  ! 23.10.93
      W=exp(dble(ETA))         ! 23.10.93
   11 R=FMH(W)/V
      A=0.
      B=0.
      do I=1,NION
         A=A+I*RION(I)
         B=B+I**2*RION(I)
      end do
      X=R+B/A

      Y=.353553*W
      C=1.1799E5*X*AUNE/(AUT*AUT*AUT)
      do N=0,NFREQ
         D=C/U(N)**2
         if (D > 5.) then
            D=-2./D
            F=2.666667*(1.+D*(.7+D*(.55+.341*D)))
         else
            G=2.*D*(1+D)
            F=D*((G+D*D*D)*log(dble(D/(2.+D)))+G+2.6666667)
         end if
         DELTA=.375*R*F/X
         SF(N)=(1.-R*DELTA-Y*FUNS(W))*
     &   (1.+U(N)*(C1+U(N)*(C2+U(N)*C3)))   !SAMPSON CORRECTION
      end do

      return

  600 FORMAT(5X,'NOT CONVERGED IN LOOP 10 OF BRCKR'/5X,'T=',1P,E10.2,', FNE=',E10.2)

      end subroutine BRCKR
! **********************************************************************
      function FUNS(A)
!
      if (A <= 0.001) then
         FUNS=1.
      else if (A <= 0.01) then
         FUNS=(1.+A*(-1.0886+A*(1.06066+A*1.101193)))/
     &     (1.+A*(0.35355+A*(0.19245+A+0.125)))
      else
         FUNS=(  1./(1.+0.81230*A)**2+
     &        0.92007/(1.+0.31754*A)**2+
     &        0.05683/(1.+0.04307*A)**2 )/
     &     (  1./(1.+0.65983*A)+
     &        0.92007/(1.+0.10083*A)+
     &        0.05683/(1.+0.00186*A)    )
      end if

      end function FUNS
! **********************************************************************
      function FMH(W)
!
!  CALCULATES FD INTEGRAL I_(-1/2)(ETA). INCLUDES FACTOR 1/GAMMA(1/2).
!  ETA=log(W)
!
      if (W <= 2.718282) then
         FMH=W*(1+W*(-.7070545+W*(-.3394862-W*6.923481E-4))/(1.+W*(1.2958546+W*.35469431)))
      else if (W <= 54.59815) then
         X=log(dble(W))
         FMH=(.6652309+X*(.7528360+X*.6494319))/(1.+X*(.8975007+X*.1153824))
      else
         X=log(dble(W))
         Y=1./X**2
         FMH=sqrt(X)*(1.1283792+(Y*(-.4597911+Y*(2.286168-Y*183.6074)))/(1.+Y*(-10.867628+Y*384.61501)))
      end if

      end function FMH
! **********************************************************************
      subroutine FDETA(X,ETA, ierr)
!
!  GIVEN X=N_e/P_e, CALCULATES FERMI-DIRAC ETA
!  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
!  MATHS. OF COMP., 21, 30, 1967.
!
      integer, intent(inout) :: ierr
      DIMENSION D(2:12)
      DATA D/
     &  3.5355339E-01, 5.7549910E-02, 5.7639604E-03, 4.0194942E-04,
     &  2.0981899E-05, 8.6021311E-07, 2.8647149E-08, 7.9528315E-10,
     &  1.8774422E-11, 3.8247505E-13, 6.8427624E-15/

      integer :: n, k

!      ierr = 0
      a=x*0.88622693

      if (X < 1) then
         v=x
         S=V
         U=V
         do N=2,12
             S=S*V
            SS=S*D(N)
            U=U+SS
            if (ABS(SS) < 1.E-6*U)GOTO 11
         end do
!         write(*,*) ' COMPLETED LOOP 10 IN FDETA'
         ierr = 11
         return
!         STOP
   11    ETA=log(dble(U))

      else
         if (a < 2) then
            E=log(dble(X))
         else
            e=pow(1.5d0*a,2d0/3d0)
         end if
         do k=1,10
            call FDF1F2(E,F1,F2)
            DE=(A-F2)*2./F1
            E=E+DE
            if (abs(dE) < 1.e-4*abs(E))GOTO 21
         end do
!         write(*,*) ' completed loop 20 IN FDETA'
         ierr = 12
         return
!         stop
   21    ETA=E

      end if

      end subroutine FDETA
! **********************************************************************
      subroutine FDF1F2(ETA,F1,F2)
!
!  CALCULATES FD INTEGRALS F1, F2=F(-1/2), F(+1/2)
!  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
!  MATHS. OF COMP., 21, 30, 1967.
!
      if (ETA <= 1) then
         X=exp(dble(ETA))
         F1=X*(1.772454+X*(-1.2532215+X*(-0.60172359-X*0.0012271551))/
     &      (1.+X*(1.2958546+X*0.35469431)))
         F2=X*(0.88622693+X*(-0.31329180+X*(-0.14275695-
     &      X*0.0010090890))/
     &      (1.+X*(0.99882853+X*0.19716967)))
      else if (ETA <= 4) then
         X=ETA
         F1=(1.17909+X*(1.334367+X*1.151088))/
     &      (1.+X*(0.8975007+X*0.1153824))
         F2=(0.6943274+X*(0.4918855+X*0.214556))/
     &      (1.+X*(-0.0005456214+X*0.003648789))
      else
         X=sqrt(ETA)
         Y=1./ETA**2
         F1=X*(2.+Y*(-0.81495847+Y*(4.0521266-Y*325.43565))/
     &      (1.+Y*(-10.867628+Y*384.61501)))
         F2=ETA*X*(0.666666667+Y*(0.822713535+Y*(5.27498049+
     &      Y*290.433403))/
     &      (1.+Y*(5.69335697+Y*322.149800)))
      end if

      end subroutine FDF1F2


      subroutine screen2(ft,fne,rion,epa,ntot,umin,umax,umesh,p)
      parameter(ipz=28)
      real :: umesh(:) ! (nptot)
      real :: p(:) ! (nptot)
      dimension rion(ipz),f(100)
      dimension x(3),wt(3)
      data x/0.415775,2.294280,6.289945/
      data wt/0.711093,0.278518,0.0103893/
      data twopi/6.283185/
      COMMON/CIMESH/U(100),AA(nptot),BB(nptot),IN(nptot),ITOT,NN
      save /cimesh/

      rydt=ft/157894.
      aune=1.48185e-25*fne

!     get alp2=1/(Debye)**2
      b=0
      do i=1,ipz
        b=b+rion(i)*i**2
      end do
        alp2=(5.8804e-19)*fne*b/(epa*ft)
      if (alp2/ft < 5e-8) return !!!!!!!!!!!

      c=1.7337*aune/sqrt(rydt)

      do i=1,itot
        w=u(i)*rydt
        f(i)=0.
        do k=1,ipz
          if (rion(k) <= 0.01) cycle
          crz=c*rion(k)*k**2
          ff=0
          do j=1,3
            e=x(j)*rydt
            fk=sqrt(e)
            fkp=sqrt(e+w)
            x1=1.+alp2/(fkp+fk)**2
            x2=1.+alp2/(fkp-fk)**2
            q=(1./x2-1./x1+log(dble(x1/x2)))*
     &        (fkp*(1.-exp(dble(-twopi*k/fkp))))/(fk*(1.-exp(dble(-twopi*k/fk))))
              ff=ff+wt(j)*q
          end do
          f(i)=f(i)+crz*ff
         end do
        end do

      p(1)=f(1)
      do n=2,nn
        w=umesh(n)*rydt
        p(n)=p(n)+(aa(n)*f(in(n)-1)+bb(n)*f(in(n)))/(w*w*w)
      end do
      do n=nn+1,ntot
        w=umesh(n)*rydt
        p(n)=p(n)+f(in(n))/(w*w*w)
      end do

      end subroutine screen2

      end module op_load
