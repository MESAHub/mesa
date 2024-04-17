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

c FORTRAN 90 module for calculation of radiative accelerations,
c based on the Opacity Project (OP) code "OPserver".
c See CHANGES_HU for changes made to the original code.
c
c Haili Hu 2010
c
      module op_load
      use math_lib
      use op_def
      logical :: have_loaded_op = .false.


      contains
C******************************************************************
      subroutine op_dload(path, cache_filename, ierr)
        implicit none
      character (len=*), intent(in) :: path, cache_filename
      integer, intent(out) :: ierr




      integer,parameter :: ipz=28
      real :: am,amm,delp,dpack
      integer :: ios,it,ite11,ite22,ite33,itt,itte1,itte2,itte3,izz,jne,ite
      integer :: jne1,jne22,jne33,jnn,k,n,ncount2,ncount3,ja,jn,jnw11
      integer :: jne11,jne2,nccc,ne,nfff,ntott,nn
      real :: orss,um,ux,umaxx,uminn,u
      real,dimension(nptot):: umesh, semesh
      integer :: ntotv
      real :: dv,dv1

      integer :: cache_version

       common /mesh/ ntotv,dv,dv1,umesh,semesh
!    common /atomdata/
!       common/atomdata/ ite1,ite2,ite3,jn1(91),jn2(91),jne3,umin,umax,ntot,
!      +  nc,nf,int(17),epatom(17,91,25),oplnck(17,91,25),ne1(17,91,25),
!      +  ne2(17,91,25),fion(-1:28,28,91,25),np(17,91,25),kp1(17,91,25),
!      +  kp2(17,91,25),kp3(17,91,25),npp(17,91,25),mx(33417000),
!      +  yy1(33417000),yy2(120000000),nx(19305000),yx(19305000)
!

      integer,dimension(ipe) :: ifl,iflp
      character num(0:9)*1,zlab(ipe)*3,tlab*6,zlabp(ipe)*3
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/

      integer :: kz(17)
      data kz/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28/      
      save  /mesh/    !HH: put common block in static memory
c
      integer :: nx_temp(nptot)
      real :: y_temp(nptot)
      integer :: nx_index, left_n, right_n, n_index
      real :: left_val, right_val, cross_section, slope

      if(allocated(yy2) .eqv. .false.) then
         ! yy2 actually needs 29,563 x 10,000 length
         ALLOCATE(yy2(30000*10000),nx(19305000),yx(19305000),stat=ierr)
         if(ierr/=0) return
         yy2=0.0
         nx=0.0
         yx=0.0
   !     write(*,*) "ierr",ierr
      end if


      ierr=0
      if (have_loaded_op) return

!$omp critical (critial_do_op_dload)

      if (have_loaded_op) goto 1001


      !path = '../OP4STARS_1.3'
      !call getenv("oppath", path)
      !if (len(trim(path)) == 0) then
      !   write(6,*) 'Define environmental variable oppath (directory of OP data)' 
      !   stop
      !endif


      ios = 0
      open(1,file=trim(cache_filename),action='read',
     >         status='old',iostat=ios,form='unformatted')
      if (ios == 0) then
         write(*,*) 'reading OP cache file ' // trim(cache_filename)
         read(1,iostat=ios) cache_version, ntotv,dv,dv1,umesh,
     >      ite1,ite2,ite3,jn1,jn2,jne3,umin,umax,ntotp,nc,nf,int,epatom,oplnck, ne1p,
     >      ne2p,fionp,np,kp1,kp2,kp3,npp,yy2,nx,yx
         write(*,*) 'done reading OP cache file'
         close(1)
         if (cache_version .ne. op_cache_version) then
            write(*,*) 'wrong version of OP cache'
            write(*,*) 'cache file path is set by op_mono_data_cache_filename'
            write(*,*) 'perhaps cache is shared between different MESA versions'
            write(*,*) 'please delete cache file and try again'
            stop 'op_load'
         end if
         if (ios == 0) then
            have_loaded_op = .true.
            CALL IMESH(UMESH,NTOTV)
            goto 1001
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

C  READ INDEX FILES
C     FIRST FILE
      NN=1
c     print*,' Opening '//'./'//zlab(1)//'.index'
      OPEN(1,FILE=trim(path)//'/'//ZLAB(1)//'.index',STATUS='OLD',
     + iostat=ios)
      if (ios /= 0) then
         write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(1) // '.index'
         ierr = -1
         goto 1001
      end if
      READ(1,*)IZZ,AMM
      READ(1,*)ITTE1,ITTE2,ITTE3
      READ(1,*)UMIN,UMAX
      READ(1,*)NC,NF
      READ(1,*)DPACK
      CLOSE(1)
      IF(IZZ.NE.KZ(1))then
         write(6,6001)zlab(1),izz,nn,kz(1)
         ierr=1
         goto 1001
      endif
      NTOTP=NF
      IF(NTOTP.GT.nptot)then
         write(6,6002)ntotp,nptot
         ierr=2
         goto 1001
      endif
      INT(1)=1
      IF(ITTE3.NE.ITE3)then
         write(6,6077)ite3,itte3,nn
         ierr=3
         goto 1001
      endif
c      ITE1=MAX(ITE1,ITTE1)
c      ITE2=MIN(ITE2,ITTE2)
C
c  READ MESH FILES
      OPEN(1,FILE=trim(path)//'/'//ZLAB(1)//'.mesh',status='old',
     + form='unformatted',iostat=ios)
      if (ios /= 0) then
         write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(1) // '.mesh'
         ierr = -1
         goto 1001
      end if
      READ(1)DV,NTOTV,(UMESH(N),N=1,NTOTV)
      umin=umesh(1)
      umax=umesh(ntotv)
      DV1=DV
      CLOSE(1)
C
C  GET MESH FOR SCREEN
      CALL IMESH(UMESH,NTOTV)
C
C     SUBSEQUENT FILES
      DO 40 N=2,ipe
         NN=N
         OPEN(1,FILE=trim(path)//'/'//ZLAB(N)//'.index',
     +   STATUS='OLD')
         READ(1,*)IZZ,AMM
         READ(1,*)ITE11,ITE22,ITE33
         READ(1,*)UMINN,UMAXX
         READ(1,*)NC,NF
         READ(1,*)DPACK
         CLOSE(1)
         IF(ITE33.NE.ITE3)then
            write(6,6077)ite3,ite33,nn
            ierr=4
            goto 1001
         endif
c         ITE1=MAX(ITE1,ITE11)
c         ITE2=MIN(ITE2,ITE22)
         IF(IZZ.NE.KZ(N))then
            write(6,6001)zlab(n),izz,nn,kz(nn)
            ierr=5
            goto 1001
         endif
         NTOTT=NF
         IF(NTOTT.GT.NTOTP)then
            write(6,6006)nn,ntott,ntotp
            ierr=6
            goto 1001
         endif
c!!         IF(UMIN.NE.UMINN.OR.UMAX.NE.UMAXX) GOTO 1003  !!
         INT(N)=NTOTP/NTOTT
         IF(INT(N)*NTOTT.NE.NTOTP)then
            WRITE(6,6009)NN,NTOTT,NTOTP
            ierr=7
            goto 1001
         endif
         IF(INT(N).NE.1)WRITE(6,6007)N,INT(N)
c
c        READ MESH FILES
c
         OPEN(1,FILE=trim(path)//'/'//ZLAB(N)//'.mesh',
     +   status='old',form='unformatted',iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open ' // trim(path) // '/' // ZLAB(N) // '.mesh'
            ierr = -1
            goto 1001
         end if
         READ(1)DV
       IF(DV.NE.DV1)THEN
!          PRINT*,' OP: N=',N,', DV=',DV,' NOT EQUAL TO DV1=',DV1
            ierr=8
          goto 1001
       ENDIF
         CLOSE(1)
   40 CONTINUE
C
C  START TEMPERATURE LOOP
C
      ncount2=0
      ncount3=0
      do it=ite1,ite2,ite3
c
C        OPEN FILES
c
         TLAB='.'//NUM(IT/100)//NUM(IT/10-10*(IT/100))//
     +    NUM(IT-10*(IT/10))
         do n=1,ipe
c            IF(SKIP(N))GOTO 70
            NN=N
            OPEN(IFL(N),FILE=trim(path)//'/'//ZLAB(N)//TLAB,
     +      FORM='UNFORMATTED',STATUS='OLD',iostat=ios)
            if (ios /= 0) then
               write(*,*) 'failed to open ' // trim(path)//'/'//ZLAB(N)//TLAB
               ierr = -1
               goto 1001
            end if
            if(n.gt.2) then
            OPEN(IFLP(N),FILE=trim(path)//'/'//ZLABP(N)//TLAB,
     +      FORM='UNFORMATTED',STATUS='OLD',iostat=ios)
            if (ios /= 0) then
               write(*,*) 'failed to open ' // trim(path)//'/'//ZLABP(N)//TLAB
               ierr = -1
               goto 1001
            end if
            endif
         end do
C        READ HEADINGS
         NN=1
         READ(IFL(1))IZZ,ITE,AM,UM,UX,NCCC,NFFF,DelP,JNE1,JNE2,JNE3
         do n=2,ipe
c            IF(SKIP(N))GOTO 80
            NN=N
            READ(IFL(N))IZZ,ITE,AM,UM,UX,NC,NF,DelP,JNE11,JNE22,JNE33
            if(n.gt.2) read(iflp(n))
            IF(JNE33.NE.JNE3)then
              write(6,6099)jne3,jne33,nn
              ierr=9
              goto 1001
            endif
            JNE1=MAX(JNE1,JNE11)
            JNE2=MIN(JNE2,JNE22)
         end do
         itt=(it-ite1)/2+1
         jn1(itt)=jne1
         jn2(itt)=jne2
C
c         WRITE(98,9802)ITE,JNE1,JNE2,JNE3
C
C        START DENSITY LOOP
C
         do n=1,ipe
           do jn=jne1,jne2,jne3
             jnn=(jn-jne1)/2+1
C
C           START LOOP ON ELEMENTS
C
   95        READ(IFL(N))JNE,EPATOM(n,itt,jnn),OPLNCK(n,itt,jnn),ORSS,
     +         NE1P(n,itt,jnn),NE2P(n,itt,jnn),
     +         (FIONP(NE,n,itt,jnn),NE=NE1P(n,itt,jnn),NE2P(n,itt,jnn))
             read(ifl(n))np(n,itt,jnn)
             if(np(n,itt,jnn).gt.0)then
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
             endif
               if(n.gt.2) then
                 read(iflp(n))ja,npp(n,itt,jnn)
                 if(npp(n,itt,jnn).gt.0) then
           read(iflp(n))(nx(k+ncount3),yx(k+ncount3),k=1,npp(n,itt,jnn))
                   kp3(n,itt,jnn)=ncount3
                   ncount3=ncount3+npp(n,itt,jnn)
                 endif
               endif
            end do
          end do
c
c     write(6,610)it
c     write(6,*)'ncount1 = ',ncount1
c     write(6,*)'ncount2 = ',ncount2
c     write(6,*)'ncount3 = ',ncount3
c
C        CLOSE FILES
c
         DO 150 N=1,ipe
            CLOSE(IFL(N))
            close(iflp(n))
  150    CONTINUE
c
      end do

      write(*,*) 'done loading OP mono data'
      have_loaded_op = .true.

      !write(*,*)'ncount1 = ',ncount1
      !write(6,*)'ncount2 = ',ncount2
      !write(6,*)'ncount3 = ',ncount3
      ios = 0
      open(1, file=trim(cache_filename), iostat=ios,
     >         action='write', form='unformatted')
      if (ios == 0) then
         write(*,*) 'write ' // trim(cache_filename)
         write(1) op_cache_version, ntotv,dv,dv1,umesh,
     >      ite1,ite2,ite3,jn1,jn2,jne3,umin,umax,ntotp,nc,nf,int,epatom,oplnck, ne1p, 
     >      ne2p,fionp,np,kp1,kp2,kp3,npp,yy2,nx,yx
         close(1)
      end if

1001  continue

C     pre-calculate semesh
      do n = 1, nptot
         u = umesh(n)
         semesh(n) = 1.d0 - exp(dble(-u))
      end do

!$omp end critical (critial_do_op_dload)



      return
610   format(10x,'Done IT= ',i3)
1004  WRITE(6,6004)ZLAB(NN),TLAB
      STOP
6001  FORMAT(//5X,'*** OP: FILE ',A3,' GIVES IZZ=',I3,
     + 'NOT EQUAL TO IZ(',I2,')=',I2,' ***')
6002  FORMAT(//5X,'*** OP: NTOT=',I7,' GREATER THAN nptot=',
     + I7,' ***')
6003  FORMAT(//5X,'*** OP: DISCREPANCY BETWEEN DATA ON FILES ',
     + A3,' AND ',A3,' ***')
6004  FORMAT(//5X,'*** OP: ERROR OPENING FILE ',A3,A6,'  ***')
6006  FORMAT(//5X,'OP: N=',I2,', NTOTT=',I7,', GREATER THAN NTOT=',I7)
6007  FORMAT(/5X,'OP: N=',I2,', INT(N)=',I4)
6009  FORMAT(' OP: N=',I5,', NTOTT=',I10,', NTOT=',I10/
     + '   NTOT NOT MULTIPLE OF NTOTT')
c6012  FORMAT(/10X,'ERROR, SEE WRITE(6,6012)'/
c     + 10X,'IT=',I3,', JN=',I3,', N=',I3,', JNE=',I3/)
6077  FORMAT(//5X,'OP: DISCREPANCY IN ITE3'/10X,I5,' READ FROM UNIT 5'/
     + 10X,I5,' FROM INDEX FILE ELEMENT',I5)
6099  FORMAT(//5X,'OP: DISCREPANCY IN JNE3'/10X,I5,' READ FOR N=1'/
     + 10X,I5,' READ FOR N=',I5)

c8000  FORMAT(5X,I5,F10.4/5X,3I5/2E10.2/2I10/10X,E10.2)
1010  print*,' OP: ERROR OPENING FILE '//'./'//ZLAB(1)//'.index'
      stop
1011  print*,' OP: ERROR OPENING FILE '//'./'//ZLAB(1)//'.mesh'
      stop
      end subroutine op_dload

c***********************************************************************
        SUBROUTINE IMESH(UMESH,NTOT)
C
      DIMENSION UMESH(nptot)
      COMMON/CIMESH/U(100),AA(nptot),BB(nptot),IN(nptot),ITOT,NN
      save /cimesh/

      UMIN=UMESH(1)
      UMAX=UMESH(NTOT)
c
      II=100
      A=(II*UMIN-UMAX)/REAL(II-1)
      B=(UMAX-UMIN)/REAL(II-1)
      DO I=1,II
        U(I)=A+B*I
      ENDDO
c
      ib=2
      ub=u(ib)
      ua=u(ib-1)
      d=ub-ua
      ibb=0
      do n=2,ntot
        if(umesh(n).gt.ub)then
          ua=ub
          ib=ib+1
          ub=u(ib)
          d=ub-ua
          if(umesh(n).gt.ub)then
            nn=n-1
            ibb=ib-1
            goto 1
          endif
        endif
        in(n)=ib
        aa(n)=(ub-umesh(n))/d
        bb(n)=(umesh(n)-ua)/d
      end do
c
    1      ib=ibb
      do n=nn+1,ntot
        ib=ib+1
        in(n)=ib
        u(ib)=umesh(n)
      end do
      itot=ib
c
        return
      end SUBROUTINE IMESH


      subroutine msh(dv, ntot, umesh, semesh, uf, dscat)
      implicit none
      integer, intent(out) :: ntot
      real, intent(out) :: dv, uf(0:100), dscat
      real, intent(out) :: umesh(:), semesh(:) ! (nptot)
      integer :: i, k, ntotv
      real :: dvp, dv1, umin, umax, umeshp(nptot), semeshp(nptot)
      common /mesh/ ntotv, dvp, dv1, umeshp, semeshp
      save /mesh/
c
      ntot = ntotv
      dv = dvp
      do i=1,ntot
         umesh(i) = umeshp(i)
      end do
      do i=1,ntot
         semesh(i) = semeshp(i)
      end do

c
      umin = umesh(1)
      umax = umesh(ntot)
      dscat = (umax - umin)*0.01
      do i = 0, 100
         uf(i) = umin + i*dscat
      end do
c
      return
c
      end subroutine msh


      subroutine solve(u,v,z,uz,ierr)
      integer, intent(inout) :: ierr
      dimension u(4)
c
c  If  P(R) =   u(1)  u(2)  u(3)  u(4)
c  for   R  =    -3    -1    1     3
c  then a cubic fit is:
      P(R)=(
     +  27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(
     +  27*(u(3)-u(2))-(u(4)-u(1))   +R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
c  First derivative is:
      PP(R)=(
     +  27*(u(3)-u(2))-(u(4)-u(1))+ 2*R*(
     +  -3*(u(2)+u(3))+3*(u(4)+u(1)) +3*R*(
     +  -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.
c
!      ierr = 0
c  Find value of z giving P(z)=v
c  First estimate
      z=(2.*v-u(3)-u(2))/(u(3)-u(2))
c  Newton-Raphson iterations
      do k=1,10
         uz=pp(z)
         d=(v-p(z))/uz
         z=z+d
         if(abs(d).lt.1.e-4)return
      end do
c
!      print*,' Not converged after 10 iterations in SOLVE'
!      print*,' v=',v
!      DO N=1,4
!         PRINT*,' N, U(N)=',N,U(N)
!      ENDDO
      ierr = 10
      return
!      stop
c
      end subroutine solve
c***********************************************************************

      SUBROUTINE BRCKR(T,FNE,RION,NION,U,NFREQ,SF, ierr)
      integer, intent(inout) :: ierr
C
C  CODE FOR COLLECTIVE EFFECTS ON THOMSON SCATTERING.
C  METHOD OF D.B. BOERCKER, AP. J., 316, L98, 1987.
C
C  INPUT:-
C     T=TEMPERATURTE IN K
C     FNE=ELECTRON DENSITY IN CM**(-3)
C     ARRAY RION (DIMENSIONED FOR 30 IONS).
C        RION(IZ) IS NUMBER OF IONS WITH NET CHARGE IZ.
C        NORMALISATION OF RION IS OF NO CONSEQUENCE.
C     NION=NUMBER OF IONS INCLUDED.
C     ARRAY U (DIMENSIONED FOR 1000). VALUES OF (H*NU/K*T).
C     NFREQ=NUMBER OF FREQUENCY POINTS.
C
C  OUTPUT:-
C     ARRAY SF, GIVING FACTORS BY WHICH THOMSON CROSS SECTION
C     SHOULD BE MULTIPLIED TO ALLOW FOR COLLECTIVE EFFECTS.
C
C  MODIFFICATIONS:-
C     (1) REPLACE (1.-Y) BY EXP(-Y) TO AVOID NEGATIVE FACTORS FOR
C         HIGHLY-DEGENERATE CASES.
C     (2) INCLUDE RELATIVISTIC CORRECTION.
C
      PARAMETER (IPZ=28,IPNC=100)
      DIMENSION RION(IPZ),U(0:IPNC),SF(0:IPNC)
C
      AUNE=1.48185E-25*FNE
      AUT=3.16668E-6*T
      C1=-1.0650E-4*AUT
      C2=+1.4746E-8*AUT**2
      C3=-2.0084E-12*AUT*AUT*AUT
      V=7.8748*AUNE/(AUT*SQRT(AUT))
      CALL FDETA(V,ETA, ierr)  ! 23.10.93
      W=exp(dble(ETA))         ! 23.10.93
   11 R=FMH(W)/V
      A=0.
      B=0.
      DO 20 I=1,NION
         A=A+I*RION(I)
         B=B+I**2*RION(I)
   20 CONTINUE
      X=R+B/A
C
      Y=.353553*W
      C=1.1799E5*X*AUNE/(AUT*AUT*AUT)
      DO 30 N=0,NFREQ
         D=C/U(N)**2
         IF(D.GT.5.)THEN
            D=-2./D
            F=2.666667*(1.+D*(.7+D*(.55+.341*D)))
         ELSE
            G=2.*D*(1+D)
            F=D*((G+D*D*D)*LOG(dble(D/(2.+D)))+G+2.6666667)
         ENDIF
         DELTA=.375*R*F/X
         SF(N)=(1.-R*DELTA-Y*FUNS(W))*
     +   (1.+U(N)*(C1+U(N)*(C2+U(N)*C3)))   !SAMPSON CORRECTION
   30 CONTINUE
C
      RETURN
C
  600 FORMAT(5X,'NOT CONVERGED IN LOOP 10 OF BRCKR'/
     +       5X,'T=',1P,E10.2,', FNE=',E10.2)
C
      END SUBROUTINE BRCKR
C***********************************************************************
      FUNCTION FUNS(A)
C
      IF(A.LE.0.001)THEN
         FUNS=1.
      ELSEIF(A.LE.0.01)THEN
         FUNS=(1.+A*(-1.0886+A*(1.06066+A*1.101193)))/
     +     (1.+A*(0.35355+A*(0.19245+A+0.125)))
      ELSE
         FUNS=(  1./(1.+0.81230*A)**2+
     +        0.92007/(1.+0.31754*A)**2+
     +        0.05683/(1.+0.04307*A)**2 )/
     +     (  1./(1.+0.65983*A)+
     +        0.92007/(1.+0.10083*A)+
     +        0.05683/(1.+0.00186*A)    )
      ENDIF
      RETURN
      END FUNCTION FUNS
C***********************************************************************
      FUNCTION FMH(W)
C
C  CALCULATES FD INTERGAL I_(-1/2)(ETA). INCLUDES FACTOR 1/GAMMA(1/2).
C  ETA=LOG(W)
C
      IF(W.LE.2.718282)THEN
         FMH=W*(1+W*(-.7070545+W*(-.3394862-W*6.923481E-4))
     +   /(1.+W*(1.2958546+W*.35469431)))
      ELSEIF(W.LE.54.59815)THEN
         X=LOG(dble(W))
         FMH=(.6652309+X*(.7528360+X*.6494319))
     +   /(1.+X*(.8975007+X*.1153824))
      ELSE
         X=LOG(dble(W))
         Y=1./X**2
         FMH=SQRT(X)*(1.1283792+(Y*(-.4597911+Y*(2.286168-Y*183.6074)))
     +   /(1.+Y*(-10.867628+Y*384.61501)))
      ENDIF
C
      RETURN
      END FUNCTION FMH
C***********************************************************************
      SUBROUTINE FDETA(X,ETA, ierr)
C
C  GIVEN X=N_e/P_e, CALCULATES FERMI-DIRAC ETA
C  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
C  MATHS. OF COMP., 21, 30, 1967.
C
      integer, intent(inout) :: ierr
      DIMENSION D(2:12)
      DATA D/
     +  3.5355339E-01, 5.7549910E-02, 5.7639604E-03, 4.0194942E-04,
     +  2.0981899E-05, 8.6021311E-07, 2.8647149E-08, 7.9528315E-10,
     +  1.8774422E-11, 3.8247505E-13, 6.8427624E-15/
C
      integer n,k
c
!      ierr = 0
      a=x*0.88622693
c
      IF(X.LT.1)THEN
         v=x
         S=V
         U=V
         DO 10 N=2,12
             S=S*V
            SS=S*D(N)
            U=U+SS
            IF(ABS(SS).LT.1.E-6*U)GOTO 11
   10    CONTINUE
!         PRINT*,' COMPLETED LOOP 10 IN FDETA'
         ierr = 11
         return
!         STOP
   11    ETA=LOG(dble(U))
c
      ELSE
         if(a.lt.2)then
            E=LOG(dble(X))
         else
            e=pow(1.5d0*a,2d0/3d0)
         endif
         do 20 k=1,10
            CALL FDF1F2(E,F1,F2)
            DE=(A-F2)*2./F1
            E=E+DE
            if(abs(dE).lt.1.e-4*abs(E))goto 21
   20    continue
!         print*,' completed loop 20 IN FDETA'
         ierr = 12
         return
!         stop
   21    ETA=E
c
      ENDIF
C
      RETURN
      END SUBROUTINE FDETA
C***********************************************************************
      SUBROUTINE FDF1F2(ETA,F1,F2)
C
C  CALCULATES FD INTEGRALS F1, F2=F(-1/2), F(+1/2)
C  USE CHEBYSHEV FITS OF W.J. CODY AND H.C. THACHER,
C  MATHS. OF COMP., 21, 30, 1967.
C
      IF(ETA.LE.1)THEN
         X=exp(dble(ETA))
         F1=X*(1.772454+X*(-1.2532215+X*(-0.60172359-X*0.0012271551))/
     +      (1.+X*(1.2958546+X*0.35469431)))
         F2=X*(0.88622693+X*(-0.31329180+X*(-0.14275695-
     +      X*0.0010090890))/
     +      (1.+X*(0.99882853+X*0.19716967)))
      ELSEIF(ETA.LE.4)THEN
         X=ETA
         F1=(1.17909+X*(1.334367+X*1.151088))/
     +      (1.+X*(0.8975007+X*0.1153824))
         F2=(0.6943274+X*(0.4918855+X*0.214556))/
     +      (1.+X*(-0.0005456214+X*0.003648789))
      ELSE
         X=SQRT(ETA)
         Y=1./ETA**2
         F1=X*(2.+Y*(-0.81495847+Y*(4.0521266-Y*325.43565))/
     +      (1.+Y*(-10.867628+Y*384.61501)))
         F2=ETA*X*(0.666666667+Y*(0.822713535+Y*(5.27498049+
     +      Y*290.433403))/
     +      (1.+Y*(5.69335697+Y*322.149800)))
      ENDIF
C
      RETURN
      END SUBROUTINE FDF1F2


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
c
      rydt=ft/157894.
      aune=1.48185e-25*fne
c
c       get alp2=1/(Debye)**2
      b=0
      do i=1,ipz
        b=b+rion(i)*i**2
      end do
        alp2=(5.8804e-19)*fne*b/(epa*ft)
      if(alp2/ft.lt.5e-8)return !!!!!!!!!!!
c
      c=1.7337*aune/sqrt(rydt)
c
      do i=1,itot
        w=u(i)*rydt
        f(i)=0.
        do 1 k=1,ipz
          if(rion(k).le.0.01)goto 1
          crz=c*rion(k)*k**2
          ff=0
          do j=1,3
            e=x(j)*rydt
            fk=sqrt(e)
            fkp=sqrt(e+w)
            x1=1.+alp2/(fkp+fk)**2
            x2=1.+alp2/(fkp-fk)**2
            q=(1./x2-1./x1+LOG(dble(x1/x2)))*
     +        (fkp*(1.-exp(dble(-twopi*k/fkp))))/(fk*(1.-exp(dble(-twopi*k/fk))))
              ff=ff+wt(j)*q
          end do
          f(i)=f(i)+crz*ff
    1     continue
        end do
c
      p(1)=f(1)
      do n=2,nn
        w=umesh(n)*rydt
        p(n)=p(n)+(aa(n)*f(in(n)-1)+bb(n)*f(in(n)))/(w*w*w)
      end do
      do n=nn+1,ntot
        w=umesh(n)*rydt
        p(n)=p(n)+f(in(n))/(w*w*w)
      end do
c
      return
      end subroutine screen2

      end module op_load
