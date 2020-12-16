      subroutine set_rotation(x, nn, icontr)
c
c  sets angular velocity into common/comgrp/, for calculation of
c  rotational splitting. 
c  This is intended as a user modifiable routine. The present version
c  simply reads in the x and the angular velocity from a file and 
c  interpolates to the computational mesh which must be provided in
c  the argument x(n), n = 1, ..., nn.
c  The argument icontr may be used to control the action of the routine
c  if needed.
c
c  The file name is read from standard input. The file is assumed to be
c  in ASCII format, and structured as
c
c       r/R   omega
c
c  Note that in the present version the angular velocity is read in
c  only in the first call of set_rotation.
c
c  Original version: 15/2/06
c
      include 'adipls.c.d.incl'
      implicit double precision (a-h, o-z)
      character*280 omega_file
      dimension x(1)
      dimension xr(nnmax), omegar(nnmax)
      common/comgrp/ isprtp, irotcp, omgrtp(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data init_rot /0/
c
      write(istdou,'(/'' Entering set_rotation'')')
      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *  write(istdpr,'(/'' Entering set_rotation'')')
      if(init_rot.eq.0) then
        read(istdin,'(a)') omega_file
        if(istdpr.gt.0) write(istdpr,'(/'' Read rotation from '',a)')
     *    omega_file
        close(99)
        open(99,file=omega_file,status='old')
        call skpcom(99)
c
        n=1
c
   10   read(99,*,end=20) xr(n), omegar(n)
        n=n+1
        go to 10
c
   20   nnr=n-1
        close(99)
c
	init_rot=1
        if(istdpr.gt.0) write(istdpr,'(/'' Rotation rate read at '',i5,
     *    '' points in s/r set_rotation''/)') nnr
c
      end if
c
c  interpolate to mesh in x
c
      do n=1,nn
	call lir(x(n),xr,omgrtp(n),omegar,1,1,nnr,n,inter)
      end do
      write(istdou,'(/'' Exiting set_rotation'')')
      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *  write(istdpr,'(/'' Exiting set_rotation'')')
      return
      end
