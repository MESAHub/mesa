c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Chemical Akzo Nobel problem
c        ODE of dimension 6
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/chemakzo.f
c
c     This is revision
c     $Id: chemakzo.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------
      subroutine chemakzo_prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character*(*) fullnm, problm, type
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(*)
      double precision t(0:*)
      logical numjac, nummas

      integer i

      fullnm = 'Chemical Akzo Nobel problem'
      problm = 'chemakzo'
      type   = 'DAE'
      neqn   = 6
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 180d0
      numjac = .false.
      mljac  = neqn
      mujac  = neqn
      mlmas  = 0
      mumas  = 0
      do 10 i=1,neqn
         ind(i) = 0
   10 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine chemakzo_init(neqn,y,yprime,consis)
      integer neqn
      double precision y(neqn),yprime(neqn)
      logical consis

      integer ierr,ipar(0)
      double precision rpar(0)

      double precision k1,k2,k3,k4,kbig,kla,ks
      parameter (
     +   ks   =115.83d0
     +)

      y(1) = 0.444d0
      y(2) = 0.00123d0
      y(3) = 0d0
      y(4) = 0.007d0
      y(5) = 0d0
      y(6) = ks*y(1)*y(4)

      consis = .true.

      call chemakzo_feval(neqn,0d0,y,y,yprime,ierr,rpar,ipar)

      return
      end
c-----------------------------------------------------------------------
      subroutine chemakzo_feval(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      double precision k1,k2,k3,k4,kbig,kla,po2,hen,ks
      parameter (
     +   k1   = 18.7d0,
     +   k2   = 0.58d0,
     +   k3   = 0.09d0,
     +   k4   = 0.42d0,
     +   kbig = 34.4d0,
     +   kla  = 3.3d0,
     +   ks   = 115.83d0,
     +   po2  = 0.9d0,
     +   hen  = 737d0
     +)
      double precision r1,r2,r3,r4,r5,fin,sqy2
      include 'formats'

      if (y(2) .lt. 0d0) then
         !write(*,*) 'y(2)', y(2)
         ierr = -1
         return
      endif
      
      sqy2 = sqrt(y(2))
      r1  = k1*y(1)*y(1)*y(1)*y(1)*sqy2
      r2  = k2*y(3)*y(4)
      r3  = k2/kbig*y(1)*y(5)
      r4  = k3*y(1)*y(4)*y(4)
      r5  = k4*(y(6)**2)*sqy2
      fin = kla*(po2/hen-y(2))

      f(1) =   ((-2d0*r1 +r2) -r3)     -r4
      f(2) = ((-0.5d0*r1             -r4)     -0.5d0*r5) + fin
      f(3) =         (r1 -r2) +r3
      f(4) =            (-r2 +r3) -2d0*r4
      f(5) =             (r2 -r3)         +r5
      f(6) = ks*y(1)*y(4)-y(6)

      return
      end
c-----------------------------------------------------------------------
      subroutine chemakzo_jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      double precision k1,k2,k3,k4,kbig,kla,ks
      parameter (
     +   k1   =18.7d0,
     +   k2   =0.58d0,
     +   k3   =0.09d0,
     +   k4   =0.42d0,
     +   kbig =34.4d0,
     +   kla  =3.3d0,
     +   ks   =115.83d0
     +)
      integer i,j
      double precision r11,r12,r23,r24,r31,r35,r41,r44,r52,r56,fin2
      double precision y13

      if (y(2) .lt. 0d0) then
         ierr = -1
         return
      endif

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j) = 0d0
   10    continue
   20 continue
      
      y13  = y(1)*y(1)*y(1)
      r11  = 4d0*k1*y13*sqrt(y(2))
      r12  = 0.5d0*k1*y(1)*y13/sqrt(y(2))
      r23  = k2*y(4)
      r24  = k2*y(3)
      r31  = (k2/kbig)*y(5)
      r35  = (k2/kbig)*y(1)
      r41  = k3*y(4)**2
      r44  = 2d0*k3*y(1)*y(4)
      r52  = 0.5d0*k4*(y(6)**2)/sqrt(y(2))
      r56  = 2d0*k4*y(6)*sqrt(y(2))
      fin2 = -kla

      dfdy(1,1) = (-2d0*r11-r31)-r41
      dfdy(1,2) = -2d0*r12
      dfdy(1,3) = r23
      dfdy(1,4) = r24-r44
      dfdy(1,5) = -r35
      dfdy(2,1) = -0.5d0*r11-r41
      dfdy(2,2) = -0.5d0*r12-0.5d0*r52+fin2
      dfdy(2,4) = -r44
      dfdy(2,6) = -0.5d0*r56
      dfdy(3,1) = r11+r31
      dfdy(3,2) = r12
      dfdy(3,3) = -r23
      dfdy(3,4) = -r24
      dfdy(3,5) = r35
      dfdy(4,1) = r31-2d0*r41
      dfdy(4,3) = -r23
      dfdy(4,4) = -r24-2d0*r44
      dfdy(4,5) = r35
      dfdy(5,1) = -r31
      dfdy(5,2) = r52
      dfdy(5,3) = r23
      dfdy(5,4) = r24
      dfdy(5,5) = -r35
      dfdy(5,6) = r56
      dfdy(6,1) = ks*y(4)
      dfdy(6,4) = ks*y(1)
      dfdy(6,6) = -1d0

      return
      end
c-----------------------------------------------------------------------
      subroutine chemakzo_meval(ldim,neqn,t,yprime,dfddy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,yprime(neqn),dfddy(ldim,neqn),rpar(*)

      integer i
      ierr = 0

      do 10 i=1,neqn-1
         dfddy(1,i)=1d0
   10 continue

      dfddy(1,neqn)=0d0

      return
      end
c-----------------------------------------------------------------------
      subroutine chemakzo_solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)

c computed using true double precision on a Cray C90

c Solving Chemical Akzo Nobel problem using PSIDE

c User input:

c give relative error tolerance: 1d-19
c give absolute error tolerance: 1d-19

c Integration characteristics:

c    number of integration steps        5535
c    number of accepted steps           5534
c    number of f evaluations          100411
c    number of Jacobian evaluations        7
c    number of LU decompositions         368

c CPU-time used:                          30.77 sec
c

      y(  1) =  0.1150794920661702d0
      y(  2) =  0.1203831471567715d-2
      y(  3) =  0.1611562887407974d0
      y(  4) =  0.3656156421249283d-3
      y(  5) =  0.1708010885264404d-1
      y(  6) =  0.4873531310307455d-2

      return
      end
