c-----------------------------------------------------------------------
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Problem VAN DER POL
c        ODE of dimension 2
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/vdpol.f
c
c     This is revision
c     $Id: vdpol.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------
      subroutine vdpol_init(neqn,y,yprime,consis)
      integer neqn
      double precision y(neqn),yprime(neqn)
      logical consis

      y(1) = 2d0
      y(2) = 0d0
      
      return
      end
c-----------------------------------------------------------------------
      subroutine vdpol_feval(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      f(1) = y(2)
      f(2) = ((1-y(1)*y(1))*y(2)-y(1))/1.0d-3
      
      return
      end
c-----------------------------------------------------------------------
      subroutine vdpol_jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j

      dfdy(1,1) = 0d0
      dfdy(1,2) = 1d0
      dfdy(2,1) = (-2.0d0*y(1)*y(2)-1d0)/1.0d-3
      dfdy(2,2) = (1d0-y(1)*y(1))/1.0d-3

      return
      end
c-----------------------------------------------------------------------
      subroutine vdpol_solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
      ! note -- this is for stiffness param = 1d-3
      y(1) =  1.7632345401889102d+00           
      y(2) = -8.3568868191466206d-01
      return
      end
