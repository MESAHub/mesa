!-----------------------------------------------------------------------
!     This file is part of the Test Set for IVP solvers
!     http://www.dm.uniba.it/~testset/
!
!        Problem VAN DER POL
!        ODE of dimension 2
!
!     DISCLAIMER: see
!     http://www.dm.uniba.it/~testset/disclaimer.php
!
!     The most recent version of this source file can be found at
!     http://www.dm.uniba.it/~testset/src/problems/vdpol.f
!
!     This is revision
!     $Id: vdpol.F,v 1.2 2006/10/02 10:29:14 testset Exp $
!
!-----------------------------------------------------------------------

module bari_vdpol

   implicit none

contains

   subroutine vdpol_init(neqn, y, yprime, consis)
      integer :: neqn
      double precision :: y(neqn), yprime(neqn)
      logical :: consis

      y(1) = 2d0
      y(2) = 0d0

   end subroutine vdpol_init

   subroutine vdpol_feval(neqn, t, y, yprime, f, ierr, rpar, ipar)
      integer :: neqn, ierr, ipar(*)
      double precision :: t, y(neqn), yprime(neqn), f(neqn), rpar(*)

      f(1) = y(2)
      f(2) = ((1 - y(1)*y(1))*y(2) - y(1))/1.0d-3

   end subroutine vdpol_feval

   subroutine vdpol_jeval(ldim, neqn, t, y, yprime, dfdy, ierr, rpar, ipar)
      integer :: ldim, neqn, ierr, ipar(*)
      double precision :: t, y(neqn), yprime(neqn), dfdy(ldim, neqn), rpar(*)

      dfdy(1, 1) = 0d0
      dfdy(1, 2) = 1d0
      dfdy(2, 1) = (-2.0d0*y(1)*y(2) - 1d0)/1.0d-3
      dfdy(2, 2) = (1d0 - y(1)*y(1))/1.0d-3

   end subroutine vdpol_jeval

   subroutine vdpol_solut(neqn, t, y)
      integer :: neqn
      double precision, intent(in) :: t
      double precision, intent(out) :: y(neqn)

      ! note -- this is for stiffness param = 1d-3
      y(1) = 1.7632345401889102d+00
      y(2) = -8.3568868191466206d-01

   end subroutine vdpol_solut

end module bari_vdpol
