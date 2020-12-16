! ***********************************************************************
!
!   Copyright (C) 2008-2019  Bill Paxton & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module vdpol
      
      ! for information about this problem,
      ! see http://pitagora.dm.uniba.it/~testset/problems/vdpol.php

      use num_def
      use num_lib
      use mtx_lib
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none
      
      ! stiffness parameter
      real(dp), parameter :: mu = 1d-3
      
         
      
      contains
      
      
      subroutine solve_vdpol
      
         ! args for isolve -- see num_isolve.dek in num/public
      
         integer, parameter :: which_solver = ros3p_solver ! as defined in num_def.f
      
         integer, parameter :: n = 2 ! the number of variables in the "vdpol" system of ODEs

         real(dp) :: x 
            ! input: initial x value
            ! output: x value for which the solution has been computed.
         real(dp), pointer :: y(:) 
            ! input: initial values for y
            ! output: values of y for final value of x.
         real(dp) :: xend ! desired final x value (positive or negative)
         real(dp) :: h 
            ! input: initial step size guess
            ! output: predicted next step size from the last accepted step
         real(dp) :: max_step_size
         integer :: max_steps
      
         ! absolute and relative error tolerances
         real(dp) :: rtol(1), atol(1) 
         integer :: itol

         ! information about the jacobian matrix
         integer :: ijac, nzmax, isparse, mljac, mujac

         ! information about the "mass" matrix
         integer :: imas, mlmas, mumas
      
         ! switch for calling the subroutine solout or nor
         integer :: iout
      
         integer :: lrd, lid
         real(dp), pointer :: rpar_decsol(:) ! (lrd)
         integer, pointer :: ipar_decsol(:) ! (lid)

         integer :: caller_id, nvar, nz
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
         
      
         ! work arrays.
         integer :: lwork, liwork
         real(dp), pointer :: work(:) ! (lwork)
         integer, pointer :: iwork(:) ! (liwork)
         
         ! parameter arrays.
         integer, parameter :: lrpar = 1, lipar = 3
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
               
         ! io unit for warnings and errors
         integer :: lout
      
         ! result code
         integer :: idid
      
         integer :: ierr, i
         real(dp) :: yexact(n), y_min, y_max
         real(dp), target :: y_ary(n)
         
         ipar => ipar_ary
         rpar => rpar_ary
         y => y_ary

         nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
         caller_id = 0
         nvar = 0
         nz = 0
         
         x = 0
         
         y(1) = 2d0
         y(2) = 0d0

         xend = 2d0

         h = 1d-10

         max_step_size = 0 
         max_steps = 500000

         rtol(1) = 1d-8
         atol(1) = 1d-8
         itol = 0
         
         y_min = -1d199
         y_max = 1d199
         
         ijac = 1
         nzmax = 0
         isparse = 0
         mljac = n ! square matrix
         mujac = n

         imas = 0
         mlmas = 0
         mumas = 0        
         
         iout = 1
         
         lid = 0
         lrd = 0

         ipar = 0
         rpar = 0         

         lout = 6

         call lapack_work_sizes(n, lrd, lid)

         call isolve_work_sizes(n, nzmax, imas, mljac, mujac, mlmas, mumas, liwork, lwork)
         
         allocate(iwork(liwork), work(lwork), ipar_decsol(lid), rpar_decsol(lrd), stat=ierr)
         if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
      
         iwork = 0
         work = 0
         
         write(*,*)
         write(*,*) 'vdpol'
         write(*,*)
         
         call isolve( &
            which_solver, n, vdpol_derivs, x, y, xend, & 
            h, max_step_size, max_steps, & 
            rtol, atol, itol, y_min, y_max, & 
            vdpol_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, & 
            null_mas, imas, mlmas, mumas, & 
            vdpol_solout, iout, & 
            lapack_decsol, null_decsols, null_decsolblk, &
            lrd, rpar_decsol, lid, ipar_decsol, &  
            caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, & 
            null_fcn_blk_dble, null_jac_blk_dble, &
            work, lwork, iwork, liwork, & 
            lrpar, rpar, lipar, ipar, & 
            lout, idid)
            
         if (idid /= 1) ierr = -1
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         
         write(*,*)
         write(*,*) 'nsteps', iwork(16)
         
         deallocate(iwork, work, ipar_decsol, rpar_decsol)
         
         ! expected solution for stiffness param = 1d-3
         yexact(1) =  1.7632345401889102d+00           
         yexact(2) = -8.3568868191466206d-01
         
         write(*,'(/,a5,99a20)') 'i', 'calculated    ', 'reference    ', 'lg(abs rel diff)'
         do i=1, n
            write(*,'(i5,2e20.10,f20.10)') i, y(i), yexact(i), &
                  log10(abs(y(i)-yexact(i))/max(1d-299, abs(yexact(i))))
         end do
         write(*,*)
         write(*,*)

      end subroutine solve_vdpol


      subroutine vdpol_derivs(n, x, h, vars, dvars_dx, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: vars(:) ! (n)
         real(dp), intent(inout) :: dvars_dx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         dvars_dx(1) = vars(2)
         dvars_dx(2) = ((1-vars(1)*vars(1))*vars(2)-vars(1))/mu
      end subroutine vdpol_derivs


      subroutine vdpol_jacob(n, x, h, y, f, dfdy, ld_dfdy, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:), dfdy(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         integer :: nz, i, j
         ierr = 0
         dfdy(1, 1) = 0d0
         dfdy(1, 2) = 1d0
         dfdy(2, 1) = (-2.0d0*y(1)*y(2)-1d0)/mu
         dfdy(2, 2) = (1d0-y(1)*y(1))/mu
         call vdpol_derivs(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
      end subroutine vdpol_jacob


      subroutine vdpol_solout(nr, xold, x, n, y, rwork, iwork, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            real(dp) function interp_y(i, s, rwork, iwork, ierr)
               use const_def, only: dp
               integer, intent(in) :: i
               real(dp), intent(in) :: s
               real(dp), intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn
         
         real(dp) :: xout, y1, y2
         integer :: ierr
         
         ierr = 0
         irtrn = 0
         xout = rpar(1)
         if (nr.eq.1) then
            write (6, 99) x, y(1), y(2), nr-1
            xout=0.2d0
         else
            do
               if (x >= xout) then
                  y1 = interp_y(1, xout, rwork, iwork, ierr)
                  if (ierr /= 0) exit
                  y2 = interp_y(2, xout, rwork, iwork, ierr)
                  write (6, 99) xout, y1, y2, nr-1
                  if (ierr /= 0) exit
                  xout=xout+0.2d0
                  cycle
               end if
               exit
            end do
         end if
         if (ierr /= 0) then
            write(*, *) 'problem with interp_y in vdpol_solout'
            irtrn = -1
         end if
         rpar(1) = xout
  99     format(1x, 'x =', f5.2, '    y =', 2e18.10, '    nstep =', i8)
  
      end subroutine vdpol_solout


      end module vdpol
      

      program sample_ode_solver
      use vdpol
      implicit none
      
      call solve_vdpol
      
      end program sample_ode_solver
