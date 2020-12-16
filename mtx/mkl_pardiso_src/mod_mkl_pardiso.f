! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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


      module mod_mkl_pardiso
      
      implicit none

      integer, parameter :: num_pardiso_ipar_decsol = 64*3
      integer, parameter :: num_pardiso_rpar_decsol = 0
      
      contains
      
      
      logical function use_mkl_pardiso()
         use_mkl_pardiso = .true.
      end function use_mkl_pardiso

      
      subroutine do_mkl_pardiso_work_sizes(n,nzmax,lrd,lid)
         integer, intent(in) :: n,nzmax
         integer, intent(out) :: lrd,lid
         lid = num_pardiso_ipar_decsol
         lrd = num_pardiso_rpar_decsol
      end subroutine do_mkl_pardiso_work_sizes

      
      subroutine do_mkl_pardiso(iop,n,nzmax,ia,ja,a,b,lrd,rpar_decsol,lid,ipar_decsol,ierr)
         integer, intent(in) :: iop, n, nzmax, lrd, lid
         integer, intent(inout) :: ia(n+1), ja(nzmax)
         double precision, intent(inout) :: a(nzmax)
         double precision, intent(inout) :: b(n)
         double precision, intent(inout), target :: rpar_decsol(lrd)
         integer, intent(inout), target :: ipar_decsol(lid)
         integer, intent(out) :: ierr
         call do_mkl_pardiso_nrhs(iop,1,n,nzmax,ia,ja,a,b,lrd,rpar_decsol,lid,ipar_decsol,ierr)
      end subroutine do_mkl_pardiso


      
      subroutine do_mkl_pardiso_nrhs(iop,nrhs,n,nzmax,ia,ja,a,b,lrd,rpar_decsol,lid,ipar_decsol,ierr)
         integer, intent(in) :: iop, nrhs, n, nzmax, lrd, lid
         integer, intent(inout) :: ia(n+1), ja(nzmax)
         double precision, intent(inout) :: a(nzmax)
         double precision, intent(inout) :: b(n)
         double precision, intent(inout), target :: rpar_decsol(lrd)
         integer, intent(inout), target :: ipar_decsol(lid)
         integer, intent(out) :: ierr
         
         integer, pointer :: pt(:) ! internal solver memory pointers
         integer, pointer :: iparm(:) ! parameters for pardiso
         integer :: mtype, phase, error, msglvl
         integer :: i, idum
         real(dp) :: ddum
         real(dp), pointer :: x(:)
         integer, parameter :: maxfct = 1, mnum = 1
         
         logical, parameter :: dbg = .false.
         
         if (lid < 64*3) then
            ierr = -1
            write(*,*) 'lid too small for mkl_pardiso_decsols', lid
            return
         end if
         iparm => ipar_decsol(1:64)
         pt => ipar_decsol(65:64*3)

         ierr = 0
         error = 0 ! initialize error flag
         !msglvl = 1 ! print statistical information
         msglvl = 0 ! no statistical information
         mtype = 11 ! real unsymmetric
         
         if (iop == 0) then ! factor
            !  initilize the internal solver memory pointer
            pt(:) = 0
            !  set up pardiso control parameters
            iparm(:) = 0
            iparm(1) = 1 ! no solver default
            iparm(2) = 2 ! fill-in reordering from metis
            iparm(3) = 4 ! use openmp    
            iparm(4) = 0 ! no iterative-direct algorithm
            iparm(5) = 0 ! no user fill-in reducing permutation
            iparm(6) = 0 ! =0 solution on the first n compoments of x
            iparm(7) = 0 ! not in use
            iparm(8) = 9 ! numbers of iterative refinement steps
            iparm(9) = 0 ! not in use
            iparm(10) = 13 ! perturbe the pivot elements with 1e-13
            iparm(11) = 1 ! use nonsymmetric permutation and scaling mps
            iparm(12) = 0 ! not in use
            iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
            iparm(14) = 0 ! output: number of perturbed pivots
            iparm(15) = 0 ! not in use
            iparm(16) = 0 ! not in use
            iparm(17) = 0 ! not in use
            iparm(18) = -1 ! output: number of nonzeros in the factor lu
            iparm(19) = -1 ! output: mflops for lu factorization
            iparm(20) = 0 ! output: numbers of cg iterations
            !  reordering and symbolic factorization, this step also allocates
            !  all memory that is necessary for the factorization
            phase = 11 ! only reordering and symbolic factorization
            if (dbg) write(*,*) 'call pardiso: phase', phase
            call pardiso(pt, maxfct, mnum, mtype, phase, 
     >            n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
            if (error .ne. 0) then
               write(*,*) 'the following error was detected: ', error
               call mesa_error(__FILE__,__LINE__)
            end if
            if (dbg) write(*,*) 'number of nonzeros in factors = ',iparm(18)
            if (dbg) write(*,*) 'fill fraction = ',dble(iparm(18))/dble(n)
         
            !  factorization.
            phase = 22 ! only factorization
            if (dbg) write(*,*) 'call pardiso: phase', phase
            call pardiso(pt, maxfct, mnum, mtype, phase, 
     >            n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
            if (error .ne. 0) then
               write(*,*) 'the following error was detected: ', error
               call mesa_error(__FILE__,__LINE__)
            endif
         
         else if (iop == 1) then ! solve
         
            !  back substitution and iterative refinement
            iparm(8) = 2 ! max numbers of iterative refinement steps
            phase = 33 ! only solve
            allocate(x(n))
            if (dbg) write(*,*) 'call pardiso: phase, n', phase, n, nrhs
            call pardiso(pt, maxfct, mnum, mtype, phase, 
     >            n, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error)
            do i=1,n 
               b(i) = x(i)
            end do
            deallocate(x)
            if (error .ne. 0) then
               write(*,*) 'the following error was detected: ', error
               call mesa_error(__FILE__,__LINE__)
            endif
         
         else if (iop == 2) then ! deallocate
            phase = -1 ! release internal memory
            if (dbg) write(*,*) 'call pardiso: phase', phase
            call pardiso(pt, maxfct, mnum, mtype, phase, 
     >            n, ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
            if (error .ne. 0) then
               write(*,*) 'the following error was detected: ', error
               call mesa_error(__FILE__,__LINE__)
            endif
         
         else 
         
            write(*,*) 'invalid iop for do_mkl_pardiso', iop
            call mesa_error(__FILE__,__LINE__)
            
         end if
         
         if (dbg) write(*,*) 'finished pardiso'
         
      end subroutine do_mkl_pardiso_nrhs


      end module mod_mkl_pardiso


