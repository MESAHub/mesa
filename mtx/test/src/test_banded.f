! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not,write to the Free Software
!   Foundation,Inc.,59 Temple Place,Suite 330,Boston,MA 02111-1307 USA
!
! ***********************************************************************


      module test_banded
      
      use mtx_lib
      use mtx_def
      use utils_lib, only: mesa_error      

      implicit none
      
      contains
      
      
      subroutine do_test_banded
         call test_diffusion      
      end subroutine do_test_banded
      
      
      subroutine test_med
         integer,parameter :: n = 400,iounit = 33,KU = 2,KL = 2,nzmax = 500
         integer,parameter :: ldAB = 2*KL+KU+1
         real(dp) :: AB(ldAB,n),work(3*n),rcond, values(nzmax), b(n), rhs(n), ref(n)
         integer :: nz,ldb,nrhs,ierr,nn,i,j,ipiv(n),iwork(n),ios,info,colptr(n+1),rowind(nzmax),norm
         
         include 'formats'
         
         open(unit=iounit,file='med',action='read',status='old',iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open file for test_banded'
            call mesa_error(__FILE__,__LINE__)
         end if
         ! read compressed column format
         read(iounit,*) nn, nz
         if (nn /= n) stop 'test_med'
         do i=1,n+1
            read(iounit,*) colptr(i)
         end do
         do i=1,nz
            read(iounit,*) rowind(i)
         end do
         do i=1,nz
            read(iounit,'(1pe26.16)') values(i)
         end do
         do i=1,n
            read(iounit,'(1pe26.16)') rhs(i)
         end do
         do i=1,n
            read(iounit,'(1pe26.16)') ref(i)
         end do
         close(iounit)
         call column_sparse_to_band(n,KL,KU,ab,ldab,nz,colptr,rowind,values,ierr)
         if (ierr /=0) stop 'test_med column_sparse_to_band'
         
         call mtx_rcond_banded('N', n, n, KU, KL, AB, ldAB, ipiv, rcond, work, iwork, info)
         if (ierr /=0) stop 'test_med mtx_rcond_banded'
         write(*,*) 'mtx_rcond_banded rcond', rcond
         
         b = rhs
         ldb = n
         nrhs = 1
         CALL DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
         write(*,*) 'info', info
         if (info /= 0) stop 'test_med'
         
         do i=1,n
            write(*,2) 'b ref', i, b(i), ref(i)
         end do
         

      end subroutine test_med


      subroutine test_diffusion
         integer, parameter :: n = 10
         real(dp), dimension(n) :: x, sig, dq, bb
         real(dp) :: mstar, dt, init_xtotal, banded_xtotal
         integer :: lrd, lid
         
         include 'formats.dek'
         
         write(*,*) 'test_diffusion'
         
         mstar =    9.5216384650402128D+33
         dt =    8.2230155734773409D+08

         x( 1  )=    3.3853553765573918D-01
         x( 2  )=    3.3853553765572147D-01
         x( 3  )=    4.5462481641995478D-01
         x( 4  )=    4.9630359938855395D-01
         x( 5  )=    4.9630359938855390D-01
         x( 6  )=    4.9630359938855384D-01
         x( 7  )=    4.9630359938855390D-01
         x( 8  )=    5.8265967459668622D-01
         x( 9  )=    6.1370782322384976D-01
         x(10  )=    6.1370782322384976D-01

         sig( 1  )=    0.0000000000000000D+00
         sig( 2  )=    7.7006096541650005D+38
         sig( 3  )=    1.4638096218318065D+39
         sig( 4  )=    1.2767160238768072D+39
         sig( 5  )=    1.1314657565143642D+39
         sig( 6  )=    1.0969967768680467D+39
         sig( 7  )=    1.0969967748843832D+39
         sig( 8  )=    9.8367595476496157D+38
         sig( 9  )=    7.5703431242949176D+38
         sig(10  )=    3.2185674473372967D+38

         dq( 1  )=    3.2180704472120336D-08
         dq( 2  )=    1.7819062402685692D-08
         dq( 3  )=    1.5396355905853213D-08
         dq( 4  )=    1.3290183455525239D-08
         dq( 5  )=    1.5048420108580579D-08
         dq( 6  )=    1.6727579700091661D-08
         dq( 7  )=    1.8662871023197620D-08
         dq( 8  )=    2.1239878600194412D-08
         dq( 9  )=    2.2654435507986098D-08
         dq(10  )=    2.3822242205812923D-08
         
         init_xtotal = mstar*dot_product(x(:), dq(:))
         
         call lapack_work_sizes(2*n,lrd,lid)
         call do_banded(lrd,lid)
         banded_xtotal = mstar*dot_product(bb(:), dq(:))
         
         write(*,1) 'init_xtotal', init_xtotal
         write(*,1) 'banded_xtotal', banded_xtotal
         write(*,*)
         
         contains
         
         
         subroutine do_banded(lrd,lid)
            integer, intent(in) :: lrd, lid
            integer, parameter :: &
               nvar=2, nz=n, neq=nvar*nz, i_x=1, i_dx=2, equx=1, equdx=2, &
               mu = 2*nvar-1, ml = 2*nvar-1, lda = 2*ml+mu+1, idiag = ml+mu+1
            real(dp), dimension(nvar,nvar,nz) :: em1, e00, ep1
            real(dp), target :: a1_ary(lda*neq), dx_ary(n), b1_ary(nvar*nz)
            integer :: ierr, i, k
            integer, target :: ip_ary(neq)
            real(dp), target :: rpar_decsol_ary(lrd)
            integer, target :: ipar_decsol_ary(lid)
            real(dp), pointer :: rpar_decsol(:), a1(:), dx(:), b1(:)
            integer, pointer :: ipar_decsol(:)
            
            real(dp), pointer :: a(:,:), b(:,:)
            integer, pointer :: ip(:)

            include 'formats.dek'
            rpar_decsol => rpar_decsol_ary
            ipar_decsol => ipar_decsol_ary
            a1(1:lda*neq) => a1_ary
            a(1:lda,1:neq) => a1(1:lda*neq)
            b1(1:nvar*nz) => b1_ary
            b(1:nvar,1:nz) => b1(1:nvar*nz)
            dx(1:n) => dx_ary
            ip(1:neq) => ip_ary
            
            dx(1) = 0
            do k=2,n
               dx(k) = x(k-1)-x(k)
            end do
            
            ! store partials
            em1 = 0; e00 = 0; ep1 = 0
            do k=1, n
               ! dx(k) - x(k-1) + x(k) = 0
               e00(equdx,i_dx,k) = 1
               if (k > 1) then
                  e00(equdx,i_x,k) = 1
                  em1(equdx,i_x,k) = -1
               end if
               ! x(k) - (sig(k)*dx(k) - sig(k+1)*dx(k+1))*dt/(mstar*dq(k)) = x0(k)
               e00(equx,i_x,k) = 1
               if (k > 1) e00(equx,i_dx,k) = -sig(k)*dt/(mstar*dq(k))
               if (k < n) ep1(equx,i_dx,k) = sig(k+1)*dt/(mstar*dq(k))
            end do
            
            ! store matrix
            call copy_all_to_3point_jacobian(nvar, nz, ldA, A, idiag, nz, em1, e00, ep1)
            
            ! factor
            call lapack_decsol(0,neq,ldA,A1,ml,mu,b1,ip,lrd,rpar_decsol,lid,ipar_decsol,ierr)
            if (ierr /= 0) then
               write(*,*) 'lapack_decsol 0', ierr
               call mesa_error(__FILE__,__LINE__)
            end if
            
            ! set rhs of equation
            b(equdx,:) = 0
            b(equx,:) = x(:)
            
            ! solve
            call lapack_decsol(1,neq,ldA,A1,ml,mu,b1,ip,lrd,rpar_decsol,lid,ipar_decsol,ierr)
            if (ierr /= 0) then
               write(*,*) 'lapack_decsol 1', ierr
               call mesa_error(__FILE__,__LINE__)
            end if
            
            bb(:) = b(equx,:)
            
            do k=1,n
               !write(*,2) 'banded: init x, soln x(k), dx(k)', k, x(k), b(equx,k), b(equdx,k)
               write(*,2) 'banded: init x, soln x(k)', k, x(k), b(equx,k)
            end do
            write(*,*)
         
         end subroutine do_banded
         
      
      end subroutine test_diffusion

      end module test_banded
