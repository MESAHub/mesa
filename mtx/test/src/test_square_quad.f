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


      module test_square_quad
      
      use mtx_lib
      use mtx_def
      use utils_lib, only: mesa_error
      
      implicit none
      
      contains
      
      
      
      subroutine do_test_square_quad
         call test_square_quad1
         call test_square_quad2
         call test_square_quad_inv
      end subroutine do_test_square_quad
      
      
      subroutine test_square_quad_inv
      
         integer,parameter :: n=3, nrhs=1
         integer :: i,info,ipiv(n),icommon(n)
         real(16) :: A1(n,n),B1(n),A2(n,n),B2(n),work(4*n),rcond
         real(16) :: A1_init(n,n),A2_init(n,n),X(n),prod(n)
         
         include 'formats'
         
         write(*,*) 'test_square_quad_inv'
         write(*,*)
      
         A1(1,1:n) = (/ 3.14_16,7.5_16, 0.00_16 /)
         A1(2,1:n) = (/ 4.1_16,3.2_16,0.3_16 /)
         A1(3,1:n) = (/ 0.00_16,1.0_16,4.1_16 /)
         A1_init = A1
      
         A2(1,1:n) = (/ 0.0_16,3.1_16,0.0_16 /)
         A2(2,1:n) = (/ 4.7_16,6.2_16,0.0_16 /)
         A2(3,1:n) = (/ 3.2_16,0.0_16,0.31_16 /)
         A2_init = A2
         
         B1(1:n) = (/ 1.0_16,2.0_16,3.0_16 /)
         B2(1:n) = (/ 1.1_16,2.1_16,3.1_16 /)
         
         info = 0
         call QGETRF(n,n,A1,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         X = B1
         ! solve A1*X = B1
         call QGETRS('N',n,nrhs,A1,n,ipiv,X,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)
         
         write(*,1) 'B1', B1(1:n)
         ! prod = A1_init*X; should get prod == B1
         
         write(*,*) 'ipiv1', ipiv(1:n)
         write(*,*)
               
         call QGETRF(n,n,A2,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         X = B2
         ! solve A2*X = B2
         call QGETRS('N',n,nrhs,A2,n,ipiv,X,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B2', B2(1:n)
      
      end subroutine test_square_quad_inv
      
      
      subroutine test_square_quad2
      
         integer,parameter :: n=3,nrhs=1
         integer :: i,info,ipiv(n),icommon(n)
         real(16) :: A1(n,n),B1(n,nrhs),A2(n,n),B2(n,nrhs),work(4*n),rcond
         
         include 'formats'
         
         write(*,*) 'test test_square_quad2'
         write(*,*)
      
         A1(1,1:n) = (/ 3.14_16,7.5_16, 0.0_16 /)
         A1(2,1:n) = (/ 4.1_16,3.2_16,0.3_16 /)
         A1(3,1:n) = (/ 0.00_16,1.0_16,4.1_16 /)
      
         A2(1,1:n) = (/ 4.7_16,6.2_16,0.0_16 /)
         A2(2,1:n) = (/ 3.2_16,0.0_16,0.31_16 /)
         A2(3,1:n) = (/ 0.0_16,3.1_16,0.0_16 /)
         
         B1(1:n,1) = (/ 1.0_16,2.0_16,3.0_16 /)
         B2(1:n,1) = (/ 1.1_16,2.1_16,3.1_16 /)
               
         info = 0
         call QGETRF(n,n,A1,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         call QGETRS('N',n,nrhs,A1,n,ipiv,B1,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B1', B1(1:n,1)
               
         call QGETRF(n,n,A2,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix'
            call mesa_error(__FILE__,__LINE__)
         end if
         call QGETRS('N',n,nrhs,A2,n,ipiv,B2,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B2', B2(1:n,1)
         write(*,*)
      
      end subroutine test_square_quad2
      
      
      subroutine test_square_quad1
      
         integer,parameter :: n=4,nrhs=1
         integer :: i,info,ipiv(n),iwork(n)
         real(16) :: A(n,n),B(n,nrhs),A2(n,n),work(4*n),rcond
         
         include 'formats'
      
         A(1,1:n) = (/ 1.80_16,  2.88_16,  2.05_16,  0.00_16 /)
         A(2,1:n) = (/ 5.25_16, -2.95_16, -0.95_16, -3.80_16 /)
         A(3,1:n) = (/ 0.00_16,  0.00_16, -2.90_16, -1.04_16 /)
         A(4,1:n) = (/-1.11_16,  0.00_16, -0.59_16,  0.80_16 /)
         B(1:n,1) = (/ 4.35_16,  5.05_16,  3.04_16, -2.05_16 /)
      
         A2 = A
         
         write(*,*) ' test_square_quad1'
         write(*,*)
         
         info = 0
         call QGETRF(n,n,A,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call QGETRS('N',n,nrhs,A,n,ipiv,B,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B', B(1:n,1)
         write(*,*)
      
      end subroutine test_square_quad1


      end module test_square_quad
