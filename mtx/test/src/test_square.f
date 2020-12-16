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


      module test_square
      
      use mtx_lib
      use mtx_def
      use utils_lib, only: mesa_error
      
      implicit none
      
      contains
      
      
      
      subroutine do_test_square
         call test_square1
         call test_square2
         call test_square_inv
      end subroutine do_test_square
      
      
      subroutine test_square_inv
      
      
         integer,parameter :: n=3, nrhs=1
         integer :: i,info,ipiv(n),icommon(n)
         double precision :: A1(n,n),B1(n),A2(n,n),B2(n),work(4*n),rcond
         double precision :: A1_init(n,n),A2_init(n,n),X(n),prod(n)
         
         include 'formats'
         
         write(*,*) 'test_square_inv'
         write(*,*)
      
         A1(1,1:n) = (/ 3.14d0,7.5d0, 0.00d0 /)
         A1(2,1:n) = (/ 4.1d0,3.2d0,0.3d0 /)
         A1(3,1:n) = (/ 0.00d0,1d0,4.1d0 /)
         A1_init = A1
      
         A2(1,1:n) = (/ 0d0,3.1d0,0d0 /)
         A2(2,1:n) = (/ 4.7d0,6.2d0,0d0 /)
         A2(3,1:n) = (/ 3.2d0,0d0,0.31d0 /)
         A2_init = A2
         
         B1(1:n) = (/ 1.0d0,2.0d0,3.0d0 /)
         B2(1:n) = (/ 1.1d0,2.1d0,3.1d0 /)
         
         
               
         call DGETRF(n,n,A1,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         X = B1
         ! solve A1*X = B1
         call DGETRS('N',n,nrhs,A1,n,ipiv,X,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)
         
         write(*,1) 'B1', B1(1:n)
         ! prod = A1_init*X; should get prod == B1
         call dgemv('N',n,n,1d0,A1_init,n,X,1,0d0,prod,1)
         write(*,1) 'A1_init*X', prod(1:n)
         ! prod = A1*X; should get prod == B1

         call DGETRF(n,n,A2,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         X = B2
         ! solve A2*X = B2
         call DGETRS('N',n,nrhs,A2,n,ipiv,X,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B2', B2(1:n)
         ! prod = A2_init*X; should get prod == B2
         call dgemv('N',n,n,1d0,A2_init,n,X,1,0d0,prod,1)
         write(*,1) 'A2_init*X', prod(1:n)
         ! prod = A2*X; should get prod == B2
     
      end subroutine test_square_inv
      
      
      subroutine test_square2
      
         integer,parameter :: n=3,nrhs=1
         integer :: i,info,ipiv(n),icommon(n)
         double precision :: A1(n,n),B1(n,nrhs),A2(n,n),B2(n,nrhs),work(4*n),rcond
         
         include 'formats'
         
         write(*,*) 'test_square2'
         write(*,*)
      
         A1(1,1:n) = (/ 3.14d0,7.5d0, 0.00d0 /)
         A1(2,1:n) = (/ 4.1d0,3.2d0,0.3d0 /)
         A1(3,1:n) = (/ 0.00d0,1d0,4.1d0 /)
      
         A2(1,1:n) = (/ 4.7d0,6.2d0,0d0 /)
         A2(2,1:n) = (/ 3.2d0,0d0,0.31d0 /)
         A2(3,1:n) = (/ 0d0,3.1d0,0d0 /)
         
         B1(1:n,1) = (/ 1.0d0,2.0d0,3.0d0 /)
         B2(1:n,1) = (/ 1.1d0,2.1d0,3.1d0 /)
               
         call DGETRF(n,n,A1,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         call DGETRS('N',n,nrhs,A1,n,ipiv,B1,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B1', B1(1:n,1)
               
         call DGETRF(n,n,A2,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         call DGETRS('N',n,nrhs,A2,n,ipiv,B2,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B2', B2(1:n,1)
         write(*,*)
      
      end subroutine test_square2
      
      
      subroutine test_square1
      
         integer,parameter :: n=4,nrhs=1
         integer :: i,info,ipiv(n),iwork(n)
         double precision :: A(n,n),B(n,nrhs),A2(n,n),work(4*n),rcond
         
         include 'formats'
      
         A(1,1:n) = (/ 1.80d0,  2.88d0,  2.05d0, 0.00d0 /)
         A(2,1:n) = (/ 5.25d0, -2.95d0, -0.95d0, -3.80d0 /)
         A(3,1:n) = (/ 0.00d0, 0.00d0, -2.90d0, -1.04d0 /)
         A(4,1:n) = (/-1.11d0, 0.00d0, -0.59d0,  0.80d0 /)
         B(1:n,1) = (/ 4.35d0,5.05d0,3.04d0,-2.05d0 /)
      
         A2 = A
         
         write(*,*) ' test_square1'
         write(*,*)
         
         call DGETRF(n,n,A,n,ipiv,info)
         if (info /= 0) then
            write(*,*) 'singular matrix?', info
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call DGETRS('N',n,nrhs,A,n,ipiv,B,n,info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,1) 'B', B(1:n,1)
         write(*,*)
      
      end subroutine test_square1


      end module test_square
