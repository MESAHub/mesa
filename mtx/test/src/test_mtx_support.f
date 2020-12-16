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


      module test_mtx_support
      
      use mtx_lib
      use mtx_def
      use utils_lib, only: mesa_error 
     
      implicit none


      
      contains
      
      
      subroutine test_format_conversion
         use mtx_def
         integer,parameter :: n=6
         integer,parameter :: nzmax = n*n,nrow = n,ncol = n,ndns = n,ndim = n
         integer,parameter :: iwk = nzmax,im = 10
         
         real(dp) :: a(ndim,n),a2(ndim,n),values(nzmax)
         integer,parameter :: ml = 1,mu = 2, ldbb = 2*ml+mu+1
         real(dp) :: b(ndim,n),b2(ndim,n),bb(ldbb,n),bb2(ldbb,n)
         integer :: ierr,nz,iptr(n+1),jind(nzmax),i,j,k,kk,hint
         
         write(*,*) 'test_format_conversion'
         
         a(1,1:n) = (/ 10d0, 0d0, 0d0, 0d0,  0d0, 0d0 /)
         a(2,1:n) = (/  0d0,12d0,-3d0,-1d0,  0d0, 0d0 /)
         a(3,1:n) = (/  0d0, 0d0,15d0, 0d0,  0d0, 0d0 /)
         a(4,1:n) = (/ -2d0, 0d0, 0d0,10d0, -1d0, 0d0 /)
         a(5,1:n) = (/ -1d0, 0d0, 0d0,-5d0,  1d0,-1d0 /)
         a(6,1:n) = (/ -1d0,-2d0, 0d0, 0d0,  0d0, 6d0 /)
         
         b(1,1:n) = (/ 10d0, 0d0, 0d0, 0d0,  0d0, 0d0 /)
         b(2,1:n) = (/ -2d0,12d0,-3d0,-1d0,  0d0, 0d0 /)
         b(3,1:n) = (/  0d0, 1d0,15d0, 0d0,  0d0, 0d0 /)
         b(4,1:n) = (/  0d0, 0d0, 0d0,10d0, -1d0, 0d0 /)
         b(5,1:n) = (/  0d0, 0d0, 0d0,-5d0,  1d0,-1d0 /)
         b(6,1:n) = (/  0d0, 0d0, 0d0, 0d0,  0d0, 6d0 /)
         
         ierr = 0
         
         write(*,*) 'dense_to_row_sparse'
         call dense_to_row_sparse(n,ndim,a,nzmax,nz,iptr,jind,values,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         a2 = -1
         write(*,*) 'find_loc_in_row_sparse'
         do i=1,n
            hint = 0
            do k=iptr(i),iptr(i+1)-1
               j = jind(k)
               call find_loc_in_sparse(compressed_row_sparse,n,nzmax,iptr,jind,i,j,hint,kk,ierr)
               if (kk /= k .or. ierr /= 0) then
                  write(*,*) 'failure in find_loc_in_row_sparse', i, j, k, kk
                  call mesa_error(__FILE__,__LINE__)
               end if
               hint = k
            end do
         end do
              
         write(*,*) 'dense_to_band'
         call dense_to_band(n,ndim,b,ml,mu,bb,ldbb,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,*) 'band_to_dense'
         a2 = -1
         call band_to_dense(n,ml,mu,bb,ldbb,ndim,a2,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         if (any(b /= a2)) call mesa_error(__FILE__,__LINE__)
         
         write(*,*) 'band_to_column_sparse'
         call band_to_column_sparse(n,ml,mu,bb,ldbb,nzmax,nz,iptr,jind,values,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         write(*,*) 'column_sparse_to_band'
         bb2 = -1
         call column_sparse_to_band(n,ml,mu,bb2,ldbb,nz,iptr,jind,values,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         if (any(bb /= bb2)) call mesa_error(__FILE__,__LINE__)
         
         write(*,*) 'band_to_row_sparse'
         call band_to_row_sparse(n,ml,mu,bb,ldbb,nzmax,nz,iptr,jind,values,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,*) 'okay'
         write(*,*)
      
      end subroutine test_format_conversion
      
      
      subroutine test_quad_tridiag
         integer, parameter :: n = 5
         real(16), dimension(n) :: DL, D, DU, DU2, B
         integer, dimension(n) :: ip
         integer :: ierr, i
         
         write(*,*) 'test_quad_tridiag'

         DL = -2  ! subdiagonal
         D =   3  ! diagonal
         DU = -1  ! superdiagonal

         do i=1,n
            b(i) = i-1
         end do
         
         ierr = 0
         ! factor
         call qgttrf(n, DL, D, DU, DU2, ip, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in factoring'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         ierr = 0
         ! solve
         call qgttrs( 'N', n, 1, DL, D, DU, DU2, ip, B, n, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in solving'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         do i=1,n
            write(*,*) i, b(i)
         end do
         write(*,*)
         write(*,*)
         
      end subroutine test_quad_tridiag
      


   end module test_mtx_support
