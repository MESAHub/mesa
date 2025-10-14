   ! ***********************************************************************
!
!   copyright (c) 2012  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module my_lapack95

      use const_def, only: dp
      use utils_lib, only: mesa_error
      implicit none
      integer, parameter :: fltp = dp

      contains

      subroutine my_gemv(m,n,a,lda,x,y)  ! y = y - a*x
         integer :: lda,m,n
         real(fltp) :: a(:,:)  ! (lda,*)
         real(fltp) :: x(:), y(:)
         real(fltp) :: tmp
         real(fltp), parameter :: zero=0
         ! trans = 'n'
         ! alpha = -1
         ! beta = 1
         ! incx = 1
         ! incy = 1
         integer :: j, i
         do j = 1,n
            tmp = x(j)
            if (tmp/=zero) then
               do i = 1,m
                  y(i) = y(i) - tmp*a(i,j)
               end do
            end if
         end do
      end subroutine my_gemv


      subroutine my_gemv_mv(m,n,a,x,b,z,y)  ! y = y - a*x - b*z
         integer :: m, n
         real(fltp) :: a(:,:), b(:,:)
         real(fltp) :: x(:), z(:), y(:)
         real(fltp) :: tmp_x, tmp_z
         real(fltp), parameter :: zero=0
         integer :: j, i
         do j = 1,n
            tmp_x = x(j)
            tmp_z = z(j)
            if (tmp_x/=zero) then
               if (tmp_z /= zero) then
                  do i = 1,m
                     y(i) = y(i) - tmp_x*a(i,j) - tmp_z*b(i,j)
                  end do
               else
                  do i = 1,m
                     y(i) = y(i) - tmp_x*a(i,j)
                  end do
               end if
            else if (tmp_z /= zero) then
               do i = 1,m
                  y(i) = y(i) - tmp_z*b(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_mv


      subroutine my_gemv_p1(m,n,a,lda,x,y)  ! y = y + a*x
         integer :: lda,m,n
         real(fltp) :: a(:,:)  ! (lda,*)
         real(fltp) :: x(:), y(:)
         real(fltp) :: tmp
         real(fltp), parameter :: zero=0
         ! trans = 'n'
         ! alpha = -1
         ! beta = 1
         ! incx = 1
         ! incy = 1
         integer :: j, i
         do j = 1,n
            tmp = x(j)
            if (tmp/=zero) then
               do i = 1,m
                  y(i) = y(i) + tmp*a(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_p1


      subroutine my_gemv_p_mv(m,n,a,x,b,z,y)  ! y = y + a*x + b*z
         integer :: m, n
         real(fltp) :: a(:,:), b(:,:)
         real(fltp) :: x(:), z(:), y(:)
         real(fltp) :: tmp_x, tmp_z
         real(fltp), parameter :: zero=0
         integer :: j, i
         do j = 1,n
            tmp_x = x(j)
            tmp_z = z(j)
            if (tmp_x/=zero) then
               if (tmp_z /= zero) then
                  do i = 1,m
                     y(i) = y(i) + tmp_x*a(i,j) + tmp_z*b(i,j)
                  end do
               else
                  do i = 1,m
                     y(i) = y(i) + tmp_x*a(i,j)
                  end do
               end if
            else if (tmp_z /= zero) then
               do i = 1,m
                  y(i) = y(i) + tmp_z*b(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_p_mv


      subroutine my_gemm(m,n,k,a,lda,b,ldb,c,ldc)  ! c := c - a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(fltp), dimension(:,:) :: a, b, c  ! a(lda,*),b(ldb,*),c(ldc,*)
         real(fltp) :: tmp
         real(fltp), parameter :: zero=0
         integer :: j, i, l
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = -1
         ! beta = 1
         ! assumes other args are valid
         do j = 1,n
            do l = 1,k
               tmp = b(l,j)
               if (tmp /= zero) then
                  do i = 1,m
                     c(i,j) = c(i,j) - tmp*a(i,l)
                  end do
               end if
            end do
         end do
      end subroutine my_gemm


      subroutine my_gemm_p1(m,n,k,a,lda,b,ldb,c,ldc)  ! c := c + a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(fltp), dimension(:,:) :: a, b, c  ! a(lda,*),b(ldb,*),c(ldc,*)
         real(fltp) :: tmp
         real(fltp), parameter :: zero=0
         integer :: j, i, l
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = 1
         ! beta = 1
         ! assumes other args are valid
         do j = 1,n
            do l = 1,k
               tmp = b(l,j)
               if (tmp /= zero) then
                  do i = 1,m
                     c(i,j) = c(i,j) + tmp*a(i,l)
                  end do
               end if
            end do
         end do
      end subroutine my_gemm_p1


      subroutine my_gemm_plus_mm(m,n,k,a,b,d,e,c)  ! c := c + a*b + d*e
         integer, intent(in) :: k,m,n
         real(fltp), dimension(:,:) :: a, b, c, d, e
         real(fltp) :: tmp_b, tmp_e
         real(fltp), parameter :: zero=0
         integer :: j, i, l
         do j = 1,n
            do l = 1,k
               tmp_b = b(l,j)
               tmp_e = e(l,j)
               if (tmp_b /= zero) then
                  if (tmp_e /= zero) then
                     do i = 1,m
                        c(i,j) = c(i,j) + tmp_b*a(i,l) + tmp_e*d(i,l)
                     end do
                  else
                     do i = 1,m
                        c(i,j) = c(i,j) + tmp_b*a(i,l)
                     end do
                  end if
               else if (tmp_e /= zero) then
                  do i = 1,m
                     c(i,j) = c(i,j) + tmp_e*d(i,l)
                  end do
               end if
            end do
         end do
      end subroutine my_gemm_plus_mm


      subroutine my_gemm0(m,n,k,a,lda,b,ldb,c,ldc)
         ! c := -a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(fltp), dimension(:,:) :: a, b, c  ! a(lda,*),b(ldb,*),c(ldc,*)
         integer :: j, i
         real(fltp), parameter :: zero=0
         include 'formats'
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = -1
         ! beta = 0
         ! assumes other args are valid
         do j=1,n
            do i=1,m
               c(i,j) = zero
            end do
         end do
         call my_gemm(m,n,k,a,lda,b,ldb,c,ldc)
      end subroutine my_gemm0


      subroutine my_gemm0_p1(m,n,k,a,lda,b,ldb,c,ldc)
         ! c := -a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(fltp), dimension(:,:) :: a, b, c  ! a(lda,*),b(ldb,*),c(ldc,*)
         integer :: j, i
         real(fltp), parameter :: zero=0
         include 'formats'
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = -1
         ! beta = 0
         ! assumes other args are valid
         do j=1,n
            do i=1,m
               c(i,j) = zero
            end do
         end do
         call my_gemm_p1(m,n,k,a,lda,b,ldb,c,ldc)
      end subroutine my_gemm0_p1


      subroutine my_getf2(m, a, lda, ipiv, sfmin, info)
         integer :: info, lda, m
         integer :: ipiv(:)
         real(fltp) :: a(:,:)
         real(fltp) :: sfmin
         real(fltp), parameter :: one=1, zero=0
         integer :: i, j, jp, ii, jj, n, mm
         real(fltp) :: tmp, da
         if (m == 4) then
            call my_getf2_4_by_4(a, lda, ipiv, sfmin, info)
            return
         else if (m == 5) then
            call my_getf2_5_by_5(a, lda, ipiv, sfmin, info)
            return
         end if
         do j = 1, m
            info = 0
            jp = j - 1 + maxloc(abs(a(j:lda,j)),dim=1)
            ipiv( j ) = jp
            if( a( jp, j )/=zero ) then
               if( jp/=j ) then  ! swap a(j,:) and a(jp,:)
                  do i=1,m
                     tmp = a(j,i)
                     a(j,i) = a(jp,i)
                     a(jp,i) = tmp
                  end do
               end if
               if( j<m ) then
                  if( abs(a( j, j )) >= sfmin ) then
                     da = one / a( j, j )
                     n = m-j
                     mm = mod(n,5)
                     if (mm /= 0) then
                        do i = 1,mm
                           a(j+i,j) = da*a(j+i,j)
                        end do
                     end if
                     if (n >= 5) then
                        do i = mm + 1,n,5
                           a(j+i,j) = da*a(j+i,j)
                           a(j+i+1,j) = da*a(j+i+1,j)
                           a(j+i+2,j) = da*a(j+i+2,j)
                           a(j+i+3,j) = da*a(j+i+3,j)
                           a(j+i+4,j) = da*a(j+i+4,j)
                        end do
                     end if
                  else  ! no scale
                    do i = 1, m-j
                       a( j+i, j ) = a( j+i, j ) / a( j, j )
                    end do
                  end if
               end if
            else if( info==0 ) then
               info = j
            end if
            if( j<m ) then
               !call dger( m-j, m-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda, a( j+1, j+1 ), lda )
               do jj = j+1, m
                  do ii = j+1, m
                     a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                  end do
               end do
            end if
         end do
      end subroutine my_getf2


      subroutine my_getf2_4_by_4(a, lda, ipiv, sfmin, info)
         integer :: info, lda  !  m=4
         integer :: ipiv(:)
         real(fltp) :: a(:,:)
         real(fltp) :: sfmin
         real(fltp), parameter :: one=1, zero=0
         integer :: jp
         real(fltp) :: tmp, da
         info = 0

         jp = maxloc(abs(a(1:lda,1)),dim=1)
         ipiv( 1 ) = jp
         if( a( jp, 1 )/=zero ) then
            if( jp/=1 ) then  ! swap a(1,:) and a(jp,:)
               tmp = a(1,1)
               a(1,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(1,2)
               a(1,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(1,3)
               a(1,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(1,4)
               a(1,4) = a(jp,4)
               a(jp,4) = tmp
            end if
            if( abs(a( 1, 1 )) >= sfmin ) then
               da = one / a( 1, 1 )
               a(2,1) = da*a(2,1)
               a(3,1) = da*a(3,1)
               a(4,1) = da*a(4,1)
            else  ! no scale
               a( 2, 1 ) = a( 2, 1 ) / a( 1, 1 )
               a( 3, 1 ) = a( 3, 1 ) / a( 1, 1 )
               a( 4, 1 ) = a( 4, 1 ) / a( 1, 1 )
            end if
         else if( info==0 ) then
            info = 1
         end if
         a(2,2) = a(2,2) - a(2,1)*a(1,2)
         a(3,2) = a(3,2) - a(3,1)*a(1,2)
         a(4,2) = a(4,2) - a(4,1)*a(1,2)
         a(2,3) = a(2,3) - a(2,1)*a(1,3)
         a(3,3) = a(3,3) - a(3,1)*a(1,3)
         a(4,3) = a(4,3) - a(4,1)*a(1,3)
         a(2,4) = a(2,4) - a(2,1)*a(1,4)
         a(3,4) = a(3,4) - a(3,1)*a(1,4)
         a(4,4) = a(4,4) - a(4,1)*a(1,4)

         jp = 1 + maxloc(abs(a(2:lda,2)),dim=1)
         ipiv( 2 ) = jp
         if( a( jp, 2 )/=zero ) then
            if( jp/=2 ) then  ! swap a(2,:) and a(jp,:)
               tmp = a(2,1)
               a(2,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(2,2)
               a(2,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(2,3)
               a(2,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(2,4)
               a(2,4) = a(jp,4)
               a(jp,4) = tmp
            end if
            if( abs(a( 2, 2 )) >= sfmin ) then
               da = one / a( 2, 2 )
               a(3,2) = da*a(3,2)
               a(4,2) = da*a(4,2)
            else  ! no scale
               a( 3, 2 ) = a( 3, 2 ) / a( 2, 2 )
               a( 4, 2 ) = a( 4, 2 ) / a( 2, 2 )
            end if
         else if( info==0 ) then
            info = 2
         end if
         a(3,3) = a(3,3) - a(3,2)*a(2,3)
         a(4,3) = a(4,3) - a(4,2)*a(2,3)
         a(3,4) = a(3,4) - a(3,2)*a(2,4)
         a(4,4) = a(4,4) - a(4,2)*a(2,4)

         jp = 2 + maxloc(abs(a(3:lda,3)),dim=1)
         ipiv( 3 ) = jp
         if( a( jp, 3 )/=zero ) then
            if( jp/=3 ) then  ! swap a(3,:) and a(jp,:)
               tmp = a(3,1)
               a(3,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(3,2)
               a(3,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(3,3)
               a(3,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(3,4)
               a(3,4) = a(jp,4)
               a(jp,4) = tmp
            end if
            if( abs(a( 3, 3 )) >= sfmin ) then
               da = one / a( 3, 3 )
               a(4,3) = da*a(4,3)
            else  ! no scale
               a( 4, 3 ) = a( 4, 3 ) / a( 3, 3 )
            end if
         else if( info==0 ) then
            info = 3
         end if
         a(4,4) = a(4,4) - a(4,3)*a(3,4)

         jp = 3 + maxloc(abs(a(4:lda,4)),dim=1)
         ipiv( 4 ) = jp
         if( a( jp, 4 )/=zero ) then
            if( jp/=4 ) then  ! swap a(4,:) and a(jp,:)
               tmp = a(4,1)
               a(4,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(4,2)
               a(4,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(4,3)
               a(4,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(4,4)
               a(4,4) = a(jp,4)
               a(jp,4) = tmp
            end if
         else if( info==0 ) then
            info = 4
         end if

      end subroutine my_getf2_4_by_4


      subroutine my_getf2_5_by_5(a, lda, ipiv, sfmin, info)
         integer :: info, lda  !  m=5
         integer :: ipiv(:)
         real(fltp) :: a(:,:)
         real(fltp) :: sfmin
         real(fltp), parameter :: one=1, zero=0
         integer :: jp
         real(fltp) :: tmp, da
         info = 0
         jp = maxloc(abs(a(1:lda,1)),dim=1)
         ipiv( 1 ) = jp
         if( a( jp, 1 )/=zero ) then
            if( jp/=1 ) then  ! swap a(1,:) and a(jp,:)
               tmp = a(1,1)
               a(1,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(1,2)
               a(1,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(1,3)
               a(1,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(1,4)
               a(1,4) = a(jp,4)
               a(jp,4) = tmp
               tmp = a(1,5)
               a(1,5) = a(jp,5)
               a(jp,5) = tmp
            end if
            if( abs(a( 1, 1 )) >= sfmin ) then
               da = one / a( 1, 1 )
               a(2,1) = da*a(2,1)
               a(3,1) = da*a(3,1)
               a(4,1) = da*a(4,1)
               a(5,1) = da*a(5,1)
            else  ! no scale
               a( 2, 1 ) = a( 2, 1 ) / a( 1, 1 )
               a( 3, 1 ) = a( 3, 1 ) / a( 1, 1 )
               a( 4, 1 ) = a( 4, 1 ) / a( 1, 1 )
               a( 5, 1 ) = a( 5, 1 ) / a( 1, 1 )
            end if
         else if( info==0 ) then
            info = 1
         end if
         a(2,2) = a(2,2) - a(2,1)*a(1,2)
         a(3,2) = a(3,2) - a(3,1)*a(1,2)
         a(4,2) = a(4,2) - a(4,1)*a(1,2)
         a(5,2) = a(5,2) - a(5,1)*a(1,2)
         a(2,3) = a(2,3) - a(2,1)*a(1,3)
         a(3,3) = a(3,3) - a(3,1)*a(1,3)
         a(4,3) = a(4,3) - a(4,1)*a(1,3)
         a(5,3) = a(5,3) - a(5,1)*a(1,3)
         a(2,4) = a(2,4) - a(2,1)*a(1,4)
         a(3,4) = a(3,4) - a(3,1)*a(1,4)
         a(4,4) = a(4,4) - a(4,1)*a(1,4)
         a(5,4) = a(5,4) - a(5,1)*a(1,4)
         a(2,5) = a(2,5) - a(2,1)*a(1,5)
         a(3,5) = a(3,5) - a(3,1)*a(1,5)
         a(4,5) = a(4,5) - a(4,1)*a(1,5)
         a(5,5) = a(5,5) - a(5,1)*a(1,5)

         jp = 1 + maxloc(abs(a(2:lda,2)),dim=1)
         ipiv( 2 ) = jp
         if( a( jp, 2 )/=zero ) then
            if( jp/=2 ) then  ! swap a(2,:) and a(jp,:)
               tmp = a(2,1)
               a(2,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(2,2)
               a(2,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(2,3)
               a(2,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(2,4)
               a(2,4) = a(jp,4)
               a(jp,4) = tmp
               tmp = a(2,5)
               a(2,5) = a(jp,5)
               a(jp,5) = tmp
            end if
            if( abs(a( 2, 2 )) >= sfmin ) then
               da = one / a( 2, 2 )
               a(3,2) = da*a(3,2)
               a(4,2) = da*a(4,2)
               a(5,2) = da*a(5,2)
            else  ! no scale
               a( 3, 2 ) = a( 3, 2 ) / a( 2, 2 )
               a( 4, 2 ) = a( 4, 2 ) / a( 2, 2 )
               a( 5, 2 ) = a( 5, 2 ) / a( 2, 2 )
            end if
         else if( info==0 ) then
            info = 2
         end if
         a(3,3) = a(3,3) - a(3,2)*a(2,3)
         a(4,3) = a(4,3) - a(4,2)*a(2,3)
         a(5,3) = a(5,3) - a(5,2)*a(2,3)
         a(3,4) = a(3,4) - a(3,2)*a(2,4)
         a(4,4) = a(4,4) - a(4,2)*a(2,4)
         a(5,4) = a(5,4) - a(5,2)*a(2,4)
         a(3,5) = a(3,5) - a(3,2)*a(2,5)
         a(4,5) = a(4,5) - a(4,2)*a(2,5)
         a(5,5) = a(5,5) - a(5,2)*a(2,5)

         jp = 2 + maxloc(abs(a(3:lda,3)),dim=1)
         ipiv( 3 ) = jp
         if( a( jp, 3 )/=zero ) then
            if( jp/=3 ) then  ! swap a(3,:) and a(jp,:)
               tmp = a(3,1)
               a(3,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(3,2)
               a(3,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(3,3)
               a(3,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(3,4)
               a(3,4) = a(jp,4)
               a(jp,4) = tmp
               tmp = a(3,5)
               a(3,5) = a(jp,5)
               a(jp,5) = tmp
            end if
            if( abs(a( 3, 3 )) >= sfmin ) then
               da = one / a( 3, 3 )
               a(4,3) = da*a(4,3)
               a(5,3) = da*a(5,3)
            else  ! no scale
               a( 4, 3 ) = a( 4, 3 ) / a( 3, 3 )
               a( 4, 3 ) = a( 4, 3 ) / a( 3, 3 )
            end if
         else if( info==0 ) then
            info = 3
         end if
         a(4,4) = a(4,4) - a(4,3)*a(3,4)
         a(5,4) = a(5,4) - a(5,3)*a(3,4)
         a(4,5) = a(4,5) - a(4,3)*a(3,5)
         a(5,5) = a(5,5) - a(5,3)*a(3,5)

         jp = 3 + maxloc(abs(a(4:lda,4)),dim=1)
         ipiv( 4 ) = jp
         if( a( jp, 4 )/=zero ) then
            if( jp/=4 ) then  ! swap a(4,:) and a(jp,:)
               tmp = a(4,1)
               a(4,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(4,2)
               a(4,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(4,3)
               a(4,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(4,4)
               a(4,4) = a(jp,4)
               a(jp,4) = tmp
               tmp = a(4,5)
               a(4,5) = a(jp,5)
               a(jp,5) = tmp
            end if
            if( abs(a( 4, 4 )) >= sfmin ) then
               da = one / a( 4, 4 )
               a(5,4) = da*a(5,4)
            else  ! no scale
              a( 5, 4 ) = a( 5, 4 ) / a( 4, 4 )
            end if
         else if( info==0 ) then
            info = 4
         end if
         a(5,5) = a(5,5) - a(5,4)*a(4,5)

         jp = 4 + maxloc(abs(a(5:lda,5)),dim=1)
         ipiv( 5 ) = jp
         if( a( jp, 5 )/=zero ) then
            if( jp/=5 ) then  ! swap a(5,:) and a(jp,:)
               tmp = a(5,1)
               a(5,1) = a(jp,1)
               a(jp,1) = tmp
               tmp = a(5,2)
               a(5,2) = a(jp,2)
               a(jp,2) = tmp
               tmp = a(5,3)
               a(5,3) = a(jp,3)
               a(jp,3) = tmp
               tmp = a(5,4)
               a(5,4) = a(jp,4)
               a(jp,4) = tmp
               tmp = a(5,5)
               a(5,5) = a(jp,5)
               a(jp,5) = tmp
            end if
         else if( info==0 ) then
            info = 5
         end if

      end subroutine my_getf2_5_by_5


      subroutine my_laswp( n,   a, lda,  k1, k2, ipiv,  incx )
         integer :: incx, k1, k2, lda, n
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32
         real(fltp) :: temp
         ! interchange row i with row ipiv(i) for each of rows k1 through k2.
         if( incx>0 ) then
            ix0 = k1
            i1 = k1
            i2 = k2
            inc = 1
         else if( incx<0 ) then
            ix0 = 1 + ( 1-k2 )*incx
            i1 = k2
            i2 = k1
            inc = -1
         else
            return
         end if
         n32 = ( n / 32 )*32
         if( n32/=0 ) then
            do j = 1, n32, 32
               ix = ix0
               do i = i1, i2, inc
                  ip = ipiv( ix )
                  if( ip/=i ) then
                     do k = j, j + 31
                        temp = a( i, k )
                        a( i, k ) = a( ip, k )
                        a( ip, k ) = temp
                     end do
                  end if
                  ix = ix + incx
               end do
            end do
         end if
         if( n32/=n ) then
            n32 = n32 + 1
            ix = ix0
            do i = i1, i2, inc
               ip = ipiv( ix )
               if( ip/=i ) then
                  do k = n32, n
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
                  end do
               end if
               ix = ix + incx
            end do
         end if
      end subroutine my_laswp


      subroutine my_laswp_4_by_1( a, lda, ipiv )
         ! n == 1, k1 == 1, k2 == 4, incx == 1
         integer :: lda
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         integer :: ip
         real(fltp) :: temp
         ! interchange row i with row ipiv(i) for each of rows k1 through k2.
         ip = ipiv( 1 )
         if( ip/=1 ) then
            temp = a( 1, 1 )
            a( 1, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 2 )
         if( ip/=2 ) then
            temp = a( 2, 1 )
            a( 2, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 3 )
         if( ip/=3 ) then
            temp = a( 3, 1 )
            a( 3, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 4 )
         if( ip/=4 ) then
            temp = a( 4, 1 )
            a( 4, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
      end subroutine my_laswp_4_by_1


      subroutine my_laswp_5_by_1( a, lda, ipiv )
         ! n == 1, k1 == 1, k2 == 5, incx == 1
         integer :: lda
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         integer :: ip
         real(fltp) :: temp
         ! interchange row i with row ipiv(i) for each of rows k1 through k2.
         ip = ipiv( 1 )
         if( ip/=1 ) then
            temp = a( 1, 1 )
            a( 1, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 2 )
         if( ip/=2 ) then
            temp = a( 2, 1 )
            a( 2, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 3 )
         if( ip/=3 ) then
            temp = a( 3, 1 )
            a( 3, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 4 )
         if( ip/=4 ) then
            temp = a( 4, 1 )
            a( 4, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
         ip = ipiv( 5 )
         if( ip/=5 ) then
            temp = a( 5, 1 )
            a( 5, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
         end if
      end subroutine my_laswp_5_by_1


      subroutine my_laswp_4_by_4( a, lda, ipiv )
         integer :: lda
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         real(fltp) :: temp
         integer :: ip
         ! interchange row i with row ipiv(i) for each of rows 1 through 4.
         ip = ipiv( 1 )
         if( ip/=1 ) then
            temp = a( 1, 1 )
            a( 1, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 1, 2 )
            a( 1, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 1, 3 )
            a( 1, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 1, 4 )
            a( 1, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
         end if
         ip = ipiv( 2 )
         if( ip/=2 ) then
            temp = a( 2, 1 )
            a( 2, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 2, 2 )
            a( 2, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 2, 3 )
            a( 2, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 2, 4 )
            a( 2, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
         end if
         ip = ipiv( 3 )
         if( ip/=3 ) then
            temp = a( 3, 1 )
            a( 3, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 3, 2 )
            a( 3, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 3, 3 )
            a( 3, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 3, 4 )
            a( 3, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
         end if
         ip = ipiv( 4 )
         if( ip/=4 ) then
            temp = a( 4, 1 )
            a( 4, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 4, 2 )
            a( 4, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 4, 3 )
            a( 4, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 4, 4 )
            a( 4, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
         end if
      end subroutine my_laswp_4_by_4


      subroutine my_laswp_5_by_5( a, lda, ipiv )
         integer :: lda
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         real(fltp) :: temp
         integer :: ip
         ! interchange row i with row ipiv(i) for each of rows 1 through 5.
         ip = ipiv( 1 )
         if( ip/=1 ) then
            temp = a( 1, 1 )
            a( 1, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 1, 2 )
            a( 1, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 1, 3 )
            a( 1, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 1, 4 )
            a( 1, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
            temp = a( 1, 5 )
            a( 1, 5 ) = a( ip, 5 )
            a( ip, 5 ) = temp
         end if
         ip = ipiv( 2 )
         if( ip/=2 ) then
            temp = a( 2, 1 )
            a( 2, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 2, 2 )
            a( 2, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 2, 3 )
            a( 2, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 2, 4 )
            a( 2, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
            temp = a( 2, 5 )
            a( 2, 5 ) = a( ip, 5 )
            a( ip, 5 ) = temp
         end if
         ip = ipiv( 3 )
         if( ip/=3 ) then
            temp = a( 3, 1 )
            a( 3, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 3, 2 )
            a( 3, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 3, 3 )
            a( 3, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 3, 4 )
            a( 3, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
            temp = a( 3, 5 )
            a( 3, 5 ) = a( ip, 5 )
            a( ip, 5 ) = temp
         end if
         ip = ipiv( 4 )
         if( ip/=4 ) then
            temp = a( 4, 1 )
            a( 4, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 4, 2 )
            a( 4, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 4, 3 )
            a( 4, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 4, 4 )
            a( 4, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
            temp = a( 4, 5 )
            a( 4, 5 ) = a( ip, 5 )
            a( ip, 5 ) = temp
         end if
         ip = ipiv( 5 )
         if( ip/=5 ) then
            temp = a( 5, 1 )
            a( 5, 1 ) = a( ip, 1 )
            a( ip, 1 ) = temp
            temp = a( 5, 2 )
            a( 5, 2 ) = a( ip, 2 )
            a( ip, 2 ) = temp
            temp = a( 5, 3 )
            a( 5, 3 ) = a( ip, 3 )
            a( ip, 3 ) = temp
            temp = a( 5, 4 )
            a( 5, 4 ) = a( ip, 4 )
            a( ip, 4 ) = temp
            temp = a( 5, 5 )
            a( 5, 5 ) = a( ip, 5 )
            a( ip, 5 ) = temp
         end if
      end subroutine my_laswp_5_by_5


      subroutine my_getrs( n, nrhs, a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb, n, nrhs
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: one=1, zero=0
         integer :: i, j, k
         info = 0

         if (nrhs == 1) then
            if (n == 4) then
               call my_getrs_4_by_1( a, lda, ipiv, b, ldb, info )
               return
            end if
            if (n == 5) then
               call my_getrs_5_by_1( a, lda, ipiv, b, ldb, info )
               return
            end if
         else if (nrhs == 4 .and. n == 4) then
            call my_getrs_4_by_4( a, lda, ipiv, b, ldb, info )
            return
         else if (nrhs == 5 .and. n == 5) then
            call my_getrs_5_by_5( a, lda, ipiv, b, ldb, info )
            return
         end if

         call my_laswp(nrhs, b, ldb, 1, n, ipiv, 1 )
         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = 1,n
               if (b(k,j)/=zero) then
                  do i = k + 1,n
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = n,1,-1
               if (b(k,j)/=zero) then
                  b(k,j) = b(k,j)/a(k,k)
                  do i = 1,k - 1
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do

      end subroutine my_getrs


      subroutine my_getrs_5_by_5( a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb  ! , n=5, nrhs=5
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: zero=0

         info = 0

         !call my_laswp(5, b, ldb, 1, 5, ipiv, 1 )
         call my_laswp_5_by_5( b, ldb, ipiv )

         b(2,1) = b(2,1) - b(1,1)*a(2,1)
         b(3,1) = b(3,1) - b(1,1)*a(3,1)
         b(4,1) = b(4,1) - b(1,1)*a(4,1)
         b(5,1) = b(5,1) - b(1,1)*a(5,1)
         b(3,1) = b(3,1) - b(2,1)*a(3,2)
         b(4,1) = b(4,1) - b(2,1)*a(4,2)
         b(5,1) = b(5,1) - b(2,1)*a(5,2)
         b(4,1) = b(4,1) - b(3,1)*a(4,3)
         b(5,1) = b(5,1) - b(3,1)*a(5,3)
         b(5,1) = b(5,1) - b(4,1)*a(5,4)

         b(2,2) = b(2,2) - b(1,2)*a(2,1)
         b(3,2) = b(3,2) - b(1,2)*a(3,1)
         b(4,2) = b(4,2) - b(1,2)*a(4,1)
         b(5,2) = b(5,2) - b(1,2)*a(5,1)
         b(3,2) = b(3,2) - b(2,2)*a(3,2)
         b(4,2) = b(4,2) - b(2,2)*a(4,2)
         b(5,2) = b(5,2) - b(2,2)*a(5,2)
         b(4,2) = b(4,2) - b(3,2)*a(4,3)
         b(5,2) = b(5,2) - b(3,2)*a(5,3)
         b(5,2) = b(5,2) - b(4,2)*a(5,4)

         b(2,3) = b(2,3) - b(1,3)*a(2,1)
         b(3,3) = b(3,3) - b(1,3)*a(3,1)
         b(4,3) = b(4,3) - b(1,3)*a(4,1)
         b(5,3) = b(5,3) - b(1,3)*a(5,1)
         b(3,3) = b(3,3) - b(2,3)*a(3,2)
         b(4,3) = b(4,3) - b(2,3)*a(4,2)
         b(5,3) = b(5,3) - b(2,3)*a(5,2)
         b(4,3) = b(4,3) - b(3,3)*a(4,3)
         b(5,3) = b(5,3) - b(3,3)*a(5,3)
         b(5,3) = b(5,3) - b(4,3)*a(5,4)

         b(2,4) = b(2,4) - b(1,4)*a(2,1)
         b(3,4) = b(3,4) - b(1,4)*a(3,1)
         b(4,4) = b(4,4) - b(1,4)*a(4,1)
         b(5,4) = b(5,4) - b(1,4)*a(5,1)
         b(3,4) = b(3,4) - b(2,4)*a(3,2)
         b(4,4) = b(4,4) - b(2,4)*a(4,2)
         b(5,4) = b(5,4) - b(2,4)*a(5,2)
         b(4,4) = b(4,4) - b(3,4)*a(4,3)
         b(5,4) = b(5,4) - b(3,4)*a(5,3)
         b(5,4) = b(5,4) - b(4,4)*a(5,4)

         b(2,5) = b(2,5) - b(1,5)*a(2,1)
         b(3,5) = b(3,5) - b(1,5)*a(3,1)
         b(4,5) = b(4,5) - b(1,5)*a(4,1)
         b(5,5) = b(5,5) - b(1,5)*a(5,1)
         b(3,5) = b(3,5) - b(2,5)*a(3,2)
         b(4,5) = b(4,5) - b(2,5)*a(4,2)
         b(5,5) = b(5,5) - b(2,5)*a(5,2)
         b(4,5) = b(4,5) - b(3,5)*a(4,3)
         b(5,5) = b(5,5) - b(3,5)*a(5,3)
         b(5,5) = b(5,5) - b(4,5)*a(5,4)

         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         b(5,1) = b(5,1)/a(5,5)
         b(1,1) = b(1,1) - b(5,1)*a(1,5)
         b(2,1) = b(2,1) - b(5,1)*a(2,5)
         b(3,1) = b(3,1) - b(5,1)*a(3,5)
         b(4,1) = b(4,1) - b(5,1)*a(4,5)
         b(4,1) = b(4,1)/a(4,4)
         b(1,1) = b(1,1) - b(4,1)*a(1,4)
         b(2,1) = b(2,1) - b(4,1)*a(2,4)
         b(3,1) = b(3,1) - b(4,1)*a(3,4)
         b(3,1) = b(3,1)/a(3,3)
         b(1,1) = b(1,1) - b(3,1)*a(1,3)
         b(2,1) = b(2,1) - b(3,1)*a(2,3)
         b(2,1) = b(2,1)/a(2,2)
         b(1,1) = b(1,1) - b(2,1)*a(1,2)
         b(1,1) = b(1,1)/a(1,1)

         b(5,2) = b(5,2)/a(5,5)
         b(1,2) = b(1,2) - b(5,2)*a(1,5)
         b(2,2) = b(2,2) - b(5,2)*a(2,5)
         b(3,2) = b(3,2) - b(5,2)*a(3,5)
         b(4,2) = b(4,2) - b(5,2)*a(4,5)
         b(4,2) = b(4,2)/a(4,4)
         b(1,2) = b(1,2) - b(4,2)*a(1,4)
         b(2,2) = b(2,2) - b(4,2)*a(2,4)
         b(3,2) = b(3,2) - b(4,2)*a(3,4)
         b(3,2) = b(3,2)/a(3,3)
         b(1,2) = b(1,2) - b(3,2)*a(1,3)
         b(2,2) = b(2,2) - b(3,2)*a(2,3)
         b(2,2) = b(2,2)/a(2,2)
         b(1,2) = b(1,2) - b(2,2)*a(1,2)
         b(1,2) = b(1,2)/a(1,1)

         b(5,3) = b(5,3)/a(5,5)
         b(1,3) = b(1,3) - b(5,3)*a(1,5)
         b(2,3) = b(2,3) - b(5,3)*a(2,5)
         b(3,3) = b(3,3) - b(5,3)*a(3,5)
         b(4,3) = b(4,3) - b(5,3)*a(4,5)
         b(4,3) = b(4,3)/a(4,4)
         b(1,3) = b(1,3) - b(4,3)*a(1,4)
         b(2,3) = b(2,3) - b(4,3)*a(2,4)
         b(3,3) = b(3,3) - b(4,3)*a(3,4)
         b(3,3) = b(3,3)/a(3,3)
         b(1,3) = b(1,3) - b(3,3)*a(1,3)
         b(2,3) = b(2,3) - b(3,3)*a(2,3)
         b(2,3) = b(2,3)/a(2,2)
         b(1,3) = b(1,3) - b(2,3)*a(1,2)
         b(1,3) = b(1,3)/a(1,1)

         b(5,4) = b(5,4)/a(5,5)
         b(1,4) = b(1,4) - b(5,4)*a(1,5)
         b(2,4) = b(2,4) - b(5,4)*a(2,5)
         b(3,4) = b(3,4) - b(5,4)*a(3,5)
         b(4,4) = b(4,4) - b(5,4)*a(4,5)
         b(4,4) = b(4,4)/a(4,4)
         b(1,4) = b(1,4) - b(4,4)*a(1,4)
         b(2,4) = b(2,4) - b(4,4)*a(2,4)
         b(3,4) = b(3,4) - b(4,4)*a(3,4)
         b(3,4) = b(3,4)/a(3,3)
         b(1,4) = b(1,4) - b(3,4)*a(1,3)
         b(2,4) = b(2,4) - b(3,4)*a(2,3)
         b(2,4) = b(2,4)/a(2,2)
         b(1,4) = b(1,4) - b(2,4)*a(1,2)
         b(1,4) = b(1,4)/a(1,1)

         b(5,5) = b(5,5)/a(5,5)
         b(1,5) = b(1,5) - b(5,5)*a(1,5)
         b(2,5) = b(2,5) - b(5,5)*a(2,5)
         b(3,5) = b(3,5) - b(5,5)*a(3,5)
         b(4,5) = b(4,5) - b(5,5)*a(4,5)
         b(4,5) = b(4,5)/a(4,4)
         b(1,5) = b(1,5) - b(4,5)*a(1,4)
         b(2,5) = b(2,5) - b(4,5)*a(2,4)
         b(3,5) = b(3,5) - b(4,5)*a(3,4)
         b(3,5) = b(3,5)/a(3,3)
         b(1,5) = b(1,5) - b(3,5)*a(1,3)
         b(2,5) = b(2,5) - b(3,5)*a(2,3)
         b(2,5) = b(2,5)/a(2,2)
         b(1,5) = b(1,5) - b(2,5)*a(1,2)
         b(1,5) = b(1,5)/a(1,1)

      end subroutine my_getrs_5_by_5


      subroutine my_getrs_5_by_1( a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb  ! , n=5, nrhs=1
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: zero=0
         info = 0
         call my_laswp_5_by_1( b, ldb, ipiv )
         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         b(2,1) = b(2,1) - b(1,1)*a(2,1)
         b(3,1) = b(3,1) - b(1,1)*a(3,1)
         b(4,1) = b(4,1) - b(1,1)*a(4,1)
         b(5,1) = b(5,1) - b(1,1)*a(5,1)
         b(3,1) = b(3,1) - b(2,1)*a(3,2)
         b(4,1) = b(4,1) - b(2,1)*a(4,2)
         b(5,1) = b(5,1) - b(2,1)*a(5,2)
         b(4,1) = b(4,1) - b(3,1)*a(4,3)
         b(5,1) = b(5,1) - b(3,1)*a(5,3)
         b(5,1) = b(5,1) - b(4,1)*a(5,4)

         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         b(5,1) = b(5,1)/a(5,5)
         b(1,1) = b(1,1) - b(5,1)*a(1,5)
         b(2,1) = b(2,1) - b(5,1)*a(2,5)
         b(3,1) = b(3,1) - b(5,1)*a(3,5)
         b(4,1) = b(4,1) - b(5,1)*a(4,5)
         b(4,1) = b(4,1)/a(4,4)
         b(1,1) = b(1,1) - b(4,1)*a(1,4)
         b(2,1) = b(2,1) - b(4,1)*a(2,4)
         b(3,1) = b(3,1) - b(4,1)*a(3,4)
         b(3,1) = b(3,1)/a(3,3)
         b(1,1) = b(1,1) - b(3,1)*a(1,3)
         b(2,1) = b(2,1) - b(3,1)*a(2,3)
         b(2,1) = b(2,1)/a(2,2)
         b(1,1) = b(1,1) - b(2,1)*a(1,2)
         b(1,1) = b(1,1)/a(1,1)

      end subroutine my_getrs_5_by_1


      subroutine my_getrs_4_by_4( a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb  ! , n=4, nrhs=4
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: zero=0
         info = 0

         call my_laswp_4_by_4( b, ldb, ipiv )

         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         b(2,1) = b(2,1) - b(1,1)*a(2,1)
         b(3,1) = b(3,1) - b(1,1)*a(3,1)
         b(4,1) = b(4,1) - b(1,1)*a(4,1)
         b(3,1) = b(3,1) - b(2,1)*a(3,2)
         b(4,1) = b(4,1) - b(2,1)*a(4,2)
         b(4,1) = b(4,1) - b(3,1)*a(4,3)

         b(2,2) = b(2,2) - b(1,2)*a(2,1)
         b(3,2) = b(3,2) - b(1,2)*a(3,1)
         b(4,2) = b(4,2) - b(1,2)*a(4,1)
         b(3,2) = b(3,2) - b(2,2)*a(3,2)
         b(4,2) = b(4,2) - b(2,2)*a(4,2)
         b(4,2) = b(4,2) - b(3,2)*a(4,3)

         b(2,3) = b(2,3) - b(1,3)*a(2,1)
         b(3,3) = b(3,3) - b(1,3)*a(3,1)
         b(4,3) = b(4,3) - b(1,3)*a(4,1)
         b(3,3) = b(3,3) - b(2,3)*a(3,2)
         b(4,3) = b(4,3) - b(2,3)*a(4,2)
         b(4,3) = b(4,3) - b(3,3)*a(4,3)

         b(2,4) = b(2,4) - b(1,4)*a(2,1)
         b(3,4) = b(3,4) - b(1,4)*a(3,1)
         b(4,4) = b(4,4) - b(1,4)*a(4,1)
         b(3,4) = b(3,4) - b(2,4)*a(3,2)
         b(4,4) = b(4,4) - b(2,4)*a(4,2)
         b(4,4) = b(4,4) - b(3,4)*a(4,3)

         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         b(4,1) = b(4,1)/a(4,4)
         b(1,1) = b(1,1) - b(4,1)*a(1,4)
         b(2,1) = b(2,1) - b(4,1)*a(2,4)
         b(3,1) = b(3,1) - b(4,1)*a(3,4)
         b(3,1) = b(3,1)/a(3,3)
         b(1,1) = b(1,1) - b(3,1)*a(1,3)
         b(2,1) = b(2,1) - b(3,1)*a(2,3)
         b(2,1) = b(2,1)/a(2,2)
         b(1,1) = b(1,1) - b(2,1)*a(1,2)
         b(1,1) = b(1,1)/a(1,1)

         b(4,2) = b(4,2)/a(4,4)
         b(1,2) = b(1,2) - b(4,2)*a(1,4)
         b(2,2) = b(2,2) - b(4,2)*a(2,4)
         b(3,2) = b(3,2) - b(4,2)*a(3,4)
         b(3,2) = b(3,2)/a(3,3)
         b(1,2) = b(1,2) - b(3,2)*a(1,3)
         b(2,2) = b(2,2) - b(3,2)*a(2,3)
         b(2,2) = b(2,2)/a(2,2)
         b(1,2) = b(1,2) - b(2,2)*a(1,2)
         b(1,2) = b(1,2)/a(1,1)

         b(4,3) = b(4,3)/a(4,4)
         b(1,3) = b(1,3) - b(4,3)*a(1,4)
         b(2,3) = b(2,3) - b(4,3)*a(2,4)
         b(3,3) = b(3,3) - b(4,3)*a(3,4)
         b(3,3) = b(3,3)/a(3,3)
         b(1,3) = b(1,3) - b(3,3)*a(1,3)
         b(2,3) = b(2,3) - b(3,3)*a(2,3)
         b(2,3) = b(2,3)/a(2,2)
         b(1,3) = b(1,3) - b(2,3)*a(1,2)
         b(1,3) = b(1,3)/a(1,1)

         b(4,4) = b(4,4)/a(4,4)
         b(1,4) = b(1,4) - b(4,4)*a(1,4)
         b(2,4) = b(2,4) - b(4,4)*a(2,4)
         b(3,4) = b(3,4) - b(4,4)*a(3,4)
         b(3,4) = b(3,4)/a(3,3)
         b(1,4) = b(1,4) - b(3,4)*a(1,3)
         b(2,4) = b(2,4) - b(3,4)*a(2,3)
         b(2,4) = b(2,4)/a(2,2)
         b(1,4) = b(1,4) - b(2,4)*a(1,2)
         b(1,4) = b(1,4)/a(1,1)

      end subroutine my_getrs_4_by_4


      subroutine my_getrs_4_by_1( a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb  ! , n=4, nrhs=1
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: zero=0

         info = 0
         call my_laswp_4_by_1( b, ldb, ipiv )

         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         b(2,1) = b(2,1) - b(1,1)*a(2,1)
         b(3,1) = b(3,1) - b(1,1)*a(3,1)
         b(4,1) = b(4,1) - b(1,1)*a(4,1)
         b(3,1) = b(3,1) - b(2,1)*a(3,2)
         b(4,1) = b(4,1) - b(2,1)*a(4,2)
         b(4,1) = b(4,1) - b(3,1)*a(4,3)

         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         b(4,1) = b(4,1)/a(4,4)
         b(1,1) = b(1,1) - b(4,1)*a(1,4)
         b(2,1) = b(2,1) - b(4,1)*a(2,4)
         b(3,1) = b(3,1) - b(4,1)*a(3,4)
         b(3,1) = b(3,1)/a(3,3)
         b(1,1) = b(1,1) - b(3,1)*a(1,3)
         b(2,1) = b(2,1) - b(3,1)*a(2,3)
         b(2,1) = b(2,1)/a(2,2)
         b(1,1) = b(1,1) - b(2,1)*a(1,2)
         b(1,1) = b(1,1)/a(1,1)

      end subroutine my_getrs_4_by_1


      subroutine my_getrs_dbg( n, nrhs, a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb, n, nrhs
         integer, pointer :: ipiv(:)
         real(fltp), pointer :: a(:,:), b(:,:)  ! a( lda, * ), b( ldb, * )
         real(fltp), parameter :: one=1, zero=0
         integer :: i, j, k
         info = 0
         call my_laswp_dbg(nrhs, b, ldb, 1, n, ipiv, 1 )
         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = 1,n
               if (b(k,j)/=zero) then
                  do i = k + 1,n
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = n,1,-1
               if (b(k,j)/=zero) then
                  b(k,j) = b(k,j)/a(k,k)
                  do i = 1,k - 1
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
      end subroutine my_getrs_dbg


      subroutine my_laswp_dbg( n,   a, lda,  k1, k2, ipiv,  incx )
         integer :: incx, k1, k2, lda, n
         integer :: ipiv(:)
         real(fltp) :: a(:,:)  ! a( lda, * )
         integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32
         real(fltp) :: temp
         include 'formats'
         write(*,2) 'n', n
         write(*,2) 'incx', incx
         write(*,2) 'k1', k1
         write(*,2) 'k2', k2
         do j = 1, n
            write(*,3) 'ipiv(j)', j, ipiv(j)
         end do
         ! interchange row i with row ipiv(i) for each of rows k1 through k2.
         if( incx>0 ) then
            ix0 = k1
            i1 = k1
            i2 = k2
            inc = 1
         else if( incx<0 ) then
            ix0 = 1 + ( 1-k2 )*incx
            i1 = k2
            i2 = k1
            inc = -1
         else
            return
         end if
         n32 = ( n / 32 )*32
         if( n32/=0 ) then
            do j = 1, n32, 32
               ix = ix0
               do i = i1, i2, inc
                  ip = ipiv( ix )
                  if( ip/=i ) then
                     do k = j, j + 31
                        temp = a( i, k )
                        a( i, k ) = a( ip, k )
                        a( ip, k ) = temp
                     end do
                  end if
                  ix = ix + incx
               end do
            end do
         end if
         if( n32/=n ) then
            n32 = n32 + 1
            ix = ix0
            do i = i1, i2, inc
               ip = ipiv( ix )


               if (ip == 0) then

                  stop 'my_lapack95  ip == 0'

               end if

               if( ip/=i ) then
                  do k = n32, n
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
                  end do
               end if
               ix = ix + incx
            end do
         end if
      end subroutine my_laswp_dbg

      end module my_lapack95
