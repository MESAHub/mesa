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

      module bipm_db
      
      use const_def, only: dp

      implicit none


      contains

      
      subroutine do_mkbipm_db(x,nx,y,ny,f1,nf2,ierr)
         use interp_1d_def
         use interp_1d_lib
         integer, intent(in) :: nx                        ! length of x vector
         integer, intent(in) :: ny                        ! length of y vector
         real(dp), intent(in), pointer :: x(:) ! (nx)            ! x vector, strict ascending
         real(dp), intent(in), pointer :: y(:) ! (ny)            ! y vector, strict ascending
         integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
         real(dp), intent(inout), pointer :: f1(:) ! (4,nf2,ny)   ! data & interpolant coefficients
         integer, intent(out) :: ierr                      ! =0 on exit if there is no error.   
         integer, parameter :: nwork = pm_work_size
         integer :: i
         real(dp), target :: work_ary(nx*nwork)
         real(dp), pointer :: work(:), f(:)
         work => work_ary
         ierr = 0
         do i=1,ny
            f(1:4*nf2) => f1(1+(i-1)*4*nf2:i*4*nf2)
            call interp_pm(x, nx, f, nwork, work, 'do_mkbipm_db', ierr)
            if (ierr /= 0) exit
         end do
      end subroutine do_mkbipm_db


      subroutine do_evbipm_db(xget,yget,x,nx,y,ny,fin1,nf2,f,ierr)
         use num_lib, only: binary_search
         use interp_1d_def
         use interp_1d_lib
         integer, intent(in) :: nx,ny
         real(dp), intent(in) :: xget,yget        ! target of this interpolation
         real(dp), intent(in), pointer :: x(:) ! (nx)            ! ordered x grid
         real(dp), intent(in), pointer :: y(:) ! (ny)            ! ordered y grid
         integer, intent(in) :: nf2
         real(dp), intent(in), pointer :: fin1(:) ! fin(4,nf2,ny)      ! function data
         real(dp), intent(out) :: f
         integer, intent(out) :: ierr                      ! error code =0 ==> no error

         integer, parameter :: nwork = pm_work_size
         real(dp) :: x0,x1,dx,y0,y1,dy,alfa,beta,ddx,f1,f2
         real(dp) :: ys(4), ynew(1), val(1)
         integer :: j, jlo, jhi, i, ix, jy, ii
         real(dp), target :: work_ary(4*nwork), ff_ary(4*4)
         real(dp), pointer :: work(:), fin(:,:,:), ff(:,:), ff1(:)
         work => work_ary
         ff1 => ff_ary
         ff(1:4,1:4) => ff_ary(1:4*4)
         fin(1:4,1:nf2,1:ny) => fin1(1:4*nf2*ny)
         
         ierr = 0
         
         ix = binary_search(nx, x, 0, xget) ! x(ix) <= xget < x(ix+1)         
         if (ix < 1 .or. ix >= nx) then
             ierr = -1
             return
         end if
         jy = binary_search(ny, y, 0, yget) ! y(jy) <= yget < y(jy+1)
         if (jy < 1 .or. jy >= ny) then
             ierr = -1
             return
         end if

         x0 = x(ix); x1 = x(ix+1)
         y0 = y(jy); y1 = y(jy+1)
         dx = xget - x0
         dy = yget - y0
         beta = dy / (y1 - y0) ! fraction of y1 result
         alfa = 1-beta ! fraction of y0 result            
         
         ynew(1) = yget
         if (jy == 1) then
            jlo = 1;; ii = 1
         else if (jy >= ny-1) then
            jlo = jy-2; ii = 3
         else
            jlo = jy-1; ii = 2
         end if
            
         do i=1,4
            j = jlo+i-1
            ys(i) = y(j)
            ff(1,i) = fin(1,ix,j) + dx*(fin(2,ix,j) + dx*(fin(3,ix,j) + dx*fin(4,ix,j)))
         end do
      
         call interp_pm(ys, 4, ff1, nwork, work, 'do_evbipm_db', ierr) 
         if (ierr /= 0) return
      
         call interp_values(ys, 4, ff1, 1, ynew, val, ierr)
         if (ierr /= 0) return

         f = val(1)

         
      end subroutine do_evbipm_db


      end module bipm_db
