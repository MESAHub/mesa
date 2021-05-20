! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module interp_1d_pm ! piecewise monotonic algorithms
      
      use const_lib, only: dp

      implicit none
      
      contains

      ! the following produce piecewise monotonic interpolants 
      ! rather than monotonicity preserving
      ! this stricter limit never introduces interpolated values exceeding the given values, 
      ! even in places where the given values are not monotonic.
      ! the downside is reduced accuracy on smooth data compared to the mp routines.
      
      
      subroutine mk_pmcub(x, nx, f1, slope_only, nwork, work1, str, ierr) 
         ! make piecewise monotonic cubic interpolant
         use interp_1d_def
         integer, intent(in) :: nx       ! length of x vector (nx >= 2)
         real(dp), intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real(dp), intent(inout), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         logical, intent(in) :: slope_only
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real(dp), intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: h, s, p
         integer :: i
         character (len=256) :: message
         logical, parameter :: dbg = .true.
         real(dp), pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         include 'formats'
         
         ierr = 0
         
         if (nx < 2) then
            return
         end if
         
         if (nx == 2) then
            call mk_pmlinear(x, f1, slope_only, nwork, work1, str, ierr)
            return
         end if
         
         if (nx == 3) then
            call mk_pmquad(x, f1, slope_only, nwork, work1, str, ierr)
            return
         end if
         
         if (nx == 4) then
            call mk_pmcub4(x, f1, slope_only, nwork, work1, str, ierr)
            return
         end if
         
         if (nwork < pm_work_size) then
            ierr = -1
            return
         end if
         
         h(1:nx) => work1(1:nx)
         s(1:nx) => work1(1+nx:2*nx)
         p(1:nx) => work1(1+2*nx:3*nx)
         
         if (dbg) then
            h(:) = 0; s(:) = 0; p(:) = 0
         end if
         
         do i=1,nx-1
            h(i) = x(i+1) - x(i) ! width of interval
         end do
         do i = 1, nx-1
            if (h(i) == 0) then
               write(*, '(a,1x,2i5,1x,a)')  &
                  'same interpolation x values at', i, i+1, 'for ' // trim(str)
               ierr = -1
               return
            end if
         end do
         
         do i=1,nx-1
            s(i) = (f(1,i+1) - f(1,i)) / h(i) ! slope across interval
         end do
         
         do i=2,nx-1 
            p(i) = (s(i-1)*h(i) + s(i)*h(i-1))/(h(i-1)+h(i)) 
            ! slope at i of parabola through i-1, i, and i+1
         end do
         do i=2,nx-1 
            f(2,i) = (dsign(1d0, s(i-1))+dsign(1d0, s(i)))* &
                        min(abs(s(i-1)), abs(s(i)), 0.5d0*abs(p(i)))
            ! "safe" slope at i to ensure monotonic -- see Steffen's paper for explanation.
         end do
         
         p(1) = s(1)*(1 + h(1) / (h(1) + h(2))) - s(2) * h(1) / (h(1) + h(2)) 
            ! slope at 1 of parabola through 1st 3 points
         if (p(1)*s(1) <= 0) then
            f(2, 1) = 0
         else if (abs(p(1)) > 2*abs(s(1))) then
            f(2, 1) = 2*s(1)
         else
            f(2, 1) = p(1)
         end if
            
         p(nx) = s(nx-1)*(1 + h(nx-1) / (h(nx-1) + h(nx-2)))  &
                     - s(nx-2)*h(nx-1) / (h(nx-1) + h(nx-2)) 
            ! slope at nx of parabola through last 3 points
         if (p(nx)*s(nx-1) <= 0) then
            f(2,nx) = 0
         else if (abs(p(nx)) > 2*abs(s(nx-1))) then
            f(2,nx) = 2*s(nx-1)
         else
            f(2,nx) = p(nx)
         end if
         
         if (slope_only) return

         do i=1,nx-1
            f(3,i) = (3*s(i) - 2*f(2,i) - f(2,i+1)) / h(i)
            f(4,i) = (f(2,i) + f(2,i+1) - 2*s(i)) / (h(i)*h(i))
         end do
         f(3,nx) = 0
         f(4,nx) = 0
               
      end subroutine mk_pmcub

      
      ! optimize special case for nx = 4
      subroutine mk_pmcub4(x, f1, slope_only, nwork, work1, str, ierr)
         use interp_1d_def
         integer, parameter :: nx = 4   ! length of x vector (nx >= 2)
         real(dp), intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real(dp), intent(inout), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         logical, intent(in) :: slope_only
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real(dp), intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: h, s, p
         integer :: i
         character (len=256) :: message
         real(dp), pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         ierr = 0
         
         if (nwork < pm_work_size) then
            ierr = -1
            return
         end if
         
         h(1:nx) => work1(1:nx)
         s(1:nx) => work1(1+nx:2*nx)
         p(1:nx) => work1(1+2*nx:3*nx)
         
         do i=1,nx-1
            h(i) = x(i+1) - x(i) ! width of interval
         end do
         do i = 1, nx-1
            if (h(i) == 0) then
               write(*, '(a,1x,2i5,1x,a)')  &
                  'same interpolation x values at', i, i+1, 'for ' // trim(str)
               ierr = -1
               return
            end if
         end do
         
         do i=1,nx-1
            s(i) = (f(1,i+1) - f(1,i)) / h(i) ! slope across interval
         end do
         do i=2,nx-1 
            p(i) = (s(i-1)*h(i) + s(i)*h(i-1))/(h(i-1)+h(i)) 
            ! slope at i of parabola through i-1, i, and i+1
         end do
         do i=2,nx-1 
            f(2,i) = (dsign(1d0, s(i-1))+dsign(1d0, s(i)))* &
                        min(abs(s(i-1)), abs(s(i)), 0.5d0*abs(p(i)))
            ! "safe" slope at i to ensure monotonic -- see Steffen's paper for explanation.
         end do

         p(1) = s(1)*(1 + h(1) / (h(1) + h(2))) - s(2) * h(1) / (h(1) + h(2)) 
            ! slope at 1 of parabola through 1st 3 points
         if (p(1)*s(1) <= 0) then
            f(2, 1) = 0
         else if (abs(p(1)) > 2*abs(s(1))) then
            f(2, 1) = 2*s(1)
         else
            f(2, 1) = p(1)
         end if
            
         p(nx) = s(nx-1)*(1 + h(nx-1) / (h(nx-1) + h(nx-2)))  &
                     - s(nx-2)*h(nx-1) / (h(nx-1) + h(nx-2)) 
            ! slope at nx of parabola through last 3 points
         if (p(nx)*s(nx-1) <= 0) then
            f(2, nx) = 0
         else if (abs(p(nx)) > 2*abs(s(nx-1))) then
            f(2, nx) = 2*s(nx-1)
         else
            f(2, nx) = p(nx)
         end if
         
         if (slope_only) return

         do i=1,nx-1
            f(3,i) = (3*s(i) - 2*f(2,i) - f(2,i+1)) / h(i)
            f(4,i) = (f(2,i) + f(2,i+1) - 2*s(i)) / (h(i)*h(i))
         end do
         f(3,nx) = 0
         f(4,nx) = 0
      
      end subroutine mk_pmcub4
      
      
      ! optimize special case for nx = 3
      subroutine mk_pmquad(x, f1, slope_only, nwork, work1, str, ierr) 
         ! make piecewise monotonic quadratic interpolant
         use interp_1d_def
         integer, parameter :: nx = 3
         real(dp), intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real(dp), intent(inout), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         logical, intent(in) :: slope_only
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real(dp), intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: h, s, p
         integer :: i
         character (len=256) :: message
         real(dp), pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         if (nwork < pm_work_size) then
            ierr = -1
            return
         end if
         ierr = 0
         
         h(1:nx) => work1(1:nx)
         s(1:nx) => work1(1+nx:2*nx)
         p(1:nx) => work1(1+2*nx:3*nx)
         
         do i=1,nx-1
            h(i) = x(i+1) - x(i) ! width of interval
         end do
         do i = 1, nx-1
            if (h(i) == 0) then
               write(*, '(a,1x,2i5,1x,a)')  &
                  'same interpolation x values at', i, i+1, 'for ' // trim(str)
               ierr = -1
               return
            end if
         end do
         
         do i=1,nx-1
            s(i) = (f(1,i+1) - f(1,i)) / h(i) ! slope across interval
         end do
         do i=2,nx-1 
            p(i) = (s(i-1)*h(i) + s(i)*h(i-1))/(h(i-1)+h(i)) 
            ! slope at i of parabola through i-1, i, and i+1
         end do
         do i=2,nx-1 
            f(2,i) = (dsign(1d0, s(i-1))+dsign(1d0, s(i)))* &
                        min(abs(s(i-1)), abs(s(i)), 0.5d0*abs(p(i)))
            ! "safe" slope at i to ensure monotonic -- see Steffen's paper for explanation.
         end do

         p(1) = s(1)*(1 + h(1) / (h(1) + h(2))) - s(2) * h(1) / (h(1) + h(2)) 
            ! slope at 1 of parabola through 1st 3 points
         if (p(1)*s(1) <= 0) then
            f(2, 1) = 0
         else if (abs(p(1)) > 2*abs(s(1))) then
            f(2, 1) = 2*s(1)
         else
            f(2, 1) = p(1)
         end if
            
         p(nx) = s(nx-1)*(1 + h(nx-1) / (h(nx-1) + h(nx-2)))  &
                     - s(nx-2)*h(nx-1) / (h(nx-1) + h(nx-2)) 
            ! slope at nx of parabola through last 3 points
         if (p(nx)*s(nx-1) <= 0) then
            f(2, nx) = 0
         else if (abs(p(nx)) > 2*abs(s(nx-1))) then
            f(2, nx) = 2*s(nx-1)
         else
            f(2, nx) = p(nx)
         end if
         
         if (slope_only) return
         
         do i=1,nx-1
            f(3,i) = (3*s(i) - 2*f(2,i) - f(2,i+1)) / h(i)
            f(4,i) = (f(2,i) + f(2,i+1) - 2*s(i)) / (h(i)*h(i))
         end do
         f(3,nx) = 0         
         f(4,nx) = 0
      
      end subroutine mk_pmquad
      
      
      ! optimize special case for nx = 2
      subroutine mk_pmlinear(x, f1, slope_only, nwork, work1, str, ierr)
         use interp_1d_def
         integer, parameter :: nx = 2
         real(dp), intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real(dp), intent(inout), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         logical, intent(in) :: slope_only
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real(dp), intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr

         real(dp) :: h, s
         character (len=256) :: message
         real(dp), pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         ierr = 0
         
         if (nwork < pm_work_size) then
            ierr = -1
            return
         end if
         
         h = x(2) - x(1) ! width of interval
         if (h == 0) then
            ierr = -1
            write(*, '(a,1x,2i5,1x,a)')  &
                  'same interpolation x values at', 1, 2, 'for ' // trim(str)
            return
         end if
         
         s = (f(1, 2) - f(1, 1)) / h ! slope across interval
         f(2, 1) = s
         f(2, 2) = 0
         
         if (slope_only) return

         f(3, 1:2) = 0
         f(4, 1:2) = 0
      
      end subroutine mk_pmlinear
      
      
      subroutine mk_pmcub_uniform(dx, n, f1, slope_only, nwork, work1, str, ierr) 
         ! make piecewise monotonic cubic interpolant on unit spaced mesh
         use interp_1d_def
         integer, intent(in) :: n     ! length of vector
         real(dp), intent(in) :: dx
         real(dp), intent(inout), pointer :: f1(:) ! =(4,n)  ! data & interpolation coefficients
         logical, intent(in) :: slope_only
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real(dp), intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: s, p
         real(dp) :: x(2)
         integer :: i
         real(dp), pointer :: f(:,:) ! (4, n)  ! data & interpolation coefficients
         f(1:4,1:n) => f1(1:4*n)
         
         ierr = 0
         
         if (n < 2) then
            return
         end if
         
         if (n == 2) then
            x(1) = 0
            x(2) = dx
            call mk_pmlinear(x, f1, slope_only, nwork, work1, str, ierr)
            return
         end if

         if (dx == 0) then
            ierr = -1
            return
         end if
         
         if (nwork < pm_work_size) then
            ierr = -1
            return
         end if
         ierr = 0
         
         s(1:n) => work1(1:n)
         p(1:n) => work1(1+n:2*n)

         ierr = 0
         
         do i=1,n-1 
            s(i) = f(1,i+1) - f(1,i) ! slope across interval
         end do
         do i=2,n-1 
            p(i) = 0.5d0*(s(i-1) + s(i)) 
            ! slope at i of parabola through i-1, i, and i+1
         end do
         do i=2,n-1 
            f(2,i) = (sign(1d0, s(i-1))+sign(1d0, s(i)))* &
               min(abs(s(i-1)), abs(s(i)), 0.5d0*abs(p(i)))
            ! "safe" slope at i to ensure monotonic -- see Steffen's paper for explanation.
         end do
         
         p(1) = 1.5d0 * s(1) - 0.5d0 * s(2)
            ! slope at 1 of parabola through 1st 3 points
         if (p(1)*s(1) <= 0) then
            f(2, 1) = 0
         else if (abs(p(1)) > 2*abs(s(1))) then
            f(2, 1) = 2*s(1)
         else
            f(2, 1) = p(1)
         end if
            
         p(n) = 1.5d0 * s(n-1) - 0.5d0 * s(n-2)
            ! slope at n of parabola through last 3 points
         if (p(n)*s(n-1) <= 0) then
            f(2, n) = 0
         else if (abs(p(n)) > 2*abs(s(n-1))) then
            f(2, n) = 2*s(n-1)
         else
            f(2, n) = p(n)
         end if
         
         f(2, 1:n) = f(2, 1:n) / dx
         
         if (slope_only) return
         
         ! 2nd and 3rd derivatives
         do i=1,n-1
            f(3,i) = (3*s(i) - 2*f(2,i) - f(2,i+1)) / dx       
            f(4,i) = (f(2,i) + f(2,i+1) - 2*s(i)) / (dx*dx)
         end do
         f(3,n) = (3*f(1, n-1) - 3*f(1, n) + (f(2, n-1) + 2*f(2, n)) * dx) / (dx*dx*dx)
         f(4,n) = (-2*f(1, n-1) + 2*f(1, n) - (f(2, n-1) + f(2, n))*dx) / (dx*dx*dx)
      
      end subroutine mk_pmcub_uniform


      end module interp_1d_pm
