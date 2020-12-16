! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module interp_1d_lib_sg
      implicit none
      
      contains

      
      ! this routine is a simply wrapper for making an interpolant and then using it.
      subroutine interpolate_vector_sg( &
               n_old, x_old, n_new, x_new, v_old, v_new, interp_vec_sg, nwork, work1, str, ierr)
         integer, intent(in) :: n_old, n_new
         real, intent(in) :: x_old(:) ! (n_old)
         real, intent(in) :: v_old(:) ! (n_old)
         real, intent(in) :: x_new(:) ! (n_new)
         real, intent(inout) :: v_new(:) ! (n_new)
         interface
            subroutine interp_vec_sg(x, nx, f1, nwork, work1, str, ierr) ! make cubic interpolant
               ! e.g., interp_pm, interp_m3a, interp_m3b, or interp_m3q
               integer, intent(in) :: nx       ! length of x vector
               real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
               real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
               integer, intent(in) :: nwork
               real, intent(inout), pointer :: work1(:)
               character (len=*) :: str ! for debugging
               integer, intent(out) :: ierr
            end subroutine interp_vec_sg
         end interface
         integer, intent(in) :: nwork
         real, intent(inout), pointer :: work1(:) ! =(:,:) ! (n_old, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         real, pointer :: f1(:), f(:,:)
         integer :: k
         ierr = 0
         allocate(f1(4*n_old), stat=ierr)
         if (ierr /= 0) return
         f(1:4,1:n_old) => f1(1:4*n_old)
         forall (k=1:n_old) f(1,k) = v_old(k)
         call interp_vec_sg(x_old, n_old, f1, nwork, work1, str, ierr) ! make interpolant
         if (ierr /= 0) then
            deallocate(f1)
            return
         end if
         call interp_values_sg(x_old, n_old, f1, n_new, x_new, v_new, ierr)
         deallocate(f1)
      end subroutine interpolate_vector_sg
      
      
      ! this routine is a simply wrapper for making an interpolant with interp_pm and then using it.
      subroutine interpolate_vector_pm_sg( &
               n_old, x_old, n_new, x_new, v_old, v_new, work1, str, ierr)
         use interp_1d_def, only: pm_work_size
         integer, intent(in) :: n_old, n_new
         real, intent(in) :: x_old(:) ! (n_old)
         real, intent(in) :: v_old(:) ! (n_old)
         real, intent(in) :: x_new(:) ! (n_new)
         real, intent(inout) :: v_new(:) ! (n_new)
         real, intent(inout), pointer :: work1(:) ! =(:,:) ! (n_old, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         real, pointer :: f1(:), f(:,:)
         integer :: k
         ierr = 0
         allocate(f1(4*n_old), stat=ierr)
         if (ierr /= 0) return
         f(1:4,1:n_old) => f1(1:4*n_old)
         forall (k=1:n_old) f(1,k) = v_old(k)
         call interp_pm_sg(x_old, n_old, f1, pm_work_size, work1, str, ierr) ! make interpolant
         if (ierr /= 0) then
            deallocate(f1)
            return
         end if
         call interp_values_sg(x_old, n_old, f1, n_new, x_new, v_new, ierr)
         deallocate(f1)
      end subroutine interpolate_vector_pm_sg
      
      
      subroutine interp_4_to_1_sg( &
               pdqm1, pdq00, pdqp1, ndq00, pfm1, pf00, pfp1, pfp2, nf00, str, ierr) 
         ! 4 points in, 1 point out
         ! piecewise monotonic cubic interpolation
         use interp_1d_def, only: pm_work_size
         real, intent(in) :: pdqm1, pdq00, pdqp1 ! spacing between input points
         real, intent(in) :: ndq00
         real, intent(in) :: pfm1, pf00, pfp1, pfp2 ! values at input points
         real, intent(out) :: nf00 ! new value at ndq00
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         integer, parameter :: n_old=4, n_new=1
         real :: x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
         real, target :: work1_ary(n_old*pm_work_size)
         real, pointer :: work1(:)
         work1 => work1_ary
         ierr = 0
         x_old(1) = 0.0
         x_old(2) = pdqm1
         x_old(3) = pdqm1+pdq00
         x_old(4) = pdqm1+pdq00+pdqp1
         v_old(1) = pfm1
         v_old(2) = pf00
         v_old(3) = pfp1
         v_old(4) = pfp2
         x_new(1) = ndq00       
         call interpolate_vector_pm_sg( &
            n_old, x_old, n_new, x_new, v_old, v_new, work1, str, ierr)        
         nf00 = v_new(1)         
      end subroutine interp_4_to_1_sg
         
         
      subroutine interp_3_to_1_sg(pdqm1, pdq00, ndqm1, pfm1, pf00, pfp1, nf00, str, ierr) 
         ! 3 points in, 1 point out
         ! piecewise monotonic quadratic interpolation
         use interp_1d_def, only: pm_work_size
         real, intent(in) :: pdqm1, pdq00 ! spacing between input points
         real, intent(in) :: ndqm1 ! new spacing to nk
         real, intent(in) :: pfm1, pf00, pfp1 ! values at pkm1, pk, pkp1
         real, intent(out) :: nf00 ! new value at ndq00
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         integer, parameter :: n_old=3, n_new=1
         real :: x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
         real, target :: work1_ary(n_old*pm_work_size)
         real, pointer :: work1(:)
         work1 => work1_ary
         ierr = 0         
         x_old(1) = 0.0
         x_old(2) = pdqm1
         x_old(3) = pdqm1+pdq00
         v_old(1) = pfm1
         v_old(2) = pf00
         v_old(3) = pfp1
         x_new(1) = ndqm1         
         call interpolate_vector_pm_sg( &
            n_old, x_old, n_new, x_new, v_old, v_new, work1, str, ierr)         
         nf00 = v_new(1)
      end subroutine interp_3_to_1_sg
         
         
      subroutine interp_3_to_2_sg(pdqm1, pdq00, ndqm1, ndq00, pfm1, pf00, pfp1, nf00, nfp1, str, ierr) 
         ! 3 points in, 2 points out
         ! piecewise monotonic quadratic interpolation
         use interp_1d_def, only: pm_work_size
         real, intent(in) :: pdqm1, pdq00 ! previous spacing to pk and pkp1
         real, intent(in) :: ndqm1, ndq00 ! new spacing to nk and nk+1
         real, intent(in) :: pfm1, pf00, pfp1 ! previous values at pkm1, pk, pkp1
         real, intent(out) :: nf00, nfp1 ! new values at nk, nk+1
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         integer, parameter :: n_old=3, n_new=2
         real :: x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
         real, target :: work1_ary(n_old*pm_work_size)
         real, pointer :: work1(:)
         work1 => work1_ary
         ierr = 0         
         x_old(1) = 0d0
         x_old(2) = pdqm1
         x_old(3) = pdqm1+pdq00
         v_old(1) = pfm1
         v_old(2) = pf00
         v_old(3) = pfp1
         x_new(1) = ndqm1       
         x_new(2) = ndqm1+ndq00        
         call interpolate_vector_pm_sg( &
            n_old, x_old, n_new, x_new, v_old, v_new, work1, str, ierr)         
         nf00 = v_new(1)
         nfp1 = v_new(2)         
      end subroutine interp_3_to_2_sg

      
      ! general routines
      
      ! these routines use previously created interpolant information (f)
      ! the interpolant can come from either the piecewise monotonic routines, or
      ! from the monotonicity preserving routines -- they use the same format for f.
      
      subroutine interp_values_sg(init_x, nx, f1, nv, x, vals, ierr)
         use interp_1d_def
         use interp_1d_misc_sg
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals vector
         real, intent(in) :: x(:) ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real, intent(inout) :: vals(:) ! (nv)
         integer, intent(out) :: ierr ! 0 means AOK
         call do_interp_values_sg(init_x, nx, f1, nv, x, vals, ierr)
      end subroutine interp_values_sg
      
      
      subroutine interp_value_sg(init_x, nx, f1, xval, val, ierr)
         use interp_1d_def
         use interp_1d_misc_sg
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         real, intent(in) :: xval  ! location where want interpolated values
         real, intent(out) :: val
         integer, intent(out) :: ierr ! 0 means AOK
         integer, parameter :: nv = 1
         real :: x(nv), vals(nv)
         x(1) = xval
         call do_interp_values_sg(init_x, nx, f1, nv, x, vals, ierr)
         val = vals(1)
      end subroutine interp_value_sg
      
      
      subroutine interp_values_and_slopes_sg(init_x, nx, f1, nv, x, vals, slopes, ierr)
         use interp_1d_def
         use interp_1d_misc_sg
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals and slopes vectors
         real, intent(in) :: x(:) ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real, intent(inout) :: vals(:) ! (nv)
         real, intent(inout) :: slopes(:) ! (nv)
         integer, intent(out) :: ierr ! 0 means AOK
         call do_interp_values_and_slopes_sg(init_x, nx, f1, nv, x, vals, slopes, ierr)
      end subroutine interp_values_and_slopes_sg
      
      
      subroutine interp_value_and_slope_sg(init_x, nx, f1, xval, val, slope, ierr)
         use interp_1d_def
         use interp_1d_misc_sg
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         real, intent(in) :: xval  ! location where want interpolated values
         real, intent(out) :: val, slope
         integer, intent(out) :: ierr ! 0 means AOK
         integer, parameter :: nv = 1
         real :: x(nv), vals(nv), slopes(nv)
         x(1) = xval
         call do_interp_values_and_slopes_sg(init_x, nx, f1, nv, x, vals, slopes, ierr)
         val = vals(1)
         slope = slopes(1)
      end subroutine interp_value_and_slope_sg

      
      subroutine integrate_values_sg(init_x, nx, f1, nv, x, vals, ierr)
         use interp_1d_def
         use interp_1d_misc_sg
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly increasing
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals vector
         real, intent(in) :: x(:) ! (nv)
            ! strictly monotonic in same way as init_x
            ! NOTE: no extrapolation allowed -- x's must be within range of init_x's
         real, intent(inout) :: vals(:) ! (nv)
            ! for i > 1, vals(i) = integral of interpolating poly from x(i-1) to x(i)
            ! vals(1) = 0
         integer, intent(out) :: ierr ! 0 means AOK

         call do_integrate_values_sg(init_x, nx, f1, nv, x, vals, ierr)
   
      end subroutine integrate_values_sg
      
      
      ! piecewise monotonic routines

      ! the following produce piecewise monotonic interpolants rather than monotonicity preserving
      ! this stricter limit never introduces interpolated values exceeding the given values, 
      ! even in places where the given values are not monotonic.
      ! the downside is reduced accuracy on smooth data compared to the mp routines.
      
      
      ! Steffen, M., "A simple method for monotonic interpolation in one dimension", 
      !        Astron. Astrophys., (239) 1990, 443-450.
      
      
      subroutine interp_pm_sg(x, nx, f1, nwork, work1, str, ierr) 
         ! make piecewise monotonic cubic interpolant
         use interp_1d_def
         use interp_1d_pm_sg
         integer, intent(in) :: nx       ! length of x vector (>= 2)
         real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx,nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr   
         call mk_pmcub_sg(x, nx, f1, .false., nwork, work1, str, ierr)         
      end subroutine interp_pm_sg
      
      
      subroutine interp_pm_slopes_only_sg(x, nx, f1, nwork, work1, str, ierr)
         ! identical to interp_pm, but only calculates slopes and stores them in f(2,:)
         ! this is a little faster for the special case in which you just want the slopes at x
         use interp_1d_def
         use interp_1d_pm_sg
         integer, intent(in) :: nx       ! length of x vector (>= 2)
         real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx,nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr   
         call mk_pmcub_sg(x, nx, f1, .true., nwork, work1, str, ierr)         
      end subroutine interp_pm_slopes_only_sg
      
      
      subroutine interp_4pt_pm_sg(x, y, a) 
         ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
         real, intent(in)    :: x(4)    ! junction points, strictly monotonic
         real, intent(in)    :: y(4)    ! data values at x's
         real, intent(inout)   :: a(3)    ! coefficients
         real :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3
         ! for x(2) <= x <= x(3) and dx = x-x(2), 
         ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
         h1 = x(2)-x(1)
         h2 = x(3)-x(2)
         h3 = x(4)-x(3)
         s1 = (y(2)-y(1))/h1
         s2 = (y(3)-y(2))/h2
         s3 = (y(4)-y(3))/h3
         p2 = (s1*h2+s2*h1)/(h1+h2)
         p3 = (s2*h3+s3*h2)/(h2+h3)
         as2 = abs(s2)
         ss2 = sign(1.0, s2)
         yp2 = (sign(1.0, s1)+ss2)*min(abs(s1), as2, 0.5*abs(p2))
         yp3 = (ss2+sign(1.0, s3))*min(as2, abs(s3), 0.5*abs(p3))
         a(1) = yp2
         a(2) = (3*s2-2*yp2-yp3)/h2
         a(3) = (yp2+yp3-2*s2)/(h2*h2)
      end subroutine interp_4pt_pm_sg
      
      
      subroutine interp_pm_on_uniform_grid_sg(dx, nx, f1, nwork, work1, str, ierr) 
         ! make piecewise monotonic cubic interpolant on uniformly spaced mesh
         use interp_1d_def
         use interp_1d_pm_sg
         real, intent(in) :: dx    ! grid spacing
         integer, intent(in) :: nx     ! length of vector (>= 2)
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= pm_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call mk_pmcub_uniform_sg(dx, nx, f1, .false., nwork, work1, str, ierr)      
      end subroutine interp_pm_on_uniform_grid_sg
      
      
      
      ! monotonicity preserving routines
      
      ! Huynh, H.T., "Accurate Monotone Cubic Interpolation", SIAM J Numer. Anal. (30) 1993, 57-100.
      
      ! Suresh, A, and H.T. Huynh, "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
      !        Time Stepping", JCP (136) 1997, 83-99.
      
      
      subroutine interp_m3_sg(x, nx, f1, which, nwork, work1, str, ierr) 
         ! make monotonicity preserving cubic interpolant on arbitrarily spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         integer, intent(in) :: nx       ! length of x vector (>= 4)
         real, intent(in) :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: which ! average, quartic, or super_bee
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_sg(x, nx, f1, which, .false., nwork, work1, str, ierr)     
      end subroutine interp_m3_sg


      subroutine interp_m3a_sg(x, nx, f1, nwork, work1, str, ierr) 
         ! make monotonicity preserving cubic interpolant on arbitrarily spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         integer, intent(in) :: nx       ! length of x vector (>= 4)
         real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_sg(x, nx, f1, average, .false., nwork, work1, str, ierr)      
      end subroutine interp_m3a_sg


      subroutine interp_m3b_sg(x, nx, f1, nwork, work1, str, ierr) 
         ! make monotonicity preserving cubic interpolant on arbitrarily spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         integer, intent(in) :: nx       ! length of x vector (>= 4)
         real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_sg(x, nx, f1, super_bee, .false., nwork, work1, str, ierr)     
      end subroutine interp_m3b_sg


      subroutine interp_m3q_sg(x, nx, f1, nwork, work1, str, ierr) 
         ! make monotonicity preserving cubic interpolant on arbitrarily spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         integer, intent(in) :: nx       ! length of x vector (>= 4)
         real, intent(in)    :: x(:) ! (nx)    ! junction points, strictly monotonic
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_sg(x, nx, f1, quartic, .false., nwork, work1, str, ierr)      
      end subroutine interp_m3q_sg
            
      
      subroutine interp_m3_on_uniform_grid_sg(dx, nx, f1, which, nwork, work1, str, ierr)
         ! make monotonicity preserving cubic interpolant on uniformly spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         real, intent(in) :: dx ! the grid spacing
         integer, intent(in) :: nx ! length of x vector (>= 4)
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: which ! average, quartic, or super_bee
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_on_uniform_grid_sg(dx, nx, f1, which, .false., nwork, work1, str, ierr)         
      end subroutine interp_m3_on_uniform_grid_sg
            
      
      subroutine interp_m3a_on_uniform_grid_sg(dx, nx, f1, nwork, work1, str, ierr)
         ! make monotonicity preserving cubic interpolant on uniformly spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         real, intent(in) :: dx ! the grid spacing
         integer, intent(in) :: nx ! length of x vector (>= 4)
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_on_uniform_grid_sg(dx, nx, f1, average, .false., nwork, work1, str, ierr)         
      end subroutine interp_m3a_on_uniform_grid_sg
            
      
      subroutine interp_m3b_on_uniform_grid_sg(dx, nx, f1, nwork, work1, str, ierr)
         ! make monotonicity preserving cubic interpolant on uniformly spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         real, intent(in) :: dx ! the grid spacing
         integer, intent(in) :: nx ! length of x vector (>= 4)
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_on_uniform_grid_sg(dx, nx, f1, super_bee, .false., nwork, work1, str, ierr)         
      end subroutine interp_m3b_on_uniform_grid_sg
            
      
      subroutine interp_m3q_on_uniform_grid_sg(dx, nx, f1, nwork, work1, str, ierr)
         ! make monotonicity preserving cubic interpolant on uniformly spaced grid
         use interp_1d_def
         use interp_1d_mp_sg
         real, intent(in) :: dx ! the grid spacing
         integer, intent(in) :: nx ! length of x vector (>= 4)
         real, intent(inout), pointer :: f1(:) ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nwork ! nwork must be >= mp_work_size (see interp_1d_def)
         real, intent(inout), pointer :: work1(:) ! =(nx, nwork)
         character (len=*) :: str ! for debugging
         integer, intent(out) :: ierr
         call m3_on_uniform_grid_sg(dx, nx, f1, quartic, .false., nwork, work1, str, ierr)        
      end subroutine interp_m3q_on_uniform_grid_sg
      

      end module interp_1d_lib_sg
