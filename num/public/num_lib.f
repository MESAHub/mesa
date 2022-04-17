! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module num_lib
      ! various numerical routines
      
      use const_def
      use math_lib
      use num_def
      use mod_integrate
      
      ! NOTE: because of copyright restrictions, 
      !       mesa doesn't use any routines from Numerical Recipes.

      
      implicit none


      contains ! the procedure interface for the library
      ! client programs should only call these routines.

      ! tests of numerical derivative
         include "num_dfridr.dek"
      
      ! safe root finding
      ! uses alternating bisection and inverse parabolic interpolation
      ! also have option to use derivative as accelerator (newton method)
      include "num_safe_root.dek"

      
      ! solvers for ODEs and DAEs.
      
      
      ! sometimes you just want a simple runge-kutta
         include "num_rk2.dek"
         include "num_rk4.dek"
         
      
      ! but there are lots of fancier options too.
      

      ! selections from the Hairer family of ODE/DAE integrators.
      ! from Ernst Hairer's website: http://www.unige.ch/~hairer/

      ! explicit ODE solvers based on methods of Dormand and Prince (for non-stiff problems)

         ! explicit Runge-Kutta ODE integrator of order 5
         ! with dense output of order 4
         include "num_dopri5.dek" !  "DOramand PRInce order 5"

         ! explicit Runge-Kutta ODE integrator of order 8
         ! with dense output of order 7
         include "num_dop853.dek" !  "DOramand Prince order 8(5, 3)"
         
         ! both integrators have automatic step size control and monitoring for stiffness.
      
         ! For a description see:
         ! Hairer, Norsett and Wanner (1993): 
         ! Solving Ordinary Differential Equations. Nonstiff Problems. 2nd edition. 
         ! Springer Series in Comput. Math., vol. 8.
         ! http://www.unige.ch/~hairer/books.html
      
      
      ! implicit solvers (for stiff problems)
      
         ! there are a bunch of implicit solvers to pick from (listed below), 
         ! but they all have pretty much the same arguments, 
         ! so I've provided a general routine, called "isolve", that let's you
         ! pass an extra argument to specify which one of the particular solvers
         ! that you want to use.
         include "num_isolve.dek"
         
         ! if possible, you should write your code to call isolve
         ! rather than calling one of the particular solvers.
         ! only call a specific solver if you need a feature it provides
         ! that isn't supported by isolve.
         
         ! you can find an example program using isolve in num/test/src/sample_ode_solver.f
         
         
      ! the implicit solver routines
      
         ! for detailed descriptions of these routines see: 
         ! Hairer and Wanner (1996): 
         ! Solving Ordinary Differential Equations. 
         ! Stiff and Differential-Algebraic Problems. 2nd edition. 
         ! Springer Series in Comput. Math., vol. 14.
         ! http://www.unige.ch/~hairer/books.html
      
         ! linearly implicit Runge-Kutta method (Rosenbrock)
         include "num_ros2.dek"    ! L-stable; 2 stages; order 2, 2 function evaluations.
            ! ros2 is suitable for use with approximate jacobians such as from numerical differences.
            ! ros2 is designed for use with Strang splitting and is reported to be able to cope
            ! with large time steps and artificial transients introduced at the beginning of split intervals.
            ! see Verwer et al, "A Second-Order Rosenbrock Method Applied to Photochemical Dispersion Problems", 
            ! SIAM J. Sci. Comput. (20), 1999, 1456-1480.
            
         include "num_rose2.dek"    ! L-stable; 3 stages; order 2, 3 function evaluations.   
            ! rose2 is unique among the implicit solvers in that the final function evaluation
            ! uses the solution vector for the step.

         include "num_rodas3.dek"  ! L-stable; 4 stages; order 3, 3 function evaluations.
         include "num_rodas4.dek"  ! L-stable; 6 stages; order 4, 6 function evaluations.     

         ! 3rd order; for parabolic equations.
         include "num_ros3p.dek"   ! A-stable; 3 stages; order 3, 2 function evaluations.
         include "num_ros3pl.dek"  ! L-stable; 4 stages; order 3, 3 function evaluations.
         ! 4th order; for parabolic equations.
         include "num_rodasp.dek"  ! L-stable; 6 stages; order 4, 6 function evaluations. 
         
         include "num_solvers_options.dek"


      ! which implicit solver should you use?
         
         ! somewhat surprisingly, in some cases the solvers
         ! that work well at high tolerances will fail with low
         ! tolerances and vice-versa.  so you need to match
         ! the solver to the problem.
      
         ! your best bet is to try them all on some typical cases.
         ! happily this isn't too hard to do since they all
         ! use the same function arguments and have (almost)
         ! identical calling sequences and options.
      
      
      ! flexible choice of linear algebra routines
      
         ! the solvers need to solve linear systems.
         ! this is typically done by first factoring the matrix A
         ! and then repeatedly using the factored form to solve
         ! A*x=b for various vectors b.
         
         ! rather than build-in a particular matrix solver, 
         ! the mesa versions of the solvers take as arguments
         ! routines to perform these tasks.  the mesa/mtx package
         ! includes several choices for implementations of the
         ! required routines.
         
      
      ! dense, banded, or sparse matrix
      
         ! All the packages allow the matrix to be in dense or banded form.
         ! the choice of sparse matrix package is not fixed by the solvers.
         ! the only constraint is the the sparse format must be either
         ! compressed row sparse or compressed column sparse.
         ! the mesa/mtx package comes with one option for a sparse
         ! package (based on SPARSKIT), and also has hooks for another (Super_LU).
         ! Since the sparse routines are passed as arguments to the solvers, 
         ! it is possible to experiment with different linear algebra
         ! packages without a great deal of effort.
         
         
      ! analytical or numerical jacobian
      
         ! to solve M*y' = f(y), the solvers need to have the jacobian matrix, df/dy.
         ! the jacobian can either be calculated analytically by a user supplied routine, 
         ! or the solver can form a numerical difference estimate by repeatedly
         ! evaluating f(y) with slightly different y's.  Such numerical jacobians
         ! are supported by all the solvers for both dense and banded matrix forms.
         ! For the sparse matrix case, only analytical jacobians are allowed.
         
         ! NOTE: for most implicit solvers, the accuracy of the jacobian influences
         ! the rate of convergence, but doesn't impact the accuracy of the solution.
         ! however, the rodas solvers are an exception to this rule.
         ! they are based on the rosenbrock method which replaces the newton iteration
         ! by formulas that directly use the jacobian in the formula for
         ! the result.  as a result, the rodas solvers depend on having
         ! accurate jacobians in order to produce accurate results.
      
      
      ! explicit or implicit ODE systems
      
         ! systems of the form y' = f(y) are called "explicit ODE systems".
         
         ! systems of the form M*y' = f(y), with M not equal to the identity matrix, 
         ! are called "implicit ODE systems".
         
         ! in addition to the usual explicit systems, 
         ! the solvers can also handle implicit ODE systems
         ! in which M is an arbitrary constant matrix, 
         ! even including the case of M singular.
         
         ! for M non-constant, see the discussion of "problems with special structure"
      
      
      ! problems with special structure
      
         ! 3 special cases can be handled easily
         
            ! case 1, second derivatives: y'' = f(t, y, y')
            ! case 2, nonconstant matrix: C(x, y)*y' = f(t, y)
            ! case 3, both of the above: C(x, y)*y'' = f(t, y, y')
            
         ! these all work by adding auxiliary variables to the problem and
         ! converting back to the standard form with a constant matrix M.
         
            ! case 1: y'' = f(t, y, y')
               ! after add auxiliary variables z, this becomes
               ! y' = z
               ! z' = f(t, y, z)
               
            ! case 2: C(x, y)*y' = f(t, y)
               ! after add auxiliary variables z, this becomes
               ! y' = z
               ! 0 = C(x, y)*z - f(t, y)
               
            ! case 3: C(x, y)*y'' = f(t, y, y')
               ! after add auxiliary variables z and u, this becomes
               ! y' = z
               ! z' = u
               ! 0 = C(x, y)*u - f(t, y, z)
         
         ! The last two cases take advantage of the ability to have M singular.
         
         ! If the matrix for df/dy is dense in these special cases, all the solvers
         ! can reduce the cost of the linear algebra operations by special treatment
         ! of the auxiliary variables.
         
      
      ! "projection" of solution to valid range of values.
      
         ! it is often the case that the n-dimensional solution
         ! is actually constrained to a subspace of full n dimensional
         ! space of numbers.  The proposed solutions at each step
         ! need to be projected back to the allowed subspace in order
         ! to maintain valid results.  The routines all provide for this
         ! option by calling a "solout" routine, supplied by the user, 
         ! after every accepted step.  The user's solout routine can modify
         ! the solution y before returning to continue the integration.
      
         
      ! "dense output"
      
         ! the routines provide estimates of the solution over entire step.
         ! useful for tabulating the solution at prescribed points
         ! or for smooth graphical presentation of the solution.
         ! also very useful for "event location" -- e.g., at what
         ! value of x do we get a solution y(x) s.t. some relation
         ! g(x, y(x))=0.  The dense output option is very helpful here.
         ! All of the solvers support dense output.
         ! BTW: there is typically a certain overhead associated with
         ! providing the option for dense output, so don't request it
         ! unless you'll really be using it.

         ! here is a special version of safe_root for use with dense output "solout" routines
         include "num_solout_root.dek"



      ! "null" implementations of routines used by the solvers
      ! are for use in cases in which you aren't actually using the routine, 
      ! but something must be provided for the required argument.


      subroutine null_fcn(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n) 
            ! okay to edit y if necessary (e.g., replace negative values by zeros)
         real(dp), intent(inout) :: f(:) ! (n) ! dy/dx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         f=0; ierr=0
      end subroutine null_fcn


      subroutine null_fcn_blk_dble(n, caller_id, nvar, nz, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: n, caller_id, nvar, nz, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout), pointer :: y(:) 
            ! (n) okay to edit y if necessary (e.g., replace negative values by zeros)
         real(dp), intent(inout), pointer :: f(:) ! (n) dy/dx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         f=0; ierr=0
      end subroutine null_fcn_blk_dble
      

      subroutine null_jac(n, x, h, y, f, dfy, ldfy, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, ldfy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n) ! dy/dx
         real(dp), intent(inout) :: dfy(:,:) ! (ldfy, n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         f=0; dfy=0; ierr=0
      end subroutine null_jac


      subroutine null_jac_blk_dble(n, caller_id, nvar, nz, x, h, y, f, lblk, dblk, ublk, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: n, caller_id, nvar, nz, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout), pointer :: y(:) ! (n)
         real(dp), intent(inout), pointer :: f(:) ! (n) dy/dx
         real(dp), dimension(:), pointer, intent(inout) :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         f=0; y=0; ierr=0
      end subroutine null_jac_blk_dble


      subroutine null_sjac(n, x, h, y, f, nzmax, ia, ja, values, lrpar, rpar, lipar, ipar, ierr)  
         ! sparse jacobian. format either compressed row or compressed column.
         integer, intent(in) :: n, nzmax, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n) ! dy/dx
         integer, intent(inout) :: ia(:) ! (n+1)
         integer, intent(inout) :: ja(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         f=0; values = 0; ia = 0; ja = 0; ierr = 0
      end subroutine null_sjac


      subroutine null_mas(n, am, lmas, lrpar, rpar, lipar, ipar)
         integer, intent(in) :: n, lmas, lrpar, lipar
         real(dp), intent(inout) :: am(:,:) ! (lmas, n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         am = 0
      end subroutine null_mas


      subroutine null_solout(nr, xold, x, n, y, rwork, iwork, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         interface
            real(dp) function interp_y(i, s, rwork, iwork, ierr)
               use const_def, only: dp
               integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
               real(dp), intent(in) :: s ! interpolation x value (between xold and x).
               real(dp), intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn
         irtrn = 0
      end subroutine null_solout
      
      
      subroutine null_dfx(n, x, y, fx, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, y(:) ! (n)
         real(dp), intent(inout) :: fx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         fx = 0
      end subroutine null_dfx
      
      
      ! Newton-Raphson iterative solver for nonlinear systems
      ! square or banded   
      ! analytic or numerical difference jacobian      
      ! where possible, reuses jacobian to improve efficiency      
      ! uses line search method to improve "global" convergence
      include "num_newton.dek"            
      


      ! minimize scalar function of many variables without using derivatives.

      ! NEWUOA
      !     This subroutine seeks the least value of a function of many variables,
      !     by applying a trust region method that forms quadratic models by
      !     interpolation. There is usually some freedom in the interpolation
      !     conditions, which is taken up by minimizing the Frobenius norm of
      !     the change to the second derivative of the model, beginning with the
      !     zero matrix.
      ! by M.J.D. Powell (mjdp@cam.ac.uk)    
      ! M.J.D. Powell, "Developments of NEWUOA for unconstrained minimization without derivatives",
      ! Department of Applied Mathematics and Theoretical Physics, Cambridge, England, report NA05, 2007.  
      include "num_newuoa.dek"


      ! BOBYQA
      !     Similar to NEWUOA, but the values of the variables are constrained 
      !     by upper and lower bounds.
      ! by M.J.D. Powell (mjdp@cam.ac.uk)      
      ! M.J.D. Powell, "The BOBYQA algorithm for bound constrained optimization without derivatives",
      ! Department of Applied Mathematics and Theoretical Physics, Cambridge, England, report NA06, 2009.  
      include "num_bobyqa.dek"
      
      
      ! Nelder-Mead Simplex Method
      !     doesn't use interpolation, so robust with noise data.
      include "num_simplex.dek"
         ! Nelder, J. A. and Mead, R.
         ! "A Simplex Method for Function Minimization."
         ! Comput. J. 7, 308-313, 1965.

      
      
      ! global or local minimum of scalar function of 1 variable
      include "num_brent.dek"


      
      ! QuickSort. ACM Algorithm 402, van Emden, 1970
      ! mesa's implementation from Joseph M. Krahn
      ! http://fortranwiki.org/fortran/show/qsort_inline
      
      subroutine qsort(index,n,vals)
         use mod_qsort, only: sortp_dp
         integer :: index(:), n
         real(dp) :: vals(:)
         call sortp_dp(n,index,vals)
      end subroutine qsort
      
      subroutine qsort_strings(index,n,strings)
         use mod_qsort, only: sortp_string
         integer :: index(:), n
         character(len=*), intent(in) :: strings(:)
         call sortp_string(n,index,strings)
      end subroutine qsort_strings
      
      subroutine qsort_string_index(index,n,string_index,strings)
         use mod_qsort, only: sortp_string_index
         integer :: index(:), n
         integer, intent(in) :: string_index(:) ! (n)
         character(len=*), intent(in) :: strings(:) ! 1..maxval(string_index)
         call sortp_string_index(n,index,string_index,strings)
      end subroutine qsort_string_index
      
      
      ! random numbers
      real(dp) function get_dp_uniform_01(seed)
         ! returns a unit pseudorandom real(dp)
         use mod_random, only: r8_uniform_01
         integer ( kind = 4 ) seed
         get_dp_uniform_01 = r8_uniform_01(seed)
      end function get_dp_uniform_01


      function get_i4_uniform(a, b, seed)
         ! The pseudorandom integer will be scaled to be uniformly distributed
         ! between a and b.
         use mod_random, only: i4_uniform
         integer ( kind = 4 ) a, b, seed, get_i4_uniform
         get_i4_uniform = i4_uniform(a, b, seed)
      end function get_i4_uniform
      
      
      subroutine get_perm_uniform ( n, base, seed, p )
         ! selects a random permutation of n integers
         use mod_random, only: perm_uniform
         integer ( kind = 4 ) n
         integer ( kind = 4 ) base
         integer ( kind = 4 ) p(n)
         integer ( kind = 4 ) seed
         call perm_uniform ( n, base, seed, p )
      end subroutine get_perm_uniform
      
      
      subroutine get_seed_for_random(seed)
         ! returns a seed for the random number generator
         use mod_random, only: get_seed
         integer ( kind = 4 ) seed
         call get_seed(seed)
      end subroutine get_seed_for_random

      
      ! binary search
      include "num_binary_search.dek"

      
      real(dp) function linear_interp(x1, y1, x2, y2, x)
         real(dp), intent(in) :: x1, y1, x2, y2, x
         if (x2 == x1) then
            linear_interp = (y1+y2)/2
         else
            linear_interp = y1 + (y2-y1)*(x-x1)/(x2-x1)
         end if
      end function linear_interp


      real(dp) function find0(xx1, yy1, xx2, yy2) result(x)
         ! find x between xx1 and xx2 s.t. linear_interp(xx1, yy1, xx2, yy2, x) == 0
         real(dp), intent(in) :: xx1, yy1, xx2, yy2
         real(dp) :: a, b
         a = (xx1*yy2)-(xx2*yy1)
         b = yy2-yy1
         if (((abs(a) .ge. abs(b)*1d99) .and. 
     >           ((yy1 .ge. 0d0 .and. yy2 .le. 0d0) .or. (yy1 .le. 0d0 .and. yy2 .ge. 0d0)))) then
            x = 0.5d0*(xx1+xx2)
         else if (b == 0d0) then
            x = 0.5d0*(xx1+xx2)
         else
            x = a/b
         end if
         if (yy1*yy2 <= 0) then ! sanity check
            if (x > max(xx1,xx2)) x = max(xx1,xx2)
            if (x < min(xx1,xx2)) x = min(xx1,xx2)
         end if
      end function find0


      real(dp) function find0_quadratic(xx1, yy1, xx2, yy2, xx3, yy3, ierr) result(x)
         ! find x between xx1 and xx3 s.t. quad_interp(xx1, yy1, xx2, yy2, xx3, yy3, x) == 0
         ! xx2 between xx1 and xx3; yy1 and yy3 different sign; yy2 between yy1 and yy3.
         real(dp), intent(in) :: xx1, yy1, xx2, yy2, xx3, yy3
         integer, intent(out) :: ierr
         real(dp) :: a, b, s2, denom
         ierr = 0; x = 0
         s2 = (xx3**2*(-yy1 + yy2) + xx2**2*(yy1 - yy3) + xx1**2*
     >        (-yy2 + yy3))**2 -  
     >        4*(xx3*(-yy1 + yy2) + xx2*(yy1 - yy3) + xx1*(-yy2 + yy3)) 
     >           *(xx1*xx3*(-xx1 + xx3)*yy2 +  
     >              xx2**2*(xx3*yy1 - xx1*yy3) + xx2*(-xx3**2*yy1 + xx1**2*yy3))
         if (s2 < 0) then
            ierr = -1
            return
         end if
         b = sqrt(s2)
         a = xx3**2*(yy1 - yy2) + xx1**2*(yy2 - yy3) + xx2**2*(-yy1 + yy3)
         denom = 2*(xx3*(yy1 - yy2) + xx1*(yy2 - yy3) + xx2*(-yy1 + yy3))
         x = (a + b)/denom
         if (x > max(xx1,xx2,xx3)) x = (a - b)/denom
         if (x < min(xx1,xx2,xx3) .or. x > max(xx1,xx2,xx3)) ierr = -1
      end function find0_quadratic


      subroutine find_max_quadratic(x1, y1, x2, y2, x3, y3, xmax, ymax, ierr)
         ! x1 < x2 < x3 or x1 > x2 > x3.   y2 > max(y1,y2)
         ! returns max location and value of quadratic fit to the three points
         real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
         real(dp), intent(out) :: xmax, ymax
         integer, intent(out) :: ierr
         real(dp) :: a, b, c, s2, denom, dx1, dx2, dxmax
         ierr = 0; xmax = 0; ymax = 0
         dx1 = x2 - x1
         dx2 = x3 - x2
         ! f(dx) = a + b dx + c dx^2; dx = x - x1
         a = y1
         b = (y2-y1)/dx1 + (y2-y1)/(dx1+dx2) + (dx1/dx2)*(y2-y3)/(dx1+dx2)
         c = (dx2*y1 - dx1*y2 - dx2*y2 + dx1*y3)/(dx1*dx2*(dx1+dx2))
         dxmax = -b/(2d0*c)
         xmax = dxmax + x1
         ymax = a + dxmax*(b + dxmax*c)
         if (ymax < y2) ierr = -1
         if (y2 < max(y1,y2)) ierr = -1
         if (.not. ((x1 < x2 .and. x2 < x3) .or. (x1 > x2 .and. x2 > x3))) ierr = -1
      end subroutine find_max_quadratic
      
      
      subroutine two_piece_linear_coeffs(x, x0, x1, x2, a0, a1, a2, ierr)
         ! interpolation value at x is a0*f(x0) + a1*f(x1) + a2*f(x2)
         real(dp), intent(in) :: x, x0, x1, x2
         real(dp), intent(out) :: a0, a1, a2
         integer, intent(out) :: ierr
         ierr = 0
         if (x0 < x1 .and. x1 < x2) then
            if (x <= x0) then
               a0 = 1; a1 = 0; a2 = 0
            else if (x >= x2) then
               a0 = 0; a1 = 0; a2 = 1
            else if (x <= x1) then
               a1 = min(1d0, max(0d0, (x - x0)/(x1 - x0)))
               a0 = 1 - a1; a2 = 0
            else if (x < x2) then
               a2 = min(1d0, max(0d0, (x - x1)/(x2 - x1))) ! a2 => 1 as x => x2
               a1 = 1 - a2; a0 = 0
            end if
         else if (x0 > x1 .and. x1 > x2) then
            if (x >= x0) then
               a0 = 1; a1 = 0; a2 = 0
            else if (x <= x2) then
               a0 = 0; a1 = 0; a2 = 1
            else if (x >= x1) then
               a1 = min(1d0, max(0d0, (x - x0)/(x1 - x0)))
               a0 = 1 - a1; a2 = 0
            else if (x > x2) then
               a2 = min(1d0, max(0d0, (x - x1)/(x2 - x1))) ! a2 => 1 as x => x2
               a1 = 1 - a2; a0 = 0
            end if
         else
            ierr = -1
         end if
      end subroutine two_piece_linear_coeffs


      real(dp) function integrate(func, minx, maxx, args, atol, rtol, max_steps, ierr)
         procedure(integrand) :: func
         real(dp),intent(in) :: minx,maxx ! Min and max values to integrate over
         real(dp), intent(in) :: args(:) ! Extra args passed to func
         real(dp), intent(in) :: atol,rtol ! Absoulate and relative tolerances
         integer, intent(in) :: max_steps ! Max number of sub-steps
         integer, intent(inout) :: ierr ! Error code

         ierr = 0

         integrate = integrator(func, minx, maxx, args, atol, rtol, max_steps, ierr)

      end function integrate

      end module num_lib

