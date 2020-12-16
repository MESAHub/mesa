   module mod_brent
   use const_def, only: dp
      
   implicit none

   contains

   ! this code was written by John Burkardt
   ! see his wonderful site for lots of numerical software.
   ! here's the link for the implementation that we use here.
   ! http://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.html


   real(dp) function eval_global_min(max_tries, a, b, c, m, machep, e, t, f, x, ierr)
      integer, intent(in) :: max_tries
      real(dp), intent(in) :: a, b, c, m, machep, e, t
      interface
         real(dp) function f(x)
            use const_def, only: dp
            real(dp), intent(in) :: x
         end function f
      end interface
      real(dp), intent(out) :: x
      integer, intent(out) :: ierr

      !  Parameters:
      !
      !    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
      !    It must be the case that A < B.
      !
      !    Input, real ( kind = 8 ) C, an initial guess for the global
      !    minimizer.  If no good guess is known, C = A or B is acceptable.
      !
      !    Input, real ( kind = 8 ) M, the bound on the second derivative.
      !
      !    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
      !    precision.
      !
      !    Input, real ( kind = 8 ) E, a positive tolerance, a bound for the
      !    absolute error in the evaluation of F(X) for any X in [A,B].
      !
      !    Input, real ( kind = 8 ) T, a positive error tolerance.
      !
      !    Input, external real ( kind = 8 ) F, the name of a user-supplied
      !    function, of the form "FUNCTION F ( X )", which evaluates the
      !    function whose global minimum is being sought.
      !
      !    Output, real ( kind = 8 ) X, the estimated value of the abscissa
      !    for which F attains its global minimum value in [A,B].
      !
      !    Output, real ( kind = 8 ) GLOMIN, the value F(X).
      !

        real    ( kind = 8 ) a0
        real    ( kind = 8 ) a2
        real    ( kind = 8 ) a3
        real    ( kind = 8 ) d0
        real    ( kind = 8 ) d1
        real    ( kind = 8 ) d2
        real    ( kind = 8 ) h
        integer ( kind = 4 ) k, iter
        real    ( kind = 8 ) m2
        real    ( kind = 8 ) p
        real    ( kind = 8 ) q
        real    ( kind = 8 ) qs
        real    ( kind = 8 ) r
        real    ( kind = 8 ) s
        real    ( kind = 8 ) sc
        real    ( kind = 8 ) y
        real    ( kind = 8 ) y0
        real    ( kind = 8 ) y1
        real    ( kind = 8 ) y2
        real    ( kind = 8 ) y3
        real    ( kind = 8 ) yb
        real    ( kind = 8 ) z0
        real    ( kind = 8 ) z1
        real    ( kind = 8 ) z2
        
        ierr = 0
        a0 = b
        x = a0
        a2 = a
        y0 = f ( b )
        yb = y0
        y2 = f ( a )
        y = y2

        if ( y0 < y ) then
          y = y0
        else
          x = a
        end if

        if ( m <= 0.0D+00 .or. b <= a ) then
          eval_global_min = y
          return
        end if

        m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

        if ( c <= a .or. b <= c ) then
          sc = 0.5D+00 * ( a + b )
        else
          sc = c
        end if

        y1 = f ( sc )
        k = 3
        d0 = a2 - sc
        h = 9.0D+00 / 11.0D+00

        if ( y1 < y ) then
          x = sc
          y = y1
        end if

        do iter = 1, max_tries+1
       
          if (iter > max_tries) then
            ierr = -1
            exit
          end if

          d1 = a2 - a0
          d2 = sc - a0
          z2 = b - a2
          z0 = y2 - y1
          z1 = y2 - y0
          r = d1 * d1 * z0 - d0 * d0 * z1
          p = r
          qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
          q = qs

          if ( k < 1000000 .or. y2 <= y ) then

            do

              if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
                z2 * m2 * r * ( z2 * q - r ) ) then
                a3 = a2 + r / q
                y3 = f ( a3 )

                if ( y3 < y ) then
                  x = a3
                  y = y3
                end if
              end if

              k = mod ( 1611 * k, 1048576 )
              q = 1.0D+00
              r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

              if ( z2 <= r ) then
                exit
              end if

            end do

          else

            k = mod ( 1611 * k, 1048576 )
            q = 1.0D+00
            r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

            do while ( r < z2 )

              if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
                z2 * m2 * r * ( z2 * q - r ) ) then
                a3 = a2 + r / q
                y3 = f ( a3 )

                if ( y3 < y ) then
                  x = a3
                  y = y3
                end if
              end if

              k = mod ( 1611 * k, 1048576 )
              q = 1.0D+00
              r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

            end do

          end if

          r = m2 * d0 * d1 * d2
          s = sqrt ( ( ( y2 - y ) + t ) / m2 )
          h = 0.5D+00 * ( 1.0D+00 + h )
          p = h * ( p + 2.0D+00 * r * s )
          q = q + 0.5D+00 * qs
          r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )

          if ( r < s .or. d0 < 0.0D+00 ) then
            r = a2 + s
          else
            r = a2 + r
          end if

          if ( 0.0D+00 < p * q ) then
            a3 = a2 + p / q
          else
            a3 = r
          end if

          do

            a3 = max ( a3, r )

            if ( b <= a3 ) then
              a3 = b
              y3 = yb
            else
              y3 = f ( a3 )
            end if

            if ( y3 < y ) then
              x = a3
              y = y3
            end if

            d0 = a3 - a2

            if ( a3 <= r ) then
              exit
            end if

            p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )

            if ( ( 1.0D+00 + 9.0D+00 * machep ) * d0 <= abs ( p ) ) then
              exit
            end if

            if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) <= &
              ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) then
              exit
            end if

            a3 = 0.5D+00 * ( a2 + a3 )
            h = 0.9D+00 * h

          end do

          if ( b <= a3 ) then
            exit
          end if

          a0 = sc
          sc = a2
          a2 = a3
          y0 = y1
          y1 = y2
          y2 = y3

        end do

        eval_global_min = y

        return
   end function eval_global_min


   real(dp) function eval_local_min(max_tries, a, b, eps, t, f, x, ierr)
      integer, intent(in) :: max_tries
      real(dp), intent(in) :: a, b, eps, t
      interface
         real(dp) function f(x)
            use const_def, only: dp
            real(dp), intent(in) :: x
         end function f
      end interface
      real(dp), intent(out) :: x
      integer, intent(out) :: ierr

   !*****************************************************************************80
   !
   !! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
   !
   !  Discussion:
   !
   !    The method used is a combination of golden section search and
   !    successive parabolic interpolation.  Convergence is never much slower
   !    than that for a Fibonacci search.  If F has a continuous second
   !    derivative which is positive at the minimum (which is not at A or
   !    B), then convergence is superlinear, and usually of the order of
   !    about 1.324....
   !
   !    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
   !    F is never evaluated at two points closer than TOL.
   !
   !    If F is a unimodal function and the computed values of F are always
   !    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
   !    LOCAL_MIN approximates the abscissa of the global minimum of F on the
   !    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
   !
   !    If F is not unimodal, then LOCAL_MIN may approximate a local, but
   !    perhaps non-global, minimum to the same accuracy.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    13 April 2008
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Richard Brent
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Richard Brent,
   !    Algorithms for Minimization Without Derivatives,
   !    Dover, 2002,
   !    ISBN: 0-486-41998-3,
   !    LC: QA402.5.B74.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
   !
   !    Input, real ( kind = 8 ) EPS, a positive relative error tolerance.
   !    EPS should be no smaller than twice the relative machine precision,
   !    and preferably not much less than the square root of the relative
   !    machine precision.
   !
   !    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
   !
   !    Input, external real ( kind = 8 ) F, the name of a user-supplied
   !    function, of the form "FUNCTION F ( X )", which evaluates the
   !    function whose local minimum is being sought.
   !
   !    Output, real ( kind = 8 ) X, the estimated value of an abscissa
   !    for which F attains a local minimum value in [A,B].
   !
   !    Output, real ( kind = 8 ) LOCAL_MIN, the value F(X).
   !

     real ( kind = 8 ) c
     real ( kind = 8 ) d
     real ( kind = 8 ) e
     real ( kind = 8 ) fu
     real ( kind = 8 ) fv
     real ( kind = 8 ) fw
     real ( kind = 8 ) fx
     real ( kind = 8 ) m
     real ( kind = 8 ) p
     real ( kind = 8 ) q
     real ( kind = 8 ) r
     real ( kind = 8 ) sa
     real ( kind = 8 ) sb
     real ( kind = 8 ) t2
     real ( kind = 8 ) tol
     real ( kind = 8 ) u
     real ( kind = 8 ) v
     real ( kind = 8 ) w
     integer ( kind = 4 ) iter

     ierr = 0
     
   !
   !  C is the square of the inverse of the golden ratio.
   !
     c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

     sa = a
     sb = b
     x = sa + c * ( sb - sa )
     fx = f ( x )
     w = x
     v = w
     e = 0.0D+00
     fw = fx
     fv = fw
     d = 0

     do iter = 1, max_tries+1

       m = 0.5D+00 * ( sa + sb )
       tol = eps * abs ( x ) + t
       t2 = 2.0D+00 * tol
   !
   !  Check the stopping criterion.
   !
       if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
         exit
       end if
       
       if (iter > max_tries) then
         ierr = -1
         exit
       end if
   !
   !  Fit a parabola.
   !
       r = 0.0D+00
       q = r
       p = q

       if ( tol < abs ( e ) ) then

         r = ( x - w ) * ( fx - fv )
         q = ( x - v ) * ( fx - fw )
         p = ( x - v ) * q - ( x - w ) * r
         q = 2.0D+00 * ( q - r )

         if ( 0.0D+00 < q ) then
           p = - p
         end if

         q = abs ( q )

         r = e
         e = d

       end if

       if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
            q * ( sa - x ) < p .and. &
            p < q * ( sb - x ) ) then
   !
   !  Take the parabolic interpolation step.
   !
         d = p / q
         u = x + d
   !
   !  F must not be evaluated too close to A or B.
   !
         if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

           if ( x < m ) then
             d = tol
           else
             d = - tol
           end if

         end if
   !
   !  A golden-section step.
   !
       else

         if ( x < m ) then
           e = sb - x
         else
           e = sa - x
         end if

         d = c * e

       end if
   !
   !  F must not be evaluated too close to X.
   !
       if ( tol <= abs ( d ) ) then
         u = x + d
       else if ( 0.0D+00 < d ) then
         u = x + tol
       else
         u = x - tol
       end if

       fu = f ( u )
   !
   !  Update A, B, V, W, and X.
   !
       if ( fu <= fx ) then

         if ( u < x ) then
           sb = x
         else
           sa = x
         end if

         v = w
         fv = fw
         w = x
         fw = fx
         x = u
         fx = fu

       else

         if ( u < x ) then
           sa = u
         else
           sb = u
         end if

         if ( fu <= fw .or. w == x ) then
           v = w
           fv = fw
           w = u
           fw = fu
         else if ( fu <= fv .or. v == x .or. v == w ) then
           v = u
           fv = fu
         end if

       end if

     end do

     eval_local_min = fx

     return
   end function eval_local_min




   real*8 function eval_brent_safe_zero ( a, b, machep, t, epsy, f, fa_in, fb_in, lrpar, rpar, lipar, ipar, ierr )

   !*****************************************************************************80
   !
   !! seeks the root of a function F(X) in an interval [A,B].
   !
   !  Discussion:
   !
   !    The interval [A,B] must be a change of sign interval for F.
   !    That is, F(A) and F(B) must be of opposite signs.  Then
   !    assuming that F is continuous implies the existence of at least
   !    one value C between A and B for which F(C) = 0.
   !
   !    The location of the zero is determined to within an accuracy
   !    of 6 * MACHEPS * abs ( C ) + 2 * T.   or if abs(f(x)) < epsy.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license. 
   !
   !  Modified:
   !
   !    12 April 2008
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Richard Brent.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Richard Brent,
   !    Algorithms for Minimization Without Derivatives,
   !    Dover, 2002,
   !    ISBN: 0-486-41998-3,
   !    LC: QA402.5.B74.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
   !    sign interval.
   !
   !    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
   !    precision.
   !
   !    Input, real ( kind = 8 ) T, a positive error tolerance.
   !
   !    Input, external real ( kind = 8 ) F, the name of a user-supplied
   !    function, of the form "FUNCTION F ( X )", which evaluates the
   !    function whose zero is being sought.
   !
   !    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
   !    the function F.
   !
     implicit none
     
     interface
#include "num_root_fcn.dek" 
     end interface
     integer, intent(in) :: lipar, lrpar
     integer, intent(inout), pointer :: ipar(:) ! (lipar)
     real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
     integer, intent(out) :: ierr

     real ( kind = 8 ) a
     real ( kind = 8 ) b
     real ( kind = 8 ) c
     real ( kind = 8 ) d
     real ( kind = 8 ) e
     real ( kind = 8 ) fa, fa_in
     real ( kind = 8 ) fb, fb_in
     real ( kind = 8 ) fc
     real ( kind = 8 ) m
     real ( kind = 8 ) machep
     real ( kind = 8 ) p
     real ( kind = 8 ) q
     real ( kind = 8 ) r
     real ( kind = 8 ) s
     real ( kind = 8 ) sa
     real ( kind = 8 ) sb
     real ( kind = 8 ) t, epsy
     real ( kind = 8 ) tol
     real ( kind = 8 ) dfdx
     
     ierr = 0
     
   !
   !  Make local copies of A and B.
   !
     sa = a
     sb = b
     fa = fa_in
     !fa = f ( sa, dfdx, lrpar, rpar, lipar, ipar, ierr )
     !if (ierr /= 0) return
     fb = fb_in
     !fb = f ( sb, dfdx, lrpar, rpar, lipar, ipar, ierr )
     !if (ierr /= 0) return

     c = sa
     fc = fa
     e = sb - sa
     d = e

     do

       if ( abs ( fc ) < abs ( fb ) ) then

         sa = sb
         sb = c
         c = sa
         fa = fb
         fb = fc
         fc = fa

       end if

       tol = 2.0D+00 * machep * abs ( sb ) + t
       m = 0.5D+00 * ( c - sb )

       if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
         exit
       end if

       if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

         e = m
         d = e

       else

         s = fb / fa

         if ( sa == c ) then

           p = 2.0D+00 * m * s
           q = 1.0D+00 - s

         else

           q = fa / fc
           r = fb / fc
           p = s * ( 2.0D+00 * m * a * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
           q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

         end if

         if ( 0.0D+00 < p ) then
           q = - q
         else
           p = - p
         end if

         s = e
         e = d

         if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
           p < abs ( 0.5D+00 * s * q ) ) then
           d = p / q
         else
           e = m
           d = e
         end if

       end if

       sa = sb
       fa = fb

       if ( tol < abs ( d ) ) then
         sb = sb + d
       else if ( 0.0D+00 < m ) then
         sb = sb + tol
       else
         sb = sb - tol
       end if

       fb = f ( sb, dfdx, lrpar, rpar, lipar, ipar, ierr )
       if (ierr /= 0) return
       if (abs(fb) < epsy) exit

       if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
            ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
         c = sa
         fc = fa
         e = sb - sa
         d = e
       end if

     end do

     eval_brent_safe_zero = sb

     return
   end function eval_brent_safe_zero



   end module mod_brent

