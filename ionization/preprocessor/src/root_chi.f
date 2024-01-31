!
! root_chi
! francois hebert, jul 18 2010
!
! subroutine to find chi = mu/kt, given ne and kt
!
! uses the brent root-finding algorithm from numerical recipes
! with an additional subroutine for finding the bracketing
! interval. giving the 'guess' parameter is optional but helpful in
! speeding up the bracketing routine
!

      module root_chi

      use def_type, only: dp
      use mod_density

      implicit none

      contains


      real (dp) function rootfind_chi(ne, tau, opt_guess)

         real (dp), intent(in) :: ne, tau
         real (dp), optional, intent(in) :: opt_guess

         integer, parameter :: max_iter = 100
         real (dp), parameter :: eps = epsilon(1.0) ! this may need re-thinking
         real (dp), parameter :: tol = 1.0e-10_dp

         integer :: i
         real (dp) :: guess = 0.0_dp
         real (dp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s
         real (dp) :: tol1, x1, x2, xm


         if (present(opt_guess)) guess = opt_guess

         rootfind_chi = 0.0_dp

         call bracket_chi(ne, tau, guess, x1, fa, x2, fb);

         if (fa*fb > 0) call alert(1, '(root_chi) root is not bracketed')

         a = x1
         b = x2
         c = x2
         d = 0.0_dp
         e = 0.0_dp
         fc = fb

         do, i=1, max_iter
            if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
               c=a
               fc=fa
               d=b-a
               e=d
            end if
            if (abs(fc) < abs(fb)) then
               a=b
               b=c
               c=a
               fa=fb
               fb=fc
               fc=fa
            end if

            tol1 = 2.0_dp*eps*abs(b) + 0.5_dp*tol

            xm = 0.5_dp * (c-b)
            if (abs(xm) <= tol1 .or. fb==0.0) then
               rootfind_chi = b
               return
            end if

            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then 
               !attempt inverse quadratic interpolation
               s = fb/fa 
               if (a == c) then 
                  p = 2.0_dp * xm * s 
                  q = 1.0_dp - s
               else
                  q = fa/fc
                  r = fb/fc
                  p = s * ( 2.0_dp*xm*q*(q-r)-(b-a) * (r-1.0_dp) )
                  q = (q-1.0_dp) * (r-1.0_dp) * (s-1.0_dp)
               end if

               !check whether in bounds
               if (p > 0.0) q = -q 
               p = abs(p)
               if (2.0_dp*p < min( 3.0_dp*xm*q-abs(tol1*q), abs(e*q) )) then
                  !accept interpolation
                  e = d
                  d = p/q 
               else
                  !interpolation failed, use bisection
                  d = xm
                  e = d 
               end if
            else
               d = xm
               e = d
            end if

            !move last best guess to a
            a = b
            fa = fb
            !evaluate new trial root
            b = b + merge(d, sign(tol1, xm), abs(d) > tol1 )

            fb = delta_ne(ne, tau, b)
         end do

         call alert(1, '(root_chi) reached max_iter steps before convergence')




         contains

         real (dp) function delta_ne(ne, tau, chi)
            real (dp), intent(in) :: ne, tau, chi
            delta_ne = ne_plasma(chi, tau) - ne
         end function delta_ne


         subroutine bracket_chi(ne, tau, guess, x1, v1, x2, v2)

            real (dp), intent(in) :: ne, tau, guess
            real (dp), intent(out) :: x1, x2, v1, v2

            real (dp) :: c1, c2, f1, f2

            c1 = guess
            f1 = delta_ne(ne, tau, c1)
            if (f1 >= 0) then
               c2 = guess - 1.0_dp
            else
               c2 = guess + 1.0_dp
            end if
            f2 = delta_ne(ne, tau, c2)

            do
               if (f1*f2 <= 0) exit
               c1 = c2
               f1 = f2
               c2 = 2*c2 - guess
               f2 = delta_ne(ne, tau, c2)
            end do

            if (c1 < c2) then
               x1 = c1
               x2 = c2
               v1 = f1
               v2 = f2
            else
               x1 = c2
               x2 = c1
               v1 = f2
               v2 = f1
            end if

         end subroutine bracket_chi

      end function rootfind_chi


      end module root_chi

   