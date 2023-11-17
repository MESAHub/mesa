module fingering_modes

    implicit none
  
    private
    public :: rfromR
    public :: lamguess
    public :: k2guess
    public :: eq1
    public :: eq2
    public :: fun
    public :: jac
    public :: gaml2max
  
  contains
  
    function rfromR(r0, tau) result(r)
      real(dp), intent(in) :: r0, tau
      real(dp) :: r
  
      r = (r0 - 1.0_dp) / (-1.0_dp + 1.0_dp/tau)
  
    end function rfromR
  
    function lamguess(pr, tau, r0) result(lam)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: lam
      real(dp) :: r
  
      r = rfromR(r0, tau)
      if (r < tau) then
         lam = sqrt(pr) - pr*sqrt(1.0_dp + tau/pr) 
      else
         if (r > 0.5_dp) then
            lam = 2.0_dp*pr*(tau/pr)*((1.0_dp/3.0_dp)*(1.0_dp - r))**1.5_dp &
                 / (1.0_dp - (1.0_dp - r)*(1.0_dp + tau/pr)/3.0_dp)
         else
            lam = sqrt(pr*tau/r) - pr*sqrt(1.0_dp + tau/pr)
         end if
      end if
  
    end function lamguess
  
    function k2guess(pr, tau, r0) result(k2)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: k2
      real(dp) :: r
  
      r = rfromR(r0, tau)
      if (r < tau) then
         k2 = (1.0_dp + tau/pr)**(-0.5_dp) - sqrt(pr) * (1.0_dp + (tau/pr)*(1.0_dp + tau/pr)**(-2.0_dp)) 
      else
         if (r > 0.5_dp) then
            k2 = sqrt((1.0_dp - r)/3.0_dp)
         else
            k2 = (1.0_dp + tau/pr)**(-0.5_dp) - 2.0_dp*sqrt(r*tau/pr) * (1.0_dp + tau/pr)**(-2.5_dp)
         end if
      end if
  
    end function k2guess
  
    function eq1(lam, k2, pr, tau, r0) result(eq)
      real(dp), intent(in) :: lam, k2, pr, tau, r0
      real(dp) :: eq
  
      real(dp) :: b2, b1, b0
  
      b2 = k2*(1.0_dp + pr + tau) 
      b1 = k2**2*(tau*pr + pr + tau) + pr*(1.0_dp - 1.0_dp/r0)
      b0 = k2**3*tau*pr + k2*pr*(tau - 1.0_dp/r0)
  
      eq = lam**3 + b2*lam**2 + b1*lam + b0
  
    end function eq1
  
    function eq2(lam, k2, pr, tau, r0) result(eq)
      real(dp), intent(in) :: lam, k2, pr, tau, r0  
      real(dp) :: eq
  
      real(dp) :: c2, c1, c0
  
      c2 = 1.0_dp + pr + tau
      c1 = 2.0_dp*k2*(tau*pr + tau + pr) 
      c0 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
  
      eq = c2*lam**2 + c1*lam + c0
  
    end function eq2
  
    function fun(x, pr, tau, r0, passk1) result(f)
      real(dp), intent(in) :: x(2), pr, tau, r0
      logical, intent(in), optional :: passk1
      real(dp) :: f(2)
  
      if (present(passk1)) then
         f(1) = eq1(x(1), x(2)**2, pr, tau, r0)  
         f(2) = eq2(x(1), x(2)**2, pr, tau, r0)
      else
         f(1) = eq1(x(1), x(2), pr, tau, r0)
         f(2) = eq2(x(1), x(2), pr, tau, r0)
      end if
  
    end function fun
  
    function jac(x, pr, tau, r0, passk1) result(dfdx)
      real(dp), intent(in) :: x(2), pr, tau, r0
      logical, intent(in), optional :: passk1
      real(dp) :: dfdx(2,2)
  
      real(dp) :: lam, k2
      real(dp) :: b2, db2dk2, b1, db1dk2, b0, db0dk2
      real(dp) :: c2, c1, dc1dk2, c0, dc0dk2
  
      lam = x(1)
      if (present(passk1)) then
         k2 = x(2)**2
      else  
         k2 = x(2)
      end if
  
      b2 = k2*(1.0_dp + pr + tau)
      db2dk2 = 1.0_dp + pr + tau
  
      b1 = k2**2*(tau*pr + pr + tau) + pr*(1.0_dp - 1.0_dp/r0) 
      db1dk2 = 2.0_dp*k2*(tau*pr + pr + tau)
  
      b0 = k2**3*tau*pr + k2*pr*(tau - 1.0_dp/r0)
      db0dk2 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
  
      dfdx(1,1) = 3.0_dp*lam**2 + 2.0_dp*b2*lam + b1
      dfdx(1,2) = lam**2*db2dk2 + lam*db1dk2 + db0dk2
      if (present(passk1)) dfdx(1,2) = dfdx(1,2)*2.0_dp*x(2)
  
      c2 = 1.0_dp + pr + tau
      c1 = 2.0_dp*k2*(tau*pr + tau + pr)
      dc1dk2 = c1/k2
      c0 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
      dc0dk2 = 6.0_dp*k2*tau*pr
  
      dfdx(2,1) = 2.0_dp*c2*lam + c1
      dfdx(2,2) = lam*dc1dk2 + dc0dk2
      if (present(passk1)) dfdx(2,2) = dfdx(2,2)*2.0_dp*x(2)
  
    end function jac
  
    function gaml2max(pr, tau, r0) result(x)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: x(2)
  
      real(dp) :: lam, k2
      real(dp) :: lamg, k2g  
      real(dp) :: f(2)
  
      lamg = lamguess(pr, tau, r0)
      k2g = k2guess(pr, tau, r0)
  
      ! First solve for lam, k2
      call hybrd1(fun, x, f, jac, 2, lamg, k2g, pr, tau, r0) 
  
      if (x(2) < 0.0_dp) then
         ! If k2 is negative, solve for lam, k instead
         lamg = lamguess(pr, tau, r0)
         k2g = sqrt(k2guess(pr, tau, r0))
         call hybrd1(fun, x, f, jac, 2, lamg, k2g, pr, tau, r0, .true.)
         ! Convert k to k2
         x(2) = x(2)**2
      end if
  
      if (any(abs(f) > 1e-12)) then
         print *, 'gaml2max did not converge!'
         print *, f
      end if
  
    contains
  
      ! Simple hybrd1 interface
      subroutine hybrd1(func, x, fvec, fjac, n, xtol, xguess, var1, var2, var3, passk1)
        external :: func, fjac
        integer, intent(in) :: n
        real(dp), intent(in) :: xtol, var1, var2, var3
        real(dp), intent(inout) :: x(n), xguess
        real(dp), intent(out) :: fvec(n)
        logical, intent(in), optional :: passk1
  
        call hybrd1_(func, n, x, fvec, fjac, xtol, maxfev, ml, mu, epsfcn, &
             diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
             wa1, wa2, wa3, wa4, xguess, var1, var2, var3, passk1)
  
      end subroutine hybrd1
  
    end function gaml2max
  
  end module fingering_modes