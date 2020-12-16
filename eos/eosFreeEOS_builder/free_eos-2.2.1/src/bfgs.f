C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 2006 Alan W. Irwin
C
C       $Id: bfgs.f 505 2007-07-08 20:08:51Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*************************************************************

C      BFGS algorithm based on 
C      R. Fletcher, "Practical Methods of Optimization" 2nd ed.
C      John Wiley and Sons, 1987.
C      N.B. in following comments this reference is referred to as PMOO.

      subroutine bfgs_constant_vector(n, constant, x)
C      set vector x(n) to constant
      implicit none
      integer n,i
      double precision constant, x(n)
      do i = 1,n
        x(i) = constant
      enddo
      end

      subroutine bfgs_set(status, name,
     &  n, x, xi, f, gradient, p)
C      status is a control string that takes the following values
C      ('with a "_set" suffix to identify for logic debugging purposes
C      exactly where the status was set):
C        'both'.  Request both function value and gradient to be supplied
C                 externally.
C        'complete'. Initialization has been completed.

C      initialize line search
      implicit none
      include 'bfgs.h'
      integer n
      double precision p(n), x(n), xi, f, gradient
      character*(*) status, name
      save
C      sanity checks:
      if(len(status).lt.12)
     &  stop 'bfgs_set: length of status string too small'
      if(len(name).lt.3)
     &  stop 'bfgs_set:length of name string too small'
      if(status(:4).ne.'both') then
        call dcopy(n, x, 1, xi, 1)
        status = 'both_set'
        name = 'set'
      else
C        initial preferred direction for the line search is in the
C        direction of the gradient, i.e, steepest ascent.  (Negative
C        sign of direction, i.e., initial steepest descent
C        figured out later.)
        call dcopy(n, gradient, 1, p, 1)
C        Set to true for Fletcher estimate of alpha1 and set to false for
C        GSL estimate of alpha1 (or Numerical Recipes) estimate of alpha1.
        fletcher_estimate = .true.
        if(fletcher_estimate) then
C         _use_large initial value of deltaf_previous to force an initial
C          step size of unity as per PMOO, 2.6.8.
          deltaf_previous = 1.d300
        else
C          set alpha_previous to 1.d0 to force an initial step size
C          of unity for the GSL (or Numerical Recipes) estimate.
          alpha_previous = 1.d0
        endif
        status = 'complete_set'
      endif
      end
      
      subroutine bfgs_take_step(n, x, p, step, lambda, x1, dx)
C      take step in the desired -lambda*p direction.
C      unchanged input:
C      x(n) is the current position vector
C      lambda*p(n) is the unit vector in the currently desired uphill
C      (not necessarily steepest ascent) direction.
C      output results:
C      dx(n) = -step*lambda*p
C      x1(n) = x + dx
      implicit none
      integer n
      double precision x(n), p(n), step, lambda, x1(n), dx(1)
      call bfgs_constant_vector(n, 0.d0, dx)
      call daxpy (n, -step*lambda, p, 1, dx, 1)
      call dcopy(n, x, 1, x1, 1)
      call daxpy (n, 1.d0, dx, 1, x1, 1)
      end

      double precision function bfgs_cubic_minimum(
     &  a, fa, fprimea, b, fb, fprimeb, x_lo, x_hi)
C     _Use_cubic interpolation or extrapolation to calculate
C      minimum value within a specified range.

C      input variables:
C      a is first defined point
C      fa and fprimea are the function and derivative values at that point.
C      b is the second defined point
C      fb and fprimeb are the function and derivative values at that point.
C      x_lo through x_hi specify the range where the minimum value
C      is to be calculated.

C      The returned function value is the independent variable x where
C      the minimum of the interpolated/extrapolated function occurs
C      within [x_lo, x_hi].

      implicit none
      double precision a, fa, fprimea, b, fb, fprimeb, x_lo, x_hi
      double precision dx, q0, q1, q2, q3, z_lo, z_hi,
     &  q_lo, q_hi, z_stationary, q_stationary, minimum,
     &  radical, quadratic_factor
      integer nsolution

      dx = b - a
C      Find coefficients of Hermite polynomial on z which ranges from
C      0 to 1 within defined interval [a, b], which agrees with fa and
C      fprimea*dx at 0 and which agrees with fb and fprimeb*dx at 1.
C      z = (x - a)/dx or x = a + z*dx
C      f(z) = q0 + z*(q1 + z*(q2 + z*q3))
C      f'(z) = q1 + 2.d0*q2*z + 3.d0*q3*z*z
      q0 = fa
      q1 = fprimea*dx
      q2 = 3.d0*(fb-fa) - (2.d0*fprimea + fprimeb)*dx
      q3 = (fprimea+fprimeb)*dx - 2.d0*(fb-fa)
C      N.B.
C      f(0) = q0 = fa
C      f(1) = q0 + q1 + q2 + q3 = fb
C      f'(0) = q1 = fprimea*dx
C      f'(1) = q1 + 2*q2 + 3*q3 = ((1-4+3)*fprimea + (-2+3)*fprimeb)*dx
C            = fprimeb*dx
C      z_lo and z_hi correspond to x_lo and x_hi.
      z_lo = (x_lo - a)/dx
      z_hi = (x_hi - a)/dx
      q_lo = q0 + z_lo*(q1 + z_lo*(q2 + z_lo*q3))
      q_hi = q0 + z_hi*(q1 + z_hi*(q2 + z_hi*q3))
      if(q_lo.lt.q_hi) then
        bfgs_cubic_minimum = x_lo
        minimum = q_lo
      else
        bfgs_cubic_minimum = x_hi
        minimum = q_hi
      endif
      if(q3.eq.0.d0) then
        nsolution = 1
        z_stationary = -q1*0.5d0/q2
      else
C        radical = b^2 - 4 a c
        radical = 4.d0*q2*q2 - 12.d0*q1*q3
        if(radical.lt.0.d0) then
C          no roots
          nsolution = 0
        elseif(radical.eq.0.d0) then
C          one root.
          nsolution = 1
          z_stationary = - q2/(3.d0*q3)
        else
C          two roots.
          nsolution = 2
          if(q2.gt.0) then
            quadratic_factor = 2.d0*q2 + sqrt(radical)
          else
            quadratic_factor = 2.d0*q2 - sqrt(radical)
          endif
          z_stationary = -quadratic_factor/(6.d0*q3)
        endif
      endif
      if(nsolution.ge.1) then
C        at least one possible solution calculated at z_stationary
        if((z_lo.le.z_stationary.and.z_stationary.le.z_hi).or.
     &      (z_hi.le.z_stationary.and.z_stationary.le.z_lo)) then
          q_stationary = q0 + z_stationary*(q1 + z_stationary*
     &      (q2 + z_stationary*q3))
          if(q_stationary.lt.minimum) then
            bfgs_cubic_minimum = a + z_stationary*dx
            minimum = q_stationary
          endif
        endif
        if(nsolution.eq.2) then
C          try alternate quadratic solution
          z_stationary = -2.d0*q1/quadratic_factor
          if((z_lo.le.z_stationary.and.z_stationary.le.z_hi).or.
     &        (z_hi.le.z_stationary.and.z_stationary.le.z_lo)) then
            q_stationary = q0 + z_stationary*(q1 + z_stationary*
     &        (q2 + z_stationary*q3))
            if(q_stationary.lt.minimum) then
              bfgs_cubic_minimum = a + z_stationary*dx
            endif
          endif
        endif
      endif
      end

      double precision function bfgs_quadratic_minimum(
     &  a, fa, fprimea, b, fb, x_lo, x_hi)
C     _Use_quadratic interpolation or extrapolation to calculate
C      minimum value within a specified range.

C      input variables:
C      a is first defined point
C      fa and fprimea are the function and derivative values at that point.
C      b is the second defined point
C      fb is the function value at that point.
C      x_lo through x_hi specify the range where the minimum value
C      is to be calculated.

C      The returned function value is the independent variable x where
C      the minimum of the interpolated/extrapolated function occurs
C      within [x_lo, x_hi].

      implicit none
      double precision a, fa, fprimea, b, fb, x_lo, x_hi
      double precision dx, q0, q1, q2, z_lo, z_hi,
     &  q_lo, q_hi, z_stationary, q_stationary, minimum
      integer nsolution

      dx = b - a
C      Find coefficients of Hermite polynomial on z which ranges from
C      0 to 1 within defined interval [a, b], which agrees with fa and
C      fprimea*dx at 0 and which agrees with fb at 1.
C      z = (x - a)/dx or x = a + z*dx
C      f(z) = q0 + z*(q1 + z*q2)
C      f'(z) = q1 + 2.d0*q2*z with stationary point z = -q1*0.5/q2
      q0 = fa
      q1 = fprimea*dx
      q2 = fb-fa-q1
C      z_lo and z_hi correspond to x_lo and x_hi.
      z_lo = (x_lo - a)/dx
      z_hi = (x_hi - a)/dx
      q_lo = q0 + z_lo*(q1 + z_lo*q2)
      q_hi = q0 + z_hi*(q1 + z_hi*q2)
      if(q_lo.lt.q_hi) then
        bfgs_quadratic_minimum = x_lo
        minimum = q_lo
      else
        bfgs_quadratic_minimum = x_hi
        minimum = q_hi
      endif
      if(q2.eq.0.d0) then
        nsolution = 0
      else
        nsolution = 1
        z_stationary = -q1*0.5d0/q2
      endif
      if(nsolution.eq.1) then
        if((z_lo.le.z_stationary.and.z_stationary.le.z_hi).or.
     &      (z_hi.le.z_stationary.and.z_stationary.le.z_lo)) then
          q_stationary = q0 + z_stationary*(q1 + z_stationary*q2)
          if(q_stationary.lt.minimum) then
            bfgs_quadratic_minimum = a + z_stationary*dx
          endif
        endif
      endif
      end

      subroutine bfgs_linesearch(status, name,
     &  fbar, epsilon, n, p, x, xi, f, gradient, alpha_ratio, dx)
C      Do a linesearch using algorithm presented by PMOO, pp. 34-39

C      status is a control string that takes the following values
C      ('function', 'gradient', 'both', and 'continue' have
C      "_searchN" suffixes to identify for logic debugging purposes
C      exactly where the staus was set):
C        'function'.  Request function value to be supplied externally
C        'gradient'.  Request gradient to be supplied externally.
C           N.B. gradient evaluation is always done at xi of prior
C           function evaluation
C        'both'.  Request both function value and gradient to be supplied
C                 externally.
C        'continue'.  bfgs update after return rather than returning
C                     immediately from bfgs_iterate.
C        'error positive'.  positive initial derivative in direction of
C                           line search.  no further progress can be made
C                           so return immediately.
C                           ('positive' currently unused and bfgs_linesearch
C                           stops instead if it runs into this condition which
C                           should be impossible.)
C        'error roundoff'.  line search did not pass roundoff error check
C                           discussed on page 38 of PMOO.  No further
C                           progress can be made so return immediately
C                           after line search.

C      name is a control string set to 'linesearch'
C      when there is a 'function', gradient', or 'both' request or
C      when there is an error.

C      input quantities:
C      fbar is a user-specified minimum possible f used as per
C        PMOO, 2.6.1 to control the bracketing step size.  If actual f
C        values are <= fbar, then this routine stops with an error so
C        be realistic in how you specify fbar making it small enough
C        to avoid the error, but large enough to provide some control
C        over the maximum size of the bracketing step.
C      epsilon is the user-specified convergence criterion for f (not
C        relative f) discussed on p. 38 of PMOO.  It is used to terminate
C        the line-search routine to work around the case when round-off
C        errors are causing convergence problems close to the minimum.
C      p(n) is the unscaled line search direction.  The PMOO scaled s vector
C        is defined by
C        s = -(dir/pnorm)*p
C        where dir = 1 if p has an acute angle with the input gradient and
C        dir = -1 otherwise.  pnorm is the norm of p.  dir and pnorm are
C        calculated internally.

C      input and output quantities:
C      x(n) is the starting point of the line search on input and on output
C        is the ending point of the line search.  N.B. x is an acceptable
C        point only if status(:8) is 'continue'.
C      N.B. both f and gradient must be precalculated on input and calculated
C        on output.  That is:
C      f is f(x(n)) on both input and output.
C      gradient(n) is the gradient(x(n)) on both input and output.

C      output quantities:
C      alpha_ratio is the ratio of the calculated step size to the initial
C      estimate, alpha1 for the step size.
C      dx(n) is the vector of differences between the initial x and final x.
C      xi(n) is the value of x that is passed back for the purposes of
C      external test calculations of f and its gradient.
      implicit none
      include 'bfgs.h'
      integer n
      double precision fbar, epsilon,
     &  p(n), x(n), xi(n), f, gradient(n), alpha_ratio, dx(n)
C      dimensionless line search constants recommended by PMOO, page. 37
C      (superseded by PMOO, page 69 for tau2 = 0.05 rather than tau2 = 0.10).
      double precision sigma, rho, tau1, tau2, tau3
C     _use_fairly accurate line search.
      parameter(sigma = 0.1d0)
C      rho must be less than or equal to sigma
      parameter(rho = 0.01d0)
C      jump size factor increase used in bracketing phase.
      parameter(tau1 = 9.d0)
C      0 < tau2 < tau3 <= 0.5 and tau2 <= sigma is advisable
      parameter(tau2 = 0.05d0)
      parameter(tau3 = 0.5d0)
      logical bracket
      double precision alpha1, pnorm, dnrm2, pg, ddot, dir,
     &  f0, fprime0, mu,
     &  alphaim, fim, fprimeim, alphai, fi, fprimei,
     &  ai, fai, fprimeai, bi, fbi, fprimebi,
     &  dalpha, alphaip, bfgs_quadratic_minimum, bfgs_cubic_minimum
      logical acceptable, fprimebi_calc, ifcubic
C     _Use_cubic interpolation/extrapolation if fprimebi_calc is true
      parameter(ifcubic=.true.)
      double precision aip, bip, fmin
      logical debug
      parameter(debug=.false.)
      character*(*) status, name
C      N.B. the *need logic replaces go to logic with typically three or four
C      if statements. However, I decided to_use_the more complicated *need
C      approach because the g77 compiler generated all sorts of warnings
C      about the go to's into the middle of control blocks leading
C      to uncertainty and doubt about whether the go to logic would actually
C      work correctly for most compilers/optimizers.

C      N.B. the various need variables and the not_need variables derived
C      from them are used to control re-entry into the routine with
C      the 'function', 'gradient', or 'both' data that are needed.
C      to understand the ordinary flow of this routine simply ignore
C      how the need and not_need variables change that flow.
C      To understand the re-entry flow of control, note that only
C      at most one of need1,...,need5 are true at any given time, and
C      for the case when one of those are true, you should be able to
C      swiftly follow the logic back to where the data are needed and
C      continue from there.
      logical need1, need2, need3, need4, need5,
     &  need12, need345,
     &  not_need12345, not_need345, not_need12, not_need2,
     &  not_need45, not_need3, not_need5
      data need1, need2, need3, need4, need5/5*.false./
      save
!      if(need1) then
!        go to 1
!      elseif(need2) then
!        go to 2
!      elseif(need3) then
!        go to 3
!      elseif(need4) then
!        go to 4
!      elseif(need5) then
!        go to 5
!      endif
      need12 = need1.or.need2
      need345 = need3.or.need4.or.need5
      not_need12345 = .not.(need12.or.need345)
      not_need345 = .not.need345
      not_need12 = .not.need12
      not_need2 = .not.need2
      not_need45 = .not.(need4.or.need5)
      not_need3 = .not.need3
      not_need5 = .not.need5
      if(not_need12345) then
        
C        sanity check
        if(n.gt.nmax_bfgs) stop 'bfgs_linesearch: n too large'
        if(len(status).lt.16)
     &    stop 'bfgs_linesearch: length of status string too small'
        if(len(name).lt.11)
     &    stop 'bfgs_linesearch: length of name string too small'
        
C        parameters of PMOO scaled s vector where s = -(dir/pnorm)*p
        pnorm = dnrm2(n, p, 1)
C        pg subsequently used to calculate fprime0
        pg = ddot(n, p, 1, gradient, 1)
C        according to documentation of sign under 'info g77',
C        dir will be 1.d0 if pg.ge.0.d0 and -1.d0 otherwise.
        dir = sign(1.d0, pg)

C        early convergence test to avoid divide by zero.
C        Note, pnorm can be zero, if p is initialized to the gradient and
C        the gradient is zero.
        if(pnorm.eq.0.d0) then
          call bfgs_constant_vector(n, 0.d0, dx)
          alpha_ratio = 0.d0
          alpha_previous = 0.d0
          deltaf_previous = 0.d0
          status = 'continue_search1'
          return
        endif
C        Initialize "0" variables:
        f0 = f
C        fprime = s dot gradient (PMOO, 1.2.6)
        fprime0 = -(dir*pg)/pnorm

        if(fprime0.ge.0.d0) then
          call bfgs_constant_vector(n, 0.d0, dx)
          alpha_ratio = 0.d0
          alpha_previous = 0.d0
          deltaf_previous = 0.d0
          if(fprime0.eq.0.d0) then
C            another early convergence test (on fprime0)
C            to avoid a divide by zero
            status = 'continue_search2'
          else
C            positive fprime0 should not be possible since dir*pg should
C            always be positive or zero and similarly for pnorm (see how
C            these quantities calculated above), and the zero cases have
C            already been intercepted.
            status = 'error positive'
C            instead of returning a status simply stop for now.
            stop
     &        'bfgs_linesearch: fprime0 > 0. should not be possible'
          endif
          return
        endif
        if(f0.le.fbar) then
          write(0,*) 'bfgs_linesearch: ERROR f0 = f(alpha=0) '//
     &      '<= fbar, the minimum possible function value.'
          write(0,*) 'bfgs_linesearch: respecify fbar and try again.'
          stop
        else
C          maximum line search range from PMOO, 2.6.1
C          Note from above test this must always be positive
C          (unless underflow zeroes it)
          mu = (fbar - f0)/(rho*fprime0)
        endif
        
        if(fletcher_estimate) then
          if(2.d0*max(deltaf_previous,10.d0*epsilon).gt.
     &        -fprime0) then
C            -2 deltaf_previous/fprime0 will generally be close to
C            norm(x).  Thus, this branch generally taken when
C            norm(x)> 1.  alpha1 = 1 implies a relative change in x
C            of 1/norm(x).
C            According to PMOO discussion this gives rapid ultimate
C            convergence for Newton-like methods.
C            PMOO, 2.6.8
            alpha1 = 1.d0
          else
C            -2 deltaf_previous/fprime0 will generally be close to
C            norm(x).  Thus, this branch generally taken when
C            norm(x) < 1.  Thus, the alpha1 extimate below (if we ignore the
C            10.d0*epsilon roundoff error safeguard) implies
C            the relative change in x is generally close to unity for
C            this branch.
C            PMOO, 2.6.7.
            alpha1 = 2.d0*max(deltaf_previous,10.d0*epsilon)/
     &        (-fprime0)
          endif
        else
C          Estimate of initial step size based on best alpha from previous
C          iteration.  This is the GSL estimate.
C          In final super-linear convergence stage this
C          should always be an overestimate, but quadratic interpolation
C          should save the day to get an immediate reasonable minimum for
C          small enough alpha_previous
          alpha1 = alpha_previous
C          Numerical Recipes estimate
C          According to PMOO discussion this gives rapid ultimate convergence
C          for Newton-like methods.
          !alpha1 = 1.d0
        endif
        if(alpha1.le.0.d0) then
          write(0,*) 'alpha1 = ', alpha1
          stop 'bfgs_linesearch: bad initial estimate of step size'
        endif
C        Initialize "i-1" iteration variables (which have "im" suffix).
        alphaim = 0.d0
        fim = f0
        fprimeim = fprime0
C        note condition that fim <= f0+0.*rho*fprime0 is
C        automatically satisfied.
        fmin = fim
        if(debug) then
          write(0,*) "alpha, f(alpha) = ", alphaim, fim
          write(0,*) "alpha, f'(alpha) = ", alphaim, fprimeim
        endif
        
C        Initialize "i" iteration variables (which have "i" suffix).
        alphai = alpha1
        
C        force at least one bracketing attempt.
        bracket = .false.
      endif  !not_need12345
      if(not_need345) then
C PMOO, 2.6.2. LINE 1
        do while (need12.or..not.bracket)
          if(not_need12) then
C            All code in this loop follows pseudo-code in PMOO, 2.6.2.
            
C PMOO, 2.6.2. LINE 2
C            Compute new trial point at alphai corresponding to
C            xi = x + alphai * s
            call bfgs_take_step(n, x, p, alphai, dir/pnorm, xi, dx)
C            EVALUATE f(alphai) = function at xi
            status = 'function_search1'
            name = 'linesearch'
            need1 = .true.
            return
          endif  !not_need12
          if(need1) then
! 1          continue
            need1 = .false.
            need12 = need1.or.need2
            not_need12345 = .not.(need12.or.need345)
            not_need12 = .not.need12
            fi = f
            if(fi.le.f0+alphai*rho*fprime0) fmin = min(fi, fmin)
            if(debug) then
              write(0,*) "alpha, f(alpha) = ", alphai, fi
            endif
C PMOO, 2.6.2. LINE 3
            if(fi.le.fbar) then
              write(0,*) 'bfgs_linesearch: ERROR fi = f(alphai) '//
     &          '<= fbar, the minimum possible function value.'
              write(0,*) 'bfgs_linesearch: respecify fbar and '//
     &          'try again.'
              stop
            endif
C            N.B. in my copy of PMOO, there is a misprint in the 2.6.2
C            pseudo-code where the rho factor in the next uncommented line is
C            ignored.  But the rho factor must be there if the conditions
C            given by 2.6.3 are to be satisfied at the end of the bracket
C            loop.
          endif  !need1
C PMOO, 2.6.2. LINE 4 (with rho factor inserted as per e-mail with Fletcher).
          if(not_need2.and.
     &        (fi.gt.f0+alphai*rho*fprime0.or.fi.ge.fim)) then
C PMOO, 2.6.2. LINE 5
            ai = alphaim
            fai = fim
            fprimeai = fprimeim
            bi = alphai
            fbi = fi
C            At this stage in bracketing iteration gradient for xi
C            is _not_ known.
            fprimebi_calc = .false.
            bracket = .true.
          else
            if(not_need2) then
C PMOO, 2.6.2. LINE 6
C              EVALUATE gradient of function at same xi as previous.
              status = 'gradient_search1'
              name = 'linesearch'
              need2 = .true.
              return
            endif  !not_need2
            if(need2) then
! 2            continue
              need2 = .false.
              need12 = need1.or.need2
              not_need12345 = .not.(need12.or.need345)
              not_need12 = .not.need12
              not_need2 = .not.need2
C              from PMOO, 1.2.6
C              f'(alphai) = s dot gradient
              fprimei = (-dir/pnorm)*ddot(n, p, 1, gradient, 1)
              if(debug) then
                write(0,*) "alpha, f'(alpha) = ", alphai, fprimei
              endif
C PMOO, 2.6.2. LINE 7
C              good converged solution to line search if xi satisfies
C              two-sided test (PMOO, equation 2.5.6)
              if(abs(fprimei).le.-sigma*fprime0) then
C                acceptable point found.
C                output x = xi
                call dcopy(n, xi, 1, x, 1)
C                output f = f(xi)
                f = fi
C                gradient at xi already stored ready for output.
                alpha_ratio = alphai/alpha1
                alpha_previous = alphai
                deltaf_previous = f0-fi
                status = 'continue_search3'
                return
              endif
C PMOO, 2.6.2. LINE 8
              if(fprimei.ge.0.d0) then
C PMOO, 2.6.2. LINE 9
                ai = alphai
                fai = fi
                fprimeai = fprimei
                bi = alphaim
C                fprimeim always known.
                fprimebi_calc = .true.
                fprimebi = fprimeim
                fbi = fim
                bracket = .true.
              else
                dalpha = alphai - alphaim
C PMOO, 2.6.2. LINE 10
                if(mu.le.alphai + dalpha) then
C PMOO, 2.6.2. LINE 11
                  alphaip = mu
                else
C PMOO, 2.6.2. LINE 12
C                  extrapolation to find possible bracket in the range
C                  [alphai + dalpha, min(mu, alphai + tau*dalpha)]
                  if(ifcubic) then
                    alphaip = bfgs_cubic_minimum(
     &                alphaim, fim, fprimeim, alphai, fi, fprimei,
     &                alphai + dalpha, min(mu, alphai + tau1*dalpha))
                  else
                    alphaip = bfgs_quadratic_minimum(
     &                alphaim, fim, fprimeim, alphai, fi,
     &                alphai + dalpha, min(mu, alphai + tau1*dalpha))
                  endif
                endif
C                prepare for the next bracket attempt
                alphaim = alphai
                fim = fi
                fprimeim = fprimei
                alphai = alphaip
              endif
            else
              stop 'bfgs_linesearch: bad logic1'
            endif  !need2
          endif
        enddo
C        found bracket where the following conditions occur:
C        (i) ai is the current best trial point (least f) that
C          satisfies PMOO, 2.5.1, i.e., fai = f(ai) <= f(0) + ai*rho*f'(0)
C        (ii) fprimeai = f'(ai) has been evaluated and satisfies
C          (bi-ai)*f'(ai) <0 but |f'(ai)| > -sigma*f'(0)
C        (iii) bi satisfies either fbi = f(bi) > f(0) + bi*rho*f'(0) or
C          fbi = f(bi) >= f(ai) or both conditions are true.
C        PMOO, Lemma 2.6.1:
C        if sigma >= rho then such [ai, bi] brackets contain an interval of
C        acceptable alpha points such that
C        (i) f(alpha) <= f(0) + alpha*rho*f'(0)  PMOO, 2.5.1.
C        (ii) |f'(alpha)'| <= -sigma*f'(0) PMOO, 2.5.6 (two sided condition).

C       _Use_sectioning (while preserving above properties) to find
C        an acceptable alpha
        acceptable = .false.
      endif  !not_need345
C PMOO, 2.6.4. LINE 1
      do while(need345.or..not.acceptable)
        if(not_need345) then
          if(debug) then
            write(0,*) 'ai, fai, fprimeai = ', ai, fai, fprimeai
            write(0,*) 'bi, fbi = ', bi, fbi
            write(0,*) 'fai.eq.fmin', fai.eq.fmin
            write(0,*) 'fai.le.f0+ai*rho*fprime0',
     &        fai.le.f0+ai*rho*fprime0
            write(0,*) '(bi-ai)*fprimeai.lt.0.d0.and.'//
     &        'abs(fprimeai).gt.-sigma*fprime0 = ',
     &        (bi-ai)*fprimeai.lt.0.d0.and.
     &        abs(fprimeai).gt.-sigma*fprime0
            write(0,*) 'fbi.gt.f0+bi*rho*fprime0.or.fbi.ge.fai = ',
     &        fbi.gt.f0+bi*rho*fprime0.or.fbi.ge.fai
          endif
C          test assertion given by PMOO, 2.6.3
          if(.not.(
     &      (fai.eq.fmin).and.
     &      (fai.le.f0+ai*rho*fprime0).and.
     &      ((bi-ai)*fprimeai.lt.0.d0).and.
     &      (abs(fprimeai).gt.-sigma*fprime0).and.
     &      (fbi.gt.f0+bi*rho*fprime0.or.fbi.ge.fai)))
     &      stop 'bfgs_linesearch: internal bracketing logic error'
          dalpha = bi - ai
C PMOO, 2.6.4. LINE 2
C          interpolation to find acceptable point in the
C          range [ai + tau2*dalpha, bi - tau3*dalpha].
          if(ifcubic.and.fprimebi_calc) then
            alphai = bfgs_cubic_minimum(
     &        ai, fai, fprimeai, bi, fbi, fprimebi,
     &        ai + tau2*dalpha, bi - tau3*dalpha)
          else
            alphai = bfgs_quadratic_minimum(
     &        ai, fai, fprimeai, bi, fbi,
     &        ai + tau2*dalpha, bi - tau3*dalpha)
          endif
        endif  !not_need345
        if(not_need45) then
C          terminate using epsilon condition on f'
C          to avoid round-off error problems (see discussion on
C          p. 38 of PMOO).  Because assertion (bi-ai)*fprimeai.lt.0.d0
C          above is tested correct, and alphai in
C          [ai + tau2*dalpha, bi - tau3*dalpha], (alphai-ai) must be
C          the same sign as (bi-ai) and the left-hand side
C          of the following condition must always be positive
C          or zero (in the case of underflows).  Thus, negative epsilon
C          should always turn off this test.
          if(need3.or.(ai-alphai)*fprimeai.le.epsilon) then
            if(not_need3) then
              if(debug) then
                write(0,*) '(bi-ai), fprimeai, (bi-ai)*fprimeai =',
     &            (bi-ai), fprimeai, (bi-ai)*fprimeai
                write(0,*)
     &            '(ai-alphai), fprimeai, (ai-alphai)*fprimeai = ',
     &            (ai-alphai), fprimeai, (ai-alphai)*fprimeai
                write(0,*) 'epsilon = ', epsilon
              endif
C              recalculate everything for alphai = ai since that is less
C              expensive than carrying copies along for the entire line
C              search (and for _every_ line search until function convergence
C              _might_ finally trigger the above if statement.
              alphai = ai
C              Compute trial point at ai corresponding to
C              xi = x + alphai * s
              call bfgs_take_step(n, x, p, alphai, dir/pnorm, xi, dx)
C              EVALUATE f(ai) = function at xi.  This should be identical to
C              fai, but we recalculate here to set up the gradient calculation
C              to satisfy the condition (quite useful for calling routines)
C              that the gradient calculation always occurs at the same xi as
C              the preceeding function evaluation. Also EVALUATE gradient of
C              function at xi.
              status = 'both_search1'
              name = 'linesearch'
              need3 = .true.
              return
            endif  !not_need3
            if(need3) then
! 3            continue
              need3 = .false.
              need345 = need3.or.need4.or.need5
              not_need12345 = .not.(need12.or.need345)
              not_need345 = .not.need345
              not_need3 = .not.need3
              fi = f
C              from PMOO, 1.2.6
C              f'(alphai) = s dot gradient
              fprimei = (-dir/pnorm)*ddot(n, p, 1, gradient, 1)
C              output x = xi
              call dcopy(n, xi, 1, x, 1)
C              output f = f(xi)
              f = fi
C              gradient at x = xi already stored ready for output.
              alpha_ratio = alphai/alpha1
              alpha_previous = alphai
              deltaf_previous = f0-fi
              status = 'error roundoff'
              return
            else
              stop 'bfgs_linesearch: bad logic2'
            endif  !need3
          endif
C PMOO, 2.6.4. LINE 3
C          Compute new trial point at alphai corresponding to
C          xi = x + alphai * s
          call bfgs_take_step(n, x, p, alphai, dir/pnorm, xi, dx)
C          EVALUATE f(alphai) = function at xi
          status = 'function_search2'
          name = 'linesearch'
          need4 = .true.
          return
        endif  !not_need45
        if(not_need5) then
! 4        continue
          need4 = .false.
          need345 = need3.or.need4.or.need5
          not_need12345 = .not.(need12.or.need345)
          not_need345 = .not.need345
          not_need45 = .not.(need4.or.need5)
          fi = f
          if(fi.le.f0+alphai*rho*fprime0) fmin = min(fi, fmin)
          if(debug) then
            write(0,*) "alpha, f(alpha) = ", alphai, fi
          endif
        endif  !not_need5
C PMOO, 2.6.4. LINE 4
        if(not_need5.and.(
     &      (fi.gt.f0 + rho*alphai*fprime0).or.(fi.ge.fai))) then
C PMOO, 2.6.4. LINE 5
          aip = ai
C          fai and fprimeai unchanged.
          bip = alphai
          fbi = fi
C          don't know gradient at xi at this point in sectioning
C          iteration.
          fprimebi_calc = .false.
        else
          if(not_need5) then
C PMOO, 2.6.4. LINE 6
C            EVALUATE gradient of function at same xi as previous.
            status = 'gradient_search2'
            name = 'linesearch'
            need5 = .true.
            return
          endif
          if(need5) then
! 5          continue
            need5 = .false.
            need345 = need3.or.need4.or.need5
            not_need12345 = .not.(need12.or.need345)
            not_need345 = .not.need345
            not_need45 = .not.(need4.or.need5)
            not_need5 = .not.need5
C            from PMOO, 1.2.6
C            f'(alphai) = s dot gradient
            fprimei = (-dir/pnorm)*ddot(n, p, 1, gradient, 1)
            if(debug) then
              write(0,*) "alpha, f'(alpha) = ", alphai, fprimei
            endif
C PMOO, 2.6.4. LINE 7
            if(abs(fprimei).le.-sigma*fprime0) then
C              acceptable point found.
C              output x = xi
              call dcopy(n, xi, 1, x, 1)
C              output f = f(xi)
              f = fi
C              gradient at xi already stored ready for output.
              alpha_ratio = alphai/alpha1
              alpha_previous = alphai
              deltaf_previous = f0-fi
              status = 'continue_search4'
              return
            endif
C PMOO, 2.6.4. LINE 8
            aip = alphai
C PMOO, 2.6.4. LINE 9
            if(dalpha*fprimei.ge.0.d0) then
              bip = ai
              fbi = fai
              fprimebi_calc = .true.
              fprimebi = fprimeai
            else
              bip = bi
C              fbi, fprimebi_calc, and fprimebi (if defined) unchanged
            endif
            fai = fi
            fprimeai = fprimei
          else
            stop 'bfgs_linesearch: bad logic 3'
          endif  !need5
        endif
        ai = aip
        bi = bip
      enddo
      end

      subroutine bfgs_iterate(status, name,
     &  fbar, epsilon, n, x, xi, f, gradient, alpha_ratio, dx)
C      do one iteration of the BFGS minimization algorithm consisting
C      of a line search plus BFGS update.

C      status is a control string that takes the following values
C      ('function', 'gradient', 'both', and 'continue' have "_set" or
C      "_searchN" suffixes to identify for logic debugging purposes
C      exactly where the status was set):
C        'start'.  Set externally to initialize bfgs iteration.
C        'function'.  Request function value to be supplied externally
C        'gradient'.  Request gradient to be supplied externally.
C           N.B. gradient evaluation is always done at xi of prior
C           function evaluation
C        'both'.  Request both function value and gradient to be supplied
C                 externally.
C        'continue'. used internally to communicate with bfgs_linesearch.
C                    this status means an acceptable point has been found
C                    by the line search and a bfgs update should be done
C                    after the call to bfgs_linesearch.
C        'complete'.  initialize or line search + bfgs update has
C                     completed successfully and the x value is the initial
C                     point (only after 'start') or an acceptable point
C                     found by the line search.
C        'error positive'.  positive initial derivative in direction of
C                           line search.  no further progress can be made
C                           so return immediately.
C                           ('positive' currently unused and bfgs_linesearch
C                           stops instead if it runs into this condition which
C                           should be impossible.)
C        'error roundoff'.  line search did not pass roundoff error check
C                           discussed on page 38 of PMOO.  No further
C                           progress can be made so return immediately
C                           after line search.
C        'error zero dxdg'.  bfgs update found zero dxdg (dot product of
C                            the change in x and the change in gradient) so
C                            was unable to do the BFGS update.
C        Note, all values of status other than 'start' are set internally by
C        the BFGS routines and passed to the calling routine to control what
C        information is supplied on the next call to bfgs_iterate.

C      name is a control string set to the name of the appropriate bfgs
C      routine (with the bfgs_ prefix removed!) when there is a 'function',
C      'gradient', or 'both' request or when there is an error.

C      input quantities:
C      fbar is a user-specified minimum possible f used as per
C        PMOO, 2.6.1 to control the bracketing step size.  If actual f
C        values are <= fbar, then this routine stops with an error so
C        be realistic in how you specify fbar making it small enough
C        to avoid the error, but large enough to provide some control
C        over the maximum size of the bracketing step.
C      epsilon is the user-specified convergence criterion for f (not
C        relative f) discussed on p. 38 of PMOO.  It is used to terminate
C        the line-search routine to work around the case when round-off
C        errors are causing convergence problems close to the minimum.

C      input and output quantities:
C      x(n) is the starting point of the line search on input and on output
C        is the ending point of the line search.
C      xi(n) is trial values of x for external evaluations.
C        for status(:8).eq.'complete' result, x and xi are identical.
C      N.B. both f and gradient must be precalculated on input and calculated
C        on output.  That is:
C      f is f(xi(n)) on both input and output.
C      gradient(n) is the gradient(xi(n)) on both input and output.

C      output quantity:
C      alpha_ratio is the ratio of the calculated step size to the initial
C      estimate, alpha1 for the step size.
C      dx(n) is the vector of differences between the initial x and final x.

C      internal quantity:
C      p(nmax_bfgs) is the unscaled line search direction.
C      The PMOO scaled s vector is defined by
C        s = -(dir/pnorm)*p
C        where dir = 1 if p has an acute angle with the input gradient and
C        dir = -1 otherwise.  pnorm is the norm of p.  dir and pnorm are
C        calculated internally.

      implicit none
      include 'bfgs.h'
      integer lnblnk
      integer n
      double precision fbar, epsilon,
     &  x(n), xi(n), f, gradient(n), alpha_ratio, dx(n)
      double precision
     &  x0(nmax_bfgs), dx0(nmax_bfgs),
     &  g0(nmax_bfgs), dg0(nmax_bfgs), p(nmax_bfgs),
     &  ddot, dnrm2, dxg, dgg, dxdg, dgnorm, b, a
      character*(*) status, name
      logical initialized, recall
      data initialized/.false./
      save
      
C      sanity check
      if(n.gt.nmax_bfgs) stop 'bfgs_iterate: n too large'
      if(len(status).lt.16)
     &  stop 'bfgs_iterate: length of status string too small'
      if(len(name).lt.10)
     &  stop 'bfgs_iterate:length of name string too small'
      if(status(:5).eq.'error')
     &  stop 'bfgs_iterate: attempt to re-enter with error condition'

      recall =
     &  status(:4).eq.'both'.or.
     &  status(:8).eq.'function'.or.
     &  status(:8).eq.'gradient'

      if(status(:5).eq.'start'.and..not.recall) then
        call bfgs_set(status, name, n, x, xi, f, gradient, p)
        if(status(:8).ne.'complete') then
          return
        else
C          problemo: status should not be complete after initial call.
          write(0,*) "status = ", status(:lnblnk(status))
          stop 'bfgs_iterate: bad logic1'
        endif
      elseif(recall.and.name(:3).eq.'set') then
        call bfgs_set(status, name, n, x, xi, f, gradient, p)
        if(status(:8).eq.'complete') then
          initialized = .true.
          return
        else
C          problemo: status should always be 'complete_set' at this stage.
          write(0,*) "status = ", status(:lnblnk(status))
          stop 'bfgs_iterate: bad logic2'
        endif
      endif

C      further sanity check:
      if(.not.initialized)
     &  stop 'bfgs_iterate: not properly initialized'

      if(.not.recall) then
C        save old values.
        call dcopy(n, x, 1, x0, 1)
        call dcopy(n, gradient, 1, g0, 1)
        call bfgs_linesearch(status, name,
     &    fbar, epsilon, n, p, x, xi, f, gradient, alpha_ratio, dx)
C        usually this is a normal return asking for function or gradient
C        evaluation, but could be error depending on status.
        if(status(:8).ne.'continue') return
      elseif(name(:10).eq.'linesearch') then
C        this is the normal call that supplies the function and gradient
C        values that have been previously requested by bfgs_linesearch.
        call bfgs_linesearch(status, name,
     &    fbar, epsilon, n, p, x, xi, f, gradient, alpha_ratio, dx)
C        this is an error return.
        if(status(:8).ne.'continue') return
      else
          write(0,*) "status = ", status(:lnblnk(status))
          stop 'bfgs_iterate: bad logic3'
      endif
      
C      This is the BFGS update (taken from the GSL code).  I have not yet
C      verified this using PMOO, 3.2.12, but it gives good results for
C      the difficult test case of the Rosenbrock function so I consider
C      it to be numerically verified.
C      B = dx.g/dx.dg
C      A = - (1+ dg.dg/dx.dg) B + dg.g/dx.dg
C      p' = g1 - A dx - B dg

C      dx0 = x - x0 (note this is the delta used in PMOO, 3.2.12).
      call dcopy(n, x, 1, dx0, 1)
      call daxpy(n, -1.d0, x0, 1, dx0, 1)

C        dg0 = gradient - g0 (note this is the gamma used in PMOO, 3.2.12).
      call dcopy(n, gradient, 1, dg0, 1)
      call daxpy(n, -1.d0, g0, 1, dg0, 1)

C      delta dot g(k+1)
      dxg = ddot(n, dx0, 1, gradient, 1)
C      gamma dot g(k+1)
      dgg = ddot(n, dg0, 1, gradient, 1)
C      delta dot gamma
      dxdg = ddot(n, dx0, 1, dg0, 1)

      dgnorm = dnrm2(n, dg0, 1)

C      Guard against either dx0 or dg0 or both equal to the
C      zero vector or else dx0 perfectly perpendicular to dg0.
      if(dxdg.eq.0.d0) then
        status = 'error zero dxdg'
        return
      endif

C      B = dx.g/dx.dg
      b = dxg/dxdg

C      A = - (1+ dg.dg/dx.dg) B + dg.g/dx.dg
      a = -(1.d0 + dgnorm*dgnorm/dxdg)*b + dgg/dxdg

C      p' = g1 - A dx - B dg
      call dcopy(n, gradient, 1, p, 1)
      call daxpy(n, -a, dx0, 1, p, 1)
      call daxpy(n, -b, dg0, 1, p, 1)
      status = 'complete'
      end
