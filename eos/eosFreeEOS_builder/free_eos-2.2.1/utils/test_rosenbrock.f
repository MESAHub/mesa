      implicit none
      integer n
      parameter(n=2)
!      parameter(n=1)
      double precision x(n), xi(n), dx(n), f, gradient(n)
      double precision fold, step_ratio, rosenbrock,
     &  deltaf_criterion, gradient_criterion, gmax
      character*100 status, name
      integer linesearch_count, function_count, gradient_count
      integer lnblnk
      double precision fbar, epsilon

C      Fletcher's value.
!      deltaf_criterion = 1.d-8
C      Instead we converge until the function no longer decreases,
C      (note the fold.gt.f+deltaf_criterion criterion below) or a zero
C      gradient is reached.
C      a zero gradient is actually reached on at least the Linux+g77
C      platform which normally causes an 'error zero dxdg' status if
C      you attempt to converge some more, but we avoid that error
C      by testing for a zero gradient.
      deltaf_criterion = 0.d0
      gradient_criterion = 0.d0
C      We now define the starting point (read from figure 1.2.2 of
C      Fletcher's book and also confirmed with a google search for
C      rosenbrock traditional starting point)

      if(n.eq.1) then
C      starting point for one-dimensional Rosenbrock with x(2) fixed at unity.
        x(1) = 0.1d0
      else
        x(1) = -1.2d0
      endif
      if(n.eq.2) then
        x(2) = 1.d0
      endif
C      Setting fbar to a large negative number essentially disables the
C      bracket-limiting test and is fine for well-behaved
C      functions.  Rosenbrock's function is zero at the minimum and
C      otherwise positive, so following fbar is fine.
      fbar = -1.d-300
C      epsilon is a delta f expected from roundoff error to add
C      some robustness in the presence of round-off error.  Turn off this
C      test by setting epsilon to a negative value.
      epsilon = -1.d-300
C      The following loop does one complete iteration of bfgs_iterate
C      consisting of optional start, linesearch + bfgs update.
      status = 'start'
      linesearch_count = 0
      function_count = 0
      gradient_count = 0
C      force at least two iterations.
      f = 1.d300
      fold = 2.d0*f
      gmax = 1.d0
      do while(linesearch_count.lt.50.and.
     &    fold.gt.f+deltaf_criterion.and.
     &    gmax.gt.gradient_criterion.and.
     &    status(:5).ne.'error')
        fold = f
        call bfgs_iterate(status, name,
     &    fbar, epsilon,
     &    n, x, xi, f, gradient, step_ratio, dx)
        do while(.not.
     &      (status(:5).eq.'error'.or.status(:8).eq.'complete'))
!          write(*,*) 'status = ',status(:lnblnk(status))
          if(status(:4).eq.'both') then
            f = rosenbrock(n, xi)
            call grosenbrock(n, xi, gradient)
            function_count = function_count + 1
            gradient_count = gradient_count + 1
!            write(*,*) 'xi = ', xi
!            write(*,*) 'f = ', f
!            write(*,*) 'gradient = ', gradient
          elseif(status(:8).eq.'function') then
            f = rosenbrock(n, xi)
            function_count = function_count + 1
!            write(*,*) 'xi = ', xi
!            write(*,*) 'f = ', f
          elseif(status(:8).eq.'gradient') then
            call grosenbrock(n, xi, gradient)
            gradient_count = gradient_count + 1
!            write(*,*) 'xi = ', xi
!            write(*,*) 'gradient = ', gradient
          else
            stop 'test_rosenbrock: bad logic'
          endif
          call bfgs_iterate(status, name,
     &      fbar, epsilon, n, x, xi, f, gradient, step_ratio, dx)
        enddo
        if(status(:8).eq.'complete') then
          write(*,*) 'linesearch_count, function_count, '//
     &      'gradient_count, '//
     &      'x-solution, dx, step_ratio, f, gradient ='
          write(*,'(3i5)') linesearch_count, function_count,
     &      gradient_count
          if(n.eq.2) then
            write(*,'(1p2e15.5)') x(1)-1.d0, x(2)-1.d0
          else
            write(*,'(1p2e15.5)') x(1)-1.d0
          endif
          if(linesearch_count.gt.0) then
            write(*,'(1p2e15.5)') dx
            write(*,'(1p2e15.5)') step_ratio
          else
            if(n.eq.2) then
              write(*,'(1p2e15.5)') 0.d0, 0.d0
            else
              write(*,'(1p2e15.5)') 0.d0
            endif
            write(*,'(1p2e15.5)') 0.d0
          endif
          write(*,'(1p2e15.5)') f
          write(*,'(1p2e15.5)') gradient
          gmax = max(abs(gradient(1)), abs(gradient(2)))
C          increment after write since linesearch = 0 corresponds to
C          status = 'start' which does initialization calculations at
C          initial x.
          linesearch_count = linesearch_count + 1
        else
          write(0,*) 'status = ', status(:lnblnk(status))
        endif
      enddo
      end

      double precision function rosenbrock(n, x)
      implicit none
      integer n
      double precision x(n), t1
      if(n.eq.1) then
        t1 = 1.d0 - x(1)**2
      else
        t1 = x(2) - x(1)**2
      endif
      rosenbrock = 100.d0*t1*t1 + (x(1) - 1.d0)**2
      end

      subroutine grosenbrock(n, x, g)
      implicit none
      integer n
      double precision x(n), g(n), t1
      if(n.eq.1) then
        t1 = 1.d0 - x(1)**2
      else
        t1 = x(2) - x(1)**2
      endif
C      Compute gradient g for the sample problem.
      g(1) = 2.d0*(x(1) - 1.d0) - 4.d0*100.d0*t1*x(1)
      if(n.gt.1) then
        g(2) = 2.d0*100.d0*t1
      endif
      end
