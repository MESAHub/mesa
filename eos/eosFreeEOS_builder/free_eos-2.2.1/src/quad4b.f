C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: quad4b.f 811 2008-06-22 22:06:57Z airwin $
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
C*******************************************************************************
C      Note:
C      The original version of this routine is in the public domain.  It was
C      written in the late 60's or early 70's by Darwin Harwood for the US
C      government, and it was brought to the attention of AWI by Forrest J.
C      Rogers. AWI relicensed the routine under the GPL and put in many
C      changes to convert to structured programming and double precision,
C      clean up the logic (especially one case where some questionable and
C      unnecessary bit fiddling was going on),_use_more precise quadrature
C      constants, and insert commentary.
C
      subroutine quad4b(ans,a,b,f,err)
c        quad4b integrates the function f from a to b with relative error
c        less than err and stores the answer in ans.
      implicit none
      integer maxnint
C      N.B. If you need more integration intervals than this, then your
C      problem is poorly posed; significance loss or integrand discontinuities
C      might be killing you, or else your integrand has a slow polynomial-type
C      change for one range, and a rapid exponential cutoff in another range.
C      In the latter case (which, e.g., happens for Fermi-Dirac integrals),
C      splitting the range into two parts solves the problem.
C
      parameter(maxnint=200000)
      double precision ans, a, b, f, err
      double precision x(maxnint),pint(maxnint),u(2),w(2)
!      data  u/.8611363115941 , .3399810435849/
!      data  w/.3478548451375 , .6521451548625/
C      calculated using gausle
      data  u/8.611363115940526d-01, 3.399810435848563d-01/
      data  w/3.478548451374538d-01, 6.521451548625462d-01/
      
C      AWI commented out this logic, and ie = ie-it below because it
C      made no sense for double precision data.  In any case, this is
C      a catch if things go wrong, and the maxnint limit does that in any case.
!      equivalence (et,ie)
!      data  it/40/
      integer ns, nint, n, ifmore, k
      double precision r, xmn, xmx, pm, dq5, s, apint, e, c, d, g, rel, et
C      The error criterion is whether the individual error in each piece
C      is less than r* total integral value as approximated at that
C      stage of the refinement.  So total relative error could theoretically
C      be maxnint*r, but certain areas of integral will dominate this sum,
C      there will be cancellation of positive and negative errors,
C      and number of intervals less than maxnint so actual error much less
C      than maxnint*r so if we adopt r = err/4 the actual error comes out
C      roughly equal to the specified err value.
      r=0.25d0*err
      xmn=a
      x(1)=xmn
      xmx=b
      x(2)=xmx
!      assign 3 to ns
      ns = 3
c        integration procedure: 4-point Gauss quadrature from xmn to xmx
   1  dq5 = 0.5d0*(xmx-xmn)
      pm=xmn+dq5
      s = dq5*(
     &  ((f(pm-u(1)*dq5)+f(pm+u(1)*dq5))*w(1)) + 
     &  ((f(pm-u(2)*dq5)+f(pm+u(2)*dq5))*w(2)))
!      go to ns
      if(ns.eq.3) then
        go to 3
      elseif(ns.eq.5) then
        go to 5
      elseif(ns.eq.6) then
        go to 6
      else
        stop 'quad4b: internal logic error'
      endif
C      Used only after initial integration over whole interval to
C      start integration refinement.
C      apint is numerical approximation to total integral
   3  apint=s
C      pint(i) is numerical approximation to integral from x(i) to x(i+1).
      pint(1)=apint
C      nint is the number of intervals.
      nint=1
C      x(n) to x(n+1) is sub-interval of integration variable we are refining.
      n=nint
C      Integrate between x(n) and half point that is half way to x(n+1)
   4  xmn=x(n)
      e=0.5d0*(x(n)+x(n+1))
      xmx=e
!      assign 5 to ns
      ns = 5
      go to 1
C      Save integral from x(n) to half-point in c and do integration from
C      half point to x(n+1).
   5  c=s
      xmn=e
      xmx=x(n+1)
!      assign 6 to ns
      ns = 6
      go to 1
   6  d=s
C      c and d are first and second half integrals from x(n) to x(n+1)
C      g is correction to previous integral.
      g=c+d-pint(n)
C      apint is refined value of _total_ integral.
      apint=apint+g
C      Set ifmore to 1 if more refinement is necessary for integral from
C      x(n) to x(n+1).
      ifmore = 0
      if(apint.ne.0.d0) then
        rel = abs(g/apint)
        if(rel-r.ge.0.d0) then
          ifmore = 1
        endif
      elseif(g.ne.0.d0) then
          ifmore = 1
      endif
      if(ifmore.eq.1) then
C        increment the number of intervals to be refined.
        nint=nint+1
        if(nint.ge.maxnint)
     &    stop 'quad4b: shouldn''t happen unless badly posed problem'
C        shift x and pint to make room
        x(nint+1)=x(nint)
        do k=nint,n+2,-1
          pint(k)=pint(k-1)
          x(k)=x(k-1)
        enddo
C        stick appropriate values into x and pint at correct indices.
C        c is first half integral from old x(n) to x(n+1),
C        d is second half integral and e is half point.
        pint(n)=c
        pint(n+1)=d
        x(n+1)=e
        et=e
C        This catch (which doesn't work for double precision in any
C        case) not needed since the maxnint limit above will stop infinite
C        loop if refinement cannot succeed.
!        ie=ie-it
        if(et-x(n).ge.0.d0) go to 4
      else
C        interval from x(n) to x(n+1) passes convergence criterion on
C        refinement so try next interval unless you have run out of them
C        and are therefore done.
        n=n+1
        if(n-nint.le.0) go to 4
      endif
      ans=apint
      end
