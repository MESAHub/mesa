C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: plsum.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine plsum(nlow, nhigh, ioffset, a, sum, dsum, dsum2)
C       calculate planck-larkin sum where
C       sum(nmin) = sum from nmin to infinity (nmin = nlow,...,nhigh) of
C       2 n^2 (exp(a/n^2) -(1+a/n^2))
C       dsum is the derivative of sum wrt to a
C       dsum2 is the second derivative of sum wrt to a
C      _use_direct summation until can_use_plsum_approx
C       if nlow too large for plsum_approx,_use_scaling rule
      implicit none
      integer nlow, nhigh, ioffset, nlow_scale, n,
     &  nstore, nstorep1, nstart
C       at this value or greater_use_scaling law
C       have proved through tests that the scaling error goes like
C       nlow_scale^3 and the relative error at this value is a maximum
C       of 2.d-5
      parameter (nlow_scale=31)
      double precision a, sum(nhigh-ioffset), dsum(nhigh-ioffset),
     &  dsum2(nhigh-ioffset), summand, summanda, summanda2
      double precision exp_max, arg, exparg, rn2, savearg, scale
C       this is a standard value for most of the approximations
C       which is used to limit their range of applicability
      parameter(exp_max = 1.d0)
      if(nlow-ioffset.lt.1) stop 'plsum: bad ioffset value'
      if(a.le.0.d0) stop 'plsum: bad a value'
      if(nlow.le.0.or.nlow.gt.nhigh) stop 'plsum: bad nlow or nhigh'
C       full sum results via scaling or by approximation.
C       minimum principal quantum number that can be done this way.
      nstart = max(2,nlow,int(sqrt(a/exp_max))+1)
      n = nstart
      do while(n.eq.nstart.or.n.le.nhigh)
C         first part of while assures sum calculated for n = nstart,
C         and next statement assures this result stored in nhigh, if
C         this do loop only done once because nstart > nhigh.
        nstore = min(n,nhigh)-ioffset
        if(n.ge.nlow_scale) then
C           scaling rule.
          scale = dble(n)/dble(nlow_scale-1)
          call plsum_approx(nlow_scale-1, a/scale/scale,
     &      sum(nstore), dsum(nstore), dsum2(nstore))
          dsum(nstore) = dsum(nstore)/(scale*scale)
          dsum2(nstore) = dsum2(nstore)/(scale*scale*scale*scale)
C           transform results to integral from nlow_scale -1 of
C           f(a/scale/scale,t) according to trapezoidal rule,
C           transform to integral from n of f(a,t) by multiplying by scale^3,
C           and finally transform to sum by trapezoidal rule
C           n.b. argprime = a/scale/scale/(n_low_scale-1)^2 = arg
          arg = a/(dble(n)*dble(n))
          if(arg.gt.0.01d0) then
C             lose a maximum of 4 significant digits
            exparg = exp(arg)
            summand = exparg - (1.d0 + arg)
            summanda = exparg - 1.d0
            summanda2 = exparg
          else
            summand = arg*arg/2.d0*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0*(
     &        1.d0 + arg/10.d0*(
     &        1.d0 + arg/11.d0)))))))))
            summanda = arg*(
     &        1.d0 + arg/2.d0*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0*(
     &        1.d0 + arg/10.d0)))))))))
            summanda2 = (
     &        1.d0 + arg*(
     &        1.d0 + arg/2.d0*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0)))))))))
          endif
          savearg = dble(n)*dble(n)
C           sum(nstore) = (sum(nstore) - savearg/scale/scale*(exparg-(1.d0+arg)))*
C      &      scale*scale*scale + savearg*(exparg-(1.d0+arg))
C           dsum(nstore) = (dsum(nstore) - 1.d0/scale/scale*(exparg-1.d0))*
C      &      scale*scale*scale + (exparg-1.d0)
C           dsum2(nstore) = (dsum2(nstore) - 1.d0/scale/scale/savearg*exparg)*
C      &      scale*scale*scale + exparg
          sum(nstore) = sum(nstore)*scale*scale*scale +
     &      savearg*summand*(1.d0-scale)
          dsum(nstore) = dsum(nstore)*scale*scale*scale +
     &      1.d0*summanda*(1.d0-scale)
          dsum2(nstore) = dsum2(nstore)*scale*scale*scale +
     &      1.d0/savearg*summanda2*(1.d0-scale)
        else
          call plsum_approx(n, a,
     &      sum(nstore), dsum(nstore), dsum2(nstore))
        endif
        n = n + 1
      enddo
C       at this point have calculated full sum results for indexes from
C       nstart to nhigh (n.b. or just nstart and stored in nhigh index
C       if nhigh smaller than nstart).
C       do explicit summation for remaining minimum principal
C       quantum numbers
      do n = nstart-1,nlow,-1
        rn2 = dble(n)*dble(n)
        arg = a/rn2
        if(arg.gt.0.01d0) then
C           lose a maximum of 4 significant digits
          exparg = exp(arg)
          summand = 2.d0*rn2*(exparg - (1.d0 + arg))
          summanda = 2.d0*(exparg - 1.d0)
          summanda2 = 2.d0*exparg/rn2
        else
          summand = rn2*arg*arg*(
     &      1.d0 + arg/3.d0*(
     &      1.d0 + arg/4.d0*(
     &      1.d0 + arg/5.d0*(
     &      1.d0 + arg/6.d0*(
     &      1.d0 + arg/7.d0*(
     &      1.d0 + arg/8.d0*(
     &      1.d0 + arg/9.d0*(
     &      1.d0 + arg/10.d0*(
     &      1.d0 + arg/11.d0)))))))))
          summanda = 2.d0*arg*(
     &      1.d0 + arg/2.d0*(
     &      1.d0 + arg/3.d0*(
     &      1.d0 + arg/4.d0*(
     &      1.d0 + arg/5.d0*(
     &      1.d0 + arg/6.d0*(
     &      1.d0 + arg/7.d0*(
     &      1.d0 + arg/8.d0*(
     &      1.d0 + arg/9.d0*(
     &      1.d0 + arg/10.d0)))))))))
          summanda2 = 2.d0*(
     &      1.d0 + arg*(
     &      1.d0 + arg/2.d0*(
     &      1.d0 + arg/3.d0*(
     &      1.d0 + arg/4.d0*(
     &      1.d0 + arg/5.d0*(
     &      1.d0 + arg/6.d0*(
     &      1.d0 + arg/7.d0*(
     &      1.d0 + arg/8.d0*(
     &      1.d0 + arg/9.d0)))))))))/rn2
        endif
        nstore = min(nhigh,n)-ioffset
        nstorep1 = min(nhigh,n+1)-ioffset
        sum(nstore) = sum(nstorep1) + summand
        dsum(nstore) = dsum(nstorep1) + summanda
        dsum2(nstore) = dsum2(nstorep1) + summanda2
      enddo
      end
