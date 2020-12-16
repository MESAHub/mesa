C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: qryd_approx.f 820 2008-06-24 19:26:56Z airwin $
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
      subroutine qryd_approx(ifpl, ifmhd, ifneutral, eps_factor,
     &  nmin, nmin_max, nmax, a, b, nb,
     &  qryd, qryda, qrydb, qryda2, qrydab, qrydb2, nmax_reached)
C      _use_explicit summations and calls to approximation routines  to:
C       calculate Rydberg partition function sums from the range
C       of minimum principal quantum numbers (given by nmin to nmin_max)
C       to nmax of 2 n^2 exp(a/n^2) times
C       planck_larkin occupation probability (ifpl true) times
C       mhd occupation probability (ifmhd true).
C       The planck-larkin occupation probability is given by
C         (1 - exp(-a)*(1+a)).
C       The *ln* of the mhd occupation probability is given by
C       -[b(1) + b(2)*n^2*(1+g(n)) + b(3)*n^4*(1+g(n))^2 +
C       b(4)*n^6*(1+g(n))^3 + b(5)*n^(15/2)*(1+h(n))].
C       g(n) = 1/(2n).
C       h(n) = (16/(3*n*K_n))^(3/2) - 1
C         = (16/(3*n))^(3/2) - 1 [for n <= 3]
C         = ((n+1)/n)^3*((n^2 + n + 1/2)/(n^2 + 7/6*n))^(3/2) - 1 [for n >=3]
C       input quantities:
C       ifpl (logical) controls whether to_use_planck-larkin
C         occupation probability.
C       ifmhd (logical) controls whether to_use_mhd occupation probability.
C       ifneutral (logical) controls when b(1) through b(4) are employed
C         in the occupation probability calculation.
C       eps_factor is a factor used to help terminate the principal
C         quantum number sum. eps_factor = exp(-c2*min(chi)/t), where
C         min(chi) = the minimum ionization potential for all species
C         with the same a value (i.e. charge).  This definition means
C         that at the limit summand contributes approximately eps (see
C         parameter below) to ln(1+q_excited/q_ground).  With exponential
C         cutoff (MHD), this insures maximum cutoff errors are of order eps.
C         if not MHD (i.e., only Planck-Larkin occupation probability) the
C         remainder of the sum is n*summand.  With nmax of order 300 000
C         this means total relative error is of order 1.d-4 which is still
C         fine.
C       nmin to nmin_max is the range of minimum principal quantum numbers
C         required.
C       nmax is the maximum principal quantum number.  In Planck-Larkin case
C         obtain relative errors of 1.d-5 for nmin of 3 if nmax = 300 000.
C       a = c2*Z^2*R/T.
C       b(nb), (nb = 5) determines the mhd occupation probability.
C       output quantities:
C       qryd(nmin_max), qryda(nmin_max), qrydb(nb, nmin_max),
C         qryda2(nmin_max), qrydab(nb, nmin_max),
C         qrydab2(nb,nb,nmin_max) (lower triangle in nb) are
C         the resulting sum plus partial derivatives wrt a and b.
C       nmax_reached is the returned maximum principal quantum number
C       that is used in the sum taking account of the convergence criteria.
      implicit none
      include 'constants.h'
      integer nmin, nmin_max, nmax, nb, nmax_reached
      double precision eps_factor, a, b(nb),
     &  qryd(nmin_max), qryda(nmin_max), qrydb(nb,nmin_max),
     &  qryda2(nmin_max), qrydab(nb,nmin_max), qrydb2(nb,nb,nmin_max)
      logical ifpl, ifmhd, ifneutral
C       internal variables:
      double precision eps
      parameter (eps=1.d-10)
      integer nb_local
      parameter(nb_local = 5)
      double precision rn2, arg, summand, summanda, summanda2,
     &  kn, hprime, gprime, occupation,
     &  lnoccupation, dlnoccupation(nb_local), dsummand(nb_local),
     &  dsummanda(nb_local), dsummand2(nb_local,nb_local),
     &  lqryd, lqryda, lqrydb(nb_local),
     &  lqryda2, lqrydab(nb_local), lqrydb2(nb_local,nb_local),
     &  lnoccupationa, dlnoccupationa(nb_local)
      integer n, nzero, nmaxa
      double precision exp_max, expmb1
C       this is a standard value for most of the approximations
C       which is used to limit their range of applicability
      parameter(exp_max = 1.d0)
C       limits on continuous interpolation between exact sum
C       and approximation
      double precision lim_sum, lim_approx
      parameter (lim_sum = -1.d-3)
      parameter (lim_approx = -1.d-4)
      double precision sarg, carg, dsarg, d2sarg,
     &  wapprox, dwapprox(nb_local), dwapprox2(nb_local,nb_local),
     &  wsum, dwsum(nb_local), dwsum2(nb_local,nb_local)
      if(nb.ne.nb_local) stop 'qryd_approx: bad nb'
      if(nmin.lt.1.or.nmin.gt.nmin_max) stop 'qryd_approx: bad nmin'
      nmax_reached = 0
      if(.not.ifmhd) then
        if(.not.ifpl)
     &    stop 'qryd_approx: bad values of ifmhd or ifpl'
        call plsum(nmin, nmin_max, 0, a, qryd, qryda, qryda2)
        if(nmax.lt.nmin_max) then
          stop 'qryd_approx: this case not programmed'
        elseif(nmax.le.300 000) then
          call plsum(nmax+1, nmax+1, nmax, a, lqryd, lqryda, lqryda2)
          do n = nmin, nmin_max
            qryd(n) = qryd(n) - lqryd
            qryda(n) = qryda(n) - lqryda
            qryda2(n) = qryda2(n) - lqryda2
          enddo
        endif
      else
        stop 'qryd_approx: excitation approximation is invalid'
C         n.b. the code below works in many cases, but there seems to
C         be bad significance loss (or else partial derivative troubles)
C         plaguing it and it is slow.  Thus, should be using truncation
C         approximation instead (see paper).
C         n.b. to reduce the size of the resulting executable, I have
C         commented out calls to approximation routines and removed them
C         from this directory.
C         adopt minimum limit since approximations
C         only work for nmin_max.ge.3 because of change in definition of
C         kn at n = 3.
        if(nmin_max.lt.3)
     &    stop 'qryd_approx: nmin_max must be 3 or greater for ifmhd'
C         n must be >= to this number in order to_use_approximations
C         max(nmin_max,... assures at least one loop for negligible a.
        nmaxa = max(nmin_max,int(sqrt(a/exp_max))+1)
        n = nmin
C         define these quantities to clear out any undefined garbage
        summand = 0.d0
        lnoccupation = 0.d0
        lnoccupationa = 0.d0
C         n.b. nmax and summand limits used below only for *very* low
C         temperatures where nmaxa can get large.  In particular nmax
C         is usually a no-op for the MHD case where you are using
C         approximations.  One could re-programme to make nmax be
C         effective for part of approximation where not in cutoff region,
C         but I have judged this is not worth it since ordinarily ignore
C         all approximations for MHD case.
        do while(.not.(n.gt.nmaxa.and.lnoccupation.gt.lim_approx).and.
     &      n.le.nmax.and.
     &      (n.eq.nmin.or.summand.gt.eps/eps_factor))
C           avoid integer overflow by doing this in double precision
          rn2 = dble(n)*dble(n)
          arg = a/rn2
          if(ifpl) then
            if(arg.gt.0.01d0) then
C               lose a maximum of 4 significant digits
              summand = 2.d0*rn2*(exp(arg) - (1.d0 + arg))
              summanda = 2.d0*(exp(arg) - 1.d0)
              summanda2 = 2.d0*exp(arg)/rn2
            else
              summand = rn2*arg*arg*(
     &          1.d0 + arg/3.d0*(
     &          1.d0 + arg/4.d0*(
     &          1.d0 + arg/5.d0*(
     &          1.d0 + arg/6.d0*(
     &          1.d0 + arg/7.d0*(
     &          1.d0 + arg/8.d0*(
     &          1.d0 + arg/9.d0*(
     &          1.d0 + arg/10.d0*(
     &          1.d0 + arg/11.d0)))))))))
              summanda = 2.d0*arg*(
     &          1.d0 + arg/2.d0*(
     &          1.d0 + arg/3.d0*(
     &          1.d0 + arg/4.d0*(
     &          1.d0 + arg/5.d0*(
     &          1.d0 + arg/6.d0*(
     &          1.d0 + arg/7.d0*(
     &          1.d0 + arg/8.d0*(
     &          1.d0 + arg/9.d0*(
     &          1.d0 + arg/10.d0)))))))))
              summanda2 = 2.d0*(
     &          1.d0 + arg*(
     &          1.d0 + arg/2.d0*(
     &          1.d0 + arg/3.d0*(
     &          1.d0 + arg/4.d0*(
     &          1.d0 + arg/5.d0*(
     &          1.d0 + arg/6.d0*(
     &          1.d0 + arg/7.d0*(
     &          1.d0 + arg/8.d0*(
     &          1.d0 + arg/9.d0)))))))))/rn2
            endif
          else
            summand = 2.d0*rn2*exp(arg)
            summanda = 2.d0*exp(arg)
            summanda2 = 2.d0*exp(arg)/rn2
          endif
C           quantum correction K_n see Hummer and Mihalas eq. 4.24
          if(n.le.3) then
            kn = 1.d0
          else
            kn = 
     &        (16.d0*rn2*(dble(n) + 7.d0/6.d0))/
     &        (dble(3*(n+1))*dble(n+1)*(rn2 + dble(n) + 0.5d0))
          endif
          hprime = (16.d0/(3.d0*dble(n)*kn))**1.5d0
C           mhd ln occupation probability for neutral-ion and
C           ion-ion interactions
          dlnoccupation(5) = -(dble(n))**7.5d0*hprime
          lnoccupation = b(5)*dlnoccupation(5)
          if(ifneutral) then
C             mhd ln occupation probability for neutral-neutral interactions
            gprime = (1.d0 + 0.5d0/dble(n))*rn2
C             n.b. b(1) applied later.
            lnoccupation = lnoccupation -
     &        (gprime*(b(2) + gprime*(b(3) + gprime*b(4))))
            dlnoccupation(2) = -gprime
            dlnoccupation(3) = dlnoccupation(2)*gprime
            dlnoccupation(4) = dlnoccupation(3)*gprime
          endif
          occupation = exp(lnoccupation)
          summand = summand*occupation
          summanda = summanda*occupation
          summanda2 = summanda2*occupation
          dsummand(5) = summand*dlnoccupation(5)
          dsummanda(5) = summanda*dlnoccupation(5)
          dsummand2(5,5) = summand*dlnoccupation(5)*dlnoccupation(5)
          if(ifneutral) then
            dsummand(2) = summand*dlnoccupation(2)
            dsummand(3) = summand*dlnoccupation(3)
            dsummand(4) = summand*dlnoccupation(4)
            dsummanda(2) = summanda*dlnoccupation(2)
            dsummanda(3) = summanda*dlnoccupation(3)
            dsummanda(4) = summanda*dlnoccupation(4)
            dsummand2(2,2) = summand*dlnoccupation(2)*dlnoccupation(2)
            dsummand2(3,2) = summand*dlnoccupation(3)*dlnoccupation(2)
            dsummand2(4,2) = summand*dlnoccupation(4)*dlnoccupation(2)
            dsummand2(5,2) = summand*dlnoccupation(5)*dlnoccupation(2)
            dsummand2(3,3) = summand*dlnoccupation(3)*dlnoccupation(3)
            dsummand2(4,3) = summand*dlnoccupation(4)*dlnoccupation(3)
            dsummand2(5,3) = summand*dlnoccupation(5)*dlnoccupation(3)
            dsummand2(4,4) = summand*dlnoccupation(4)*dlnoccupation(4)
            dsummand2(5,4) = summand*dlnoccupation(5)*dlnoccupation(4)
          endif
          if(n.eq.nmaxa) then
            lnoccupationa = lnoccupation
            dlnoccupationa(5) = dlnoccupation(5)
            if(ifneutral) then
              dlnoccupationa(2) = dlnoccupation(2)
              dlnoccupationa(3) = dlnoccupation(3)
              dlnoccupationa(4) = dlnoccupation(4)
            endif
          endif
          if(n.ge.nmaxa) then
C             n.b. by definition of nmaxa, n >= nmin_max.
C             n.n.b. qryd contains sum from nmin_max to nmaxa-1 of summand, 
C             thus, must zero it in this  special case.
            if(nmaxa.eq.nmin_max) then
              qryd(nmin_max) = 0.d0
              qryda(nmin_max) = 0.d0
              qryda2(nmin_max) = 0.d0
              qrydb(5,nmin_max) = 0.d0
              qrydab(5,nmin_max) = 0.d0
              qrydb2(5,5,nmin_max) = 0.d0
              if(ifneutral) then
                qrydb(2,nmin_max) = 0.d0
                qrydb(3,nmin_max) = 0.d0
                qrydb(4,nmin_max) = 0.d0
                qrydab(2,nmin_max) = 0.d0
                qrydab(3,nmin_max) = 0.d0
                qrydab(4,nmin_max) = 0.d0
                qrydb2(2,2,nmin_max) = 0.d0
                qrydb2(3,2,nmin_max) = 0.d0
                qrydb2(4,2,nmin_max) = 0.d0
                qrydb2(5,2,nmin_max) = 0.d0
                qrydb2(3,3,nmin_max) = 0.d0
                qrydb2(4,3,nmin_max) = 0.d0
                qrydb2(5,3,nmin_max) = 0.d0
                qrydb2(4,4,nmin_max) = 0.d0
                qrydb2(5,4,nmin_max) = 0.d0
              endif
            endif
            if(n.eq.nmaxa) then
              lqryd = summand
              lqryda = summanda
              lqryda2 = summanda2
              lqrydb(5) = dsummand(5)
              lqrydab(5) = dsummanda(5)
              lqrydb2(5,5) = dsummand2(5,5)
              if(ifneutral) then
                lqrydb(2) = dsummand(2)
                lqrydb(3) = dsummand(3)
                lqrydb(4) = dsummand(4)
                lqrydab(2) = dsummanda(2)
                lqrydab(3) = dsummanda(3)
                lqrydab(4) = dsummanda(4)
                lqrydb2(2,2) = dsummand2(2,2)
                lqrydb2(3,2) = dsummand2(3,2)
                lqrydb2(4,2) = dsummand2(4,2)
                lqrydb2(5,2) = dsummand2(5,2)
                lqrydb2(3,3) = dsummand2(3,3)
                lqrydb2(4,3) = dsummand2(4,3)
                lqrydb2(5,3) = dsummand2(5,3)
                lqrydb2(4,4) = dsummand2(4,4)
                lqrydb2(5,4) = dsummand2(5,4)
              endif
            else
              lqryd = lqryd +
     &          summand
              lqryda = lqryda +
     &          summanda
              lqryda2 = lqryda2 +
     &          summanda2
              lqrydb(5) = lqrydb(5) +
     &          dsummand(5)
              lqrydab(5) = lqrydab(5) +
     &          dsummanda(5)
              lqrydb2(5,5) = lqrydb2(5,5) +
     &          dsummand2(5,5)
              if(ifneutral) then
                lqrydb(2) = lqrydb(2) +
     &            dsummand(2)
                lqrydb(3) = lqrydb(3) +
     &            dsummand(3)
                lqrydb(4) = lqrydb(4) +
     &            dsummand(4)
                lqrydab(2) = lqrydab(2) +
     &            dsummanda(2)
                lqrydab(3) = lqrydab(3) +
     &            dsummanda(3)
                lqrydab(4) = lqrydab(4) +
     &            dsummanda(4)
                lqrydb2(2,2) = lqrydb2(2,2) +
     &            dsummand2(2,2)
                lqrydb2(3,2) = lqrydb2(3,2) +
     &            dsummand2(3,2)
                lqrydb2(4,2) = lqrydb2(4,2) +
     &            dsummand2(4,2)
                lqrydb2(5,2) = lqrydb2(5,2) +
     &            dsummand2(5,2)
                lqrydb2(3,3) = lqrydb2(3,3) +
     &            dsummand2(3,3)
                lqrydb2(4,3) = lqrydb2(4,3) +
     &            dsummand2(4,3)
                lqrydb2(5,3) = lqrydb2(5,3) +
     &            dsummand2(5,3)
                lqrydb2(4,4) = lqrydb2(4,4) +
     &            dsummand2(4,4)
                lqrydb2(5,4) = lqrydb2(5,4) +
     &            dsummand2(5,4)
              endif
            endif
          else
            if(n.le.nmin_max) then
              qryd(n) = summand
              qryda(n) = summanda
              qryda2(n) = summanda2
              qrydb(5,n) = dsummand(5)
              qrydab(5,n) = dsummanda(5)
              qrydb2(5,5,n) = dsummand2(5,5)
              if(ifneutral) then
                qrydb(2,n) = dsummand(2)
                qrydb(3,n) = dsummand(3)
                qrydb(4,n) = dsummand(4)
                qrydab(2,n) = dsummanda(2)
                qrydab(3,n) = dsummanda(3)
                qrydab(4,n) = dsummanda(4)
                qrydb2(2,2,n) = dsummand2(2,2)
                qrydb2(3,2,n) = dsummand2(3,2)
                qrydb2(4,2,n) = dsummand2(4,2)
                qrydb2(5,2,n) = dsummand2(5,2)
                qrydb2(3,3,n) = dsummand2(3,3)
                qrydb2(4,3,n) = dsummand2(4,3)
                qrydb2(5,3,n) = dsummand2(5,3)
                qrydb2(4,4,n) = dsummand2(4,4)
                qrydb2(5,4,n) = dsummand2(5,4)
              endif
            else
              qryd(nmin_max) = qryd(nmin_max) +
     &          summand
              qryda(nmin_max) = qryda(nmin_max) +
     &          summanda
              qryda2(nmin_max) = qryda2(nmin_max) +
     &          summanda2
              qrydb(5,nmin_max) = qrydb(5,nmin_max) +
     &          dsummand(5)
              qrydab(5,nmin_max) = qrydab(5,nmin_max) +
     &          dsummanda(5)
              qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +
     &          dsummand2(5,5)
              if(ifneutral) then
                qrydb(2,nmin_max) = qrydb(2,nmin_max) +
     &            dsummand(2)
                qrydb(3,nmin_max) = qrydb(3,nmin_max) +
     &            dsummand(3)
                qrydb(4,nmin_max) = qrydb(4,nmin_max) +
     &            dsummand(4)
                qrydab(2,nmin_max) = qrydab(2,nmin_max) +
     &            dsummanda(2)
                qrydab(3,nmin_max) = qrydab(3,nmin_max) +
     &            dsummanda(3)
                qrydab(4,nmin_max) = qrydab(4,nmin_max) +
     &            dsummanda(4)
                qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &            dsummand2(2,2)
                qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &            dsummand2(3,2)
                qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &            dsummand2(4,2)
                qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &            dsummand2(5,2)
                qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &            dsummand2(3,3)
                qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &            dsummand2(4,3)
                qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &            dsummand2(5,3)
                qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &            dsummand2(4,4)
                qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &            dsummand2(5,4)
              endif
            endif
          endif
          n = n + 1
        enddo
        nmax_reached = max(nmax_reached, n-1)
        if(n-1.lt.nmaxa) then
C           in this case, sum terminated by nmax or summand conditions.
C           only remaining thing to do is zero remainder of qryd if necessary.
          do nzero = n,nmin_max
            qryd(nzero) = 0.d0
            qryda(nzero) = 0.d0
            qryda2(nzero) = 0.d0
            qrydb(5,nzero) = 0.d0
            qrydab(5,nzero) = 0.d0
            qrydb2(5,5,nzero) = 0.d0
            if(ifneutral) then
              qrydb(2,nzero) = 0.d0
              qrydb(3,nzero) = 0.d0
              qrydb(4,nzero) = 0.d0
              qrydab(2,nzero) = 0.d0
              qrydab(3,nzero) = 0.d0
              qrydab(4,nzero) = 0.d0
              qrydb2(2,2,nzero) = 0.d0
              qrydb2(3,2,nzero) = 0.d0
              qrydb2(4,2,nzero) = 0.d0
              qrydb2(5,2,nzero) = 0.d0
              qrydb2(3,3,nzero) = 0.d0
              qrydb2(4,3,nzero) = 0.d0
              qrydb2(5,3,nzero) = 0.d0
              qrydb2(4,4,nzero) = 0.d0
              qrydb2(5,4,nzero) = 0.d0
            endif
          enddo
C           n.b. in all other cases n-1 reached nmaxa >= nmin_max, so
C           (1) lnoccupationa defined
C           (2) zero (when nmaxa = nmin_max) or
C             sum from nmin_max to nmaxa-1 of summand stored in qryd(nmin_max)
C           (3) sum from nmaxa to n-1 of summand stored in lqryd
        elseif(lnoccupationa.le.lim_sum) then
C           in this case sum terminated by nmax or summand conditions and
C           to finish simply add lqryd to qrd
          qryd(nmin_max) = qryd(nmin_max) + lqryd
          qryda(nmin_max) = qryda(nmin_max) + lqryda
          qryda2(nmin_max) = qryda2(nmin_max) + lqryda2
          qrydb(5,nmin_max) = qrydb(5,nmin_max) + lqrydb(5)
          qrydab(5,nmin_max) = qrydab(5,nmin_max) + lqrydab(5)
          qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) + lqrydb2(5,5)
          if(ifneutral) then
            qrydb(2,nmin_max) = qrydb(2,nmin_max) + lqrydb(2)
            qrydb(3,nmin_max) = qrydb(3,nmin_max) + lqrydb(3)
            qrydb(4,nmin_max) = qrydb(4,nmin_max) + lqrydb(4)
            qrydab(2,nmin_max) = qrydab(2,nmin_max) + lqrydab(2)
            qrydab(3,nmin_max) = qrydab(3,nmin_max) + lqrydab(3)
            qrydab(4,nmin_max) = qrydab(4,nmin_max) + lqrydab(4)
            qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &        lqrydb2(2,2)
            qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &        lqrydb2(3,2)
            qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &        lqrydb2(4,2)
            qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &        lqrydb2(5,2)
            qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &        lqrydb2(3,3)
            qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &        lqrydb2(4,3)
            qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &        lqrydb2(5,3)
            qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &        lqrydb2(4,4)
            qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &        lqrydb2(5,4)
          endif
        elseif(lnoccupationa.le.lim_approx) then
C           lim_sum < lnoccupationa <= lim_approx
C           interpolate between
C           complete approximation at lnoccupationa = lim_approx and
C           complete summation at lnoccupationa = lim_sum using interpolation
C           calculate interpolating function with zero first derivatives
C           at each limit.
          sarg = 0.5d0*sin(pi*0.5d0*(
     &      (lnoccupationa-lim_sum) -
     &      (lim_approx-lnoccupationa)
     &      )/(lim_approx-lim_sum))
          carg = 0.5d0*cos(pi*0.5d0*(
     &      (lnoccupationa-lim_sum) -
     &      (lim_approx-lnoccupationa)
     &      )/(lim_approx-lim_sum))
          dsarg = carg*pi/(lim_approx-lim_sum)
          d2sarg = -sarg*pi*pi/
     &      ((lim_approx-lim_sum)*(lim_approx-lim_sum))
          wapprox = 0.5d0 + sarg
          wsum  = 0.5d0 - sarg
          dwapprox(5) = dsarg*dlnoccupationa(5)
          dwsum(5) = -dwapprox(5)
          dwapprox2(5,5) = d2sarg*dlnoccupationa(5)*dlnoccupationa(5)
          dwsum2(5,5) = -dwapprox2(5,5)
          qryd(nmin_max) = qryd(nmin_max) + lqryd*wsum
          qryda(nmin_max) = qryda(nmin_max) + lqryda*wsum
          qryda2(nmin_max) = qryda2(nmin_max) + lqryda2*wsum
          qrydb(5,nmin_max) = qrydb(5,nmin_max) +
     &      lqrydb(5)*wsum + lqryd*dwsum(5)
          qrydab(5,nmin_max) = qrydab(5,nmin_max) +
     &      lqrydab(5)*wsum + lqryda*dwsum(5)
          qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +
     &      lqrydb2(5,5)*wsum +
     &      lqrydb(5)*dwsum(5) + lqrydb(5)*dwsum(5) +
     &      lqryd*dwsum2(5,5)
          if(ifneutral) then
            dwapprox(2) = dsarg*dlnoccupationa(2)
            dwapprox(3) = dsarg*dlnoccupationa(3)
            dwapprox(4) = dsarg*dlnoccupationa(4)
            dwsum(2) = -dwapprox(2)
            dwsum(3) = -dwapprox(3)
            dwsum(4) = -dwapprox(4)
            dwapprox2(2,2) = d2sarg*dlnoccupationa(2)*dlnoccupationa(2)
            dwapprox2(3,2) = d2sarg*dlnoccupationa(3)*dlnoccupationa(2)
            dwapprox2(4,2) = d2sarg*dlnoccupationa(4)*dlnoccupationa(2)
            dwapprox2(5,2) = d2sarg*dlnoccupationa(5)*dlnoccupationa(2)
            dwapprox2(3,3) = d2sarg*dlnoccupationa(3)*dlnoccupationa(3)
            dwapprox2(4,3) = d2sarg*dlnoccupationa(4)*dlnoccupationa(3)
            dwapprox2(5,3) = d2sarg*dlnoccupationa(5)*dlnoccupationa(3)
            dwapprox2(4,4) = d2sarg*dlnoccupationa(4)*dlnoccupationa(4)
            dwapprox2(5,4) = d2sarg*dlnoccupationa(5)*dlnoccupationa(4)
            dwsum2(2,2) = -dwapprox2(2,2)
            dwsum2(3,2) = -dwapprox2(3,2)
            dwsum2(4,2) = -dwapprox2(4,2)
            dwsum2(5,2) = -dwapprox2(5,2)
            dwsum2(3,3) = -dwapprox2(3,3)
            dwsum2(4,3) = -dwapprox2(4,3)
            dwsum2(5,3) = -dwapprox2(5,3)
            dwsum2(4,4) = -dwapprox2(4,4)
            dwsum2(5,4) = -dwapprox2(5,4)
            qrydb(2,nmin_max) = qrydb(2,nmin_max) +
     &        lqrydb(2)*wsum + lqryd*dwsum(2)
            qrydb(3,nmin_max) = qrydb(3,nmin_max) +
     &        lqrydb(3)*wsum + lqryd*dwsum(3)
            qrydb(4,nmin_max) = qrydb(4,nmin_max) +
     &        lqrydb(4)*wsum + lqryd*dwsum(4)
            qrydab(2,nmin_max) = qrydab(2,nmin_max) +
     &        lqrydab(2)*wsum + lqryda*dwsum(2)
            qrydab(3,nmin_max) = qrydab(3,nmin_max) +
     &        lqrydab(3)*wsum + lqryda*dwsum(3)
            qrydab(4,nmin_max) = qrydab(4,nmin_max) +
     &        lqrydab(4)*wsum + lqryda*dwsum(4)
            qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &        lqrydb2(2,2)*wsum +
     &        lqrydb(2)*dwsum(2) + lqrydb(2)*dwsum(2) +
     &        lqryd*dwsum2(2,2)
            qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &        lqrydb2(3,2)*wsum +
     &        lqrydb(3)*dwsum(2) + lqrydb(2)*dwsum(3) +
     &        lqryd*dwsum2(3,2)
            qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &        lqrydb2(4,2)*wsum +
     &        lqrydb(4)*dwsum(2) + lqrydb(2)*dwsum(4) +
     &        lqryd*dwsum2(4,2)
            qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &        lqrydb2(5,2)*wsum +
     &        lqrydb(5)*dwsum(2) + lqrydb(2)*dwsum(5) +
     &        lqryd*dwsum2(5,2)
            qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &        lqrydb2(3,3)*wsum +
     &        lqrydb(3)*dwsum(3) + lqrydb(3)*dwsum(3) +
     &        lqryd*dwsum2(3,3)
            qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &        lqrydb2(4,3)*wsum +
     &        lqrydb(4)*dwsum(3) + lqrydb(3)*dwsum(4) +
     &        lqryd*dwsum2(4,3)
            qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &        lqrydb2(5,3)*wsum +
     &        lqrydb(5)*dwsum(3) + lqrydb(3)*dwsum(5) +
     &        lqryd*dwsum2(5,3)
            qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &        lqrydb2(4,4)*wsum +
     &        lqrydb(4)*dwsum(4) + lqrydb(4)*dwsum(4) +
     &        lqryd*dwsum2(4,4)
            qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &        lqrydb2(5,4)*wsum +
     &        lqrydb(5)*dwsum(4) + lqrydb(4)*dwsum(5) +
     &        lqryd*dwsum2(5,4)
          endif
C           lnoccupationa > lim_sum so may calculate approximation
C           lnoccupationa evaluated for n = nmaxa which satisfies lnoccupation
C           criterion for approximation (note that both hprime and gprime
C           are always greater than unity.)
C           also note that n-1 >= to nmaxa which satisfies"a" criterion
!          call pi_totalsum_approx(ifpl, ifneutral, nmaxa,
!     &      a, b, nb, lqryd, lqryda, lqrydb,
!     &      lqryda2, lqrydab, lqrydb2)
          qryd(nmin_max) = qryd(nmin_max) + lqryd*wapprox
          qryda(nmin_max) = qryda(nmin_max) + lqryda*wapprox
          qryda2(nmin_max) = qryda2(nmin_max) + lqryda2*wapprox
          qrydb(5,nmin_max) = qrydb(5,nmin_max) +
     &      lqrydb(5)*wapprox + lqryd*dwapprox(5)
          qrydab(5,nmin_max) = qrydab(5,nmin_max) +
     &      lqrydab(5)*wapprox + lqryda*dwapprox(5)
          qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +
     &      lqrydb2(5,5)*wapprox +
     &      lqrydb(5)*dwapprox(5) + lqrydb(5)*dwapprox(5) +
     &      lqryd*dwapprox2(5,5)
          if(ifneutral) then
            qrydb(2,nmin_max) = qrydb(2,nmin_max) +
     &        lqrydb(2)*wapprox + lqryd*dwapprox(2)
            qrydb(3,nmin_max) = qrydb(3,nmin_max) +
     &        lqrydb(3)*wapprox + lqryd*dwapprox(3)
            qrydb(4,nmin_max) = qrydb(4,nmin_max) +
     &        lqrydb(4)*wapprox + lqryd*dwapprox(4)
            qrydab(2,nmin_max) = qrydab(2,nmin_max) +
     &        lqrydab(2)*wapprox + lqryda*dwapprox(2)
            qrydab(3,nmin_max) = qrydab(3,nmin_max) +
     &        lqrydab(3)*wapprox + lqryda*dwapprox(3)
            qrydab(4,nmin_max) = qrydab(4,nmin_max) +
     &        lqrydab(4)*wapprox + lqryda*dwapprox(4)
            qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &        lqrydb2(2,2)*wapprox +
     &        lqrydb(2)*dwapprox(2) + lqrydb(2)*dwapprox(2) +
     &        lqryd*dwapprox2(2,2)
            qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &        lqrydb2(3,2)*wapprox +
     &        lqrydb(3)*dwapprox(2) + lqrydb(2)*dwapprox(3) +
     &        lqryd*dwapprox2(3,2)
            qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &        lqrydb2(4,2)*wapprox +
     &        lqrydb(4)*dwapprox(2) + lqrydb(2)*dwapprox(4) +
     &        lqryd*dwapprox2(4,2)
            qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &        lqrydb2(5,2)*wapprox +
     &        lqrydb(5)*dwapprox(2) + lqrydb(2)*dwapprox(5) +
     &        lqryd*dwapprox2(5,2)
            qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &        lqrydb2(3,3)*wapprox +
     &        lqrydb(3)*dwapprox(3) + lqrydb(3)*dwapprox(3) +
     &        lqryd*dwapprox2(3,3)
            qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &        lqrydb2(4,3)*wapprox +
     &        lqrydb(4)*dwapprox(3) + lqrydb(3)*dwapprox(4) +
     &        lqryd*dwapprox2(4,3)
            qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &        lqrydb2(5,3)*wapprox +
     &        lqrydb(5)*dwapprox(3) + lqrydb(3)*dwapprox(5) +
     &        lqryd*dwapprox2(5,3)
            qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &        lqrydb2(4,4)*wapprox +
     &        lqrydb(4)*dwapprox(4) + lqrydb(4)*dwapprox(4) +
     &        lqryd*dwapprox2(4,4)
            qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &        lqrydb2(5,4)*wapprox +
     &        lqrydb(5)*dwapprox(4) + lqrydb(4)*dwapprox(5) +
     &        lqryd*dwapprox2(5,4)
          endif
        else
C           lnoccupationa > lim_approx so may calculate approximation
C           lnoccupationa evaluated for n = nmaxa which satisfies lnoccupation
C           criterion for approximation (note that both hprime and gprime
C           are always greater than unity.)
C           also note that n-1 >= to nmaxa which satisfies"a" criterion
!          call pi_totalsum_approx(ifpl, ifneutral, nmaxa,
!     &      a, b, nb, lqryd, lqryda, lqrydb,
!     &      lqryda2, lqrydab, lqrydb2)
C          _use_pure approximation
          qryd(nmin_max) = qryd(nmin_max) + lqryd
          qryda(nmin_max) = qryda(nmin_max) + lqryda
          qryda2(nmin_max) = qryda2(nmin_max) + lqryda2
          qrydb(5,nmin_max) = qrydb(5,nmin_max) + lqrydb(5)
          qrydab(5,nmin_max) = qrydab(5,nmin_max) + lqrydab(5)
          qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) + lqrydb2(5,5)
          if(ifneutral) then
            qrydb(2,nmin_max) = qrydb(2,nmin_max) + lqrydb(2)
            qrydb(3,nmin_max) = qrydb(3,nmin_max) + lqrydb(3)
            qrydb(4,nmin_max) = qrydb(4,nmin_max) + lqrydb(4)
            qrydab(2,nmin_max) = qrydab(2,nmin_max) + lqrydab(2)
            qrydab(3,nmin_max) = qrydab(3,nmin_max) + lqrydab(3)
            qrydab(4,nmin_max) = qrydab(4,nmin_max) + lqrydab(4)
            qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &        lqrydb2(2,2)
            qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &        lqrydb2(3,2)
            qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &        lqrydb2(4,2)
            qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &        lqrydb2(5,2)
            qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &        lqrydb2(3,3)
            qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &        lqrydb2(4,3)
            qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &        lqrydb2(5,3)
            qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &        lqrydb2(4,4)
            qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &        lqrydb2(5,4)
          endif
        endif
C         calculate constant occupation probability factor and apply it.
        if(ifneutral) then
          expmb1 = exp(-b(1))
          dlnoccupation(1) = -1.d0
        else
          expmb1 = 1.d0
        endif
        qryd(nmin_max) = qryd(nmin_max)*expmb1
        qryda(nmin_max) = qryda(nmin_max)*expmb1
        qryda2(nmin_max) = qryda2(nmin_max)*expmb1
        qrydb(5,nmin_max) = qrydb(5,nmin_max)*expmb1
        qrydab(5,nmin_max) = qrydab(5,nmin_max)*expmb1
        qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max)*expmb1
        if(ifneutral) then
          qrydb(1,nmin_max) = qryd(nmin_max)*dlnoccupation(1)
          qrydb(2,nmin_max) = qrydb(2,nmin_max)*expmb1
          qrydb(3,nmin_max) = qrydb(3,nmin_max)*expmb1
          qrydb(4,nmin_max) = qrydb(4,nmin_max)*expmb1
          qrydab(1,nmin_max) = qryda(nmin_max)*dlnoccupation(1)
          qrydab(2,nmin_max) = qrydab(2,nmin_max)*expmb1
          qrydab(3,nmin_max) = qrydab(3,nmin_max)*expmb1
          qrydab(4,nmin_max) = qrydab(4,nmin_max)*expmb1
          qrydb2(1,1,nmin_max) = qrydb(1,nmin_max)*dlnoccupation(1)
          qrydb2(2,1,nmin_max) = qrydb(2,nmin_max)*dlnoccupation(1)
          qrydb2(3,1,nmin_max) = qrydb(3,nmin_max)*dlnoccupation(1)
          qrydb2(4,1,nmin_max) = qrydb(4,nmin_max)*dlnoccupation(1)
          qrydb2(5,1,nmin_max) = qrydb(5,nmin_max)*dlnoccupation(1)
          qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max)*expmb1
          qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max)*expmb1
          qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max)*expmb1
          qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max)*expmb1
          qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max)*expmb1
          qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max)*expmb1
          qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max)*expmb1
          qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max)*expmb1
          qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max)*expmb1
        endif
C         convert from summand to sum for lower indices
        do n = nmin_max-1,nmin,-1
          qryd(n) = qryd(n)*expmb1
          qryda(n) = qryda(n)*expmb1
          qryda2(n) = qryda2(n)*expmb1
          qrydb(5,n) = qrydb(5,n)*expmb1
          qrydab(5,n) = qrydab(5,n)*expmb1
          qrydb2(5,5,n) = qrydb2(5,5,n)*expmb1
          if(ifneutral) then
            qrydb(1,n) = qryd(n)*dlnoccupation(1)
            qrydb(2,n) = qrydb(2,n)*expmb1
            qrydb(3,n) = qrydb(3,n)*expmb1
            qrydb(4,n) = qrydb(4,n)*expmb1
            qrydab(1,n) = qrydab(1,n+1) + qryda(n)*dlnoccupation(1)
            qrydab(2,n) = qrydab(2,n+1) + qrydab(2,n)*expmb1
            qrydab(3,n) = qrydab(3,n+1) + qrydab(3,n)*expmb1
            qrydab(4,n) = qrydab(4,n+1) + qrydab(4,n)*expmb1
            qrydb2(1,1,n) = qrydb2(1,1,n+1) +
     &        qrydb(1,n)*dlnoccupation(1)
            qrydb2(2,1,n) = qrydb2(2,1,n+1) +
     &        qrydb(2,n)*dlnoccupation(1)
            qrydb2(3,1,n) = qrydb2(3,1,n+1) +
     &        qrydb(3,n)*dlnoccupation(1)
            qrydb2(4,1,n) = qrydb2(4,1,n+1) +
     &        qrydb(4,n)*dlnoccupation(1)
            qrydb2(5,1,n) = qrydb2(5,1,n+1) +
     &        qrydb(5,n)*dlnoccupation(1)
            qrydb(1,n) = qrydb(1,n+1) + qrydb(1,n)
            qrydb(2,n) = qrydb(2,n+1) + qrydb(2,n)
            qrydb(3,n) = qrydb(3,n+1) + qrydb(3,n)
            qrydb(4,n) = qrydb(4,n+1) + qrydb(4,n)
            qrydb2(2,2,n) = qrydb2(2,2,n+1) + qrydb2(2,2,n)*expmb1
            qrydb2(3,2,n) = qrydb2(3,2,n+1) + qrydb2(3,2,n)*expmb1
            qrydb2(4,2,n) = qrydb2(4,2,n+1) + qrydb2(4,2,n)*expmb1
            qrydb2(5,2,n) = qrydb2(5,2,n+1) + qrydb2(5,2,n)*expmb1
            qrydb2(3,3,n) = qrydb2(3,3,n+1) + qrydb2(3,3,n)*expmb1
            qrydb2(4,3,n) = qrydb2(4,3,n+1) + qrydb2(4,3,n)*expmb1
            qrydb2(5,3,n) = qrydb2(5,3,n+1) + qrydb2(5,3,n)*expmb1
            qrydb2(4,4,n) = qrydb2(4,4,n+1) + qrydb2(4,4,n)*expmb1
            qrydb2(5,4,n) = qrydb2(5,4,n+1) + qrydb2(5,4,n)*expmb1
          endif
          qryd(n) = qryd(n+1) + qryd(n)
          qryda(n) = qryda(n+1) + qryda(n)
          qryda2(n) = qryda2(n+1) + qryda2(n)
          qrydb(5,n) = qrydb(5,n+1) + qrydb(5,n)
          qrydab(5,n) = qrydab(5,n+1) + qrydab(5,n)
          qrydb2(5,5,n) = qrydb2(5,5,n+1) + qrydb2(5,5,n)
        enddo
      endif
      end
