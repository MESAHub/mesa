C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: qryd_calc.f 820 2008-06-24 19:26:56Z airwin $
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
      subroutine qryd_calc(ifpl, ifmhd, ifneutral, eps_factor,
     &  nmin, nmin_max, nmax, a, b, nb,
     &  qryd, qryda, qrydb, qryda2, qrydab, qrydb2, nmax_reached)
C       calculate Rydberg partition function sums from the range
C       of minimum principal quantum numbers (given by nmin to nmin_max)
C       to nmax of the sum 2 n^2 exp(alpha/n^2) times
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
C       b(nb), where nb = 5 determines the mhd occupation probability.
C       output quantities:
C       qryd(nmin_max), qryda(nmin_max), qrydb(nb, nmin_max),
C         qryda2(nmin_max), qrydab(nb, nmin_max), qrydb2(nb, nb, nmin_max) are
C         the resulting sum plus partial derivatives wrt a and b.
C       nmax_reached is the returned maximum principal quantum number
C       that is used in the sum taking account of the convergence criteria.
      implicit none
      integer nmin, nmin_max, nmax, nb, nmax_reached
      double precision eps_factor, a, b(nb),
     &  qryd(nmin_max), qryda(nmin_max), qrydb(nb, nmin_max),
     &  qryda2(nmin_max), qrydab(nb, nmin_max),
     &  qrydb2(nb, nb, nmin_max)
      logical ifpl, ifmhd, ifneutral
C       internal variables:
      double precision eps
      parameter (eps=1.d-10)
      integer nb_local
      parameter(nb_local = 5)
      double precision rn2, arg, summand, summanda, summanda2,
     &  kn, hprime, gprime, occupation,
     &  lnoccupation, dlnoccupation(nb_local), dsummand(nb_local),
     &  dsummanda(nb_local), dsummand2(nb_local,nb_local)
      integer n, nzero
      nmax_reached = 0
C       define for first time through
      lnoccupation = 1.d0
      summand = 1.d0
C       n.b. note that returns zero for all values if nmax < nmin, which
C       is proper thing to do.
      n = nmin
      do while(n.le.nmax.and.
     &    (n.eq.nmin.or.(summand.gt.eps/eps_factor.or.
     &    (ifmhd.and.lnoccupation.ge.-1.d0))))
C         avoid integer overflow by doing this in double precision
        rn2 = dble(n)*dble(n)
        arg = a/rn2
        if(ifpl) then
          if(arg.gt.0.01d0) then
C             lose a maximum of 4 significant digits
            summand = 2.d0*rn2*(exp(arg) - (1.d0 + arg))
            summanda = 2.d0*(exp(arg) - 1.d0)
            summanda2 = 2.d0*exp(arg)/rn2
          else
            summand = rn2*arg*arg*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0*(
     &        1.d0 + arg/10.d0*(
     &        1.d0 + arg/11.d0)))))))))
            summanda = 2.d0*arg*(
     &        1.d0 + arg/2.d0*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0*(
     &        1.d0 + arg/10.d0)))))))))
            summanda2 = 2.d0*(
     &        1.d0 + arg*(
     &        1.d0 + arg/2.d0*(
     &        1.d0 + arg/3.d0*(
     &        1.d0 + arg/4.d0*(
     &        1.d0 + arg/5.d0*(
     &        1.d0 + arg/6.d0*(
     &        1.d0 + arg/7.d0*(
     &        1.d0 + arg/8.d0*(
     &        1.d0 + arg/9.d0)))))))))/rn2
          endif
        else
          summand = 2.d0*rn2*exp(arg)
          summanda = 2.d0*exp(arg)
          summanda2 = 2.d0*exp(arg)/rn2
        endif
        if(ifmhd) then
C           quantum correction K_n see Hummer and Mihalas eq. 4.24
          if(n.le.3) then
            kn = 1.d0
          else
            kn = 
     &        (16.d0*rn2*(dble(n) + 7.d0/6.d0))/
     &        (dble(3*(n+1))*dble(n+1)*(rn2 + dble(n) + 0.5d0))
          endif
          hprime = (16.d0/(3.d0*dble(n)*kn))**1.5d0
C           mhd occupation probability for neutral-ion and ion-ion interactions
          dlnoccupation(5) = -(dble(n))**7.5d0*hprime
          lnoccupation = b(5)*dlnoccupation(5)
          if(ifneutral) then
C             mhd occupation probability for neutral-neutral interactions
C             l = l_max = n-1 --> 1/2 * (3n^2 - l(l+1)) = n^2*(1+1/2n)
            gprime = (1.d0 + 0.5d0/dble(n))*rn2
C             l = l_min = 0 --> 1/2 * (3n^2 - l(l+1)) = n^2*(3/2) (test case)
!            gprime = (1.5d0)*rn2
            lnoccupation = lnoccupation -
     &        (b(1) + gprime*(b(2) + gprime*(b(3) + gprime*b(4))))
            dlnoccupation(1) = -1.d0
            dlnoccupation(2) = dlnoccupation(1)*gprime
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
            dsummand(1) = summand*dlnoccupation(1)
            dsummand(2) = summand*dlnoccupation(2)
            dsummand(3) = summand*dlnoccupation(3)
            dsummand(4) = summand*dlnoccupation(4)
            dsummanda(1) = summanda*dlnoccupation(1)
            dsummanda(2) = summanda*dlnoccupation(2)
            dsummanda(3) = summanda*dlnoccupation(3)
            dsummanda(4) = summanda*dlnoccupation(4)
            dsummand2(1,1) = summand*dlnoccupation(1)*dlnoccupation(1)
            dsummand2(2,1) = summand*dlnoccupation(2)*dlnoccupation(1)
            dsummand2(3,1) = summand*dlnoccupation(3)*dlnoccupation(1)
            dsummand2(4,1) = summand*dlnoccupation(4)*dlnoccupation(1)
            dsummand2(5,1) = summand*dlnoccupation(5)*dlnoccupation(1)
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
          if(n.le.nmin_max) then
            qrydb(5,n) = dsummand(5)
            qrydab(5,n) = dsummanda(5)
            qrydb2(5,5,n) = dsummand2(5,5)
            if(ifneutral) then
              qrydb(1,n) = dsummand(1)
              qrydb(2,n) = dsummand(2)
              qrydb(3,n) = dsummand(3)
              qrydb(4,n) = dsummand(4)
              qrydab(1,n) = dsummanda(1)
              qrydab(2,n) = dsummanda(2)
              qrydab(3,n) = dsummanda(3)
              qrydab(4,n) = dsummanda(4)
              qrydb2(1,1,n) = dsummand2(1,1)
              qrydb2(2,1,n) = dsummand2(2,1)
              qrydb2(3,1,n) = dsummand2(3,1)
              qrydb2(4,1,n) = dsummand2(4,1)
              qrydb2(5,1,n) = dsummand2(5,1)
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
            qrydb(5,nmin_max) = qrydb(5,nmin_max) +
     &        dsummand(5)
            qrydab(5,nmin_max) = qrydab(5,nmin_max) +
     &        dsummanda(5)
            qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +
     &        dsummand2(5,5)
            if(ifneutral) then
              qrydb(1,nmin_max) = qrydb(1,nmin_max) +
     &          dsummand(1)
              qrydb(2,nmin_max) = qrydb(2,nmin_max) +
     &          dsummand(2)
              qrydb(3,nmin_max) = qrydb(3,nmin_max) +
     &          dsummand(3)
              qrydb(4,nmin_max) = qrydb(4,nmin_max) +
     &          dsummand(4)
              qrydab(1,nmin_max) = qrydab(1,nmin_max) +
     &          dsummanda(1)
              qrydab(2,nmin_max) = qrydab(2,nmin_max) +
     &          dsummanda(2)
              qrydab(3,nmin_max) = qrydab(3,nmin_max) +
     &          dsummanda(3)
              qrydab(4,nmin_max) = qrydab(4,nmin_max) +
     &          dsummanda(4)
              qrydb2(1,1,nmin_max) = qrydb2(1,1,nmin_max) +
     &          dsummand2(1,1)
              qrydb2(2,1,nmin_max) = qrydb2(2,1,nmin_max) +
     &          dsummand2(2,1)
              qrydb2(3,1,nmin_max) = qrydb2(3,1,nmin_max) +
     &          dsummand2(3,1)
              qrydb2(4,1,nmin_max) = qrydb2(4,1,nmin_max) +
     &          dsummand2(4,1)
              qrydb2(5,1,nmin_max) = qrydb2(5,1,nmin_max) +
     &          dsummand2(5,1)
              qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +
     &          dsummand2(2,2)
              qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +
     &          dsummand2(3,2)
              qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +
     &          dsummand2(4,2)
              qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +
     &          dsummand2(5,2)
              qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +
     &          dsummand2(3,3)
              qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +
     &          dsummand2(4,3)
              qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +
     &          dsummand2(5,3)
              qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +
     &          dsummand2(4,4)
              qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +
     &          dsummand2(5,4)
            endif
          endif
        endif
        if(n.le.nmin_max) then
          qryd(n) = summand
          qryda(n) = summanda
          qryda2(n) = summanda2
        else
          qryd(nmin_max) = qryd(nmin_max) + summand
          qryda(nmin_max) = qryda(nmin_max) + summanda
          qryda2(nmin_max) = qryda2(nmin_max) + summanda2
        endif
        n = n + 1
      enddo
      nmax_reached = max(nmax_reached, n-1)
C       zero remainder of qryd if necessary.
      do nzero = n,nmin_max
        qryd(nzero) = 0.d0
        qryda(nzero) = 0.d0
        qryda2(nzero) = 0.d0
        if(ifmhd) then
          qrydb(5,nzero) = 0.d0
          qrydab(5,nzero) = 0.d0
          qrydb2(5,5,nzero) = 0.d0
          if(ifneutral) then
            qrydb(1,nzero) = 0.d0
            qrydb(2,nzero) = 0.d0
            qrydb(3,nzero) = 0.d0
            qrydb(4,nzero) = 0.d0
            qrydab(1,nzero) = 0.d0
            qrydab(2,nzero) = 0.d0
            qrydab(3,nzero) = 0.d0
            qrydab(4,nzero) = 0.d0
            qrydb2(1,1,nzero) = 0.d0
            qrydb2(2,1,nzero) = 0.d0
            qrydb2(3,1,nzero) = 0.d0
            qrydb2(4,1,nzero) = 0.d0
            qrydb2(5,1,nzero) = 0.d0
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
        endif
      enddo
C       convert from summand to sum for lower indices
      do n = nmin_max-1,nmin,-1
        qryd(n) = qryd(n) + qryd(n+1)
        qryda(n) = qryda(n) + qryda(n+1)
        qryda2(n) = qryda2(n) + qryda2(n+1)
        if(ifmhd) then
          qrydb(5,n) = qrydb(5,n) + qrydb(5,n+1)
          qrydab(5,n) = qrydab(5,n) + qrydab(5,n+1)
          qrydb2(5,5,n) = qrydb2(5,5,n) + qrydb2(5,5,n+1)
          if(ifneutral) then
            qrydb(1,n) = qrydb(1,n) + qrydb(1,n+1)
            qrydb(2,n) = qrydb(2,n) + qrydb(2,n+1)
            qrydb(3,n) = qrydb(3,n) + qrydb(3,n+1)
            qrydb(4,n) = qrydb(4,n) + qrydb(4,n+1)
            qrydab(1,n) = qrydab(1,n) + qrydab(1,n+1)
            qrydab(2,n) = qrydab(2,n) + qrydab(2,n+1)
            qrydab(3,n) = qrydab(3,n) + qrydab(3,n+1)
            qrydab(4,n) = qrydab(4,n) + qrydab(4,n+1)
            qrydb2(1,1,n) = qrydb2(1,1,n) + qrydb2(1,1,n+1)
            qrydb2(2,1,n) = qrydb2(2,1,n) + qrydb2(2,1,n+1)
            qrydb2(3,1,n) = qrydb2(3,1,n) + qrydb2(3,1,n+1)
            qrydb2(4,1,n) = qrydb2(4,1,n) + qrydb2(4,1,n+1)
            qrydb2(5,1,n) = qrydb2(5,1,n) + qrydb2(5,1,n+1)
            qrydb2(2,2,n) = qrydb2(2,2,n) + qrydb2(2,2,n+1)
            qrydb2(3,2,n) = qrydb2(3,2,n) + qrydb2(3,2,n+1)
            qrydb2(4,2,n) = qrydb2(4,2,n) + qrydb2(4,2,n+1)
            qrydb2(5,2,n) = qrydb2(5,2,n) + qrydb2(5,2,n+1)
            qrydb2(3,3,n) = qrydb2(3,3,n) + qrydb2(3,3,n+1)
            qrydb2(4,3,n) = qrydb2(4,3,n) + qrydb2(4,3,n+1)
            qrydb2(5,3,n) = qrydb2(5,3,n) + qrydb2(5,3,n+1)
            qrydb2(4,4,n) = qrydb2(4,4,n) + qrydb2(4,4,n+1)
            qrydb2(5,4,n) = qrydb2(5,4,n) + qrydb2(5,4,n+1)
          endif
        endif
      enddo
      end
