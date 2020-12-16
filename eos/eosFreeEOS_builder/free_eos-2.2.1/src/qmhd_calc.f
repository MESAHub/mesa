C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: qmhd_calc.f 627 2007-07-19 01:25:19Z airwin $
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
      subroutine qmhd_calc(ifpi_fit, weight_ion,
     &  ifpl, ifapprox, eps_factor, nmin_rydberg, nmax,
     &  tl, iz, rfactor, nlevel, weight, neff, ell, nx, x,
     &  qmhd, qmhdt, qmhdx, qmhdt2, qmhdtx, qmhdx2)
C       calculate mhd partition function sums for nlevel terms
C       summand is w_pl*w_mhd*weight*exp(a), where a = c_2*R/neff^2/T
C       w_pl is the planck_larkin occupation probability (ifpl true)
C       w_mhd is the mhd occupation probability.
C       w_pl is given by (1 - exp(-a)*(1+a)).
C       The *ln* of the mhd occupation probability is given by
C       -4pi/3*[x(1) + x(2)*r + x(3)*r^2 +
C       x(4)*r^3+ x(5)*rion^3], where
C       r = bohr/2*(3neff^2-ell(ell-1))
C       rion^3 = 16*(iz^0.5 e^2 neff^2)^3/(K^0.5 hc Rydberg)
C       K(neff) is a quantum correction term
C       input quantities:
C       ifpi_fit = 2,_use_best fit to Saumon table
C       ifpi_fit = 1,_use_best fit to opal table + extensions
C       ifpi_fit = 0,_use_best fit to original MDH table.
C       ifpl (logical) controls whether to_use_planck-larkin
C         occupation probability.
C       weight_ion = statistical weight of next higher ion =
C         the factor that multiplies Rydberg part of partition function
C       ifapprox controls whether to_use_approximation or exact sum for
C         rydberg states.
C       eps_factor = exp(-c2 bion(ion)/t) is used to terminate rydberg
C         state sum.
C       nmin_rydberg is the lowest principal quantum number treated as
C         a rydberg level.
C       nmax is maximum principal number cutoff (only an issue for this
C         call to qstar_calc when ifapprox = .false.)  Also used here
C         to cutoff exact non-Rydberg He calculation.
C       tl = ln T
C       iz = core charge = 1 for neutrals, 2 for first ions, etc.
C       rfactor = finite mass factor for Rydberg = 1/(1+m_e/M)
C       nlevel = number of levels included in sum before_use_rydberg formulas
C         (other routines) to finish off to principal quantum number
C         of infinity.
C       weight(nlevel) is the statistical weight
C       neff(nlevel) is the effective principal quantum number
C       ell(nlevel) is the azimuthal quantum number
C       x(nx=5) is ordered as extrasum(4), extrasum(3), extrasum(2), extrasum(1)
C         extrasum(nextrasum-1)
C       output quantities:
C       qmhd, qmhdt, qmhdx(nx), qmhdt2, qmhdtx(nx), qmhdx2(nx, nx) are
C         the resulting sum plus partial derivatives wrt tl and x.
      implicit none
      include 'constants.h'
      include 'pi_fit.h'
      logical ifpl, ifneutral, ifapprox
      integer ifpi_fit, iz, nlevel, nx
      double precision weight_ion, eps_factor, tl, rfactor,
     &  weight(nlevel), neff(nlevel), ell(nlevel),
     &  x(nx),
     &  qmhd, qmhdt, qmhdx(nx), qmhdt2, qmhdtx(nx), qmhdx2(nx, nx)
C       internal variables
      integer nx_local
      parameter(nx_local = 5)
      double precision t, rn2, energy, arg, summand, summandt, summandt2,
     &  kn, occupation,
     &  lnoccupation, dlnoccupation(nx_local), dsummand(nx_local),
     &  dsummandt(nx_local), dsummand2(nx_local,nx_local),
     &  const, rionconst3, r, rnapprox, ellapprox
      parameter (const = 4.d0*pi/3.d0)
      parameter (rionconst3 = const*16.d0*echarge*echarge*
     &  echarge*echarge*echarge*echarge)
      integer ilevel, nmin_rydberg, nmax, nmin_max
      integer min_nmin_max
C       minimum value of principal quantum number where approximations are
C         used
      parameter(min_nmin_max = 3)
C       maximum value of mininum principal quantum number
      parameter (nmin_max = 11)
      double precision
     &  qstar(nmin_max), qstart(nmin_max), qstarx(nx_local,nmin_max),
     &  qstart2(nmin_max), qstartx(nx_local,nmin_max),
     &  qstarx2(nx_local,nx_local,nmin_max)

      double precision neutral_factor, ionized_factor
      if(nmin_rydberg.gt.nmin_max)
     &  stop 'qmhd_calc: nmin_rydberg too large'
      if(nx_local.ne.nx) stop 'qreal_calc: invalid nx'
      ifneutral = iz.eq.1
      if(nmax.ge.nmin_rydberg) then
C         n.b max(...) in argument list protects against using approximation
C         inappropriately, (i.e., for n < 3)
        call qstar_calc(ifpi_fit,
     &    ifpl, .true. , ifneutral, ifapprox,
     &    eps_factor, nmin_rydberg, max(min_nmin_max,nmin_rydberg),
     &    nmax,
     &    iz, tl, x, nx,
     &    qstar, qstart, qstarx, qstart2, qstartx, qstarx2)
        qmhd = weight_ion*qstar(nmin_rydberg)
        qmhdt = weight_ion*qstart(nmin_rydberg)
        qmhdt2 = weight_ion*qstart2(nmin_rydberg)
        qmhdx(5) = weight_ion*qstarx(5,nmin_rydberg)
        qmhdtx(5) = weight_ion*qstartx(5,nmin_rydberg)
        qmhdx2(5,5) = weight_ion*qstarx2(5,5,nmin_rydberg)
        if(iz.eq.1) then
          qmhdx(1) = weight_ion*qstarx(1,nmin_rydberg)
          qmhdx(2) = weight_ion*qstarx(2,nmin_rydberg)
          qmhdx(3) = weight_ion*qstarx(3,nmin_rydberg)
          qmhdx(4) = weight_ion*qstarx(4,nmin_rydberg)
          qmhdtx(1) = weight_ion*qstartx(1,nmin_rydberg)
          qmhdtx(2) = weight_ion*qstartx(2,nmin_rydberg)
          qmhdtx(3) = weight_ion*qstartx(3,nmin_rydberg)
          qmhdtx(4) = weight_ion*qstartx(4,nmin_rydberg)
          qmhdx2(1,1) = weight_ion*qstarx2(1,1,nmin_rydberg)
          qmhdx2(2,1) = weight_ion*qstarx2(2,1,nmin_rydberg)
          qmhdx2(3,1) = weight_ion*qstarx2(3,1,nmin_rydberg)
          qmhdx2(4,1) = weight_ion*qstarx2(4,1,nmin_rydberg)
          qmhdx2(5,1) = weight_ion*qstarx2(5,1,nmin_rydberg)
          qmhdx2(2,2) = weight_ion*qstarx2(2,2,nmin_rydberg)
          qmhdx2(3,2) = weight_ion*qstarx2(3,2,nmin_rydberg)
          qmhdx2(4,2) = weight_ion*qstarx2(4,2,nmin_rydberg)
          qmhdx2(5,2) = weight_ion*qstarx2(5,2,nmin_rydberg)
          qmhdx2(3,3) = weight_ion*qstarx2(3,3,nmin_rydberg)
          qmhdx2(4,3) = weight_ion*qstarx2(4,3,nmin_rydberg)
          qmhdx2(5,3) = weight_ion*qstarx2(5,3,nmin_rydberg)
          qmhdx2(4,4) = weight_ion*qstarx2(4,4,nmin_rydberg)
          qmhdx2(5,4) = weight_ion*qstarx2(5,4,nmin_rydberg)
        endif
      else
C         nmax cutoff allows no rydberg part of helium partition function
        qmhd = 0.d0
        qmhdt = 0.d0
        qmhdt2 = 0.d0
        qmhdx(5) = 0.d0
        qmhdtx(5) = 0.d0
        qmhdx2(5,5) = 0.d0
        if(iz.eq.1) then
          qmhdx(1) = 0.d0
          qmhdx(2) = 0.d0
          qmhdx(3) = 0.d0
          qmhdx(4) = 0.d0
          qmhdtx(1) = 0.d0
          qmhdtx(2) = 0.d0
          qmhdtx(3) = 0.d0
          qmhdtx(4) = 0.d0
          qmhdx2(1,1) = 0.d0
          qmhdx2(2,1) = 0.d0
          qmhdx2(3,1) = 0.d0
          qmhdx2(4,1) = 0.d0
          qmhdx2(5,1) = 0.d0
          qmhdx2(2,2) = 0.d0
          qmhdx2(3,2) = 0.d0
          qmhdx2(4,2) = 0.d0
          qmhdx2(5,2) = 0.d0
          qmhdx2(3,3) = 0.d0
          qmhdx2(4,3) = 0.d0
          qmhdx2(5,3) = 0.d0
          qmhdx2(4,4) = 0.d0
          qmhdx2(5,4) = 0.d0
        endif
      endif
      if(ifpi_fit.eq.0) then
        neutral_factor = exp(pi_fitx_neutral_ln_original/3.d0)
        ionized_factor = exp(pi_fitx_ion_ln_original)
      elseif(ifpi_fit.eq.1) then
        neutral_factor = exp(pi_fitx_neutral_ln/3.d0)
        ionized_factor = exp(pi_fitx_ion_ln)
      elseif(ifpi_fit.eq.2) then
        neutral_factor = exp(pi_fitx_neutral_ln_saumon/3.d0)
        ionized_factor = exp(pi_fitx_ion_ln_saumon)
      else
        stop 'qmhd_calc: bad ifpi_fit value'
      endif
      t = exp(tl)
      do ilevel = 1, nlevel
C        n.b. checked that for helium energy levels, nint(neff) = n for
C        *all* excited states. int(neff)+1 doesn't cut it because quantum
C        defect sometimes slightly negative (level slightly too high).
        if(nint(neff(ilevel)).le.min(nmin_rydberg-1,nmax)) then
          rn2 = neff(ilevel)*neff(ilevel)
C          energy below ionization continuum
          energy = ergspercmm1*rydberg*rfactor*dble(iz*iz)/rn2
          arg = energy/(boltzmann*t)
          if(ifpl) then
            if(arg.gt.0.01d0) then
C             lose a maximum of 4 significant digits
              summand = weight(ilevel)*(exp(arg) - (1.d0 + arg))
              summandt = -arg*weight(ilevel)*(exp(arg) - 1.d0)
              summandt2 = arg*weight(ilevel)*
     &          (exp(arg)*(1.d0+arg) - 1.d0)
            else
              summand = 0.5d0*weight(ilevel)*arg*arg*(
     &          1.d0 + arg/3.d0*(
     &          1.d0 + arg/4.d0*(
     &          1.d0 + arg/5.d0*(
     &          1.d0 + arg/6.d0*(
     &          1.d0 + arg/7.d0*(
     &          1.d0 + arg/8.d0*(
     &          1.d0 + arg/9.d0*(
     &          1.d0 + arg/10.d0*(
     &          1.d0 + arg/11.d0)))))))))
              summandt = -weight(ilevel)*arg*arg*(
     &          1.d0 + arg/2.d0*(
     &          1.d0 + arg/3.d0*(
     &          1.d0 + arg/4.d0*(
     &          1.d0 + arg/5.d0*(
     &          1.d0 + arg/6.d0*(
     &          1.d0 + arg/7.d0*(
     &          1.d0 + arg/8.d0*(
     &          1.d0 + arg/9.d0*(
     &          1.d0 + arg/10.d0)))))))))
              summandt2 = weight(ilevel)*arg*arg*(
     &          2.d0 + arg/2.d0*(
     &          3.d0 + arg/3.d0*(
     &          4.d0 + arg/4.d0*(
     &          5.d0 + arg/5.d0*(
     &          6.d0 + arg/6.d0*(
     &          7.d0 + arg/7.d0*(
     &          8.d0 + arg/8.d0*(
     &          9.d0 + arg/9.d0*(
     &          10.d0 + arg*11.d0/10.d0)))))))))
            endif
          else
            summand = weight(ilevel)*exp(arg)
            summandt = -weight(ilevel)*arg*exp(arg)
            summandt2 = weight(ilevel)*arg*exp(arg)*(1.d0+arg)
          endif
C          quantum correction K_n see Hummer and Mihalas eq. 4.24
          if(neff(ilevel).le.3) then
            kn = 1.d0
          else
            kn = 
     &        (16.d0*rn2*(neff(ilevel) + 7.d0/6.d0))/
     &        (3.d0*(neff(ilevel)+1.d0)*(neff(ilevel)+1.d0)*
     &        (rn2 + neff(ilevel) + 0.5d0))
          endif
C          mhd occupation probability for neutral-ion and ion-ion interactions
          dlnoccupation(5) = -ionized_factor*rionconst3*
     &      (sqrt(dble(iz)/kn)/energy)**3
          lnoccupation = x(5)*dlnoccupation(5)
          if(iz.eq.1) then
C            first expression is "exact", but details of this part of
C            occupation probability don't matter (first part much larger) so
C            simplify to be consistent with simplified presentation
C            r = neutral_factor*0.5d0*bohr*
C     &        (3.d0*rn2 - ell(ilevel)*(ell(ilevel)+1.d0))
            rnapprox = dble(nint(neff(ilevel)))
            ellapprox = rnapprox - 1.d0
            r = neutral_factor*0.5d0*bohr*
     &        (3.d0*rnapprox*rnapprox - ellapprox*(ellapprox+1.d0))
            dlnoccupation(1) = -const
            dlnoccupation(2) = dlnoccupation(1)*r*3.d0
            dlnoccupation(3) = dlnoccupation(2)*r
            dlnoccupation(4) = dlnoccupation(3)*r/3.d0
            lnoccupation = lnoccupation +
     &        dlnoccupation(1)*x(1) +
     &        dlnoccupation(2)*x(2) +
     &        dlnoccupation(3)*x(3) +
     &        dlnoccupation(4)*x(4)
          endif
          occupation = exp(lnoccupation)
          summand = summand*occupation
          summandt = summandt*occupation
          summandt2 = summandt2*occupation
          dsummand(5) = summand*dlnoccupation(5)
          dsummandt(5) = summandt*dlnoccupation(5)
          dsummand2(5,5) = summand*dlnoccupation(5)*dlnoccupation(5)
          if(iz.eq.1) then
            dsummand(1) = summand*dlnoccupation(1)
            dsummand(2) = summand*dlnoccupation(2)
            dsummand(3) = summand*dlnoccupation(3)
            dsummand(4) = summand*dlnoccupation(4)
            dsummandt(1) = summandt*dlnoccupation(1)
            dsummandt(2) = summandt*dlnoccupation(2)
            dsummandt(3) = summandt*dlnoccupation(3)
            dsummandt(4) = summandt*dlnoccupation(4)
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
          qmhd = qmhd + summand
          qmhdt = qmhdt + summandt
          qmhdt2 = qmhdt2 + summandt2
          qmhdx(5) = qmhdx(5) + dsummand(5)
          qmhdtx(5) = qmhdtx(5) + dsummandt(5)
          qmhdx2(5,5) = qmhdx2(5,5) + dsummand2(5,5)
          if(iz.eq.1) then
            qmhdx(1) = qmhdx(1) + dsummand(1)
            qmhdx(2) = qmhdx(2) + dsummand(2)
            qmhdx(3) = qmhdx(3) + dsummand(3)
            qmhdx(4) = qmhdx(4) + dsummand(4)
            qmhdtx(1) = qmhdtx(1) + dsummandt(1)
            qmhdtx(2) = qmhdtx(2) + dsummandt(2)
            qmhdtx(3) = qmhdtx(3) + dsummandt(3)
            qmhdtx(4) = qmhdtx(4) + dsummandt(4)
            qmhdx2(1,1) = qmhdx2(1,1) + dsummand2(1,1)
            qmhdx2(2,1) = qmhdx2(2,1) + dsummand2(2,1)
            qmhdx2(3,1) = qmhdx2(3,1) + dsummand2(3,1)
            qmhdx2(4,1) = qmhdx2(4,1) + dsummand2(4,1)
            qmhdx2(5,1) = qmhdx2(5,1) + dsummand2(5,1)
            qmhdx2(2,2) = qmhdx2(2,2) + dsummand2(2,2)
            qmhdx2(3,2) = qmhdx2(3,2) + dsummand2(3,2)
            qmhdx2(4,2) = qmhdx2(4,2) + dsummand2(4,2)
            qmhdx2(5,2) = qmhdx2(5,2) + dsummand2(5,2)
            qmhdx2(3,3) = qmhdx2(3,3) + dsummand2(3,3)
            qmhdx2(4,3) = qmhdx2(4,3) + dsummand2(4,3)
            qmhdx2(5,3) = qmhdx2(5,3) + dsummand2(5,3)
            qmhdx2(4,4) = qmhdx2(4,4) + dsummand2(4,4)
            qmhdx2(5,4) = qmhdx2(5,4) + dsummand2(5,4)
          endif
        endif
      enddo
      end
