C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fjs_pi.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine fjs_pi(ifnr, maxfjs_aux, ifmodified,
     &  t, nux, nuxf, nuxt, nuy, nuyf, nuyt, nuz, nuzf, nuzt,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, ne, nef, net,
     &  dvh, dvhf, dvht, dvh_aux,
     &  dvhe1, dvhe1f, dvhe1t, dvhe1_aux,
     &  dvhe2, dvhe2f, dvhe2t, dvhe2_aux,
     &  dve, dvef, dvet, dve_aux)
C       Fritz Swenson pressure ionization term (generalization of EFF).
C       the call to fjs_pi yields the ln change to the equilibrium
C       constant for H+, He+, He++ (relative to He+), and generic metal ion 
C       relative to next lower ion.
C       the later entry point fjs_pi_end yields the change
C       to the pressure, entropy, and internal energy.
C       n.b. the functional form for this free energy is such that the
C         free energy per unit volume is exactly dpi.
C         all other values are derived from the expression for dpi below.
C       ifmodified <= 0 (use original FJS expression)
C       ifmodified > 0 (use expression modified for best fit to opal
C         and extension.)
C       nux, nuy, and nuz and derivatives are number/volume of maximum 
C       possible ionization electrons for hydrogen, helium, and metals.
C       h, hd, he, ne and derivatives are H+, He+, He++, e- number/volume.
C       n.b. all these values are number densities, *not* number densities
C       divided by rho*avogadro as in old code (except for final call to
C         fjs_pi_end).
C       input quantities:
C       ifnr = 0, calculate fl and tl derivatives of output using chain rule
C         and auxf and auxt
C       ifnr = 1, calculate aux derivatives of output with fl, tl fixed
C       ifnr = 2, calculate fl and tl derivatives of output with aux fixed.
C       ifnr = 3 is combination of ifnr = 1 and ifnr = 2.
      implicit none
      include 'constants.h'
      integer ifnr, maxfjs_aux
      double precision t, rho, rf, rt,
     &  nux, nuxf, nuxt, nuy, nuyf, nuyt, nuz, nuzf, nuzt,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, ne, nef, net,
     &  api1, api2, api3, api4, chipi1, chipi2, chipi3, chipi4,
     &  chpi1, chpi2, chpi3, chpi4,
     &  omht, omh, omhe1, omhe1t, omhe2, omhe2t, ome, omet,
     &  ppi1, ppi2, ppi3, ppi4,
     &  ppi1f, ppi2f, ppi3f, ppi4f,
     &  ppi1t, ppi2t, ppi3t, ppi4t,
     &  ppi11, ppi12, ppi13, ppi14,
     &  ppi21, ppi22, ppi23, ppi24,
     &  ppi31, ppi32, ppi33, ppi34,
     &  ppi41, ppi42, ppi43, ppi44,
     &  ppi, ppif, ppit, spi, spif, spit, upi,
     &  ne0, ne0f, ne0t, onedkt,
     &  dvh, dvhf, dvht, dvh_aux(maxfjs_aux),
     &  dvhe1, dvhe1f, dvhe1t, dvhe1_aux(maxfjs_aux),
     &  dvhe2, dvhe2f, dvhe2t, dvhe2_aux(maxfjs_aux),
     &  dve, dvef, dvet, dve_aux(maxfjs_aux),
     &  api1_orig, chipi1_orig, api1_mod,
     &  api2_orig, chipi2_orig, api2_mod,
     &  api3_orig, chipi3_orig, api3_mod,
     &  api4_orig, chipi4_orig, api4_mod
      integer ifmodified, ifmodifiedold
C      need invalid value
      data ifmodifiedold/-1000/
C      original fjs set.
      data api1_orig/1.2d-8/,api2_orig/0.45d-8/,
     &  api3_orig/0.3d-8/,api4_orig/0.57d-8/,
     &  chipi1_orig/0.5d0/,chipi2_orig/36.0d0/,
     &  chipi3_orig/4.0d0/,chipi4_orig/8.0d0/
C      working set + fudge that gives reasonable fit to opal + extension
      data api1_mod/1.2d-8/,api2_mod/0.45d-8/,
     &  api3_mod/0.3d-8/,api4_mod/0.57d-8/
C      fudge to best fit extended m* results.
      double precision fudge_ln(4)
C      n.b. (3) must be sufficiently larger than (2) otherwise, He+ is
C      never pressure ionized to He++.
!old  data fudge_ln /1.d0, 4.5d0, 5.5d0, -15.0d0/
C      small adjustment made to fit pure He opal tables.  Didn't degrade
C      rest of table fit or extension fit.
C      data fudge_ln /1.d0, 5.0d0, 5.5d0, -5.0d0/
C      now adjust again against best eos rather than opal especially for
C      zero X case
      data fudge_ln /1.d0, 2.5d0, 3.0d0, -5.0d0/
      double precision free_pi, free_pif, free_pit, free_pi2_daux(4)
      logical ifnr13, ifnr23
      save
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      ifnr23 = ifnr.eq.2.or.ifnr.eq.3
C      set up some oft used variables
      if(ifmodified.ne.ifmodifiedold) then
        ifmodifiedold = ifmodified
        if(ifmodified.le.0) then
          api1 = api1_orig
          api2 = api2_orig
          api3 = api3_orig
          api4 = api4_orig
          chipi1 = chipi1_orig
          chipi2 = chipi2_orig
          chipi3 = chipi3_orig
          chipi4 = chipi4_orig
        else
          api1 = api1_mod*exp(fudge_ln(1)/3.d0)
          api2 = api2_mod*exp(fudge_ln(2)/3.d0)
          api3 = api3_mod*exp(fudge_ln(3)/3.d0)
          api4 = api4_mod*exp(fudge_ln(4)/3.d0)
          chipi1 = 0.d0
          chipi2 = 0.d0
          chipi3 = 0.d0
          chipi4 = 0.d0
        endif
C        derived constants.
        omht = api1**3*boltzmann
        omhe1t = api2**3*boltzmann
        omhe2t = api3**3*boltzmann
        omet = api4**3*boltzmann
        chpi1 = api1**3*8.d0*pi*ergsperev*chipi1
        chpi2 = api2**3*8.d0*pi*ergsperev*chipi2
        chpi3 = api3**3*8.d0*pi*ergsperev*chipi3
        chpi4 = api4**3*8.d0*pi*ergsperev*chipi4
      endif
      onedkt = 1.d0/(boltzmann*t)
      omh = omht*t+chpi1
      omhe1 = omhe1t*t+chpi2
      omhe2 = omhe2t*t+chpi3
      ome = omet*t+chpi4
C      compute chemical potential of lower ion - ion - ne 
C      (all divided by kt)
      dvh = ((ne+h+nuy-hd-2.d0*he)*omh+2.d0*ne*ome)*onedkt
      dvhe1 = ((ne+nux-h+hd)*omhe1+2.d0*ne*ome)*onedkt
      dvhe2 = ((hd-ne-nux+h)*omhe1+(nuy+2.d0*(ne+nux-h-hd))*omhe2
     &  + 2.d0*ne*ome)*onedkt
      dve = (h*omh+hd*omhe1+2.d0*he*omhe2+2.d0*ne*ome)*onedkt
      if(ifnr.eq.0) then
C        partial derivatives assuming all number densities are a
C        function of fl and tl.
        dvhf = ((nef+hf+nuyf-hdf-2.d0*hef)*omh+2.d0*nef*ome)
     &    *onedkt
        dvht = -dvh +((net+ht+nuyt-hdt-2.d0*het)
     &    *omh+omht*t*(ne+h+nuy-hd-2.d0*he)+2.d0*omet*t*ne+2.d0*ome*net)
     &    *onedkt
        dvhe1f = ((nef+nuxf-hf+hdf)*omhe1+2.d0*nef*ome)
     &    *onedkt
        dvhe1t = -dvhe1+((net+nuxt-ht+hdt)*omhe1+omhe1t*t*(ne+nux-h+hd)
     &    +2.d0*omet*t*ne+2.d0*ome*net)*onedkt
        dvhe2f = (omhe1*(hdf-nef-nuxf+hf)+
     &    omhe2*(nuyf + 2.d0*(nef+nuxf-hf-hdf))
     &    +2.d0*ome*nef)*onedkt
        dvhe2t = -dvhe2+((hd-ne-nux+h)*omhe1t*t
     &    +(hdt-net-nuxt+ht)*omhe1+(nuy+2.d0*(ne+nux-h-hd))*omhe2t*t
     &    +(nuyt+2.d0*(net+nuxt-ht-hdt))*omhe2+
     &    2.d0*omet*t*ne+2.d0*ome*net)
     &    *onedkt
        dvef = (omh*hf+omhe1*hdf+2.d0*omhe2*hef
     &    +2.d0*ome*nef)*onedkt
        dvet = -dve+(omht*t*h+omhe1t*t*hd+2.d0*omhe2t*t*he+
     &    2.d0*omet*t*ne
     &    +omh*ht+omhe1*hdt+2.d0*omhe2*het+2.d0*ome*net)*onedkt
      elseif(ifnr13) then
C        derivative (1) is wrt to h
C        derivative (2) is wrt to hd
C        derivative (3) is wrt to he
C        nux, nuy, nuz are proportional to rho.
C        derivative (4) is wrt ln rho.
        dvh_aux(1) = omh*onedkt
        dvh_aux(2) = -omh*onedkt
        dvh_aux(3) = -2.d0*omh*onedkt
        dvh_aux(4) = nuy*omh*onedkt
        dvhe1_aux(1) = -omhe1*onedkt
        dvhe1_aux(2) = omhe1*onedkt
        dvhe1_aux(3) = 0.d0
        dvhe1_aux(4) = nux*omhe1*onedkt
        dvhe2_aux(1) = (omhe1-2.d0*omhe2)*onedkt
        dvhe2_aux(2) = (omhe1-2.d0*omhe2)*onedkt
        dvhe2_aux(3) = 0.d0
        dvhe2_aux(4) = (-nux*omhe1+(nuy+2.d0*nux)*omhe2)*onedkt
        dve_aux(1) = omh*onedkt
        dve_aux(2) = omhe1*onedkt
        dve_aux(3) = 2.d0*omhe2*onedkt
        dve_aux(4) = 0.d0
      endif
      if(ifnr23) then
C        partial derivatives assuming all number densities are constant
C        except ne which is a function of fl and tl.
        dvhf = ((nef)*omh+2.d0*nef*ome)
     &    *onedkt
        dvht = -dvh +((net)
     &    *omh+omht*t*(ne+h+nuy-hd-2.d0*he)+2.d0*omet*t*ne+2.d0*ome*net)
     &    *onedkt
        dvhe1f = ((nef)*omhe1+2.d0*nef*ome)
     &    *onedkt
        dvhe1t = -dvhe1+((net)*omhe1+omhe1t*t*(ne+nux-h+hd)
     &    +2.d0*omet*t*ne+2.d0*ome*net)*onedkt
        dvhe2f = (omhe1*(-nef)+
     &    omhe2*(2.d0*(nef))
     &    +2.d0*ome*nef)*onedkt
        dvhe2t = -dvhe2+((hd-ne-nux+h)*omhe1t*t
     &    +(-net)*omhe1+(nuy+2.d0*(ne+nux-h-hd))*omhe2t*t
     &    +(2.d0*(net))*omhe2+2.d0*omet*t*ne+2.d0*ome*net)
     &    *onedkt
        dvef = (
     &    +2.d0*ome*nef)*onedkt
        dvet = -dve+(omht*t*h+omhe1t*t*hd+2.d0*omhe2t*t*he+
     &    2.d0*omet*t*ne
     &    +2.d0*ome*net)*onedkt
      endif
      return
      entry fjs_pi_free(nux, nuy, nuz, h, hd, he,
     &  t, ne, nef, net,
     &  free_pi, free_pif, free_pit, free_pi2_daux)
C      Calculate pressure and free energy components due to fjs
C      form of pressure ionization.
C      n.b. ppi and all of its derivatives are exactly equivalent to
C      free_pi and all of its derivatives so don't even bother returning
C      ppi and derivatives.
C      n.b. the following input parameters are all in "n" form
C      rather than "nu" form.
C      h is the first auxiliary variable,
C      hd is the second auxiliary variable
C      he is the third auxiliary variable
C      ne0, nux, nuy, and nuz are proportional to rho, where ln rho
C      is the 4th auxiliary variable.
C      ne is actually the number density form (n_e outside) and the only
C      input parameter that provides an fl dependence for fixed auxiliary
C      variables since n_e is a function of fl, tl.
      ne0 = nux+nuy+nuz
      ppi1 = nux*ne0-h*(ne+nuy-2.d0*he-hd)
      ppi1f = -h*nef
      ppi1t = -h*net
      ppi11 = -(ne+nuy-2.d0*he-hd)
      ppi12 = h
      ppi13 = 2.d0*h
      ppi14 = 2.d0*nux*ne0 - h*nuy
      ppi2 = hd*(h-ne-nux)
      ppi2f = -hd*nef
      ppi2t = -hd*net
      ppi21 = hd
      ppi22 = (h-ne-nux)
      ppi23 = 0.d0
      ppi24 = -hd*nux
      ppi3 = nuy*ne0-he*(nuy+2.d0*(ne+nux-h-he-hd))
      ppi3f = -2.d0*he*nef
      ppi3t = -2.d0*he*net
      ppi31 = 2.d0*he
      ppi32 = 2.d0*he
      ppi33 = -(nuy+2.d0*(ne+nux-h-he-hd)) + 2.d0*he
      ppi34 = 2.d0*nuy*ne0 - he*(nuy+2.d0*nux)
      ppi4 = (ne0-ne)*(ne0+ne)
      ppi4f = -2.d0*ne*nef
      ppi4t = -2.d0*ne*net
      ppi41 = 0.d0
      ppi42 = 0.d0
      ppi43 = 0.d0
      ppi44 = 2.d0*ne0*ne0
C      n.b. from functional form we have free energy/unit volume = ppi
      free_pi = (omh*ppi1+omhe1*ppi2+omhe2*ppi3+ome*ppi4)
C      partial wrt fl with tl and auxiliary variables fixed.
      free_pif = (omh*ppi1f+omhe1*ppi2f+omhe2*ppi3f+ome*ppi4f)
C      partial wrt tl with fl and auxiliary variables fixed.
      free_pit = (omh*ppi1t+omhe1*ppi2t+omhe2*ppi3t+ome*ppi4t) +
     &  (omht*ppi1+omhe1t*ppi2+omhe2t*ppi3+omet*ppi4)*t
C      partial wrt 4 auxiliary variables with fl, tl fixed.
      free_pi2_daux(1) = (omh*ppi11+omhe1*ppi21+omhe2*ppi31+ome*ppi41)
      free_pi2_daux(2) = (omh*ppi12+omhe1*ppi22+omhe2*ppi32+ome*ppi42)
      free_pi2_daux(3) = (omh*ppi13+omhe1*ppi23+omhe2*ppi33+ome*ppi43)
      free_pi2_daux(4) = (omh*ppi14+omhe1*ppi24+omhe2*ppi34+ome*ppi44)
      return
      entry fjs_pi_end(t, rho, rf, rt,
     &  nux, nuxf, nuxt, nuy, nuyf, nuyt, nuz, nuzf, nuzt,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, ne, nef, net,
     &  ppi, ppif, ppit, spi, spif, spit, upi)
C      compute remaining quantities having found ionization balance
C      n.b. all input quantities are in *nu* (= n/(rho*avogadro) form
C      for the call to fjs_pi_end unlike the input quantities to fjs_pi
C      (which are in *n* form).
      ne0 = nux+nuy+nuz
      ne0f = nuxf+nuyf+nuzf
      ne0t = nuxt+nuyt+nuzt
      ppi1 = nux*ne0-h*(ne+nuy-2.d0*he-hd)
      ppi1f = nuxf*ne0+nux*ne0f-hf*(ne+nuy-2.d0*he-hd)
     &  -h*(nef+nuyf-2.d0*hef-hdf)
      ppi1t = nuxt*ne0+nux*ne0t-ht*(ne+nuy-2.d0*he-hd)
     &  -h*(net+nuyt-2.d0*het-hdt)
      ppi2 = hd*(h-ne-nux)
      ppi2f = hdf*(h-ne-nux)+hd*(hf-nef-nuxf)
      ppi2t = hdt*(h-ne-nux)+hd*(ht-net-nuxt)
      ppi3 = nuy*ne0-he*(nuy+2.d0*(ne+nux-h-he-hd))
      ppi3f = nuyf*ne0+nuy*ne0f
     &  -hef*(nuy+2.d0*(ne+nux-h-he-hd))
     &  -he*(nuyf+2.d0*(nef+nuxf-hf-hef-hdf))
      ppi3t = nuyt*ne0+nuy*ne0t
     &  -het*(nuy+2.d0*(ne+nux-h-he-hd))
     &  -he*(nuyt+2.d0*(net+nuxt-ht-het-hdt))
      ppi4 = (ne0-ne)*(ne0+ne)
      ppi4f = 2.d0*(ne0*ne0f -ne*nef)
      ppi4t = 2.d0*(ne0*ne0t -ne*net)
C      n.b. from functional form we have free energy/unit volume = ppi
C      n.b. multiply by (rho*avogadro)^2 to convert to number density ^2.
      ppi = (omh*ppi1+omhe1*ppi2+omhe2*ppi3+ome*ppi4)*
     &  (rho*avogadro)*(rho*avogadro)
      ppif = (omh*ppi1f+omhe1*ppi2f+omhe2*ppi3f+ome*ppi4f)*
     &  (rho*avogadro)*(rho*avogadro) + 2.d0*ppi*rf
      ppit = ((omh*ppi1t+omhe1*ppi2t+omhe2*ppi3t+ome*ppi4t) +
     &  (omht*ppi1+omhe1t*ppi2+omhe2t*ppi3+omet*ppi4)*t)*
     &  (rho*avogadro)*(rho*avogadro) + 2.d0*ppi*rt
C      from functional form fpi = free energy/unit volume = ppi
C      spi = entropy/unit mass = -partial ppi(t, V, n)/partial t/rho
      spi = -(omht*ppi1+omhe1t*ppi2+omhe2t*ppi3+omet*ppi4)*
     &  (avogadro)*(rho*avogadro)
      spif = -(omht*ppi1f+omhe1t*ppi2f+omhe2t*ppi3f+omet*ppi4f)*
     &  (avogadro)*(rho*avogadro) + spi*rf
      spit = -(omht*ppi1t+omhe1t*ppi2t+omhe2t*ppi3t+omet*ppi4t)*
     &  (avogadro)*(rho*avogadro) + spi*rt
C      internal energy/unit mass = fpi*V/m + t*spi = ppi/rho + t*spi
      upi = (chpi1*ppi1+chpi2*ppi2+chpi3*ppi3+chpi4*ppi4)*
     &  (avogadro)*(rho*avogadro)
      end
