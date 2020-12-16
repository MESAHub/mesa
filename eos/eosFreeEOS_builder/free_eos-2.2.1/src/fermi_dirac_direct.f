C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_direct.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine fermi_dirac_direct(
     &  akindex_local, eta_local, beta_local, ifderivative, answer)
C      Calculate Fermi-Dirac integral (C.G. eq. 24,97) (answer(1)) and
C      its partial derivative with respect to e=eta (answer(2)),
C      b=beta (answer(3)), ee (answer(4)), eb (answer(5)), bb (answer(6)),
C      eee (answer(7)), eeb (answer(8)), and ebb (answer(9)).
C      N.B.
C      only answer(1) through answer(min(1, max(9,ifderivative+1)))
C      is returned.
      implicit none
      include 'quad4block.h'
      double precision akindex_local, eta_local, beta_local
      double precision fermi_fun, dfermi_fun_de, dfermi_fun_db,
     &  dfermi_fun_dee, dfermi_fun_deb, dfermi_fun_dbb,
     &  dfermi_fun_deee, dfermi_fun_deeb, dfermi_fun_debb,
     &  fermi_fun0, dfermi_fun0_db, dfermi_fun0_dbb
      external fermi_fun, dfermi_fun_de, dfermi_fun_db,
     &  dfermi_fun_dee, dfermi_fun_deb, dfermi_fun_dbb,
     &  dfermi_fun_deee, dfermi_fun_deeb, dfermi_fun_debb,
     &  fermi_fun0, dfermi_fun0_db, dfermi_fun0_dbb
      double precision q4err, amp, ampmid
      double precision answer(9), anslow, ansmid
      integer ifderivative
      data q4err /1.d-09/           ! accuracy of integration
C      copy argument list to common block
      akindex = akindex_local
      eta = eta_local
      beta = beta_local

C      range of integration.
C      1/exp(40) ~ 4.d-18
      amp= 40.d0+abs(eta)
      if(eta.le.0.d0.or.amp-80.d0.le.0.d0) then
C        divide integral into two ranges:
C        first range where exp factor is variable,
C        but not overwhelmingly so, and upper range where exponential
C        cutoff factor is most active.
        if(eta.le.0.d0) then
C          Most of the integral determined below this value.
          ampmid = 1.d0
        else
C        N.B. eta > 0 ==> amp = eta + 40 <= 80 ==> eta <= 40 ==> amp - 80 < eta
          ampmid = 1.d0 + eta
        endif
        call quad4b(ansmid,0.d0,ampmid,fermi_fun,q4err)
        call quad4b(answer(1),ampmid,amp,fermi_fun,q4err)
        answer(1) = ansmid + answer(1)
        if(ifderivative.gt.0) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_de,q4err)
          call quad4b(answer(2),ampmid,amp,dfermi_fun_de,q4err)
          answer(2) = ansmid + answer(2)
        endif
        if(ifderivative.gt.1) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_db,q4err)
          call quad4b(answer(3),ampmid,amp,dfermi_fun_db,q4err)
          answer(3) = ansmid + answer(3)
        endif
        if(ifderivative.gt.2) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_dee,q4err)
          call quad4b(answer(4),ampmid,amp,dfermi_fun_dee,q4err)
          answer(4) = ansmid + answer(4)
        endif
        if(ifderivative.gt.3) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_deb,q4err)
          call quad4b(answer(5),ampmid,amp,dfermi_fun_deb,q4err)
          answer(5) = ansmid + answer(5)
        endif
        if(ifderivative.gt.4) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_dbb,q4err)
          call quad4b(answer(6),ampmid,amp,dfermi_fun_dbb,q4err)
          answer(6) = ansmid + answer(6)
        endif
        if(ifderivative.gt.5) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_deee,q4err)
          call quad4b(answer(7),ampmid,amp,dfermi_fun_deee,q4err)
          answer(7) = ansmid + answer(7)
        endif
        if(ifderivative.gt.6) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_deeb,q4err)
          call quad4b(answer(8),ampmid,amp,dfermi_fun_deeb,q4err)
          answer(8) = ansmid + answer(8)
        endif
        if(ifderivative.gt.7) then
          call quad4b(ansmid,0.d0,ampmid,dfermi_fun_debb,q4err)
          call quad4b(answer(9),ampmid,amp,dfermi_fun_debb,q4err)
          answer(9) = ansmid + answer(9)
        endif
      else
C        divide integral into three ranges: lower range where cutoff factor
C        1/(exp(x-eta)+1) is unity, middle range where exp factor is variable,
C        but not overwhelmingly so,and upper range where exponential
C        cutoff factor is most active.
C        N.B. eta > 0 ==> amp = eta + 40 > 80 ==> eta > 40 ==> amp - 80 < eta
        call quad4b(anslow,0.d0,amp-80.d0,fermi_fun0,q4err)
        call quad4b(ansmid,amp-80.d0,eta,fermi_fun,q4err)
        call quad4b(answer(1),eta,amp,fermi_fun,q4err)
        answer(1) = anslow + ansmid + answer(1)
        if(ifderivative.gt.0) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_de,q4err)
          call quad4b(answer(2),eta,amp,dfermi_fun_de,q4err)
          answer(2) = ansmid + answer(2)
        endif
        if(ifderivative.gt.1) then
          call quad4b(anslow,0.d0,amp-80.d0,dfermi_fun0_db,q4err)
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_db,q4err)
          call quad4b(answer(3),eta,amp,dfermi_fun_db,q4err)
          answer(3) = anslow + ansmid + answer(3)
        endif
        if(ifderivative.gt.2) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_dee,q4err)
          call quad4b(answer(4),eta,amp,dfermi_fun_dee,q4err)
          answer(4) = ansmid + answer(4)
        endif
        if(ifderivative.gt.3) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_deb,q4err)
          call quad4b(answer(5),eta,amp,dfermi_fun_deb,q4err)
          answer(5) = ansmid + answer(5)
        endif
        if(ifderivative.gt.4) then
          call quad4b(anslow,0.d0,amp-80.d0,dfermi_fun0_dbb,q4err)
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_dbb,q4err)
          call quad4b(answer(6),eta,amp,dfermi_fun_dbb,q4err)
          answer(6) = anslow + ansmid + answer(6)
        endif
        if(ifderivative.gt.5) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_deee,q4err)
          call quad4b(answer(7),eta,amp,dfermi_fun_deee,q4err)
          answer(7) = ansmid + answer(7)
        endif
        if(ifderivative.gt.6) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_deeb,q4err)
          call quad4b(answer(8),eta,amp,dfermi_fun_deeb,q4err)
          answer(8) = ansmid + answer(8)
        endif
        if(ifderivative.gt.7) then
C          note: lower range integral is zero so just have to calculate
C          middle and upper range.
          call quad4b(ansmid,amp-80.d0,eta,dfermi_fun_debb,q4err)
          call quad4b(answer(9),eta,amp,dfermi_fun_debb,q4err)
          answer(9) = ansmid + answer(9)
        endif
      endif
      end

      double precision function fermi_fun(x)
C      Integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2
      fermi_fun=x**akindex
      a1=sqrt(1.d0+0.5d0*beta*x)
      a2=exp(x-eta) + 1.d0
      fermi_fun=fermi_fun*a1/a2
      end

      double precision function dfermi_fun_de(x)
C      Derivative wrt eta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_de = x**akindex
      a1=sqrt(1.d0+0.5d0*beta*x)
C      1.d300 is roughly exp(690).  So guard against under and overflow.
      if(x-eta.le.-690.d0) then
C        underflow. a2 ~ 1, a4 gets very small
        dfermi_fun_de = 0.d0
      elseif(x-eta.le.690.d0) then
        a3 = exp(x-eta)
        a2 = a3 + 1.d0
        a4 = a3/a2
        dfermi_fun_de = dfermi_fun_de*a1*a4/a2
      else
C        overflow.  a4 ~ 1, a2 gets very large
        dfermi_fun_de = 0.d0
      endif
      end

      double precision function dfermi_fun_db(x)
C      Derivative wrt beta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2
      dfermi_fun_db = 0.25d0*x**(akindex+1.d0)
      a1=sqrt(1.d0+0.5d0*beta*x)
      a2=exp(x-eta) + 1.d0
      dfermi_fun_db = dfermi_fun_db/(a1*a2)
      end

      double precision function dfermi_fun_dee(x)
C      Derivative wrt eta, eta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_dee = x**akindex
      a1=sqrt(1.d0+0.5d0*beta*x)
C      1.d300 is roughly exp(690).  So guard against under and overflow.
      if(x-eta.le.-690.d0) then
C        underflow. a2 ~ 1, a4 gets very small
        dfermi_fun_dee = 0.d0
      elseif(x-eta.le.690.d0) then
        a3 = exp(x-eta)
        a2 = a3 + 1.d0
        a4 = a3/a2
        dfermi_fun_dee = dfermi_fun_dee*a1*
     &    (2.d0*a4-1.d0)*a4/a2
      else
C        overflow.  a4 ~ 1, a2 gets very large
        dfermi_fun_dee = 0.d0
      endif
      end

      double precision function dfermi_fun_deb(x)
C      Derivative wrt eta, beta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_deb = 0.25d0*x**(akindex+1.d0)
      a1=sqrt(1.d0+0.5d0*beta*x)
C      1.d300 is roughly exp(690).  So guard against under and overflow.
      if(x-eta.le.-690.d0) then
C        underflow. a2 ~ 1, a4 gets very small
        dfermi_fun_deb = 0.d0
      elseif(x-eta.le.690.d0) then
        a3 = exp(x-eta)
        a2 = a3 + 1.d0
        a4 = a3/a2
        dfermi_fun_deb = dfermi_fun_deb*a4/(a1*a2)
      else
C        overflow.  a4 ~ 1, a2 gets very large
        dfermi_fun_deb = 0.d0
      endif
      end

      double precision function dfermi_fun_dbb(x)
C      Derivative wrt beta, beta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2
      dfermi_fun_dbb = -0.0625d0*x**(akindex+2.d0)
      a1=sqrt(1.d0+0.5d0*beta*x)*(1.d0+0.5d0*beta*x)
      a2=exp(x-eta) + 1.d0
      dfermi_fun_dbb = dfermi_fun_dbb/(a1*a2)
      end

      double precision function dfermi_fun_deee(x)
C      Derivative wrt eta, eta, eta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_deee = x**akindex
      a1=sqrt(1.d0+0.5d0*beta*x)
C      1.d300 is roughly exp(690).  So guard against under and overflow.
      if(x-eta.le.-690.d0) then
C        underflow. a2 ~ 1, a4 gets very small
        dfermi_fun_deee = 0.d0
      elseif(x-eta.le.690.d0) then
        a3 = exp(x-eta)
        a2 = a3 + 1.d0
        a4 = a3/a2
C        1 - a4 = (a2 - a3)/a2 = 1/a2
        dfermi_fun_deee = dfermi_fun_deee*a1*
     &    (1.d0 - 6.d0*a4/a2)*a4/a2
      else
C        overflow.  a4 ~ 1, a2 gets very large
        dfermi_fun_deee = 0.d0
      endif
      end

      double precision function dfermi_fun_deeb(x)
C      Derivative wrt eta, eta, beta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_deeb = 0.25d0*x**(akindex+1.d0)
      a1=sqrt(1.d0+0.5d0*beta*x)
C      1.d300 is roughly exp(690).  So guard against under and overflow.
      if(x-eta.le.-690.d0) then
C        underflow. a2 ~ 1, a4 gets very small
        dfermi_fun_deeb = 0.d0
      elseif(x-eta.le.690.d0) then
        a3 = exp(x-eta)
        a2 = a3 + 1.d0
        a4 = a3/a2
        dfermi_fun_deeb = dfermi_fun_deeb*
     &    (2.d0*a4-1.d0)*a4/(a1*a2)
      else
C        overflow.  a4 ~ 1, a2 gets very large
        dfermi_fun_deeb = 0.d0
      endif
      end

      double precision function dfermi_fun_debb(x)
C      Derivative wrt eta, beta, beta of integrand of Fermi-Dirac integral.
      implicit none
      include 'quad4block.h'
      double precision x, a1, a2, a3, a4
      dfermi_fun_debb = -0.0625d0*x**(akindex+2.d0)
      a1=sqrt(1.d0+0.5d0*beta*x)*(1.d0+0.5d0*beta*x)
      a3 = exp(x-eta)
      a2 = a3  + 1.d0
      a4 = a3/a2
      dfermi_fun_debb = dfermi_fun_debb*a4/(a1*a2)
      end

      double precision function fermi_fun0(x)
C      Integrand of Fermi-Dirac integral
C      in x-eta >> 1 limit where denominator = 1.
      implicit none
      include 'quad4block.h'
      double precision x
      fermi_fun0=(x**akindex)*sqrt(1.d0+0.5d0*beta*x)
      end

      double precision function dfermi_fun0_db(x)
C      Derivative wrt beta of integrand of Fermi-Dirac integral
C      in x-eta >> 1 limit where denominator = 1.
      implicit none
      include 'quad4block.h'
      double precision x
      dfermi_fun0_db = (0.25d0*x**(akindex+1.d0))/
     &  sqrt(1.d0+0.5d0*beta*x)
      end
      double precision function dfermi_fun0_dbb(x)
C      Derivative wrt beta, beta of integrand of Fermi-Dirac integral
C      in x-eta >> 1 limit where denominator = 1.
      implicit none
      include 'quad4block.h'
      double precision x
      dfermi_fun0_dbb = (-0.0625d0*x**(akindex+2.d0))/
     &  (sqrt(1.d0+0.5d0*beta*x)*(1.d0+0.5d0*beta*x))
      end
