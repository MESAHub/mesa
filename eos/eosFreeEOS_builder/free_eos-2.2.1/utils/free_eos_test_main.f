C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Irwin stellar interior equation of state code
C       Copyright (C) 1996 by Alan W. Irwin
C
C       $Id: free_eos_test_main.f 486 2007-06-22 17:40:07Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@uvastro.phys.uvic.ca.
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
C       of the Irwin stellar interior equation of state code
C*******************************************************************************
      program free_eos_test_main
C       main programme for standalone testing of free_eos
      implicit none
      integer narg, iarg
      parameter(narg=2)  !special for free_eos
      integer nfunc, ifunc
      parameter(nfunc=6)  !special for free_eos
      integer jarg1, jarg2  !special for free_eos
      integer marg(narg), nsig, id, kif
      double precision arglo(narg), arghi(narg), darg(narg), 
     &  arg(narg), argsave(narg), deltaarg(narg), deltaarg_first(narg), 
     &  funcsave(1+narg, nfunc), denslim, tlim,
     &  funcp(1+narg, nfunc), funcm(1+narg, nfunc), err
      character*25 version, version_string
      integer lnblnk
!  integer ieee_handler, error_handler, ieeer
!c  Statements that go in "main" to set up error handling.
!c  'all' for all floating point errors (inexact, under, over, divide, invalid.
!c  'common' for most important floating point errors (over, divide, invalid).
!  ieeer = ieee_handler ('clear','all',error_handler)
!  ieeer = ieee_handler ('set','common',error_handler)
!  ieeer = ieee_handler ('set','under',error_handler)
!c  (possibly) turn off ieee special error_handling until required.
!c  ieeer = ieee_handler ('clear','all',error_handler)
!  if(ieeer.ne.0) then
!    write (0,*) 'Failure - error handler not set!'
!    stop
!  endif
      version_string = version()
      write(*,'(a)') 'FreeEOS version = '//
     &  version_string(:lnblnk(version_string))
      write(*,'(a)')
     &  ' marg, arglo, arghi, deltaarg_first, denslim, tlim, kif'
      read(*,*)
     &  marg, arglo, arghi, deltaarg_first,denslim, tlim, kif
      write(*,'(2i10,/,(1p2d25.16))')
     &  marg, arglo, arghi, deltaarg_first, denslim, tlim
      write(*,*) 'kif = ', kif
C       deltaarg_first is natural log, all other (real) input is
C      log base 10 (unless kif = 0) which is now converted to natural log.
      if(kif.ne.0) then
        do iarg = 1, narg
          arglo(iarg) = log(10.d0)*arglo(iarg)
          arghi(iarg) = log(10.d0)*arghi(iarg)
        enddo
      else
        do iarg = 2, narg
          arglo(iarg) = log(10.d0)*arglo(iarg)
          arghi(iarg) = log(10.d0)*arghi(iarg)
        enddo
      endif
      denslim = log(10.d0)*denslim
      tlim = log(10.d0)*tlim
      do iarg = 1,narg
        darg(iarg) = (arghi(iarg)-arglo(iarg))/dble(max(1,marg(iarg)-1))
      enddo
C      special for free_eos, reverse outer and inner do loops
C      so that you end up doing isotherms, not isobars.
      do jarg2 = 1, marg(2)
        arg(2) = arglo(2) + dble(jarg2-1)*darg(2)
        argsave(2) = arg(2)
C       special for free_eos, stop with 2 arguments.
      do jarg1 = 1, marg(1)
        arg(1) = arglo(1) + dble(jarg1-1)*darg(1)
        argsave(1) = arg(1)
C       if kif.eq.2, then first argument is log rho second argument is
C       log T, otherwise it is the user's responsibility to make tlim = 0
C       so this test always passes.
      if(arg(1).le.1.5d0*(arg(2)-log(1.d5)) + denslim.or.
     &  arg(2).ge.tlim) then
        do iarg = 1,narg
          deltaarg(iarg) = 10.d0*deltaarg_first(iarg)
        enddo
        call deriv(argsave, narg, funcsave, nfunc,1)
        write(*,'(a/(1p2d25.16))') ' argsave = ',argsave
        write(*,'(a/(1p2d25.16))') ' funcsave = ', funcsave
        do id = 1,8
          write(*,'(a,i10)') ' id = ', id
          do iarg=1,narg
            deltaarg(iarg) = 0.1d0*deltaarg(iarg)
            arg(iarg) = argsave(iarg) + deltaarg(iarg)
            call deriv(arg, narg, funcp, nfunc,0)
            arg(iarg) = argsave(iarg) - deltaarg(iarg)
            call deriv(arg, narg, funcm, nfunc,0)
            arg(iarg) = argsave(iarg)
            do ifunc = 1, nfunc
C               calculate number of significant digits in calculation.
              if(funcsave(1,ifunc).ne.0.d0) then
                nsig = 17 + nint(log10(max(1.d-18,
     &            abs((funcp(1,ifunc)-funcm(1,ifunc))/
     &            funcsave(1,ifunc)))))
              else
                nsig = +30
              endif
              if(funcsave(1+iarg,ifunc).ne.0.d0) then
                err = (funcp(1,ifunc)-funcm(1,ifunc))/
     &            (funcsave(1+iarg,ifunc)*2.d0*deltaarg(iarg)) - 1.d0
              else
                err = 1.d30
              endif
              write(*,'(2i3, 1pd10.2,i3)') iarg, ifunc, err,nsig
            enddo
          enddo
        enddo
      endif
      enddo
      enddo  !special for free_eos, stop with 2 arguments
      end
      subroutine deriv(arg, narg, func, nfunc, ifpart)
C       interface between outside routine, and routine with derivatives that
C       is actually being tested.
C       arg(narg) = argument variables.
C       func(1+narg, nfunc) = nfunc func values with narg derivatives each.
      implicit none
      integer narg, nfunc, ifpart
      double precision arg(narg), func(1+narg,nfunc),
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e
      include 'constants.h'
      integer neps
      parameter (neps = 20)
      double precision eps(neps), match_variable, tl, fl,
     &  t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,
     &     rmue, fh2, fhe2, fhe3, xmu1, xmu3, psi

      double precision ln10

      
      integer ifoption, ifmodified, ifmtrace, kif,
     &  iteration_count
      integer iffirst
      data iffirst/1/
      save

      ln10 = log(10.0d0)
      
      if(iffirst.eq.1) then
        iffirst = 0
        call abund_process(eps, neps)
        write(*,'(a)')
     &    ' ifoption, ifmodified, ifmtrace, kif?'
        read(*,*)
     &    ifoption, ifmodified, ifmtrace, kif
        write(*,'(4i5)')
     &    ifoption, ifmodified, ifmtrace, kif
      endif
      if(narg.ne.2.or.nfunc.ne.6)
     &  stop 'deriv call not right for free_eos'
      match_variable = arg(1)
      tl = arg(2)


      ifoption = 1
      ifmodified = 101 
      ifmtrace = 0
      
      match_variable = 3.5d0 * ln10
      tl = 6.0d0 * ln10

      eps(1) = 0.11904761904761904d0
      eps(2) = 0.21236196472293009d0
      eps(3) = 4.3287958638259602d-4
      eps(4) = 1.1389620681530981d-4
      eps(5) = 9.0453016988137054d-4
      eps(6) = 1.4669953912483275d-4
      eps(7) = 0.0000000000000000d0
      eps(8) = 2.3752890351779740d-4
      eps(9:20) = 0.0000000000000000d0

        write(*,'(4i5)') ifoption, ifmodified, ifmtrace, kif

      
      call free_eos(ifoption, ifmodified,
     &  ifmtrace, kif, eps, neps, match_variable, tl, fl,
     &  t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,
     &  rmue, fh2, fhe2, fhe3, xmu1, xmu3, psi,
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,
     &  func(1,1), func(1,2), func(1,3),
     &  func(1,4), func(1,5), func(1,6),
     &  iteration_count)
      if(iteration_count.lt.0)
     &  stop 'deriv: bad status for free_eos'
      if(ifpart.eq.1) then
         continue
      endif
        write(*,'(a,/,1p3d25.15)') ' psi,exp(tl),exp(pl) = ',
     &    psi,exp(tl),exp(pl)
        write(*,'(a,/,1p3d25.15)') ' degeneracy = ',
     &    func(1,1), func(2,1), func(3,1)
        write(*,'(a,/,1p3d25.15)') ' pressure = ',
     &    func(1,2), func(2,2), func(3,2)
        write(*,'(a,/,1p3d25.15)') ' density = ',
     &    func(1,3), func(2,3), func(3,3)
        write(*,'(a,/,1p3d25.15)') ' energy = ',
     &    func(1,4), func(2,4), func(3,4)
        write(*,'(a,/,1p3d25.15)') ' enthalpy = ',
     &    func(1,5), func(2,5), func(3,5)
        write(*,'(a,/,1p3d25.15)') ' entropy = ',
     &    func(1,6), func(2,6), func(3,6)
        write(*,'(a,/,1p3d25.15)') ' gamma1,gamma2,gamma3 = ',
     &    gamma1,gamma2,gamma3
        write(*,'(a,/,1p2d25.15)') ' h2rat,h2plusrat = ',
     &    h2rat,h2plusrat
        write(*,'(a,/,1p3d25.15)') ' fh2,fhe2,fhe3 = ',
     &    fh2,fhe2,fhe3
        write(*,'(a,/,1p2d25.15)') ' lambda,gamma = ',
     &    lambda,lambda**(2.d0/3.d0)/3.d0**(1.d0/3.d0)
        write(*,'(a, i10)') ' iteration count = ', iteration_count

        write(*,*) ' p, e, s ', p, func(1,4), func(1,6)
        
!endif

        stop
      end
