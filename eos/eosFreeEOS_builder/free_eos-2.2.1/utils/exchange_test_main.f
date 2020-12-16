C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Irwin stellar interior equation of state code
C       Copyright (C) 1996 by Alan W. Irwin
C
C       $Id: exchange_test_main.f,v 1.3 2007/07/24 00:04:45 irwin Exp $
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
      program exchange_test_main
C       main programme for standalone testing of exchange_gcpf
      implicit none
      include 'constants.h' !special for exchange_gcpf to transform gl to tl
      integer narg, iarg
      parameter(narg=2)  !special for exchange_gcpf
      integer nfunc, ifunc
      parameter(nfunc=8)  !special for exchange_gcpf
      integer jarg1, jarg2  !special for exchange_gcpf
C      special for exchange_gcpf.  transform arg(2) from gl to tl
      double precision f, dpsidf, g, t
      integer marg(narg), nsig, id
      double precision arglo(narg), arghi(narg), darg(narg), 
     &  arg(narg), argsave(narg), deltaarg(narg),
     &  funcsave(1+narg, nfunc),
     &  funcp(1+narg, nfunc), funcm(1+narg, nfunc), err
!  integer ieee_handler, error_handler, ieeer
!c  Statements that go in "main" to set up error handling.
!c  'all' for all floating point errors (inexact, under, over, divide, invalid.
!c  'common' for most important floating point errors (over, divide, invalid).
!  ieeer = ieee_handler ('clear','all',error_handler)
!  ieeer = ieee_handler ('set','common',error_handler)
!  ieeer = ieee_handler ('set','under',error_handler)
!c  turn off ieee special error_handling until required.
!c  ieeer = ieee_handler ('clear','all',error_handler)
!  if(ieeer.ne.0) then
!    write (0,*) 'Failure - error handler not set!'
!    stop
!  endif
      write(*,'(a)') 'marg, arglo, arghi'
      read(*,*) marg, arglo, arghi
      write(*,'(2i10,/,(1p2d25.16))') marg, arglo, arghi
      do iarg = 1,narg
        darg(iarg) = (arghi(iarg)-arglo(iarg))/
     &    dble(max(1,marg(iarg)-1))
      enddo
      do jarg1 = 1, marg(1)
        arg(1) = arglo(1) + dble(jarg1-1)*darg(1)
      do jarg2 = 1, marg(2)
C         special for exchange_gcpf, stop with 2 arguments.
        arg(2) = arglo(2) + dble(jarg2-1)*darg(2)
C        special for exchange_gcpf.  transform arg(2) from gl to tl
        f = exp(arg(1))
        dpsidf = sqrt(1.d0+f)
        g = exp(arg(2))
        t = g/(dpsidf*ct)
        arg(2) = log(t)
C        end of special for exchange_gcpf.  transform arg(2) from gl to tl
        do iarg = 1,narg
          argsave(iarg) = arg(iarg)
          deltaarg(iarg) = 1.d0  !depends on exchange_gcpf arguments
        enddo
        call deriv(argsave, narg, funcsave, nfunc)
        write(*,'(a/(1p2d25.16))') 'argsave = ',argsave
        write(*,'(a/(1p2d25.16))') 'funcsave = ', funcsave
        do id = 1,8
          write(*,*) 'id = ', id
          do iarg=1,narg
            deltaarg(iarg) = 0.1d0*deltaarg(iarg)
            arg(iarg) = argsave(iarg) + deltaarg(iarg)
            call deriv(arg, narg, funcp, nfunc)
            arg(iarg) = argsave(iarg) - deltaarg(iarg)
            call deriv(arg, narg, funcm, nfunc)
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
              write(*,'(2i3, 1pe10.2,i3)') iarg, ifunc, err,nsig
            enddo
          enddo
        enddo
      enddo
      enddo
      end
      subroutine deriv(arg, narg, func, nfunc)
C       interface between outside routine, and routine with derivatives that
C       is actually being tested.
C       arg(narg) = argument variables.
C       func(1+narg, nfunc) = nfunc func values with narg derivatives each.
      implicit none
      integer narg, nfunc
      double precision arg(narg), func(1+narg,nfunc)
      integer ifexchange, iffirst
      data iffirst/1/
      save
      if(iffirst.eq.1) then
        iffirst = 0
        write(*,'(a)') 'ifexchange?'
        read(*,*) ifexchange
        write(*,'(a,i10)') 'ifexchange = ', ifexchange
      endif
      if(narg.ne.2.or.nfunc.ne.8)
     &  stop 'deriv call not right for exchange_gcpf'
!      call exchange_gcpf(ifexchange, fl, tl,
!     &  psi, dpsidf, dpsidf2,
!     &  fex, fexf, fext,
!     &  fexf2, fexft, fext2,
!     &  fex_psi, fex_psif, fex_psit,
!     &  fex_psif2, fex_psift, fex_psit2)
      call exchange_gcpf(ifexchange, arg(1), arg(2),
     &  func(1,1), func(1,2), func(2,2),
     &  func(1,3), func(2,3), func(3,3),
     &  func(2,4), func(3,4), func(3,5),
     &  func(1,6), func(2,6), func(3,6),
     &  func(2,7), func(3,7), func(3,8))
      func(2,1) = func(1,2)
      func(3,1) = 0.d0
      func(3,2) = 0.d0
      func(1,4) = func(2,3)
      func(1,5) = func(3,3)
      func(2,5) = func(3,4)
      func(1,7) = func(2,6)
      func(1,8) = func(3,6)
      func(2,8) = func(3,7)
      if(abs(log(func(1,6)*func(1,2)/func(2,3))).gt.1.d-14) then
        write(0,*) abs(log(func(1,6)*func(1,2)/func(2,3))),
     &    ' = abs(log(fex_psi*dpsidf/fexf)) > 1.d-14'
        stop 'fex_psi*dpsidf != fexf'
      endif
      end
