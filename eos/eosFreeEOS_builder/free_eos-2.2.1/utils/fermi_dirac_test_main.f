C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Irwin stellar interior equation of state code
C       Copyright (C) 1996 by Alan W. Irwin
C
C       $Id: fermi_dirac_test_main.f,v 1.3 2005/07/03 20:47:07 irwin Exp $
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
      program fermi_dirac_test_main
C       main programme for standalone testing of fermi_dirac
      implicit none
      integer narg, iarg
      parameter(narg=2)  !special for fermi_dirac
      integer nfunc, ifunc
      parameter(nfunc=12)  !special for fermi_dirac
C      special for fermi_dirac.  transform arg(2) from gl to tcl = ln beta.
      double precision f, dpsidf, g, beta
      integer jarg1, jarg2  !special for fermi_dirac
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
C         special for fermi_dirac, stop with 2 arguments.
        arg(2) = arglo(2) + dble(jarg2-1)*darg(2)
C        special for fermi_dirac.  transform arg(2) from gl to tcl = ln beta
        f = exp(arg(1))
        dpsidf = sqrt(1.d0+f)
        g = exp(arg(2))
        beta = g/dpsidf
        arg(2) = log(beta)
        do iarg = 1,narg
          argsave(iarg) = arg(iarg)
          deltaarg(iarg) = 1.d0  !depends on fermi_dirac arguments
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
      double precision arg(narg), func(1+narg,nfunc),
     &  funcfd(9,4)
      integer morder, iffirst, index, ioffset
      data iffirst/1/
      save
      if(iffirst.eq.1) then
        iffirst = 0
        write(*,'(a)') 'morder?'
        read(*,*) morder
        write(*,'(a,2i10)') 'morder = ', morder
      endif
      if(narg.ne.2.or.nfunc.ne.12)
     &  stop 'deriv call not right for fermi_dirac'
      call fermi_dirac(arg(1),arg(2),
     &  funcfd(1,1), funcfd(1,2), funcfd(1,3), funcfd(1,4), 3, morder)
      ioffset = 0
      do index = 1,4
        func(1,1+ioffset) = log(funcfd(1,index))
        func(2,1+ioffset) = funcfd(2,index)
        func(3,1+ioffset) = funcfd(3,index)
        if(index.lt.3) then
          func(1,2+ioffset) = funcfd(2,index)
          func(2,2+ioffset) = funcfd(4,index)
          func(3,2+ioffset) = funcfd(5,index)
          func(1,3+ioffset) = funcfd(3,index)
          func(2,3+ioffset) = funcfd(5,index)
          func(3,3+ioffset) = funcfd(6,index)
          func(1,4+ioffset) = funcfd(4,index)
          func(2,4+ioffset) = funcfd(7,index)
          func(3,4+ioffset) = funcfd(8,index)
          func(1,5+ioffset) = funcfd(5,index)
          func(2,5+ioffset) = funcfd(8,index)
          func(3,5+ioffset) = funcfd(9,index)
          ioffset = ioffset + 5
        else
          ioffset = ioffset + 1
        endif
      enddo
      end
