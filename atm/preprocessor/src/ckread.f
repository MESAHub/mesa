! ***********************************************************************
!
!   Copyright (C) 2010  Aaron Dotter
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
      program ckread
      use const_def, only: dp
      implicit none
! this program reads in a Castelli & Kurucz 2003 atmosphere structure file, input as 
! the first command line argument, and writes out a formatted table of Pgas @ T=Teff
! for the full range of Teff and logg covered by the model grid.
! the program should be executed with both input and output filenames specified on 
! the command line, followed by tau_base ('phot' or '100').  for example: 
!
!              ./ckread ap00k2odfnew.dat ap00k2odfnew.tbl
!                      (structures file)  (output table) phot
!
! the a*k2odfnew.dat atm structure files are provided for convenience in
! atm_input_data/ck03/ck03_structures.tar.xz 
!
! the output files can be used as input to the table_merge program
      integer :: i, j, k, j1, it, npi
      integer, parameter :: nint=4, nt=76, np=73, maxg=11
      integer :: ng(nt)
      real(dp) :: teff(nt),logg(nt,maxg),pres(nt,maxg),temp(nt,maxg),t(np),p(np),ft(nint),qt(nint)
      character(len=200) :: input_file, output_file, ctau_base, fname
      logical :: do_phot
   
      ng(1:11) = 11
      ng(12:17)= 10
      ng(18:20)=  9
      ng(21:23)=  8
      ng(24:34)=  7
      ng(35:39)=  6
      ng(40)   =  7
      ng(41:45)=  6
      ng(46:52)=  5
      ng(53:57)=  4
      ng(58:65)=  3
      ng(66:75)=  2
      ng(76)   =  1
   
      if ( COMMAND_ARGUMENT_COUNT() /= 2) stop 'usage: ./ckread [input file] [phot or 100]'

      call GET_COMMAND_ARGUMENT(1,input_file)
      call GET_COMMAND_ARGUMENT(2,ctau_base)
      if (trim(adjustl(ctau_base)) == 'phot') then
         do_phot = .true.
      else if (trim(adjustl(ctau_base)) == '100') then
         do_phot = .false.
      else
         write(*,*) 'third arg is currently limited to phot or 100'
         error stop 1
      end if
      output_file = trim(adjustl(ctau_base)) // '_' // trim(input_file) // '.tab'
      fname = 'atm_input_data/ck03/ck03_structures/' // trim(input_file)
      write(*,*) 'read ' // trim(fname)
      open(1,file=trim(fname))
      do i=1,nt
         j1=maxg-ng(i)+1
         if(j1>0)then
            do j=1,j1-1
               pres(i,j)=0.0d0
               temp(i,j)=0.0d0
            enddo
         endif
         do j=j1,maxg
            read(1,'(4x,f8.0,9x,f8.5)') teff(i),logg(i,j)
            do k=1,21
               read(1,*)
            enddo
            read(1,'(10x,i3)') npi
            do k=1,npi
               read(1,'(15x,f9.1,e10.3)') t(k),p(k)
            enddo
            if (do_phot) then
               call locate(t,npi,teff(i),it)
               do k=1,nint
                  ft(k)=t(it+k-nint/2)
               enddo
               call interp(ft,qt,teff(i),nint)
               pres(i,j) = 0.0d0
               do k=1,nint
                  pres(i,j) = pres(i,j) + qt(k)*p(it+k-nint/2)
               enddo
            else ! the last row of data is for tau=100
               temp(i,j) = t(npi)
               pres(i,j) = p(npi)
            end if
            read(1,*)
            read(1,*)
         enddo
      enddo

      open(2,file=trim(output_file))
      write(2,'("#Teff(K)| Pgas@",13("  log g =",f5.2," "))') (-0.5+0.5*(i-1),i=1,maxg+2)
      do i=1,nt
         write(2,'(1p,20e15.7)') teff(i),0d0,(pres(i,j),j=1,maxg),0d0
      enddo
      if (.not. do_phot) then
         write(2,'("#Teff(K)|    T@",13("  log g =",f5.2," "))') (-0.5+0.5*(i-1),i=1,maxg+2)
         do i=1,nt
            write(2,'(1p,20e15.7)') teff(i),0d0,(temp(i,j),j=1,maxg),0d0
         enddo
      end if

      close(1)
      close(2)
      
      contains

      subroutine interp(a,b,x,n)
      ! {a} are the tabulated values for use in interpolation, must be sorted
      ! {b} are coefficients of the interpolating polynomial returned by interp
      !  x  is the abscissa to be interpolated
      !  n  is the number of points to be used
      !  the resulting interpolating polynomial has order n-1 
      implicit none
      integer :: i,j,n
      real(dp) :: a(n),b(n),x
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      return
      end subroutine interp

      subroutine locate(xx,n,x,j)
      implicit none
      integer :: n, j, jl, ju, jm
      real(dp) :: xx(n), x
      jl=0
      ju=n+1
      do while(ju-jl > 1)
         jm=(ju+jl)/2
         if((xx(n)>xx(1)).eqv.(x>xx(jm)))then
            jl=jm
         else
            ju=jm
         endif
      enddo
      j=jl
      end subroutine locate

      end program ckread
