! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
      module mod_wd_tau_25
      use const_def
      use utils_lib, only: mesa_error

      implicit none
      
      integer, parameter :: nt = 381, ng = 41
      real(dp), parameter :: min_logg = 5.5d0, max_logg = 9.5d0, dlogg = 0.1d0
      real(dp), parameter :: min_Teff = 2d3, max_Teff = 40d3, dTeff = 100d0
      real(dp), dimension(ng,nt) :: lnT6, lnP15, lnq, lnX
      real(dp) :: Teff(nt), logg(ng), T(ng,nt), Pgas(ng,nt)
      
      contains
      
      
      subroutine build_wd_tau_25_tables
         use const_lib
         character (len=256) :: fname_in, fname_out, line, my_mesa_dir
         integer :: i, j, iounit, ierr, iounit_out, k
         real(dp) :: Prad, Teff_in, logg_in
         
         include 'formats'

         my_mesa_dir = '../..'
         call const_init(my_mesa_dir, ierr)      
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
         
         do i=1,ng
            logg(i) = min_logg + (i-1)*dlogg
         end do
         if (abs(max_logg - logg(ng)) > 1d-6) then
            do i=1,ng
               write(*,2) 'logg(i)', i, logg(i)
            end do
            stop 'bad logg range for wd tau25'
         end if
         
         do i=1,nT
            Teff(i) = min_Teff + (i-1)*dTeff
         end do
         if (abs(max_Teff - Teff(nT)) > 1d-6) then
            stop 'bad Teff range for wd tau25'
         end if
         
         
         ierr = 0
         iounit = 99
         iounit_out = 98
         
         fname_in = 'atm_input_data/wd_atm_tables/boundary_conditions.dat'
         open(iounit,file=trim(fname_in),action='read',status='old',iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname_in)
            call mesa_error(__FILE__,__LINE__)
         end if
         do j= 1, 18
            read(iounit,*,iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in reading header ' // trim(fname_in)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         k = 19
         do j = 1,ng
            do i=nt,1,-1 ! input has Teff in decreasing order
               read(iounit,*,iostat=ierr) &
                  Teff_in, logg_in, lnT6(j,i), lnP15(j,i), lnq(j,i), lnX(j,i)
               if (ierr /= 0 .or. &
                     abs(Teff_in - Teff(i)) > 1d-6 .or. &
                     abs(logg_in - logg(j)) > 1d-6) then
                  write(*,*) 'failed while reading ' // trim(fname_in), i, j
                  write(*,2) 'line number', k, Teff(i), logg(j)
                  call mesa_error(__FILE__,__LINE__)
               end if
               k = k+1
            end do
         end do
         close(iounit)

         
         do i=1,nt
            do j=1,ng
               T(j,i) = exp(lnT6(j,i))*1d6
               Prad = (crad*T(j,i)**4)/3
               Pgas(j,i) = exp(lnP15(j,i))*1d15 - Prad
               if (Pgas(j,i) <= 0) then
                  write(*,*) 'bad value for Pgas', j, i, Pgas(j,i), Prad, T(j,i)
                  call mesa_error(__FILE__,__LINE__)
               end if
            end do
         end do
         
         
         fname_out = 'atm_data/wd_25.tbl'
         open(iounit_out,file=trim(fname_out),action='write',iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname_out)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         fname_in = 'atm_data/wd_tau25_header.txt'
         open(iounit,file=trim(fname_in),action='read',status='old',iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname_in)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         do i=1,2
            read(iounit,'(a)',iostat=ierr) line
            if (ierr /= 0) then
               write(*,*) 'failed while reading ' // trim(fname_in), nt
               call mesa_error(__FILE__,__LINE__)
            end if
            write(iounit_out,'(a)') line
         end do

         write(iounit_out,'(a15)',advance='no') "#Teff(K)| Pgas@" 
         do i = 1, ng
            write(iounit_out,fmt='("  log g =",f5.2," ")',advance='no') logg(i)
         end do
         write(iounit_out,*)
         do j = 1, nT
            write(iounit_out,fmt='(e15.7)',advance='no') Teff(j)
            do i = 1, ng
               write(iounit_out,fmt='(e15.7)',advance='no') Pgas(i,j)
            end do
            write(iounit_out,*)
         end do

         write(iounit_out,'(a15)',advance='no') "#Teff(K)|    T@" 
         do i = 1, ng
            write(iounit_out,fmt='("  log g =",f5.2," ")',advance='no') logg(i)
         end do
         write(iounit_out,*)
         do j = 1, nT
            write(iounit_out,fmt='(e15.7)',advance='no') Teff(j)
            do i = 1, ng
               write(iounit_out,fmt='(e15.7)',advance='no') T(i,j)
            end do
            write(iounit_out,*)
         end do
         
         close(iounit_out)

      end subroutine build_wd_tau_25_tables
   
      
      end module mod_wd_tau_25
      

      program create_wd_tau_25
      use mod_wd_tau_25
      call build_wd_tau_25_tables
      end program create_wd_tau_25
