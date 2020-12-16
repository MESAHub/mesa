
! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
      module nextgen_support

      implicit none

      integer, parameter :: table100_version = 2

      character(len=64) :: f_output, head_phx

      contains
      
      
      subroutine make_nextgen

         character (len=256) :: filename
         
         integer, parameter :: ng = 5, nT = 46

         integer :: i, iT, j, io, ierr
         integer :: T2s(nT), ibound(ng)
         real(dp) :: loggs(ng)
         integer, parameter :: n = 50
         integer :: layers, k
         real(dp) :: vals(n)
         real(dp), dimension(ng,nT) :: &
            T, Pgas, Pe, logT, logPgas, logPe, Psum, logPsum
         
         logical, parameter :: write_plot_files = .false.

         logical, parameter :: dbg = .false.
         
         head_phx = '#[Z/Z_SOLAR]= 0.00 [A/Fe]= 0.0 GN93'  ! for now
         
         i = 0
         do iT = 21, 40
            i = i+1
            T2s(i) = iT
         end do
         do iT = 42, 86, 2
            i = i+1
            T2s(i) = iT
         end do
         do iT = 90, 98, 4
            i = i+1
            T2s(i) = iT
         end do
         if (i /= nT) then
            write(*,*) 'bad i for Ts', i
            call mesa_error(__FILE__,__LINE__)
         end if
         
         loggs(:) = (/ 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0 /)
         ibound(:) = nT
         
         io = 33
         ierr = 0
               
         do i = 1, ng
            do j = 1, nT
               write(filename,'(a,i2,a,f3.1,a)') &
                  'atm_input_data/NextGen_solar/lte', T2s(j), '-', loggs(i), '-0.0.NextGen.20'
               write(*,*) trim(filename)
               open(io,file=trim(filename),action='read',status='old',iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open ' // trim(filename)
                  call mesa_error(__FILE__,__LINE__)
               end if
               call read_nextgen(i,j)
               close(io)
            end do
            write(*,*)
         end do
         
         if (write_plot_files) then
         
            open(io,file='plot_nextgen/logT.data',action='write')
            write(io,'(e20.10)') logT(:,:)
            close(io)
         
            open(io,file='plot_nextgen/logPgas.data',action='write')
            write(io,'(e20.10)') logPgas(:,:)
            close(io)
         
            open(io,file='plot_nextgen/logPe.data',action='write')
            write(io,'(e20.10)') logPe(:,:)
            close(io)
         
            open(io,file='plot_nextgen/logPsum.data',action='write')
            write(io,'(e20.10)') logPsum(:,:)
            close(io)
         
            open(io,file='plot_nextgen/Teff.data',action='write')
            write(io,'(e20.10)') T2s(:)*1d2
            close(io)
         
            open(io,file='plot_nextgen/logg.data',action='write')
            write(io,'(e20.10)') loggs(:)
            close(io)
            
         end if
         
         
         !write out combined table
         f_output = 'atm_data/100_Zp00.tbl'
         open(io,file=f_output)
         write(io,'("Table Version",i4)') table100_version
         write(io,'(a40,a15,20i4)') head_phx,'| VALID RANGE: ', ibound
         write(io,'("#Teff(K)| Pgas@")',advance='no')
         call write_logg_header
         do i=1,nT
            write(io,'(1p,20e15.7)') 1d2*T2s(i), Psum(:,i)
         enddo
         write(io,'("#Teff(K)|    T@")',advance='no')
         call write_logg_header
         do i=1,nT
            write(io,'(1p,20e15.7)') 1d2*T2s(i), T(:,i)
         enddo
         close(io)
         
         
         
         contains
         
         subroutine write_logg_header
            do i=1,ng
               write(io,'("  log g =",f5.2," ")',advance='no') loggs(i)
            end do
            write(io,*)
         end subroutine write_logg_header
         
         subroutine read_nextgen(i,j)
            integer, intent(in) :: i, j            
            read(io,*)
            read(io,'(i5)') layers
            if (layers /= n) then
               write(*,*) 'layers /= n', i, j, layers
               call mesa_error(__FILE__,__LINE__)
            end if
            call read_vals(.false.) ! tau
            call read_vals(.false.) ! temperature
            T(i,j) = vals(n)
            logT(i,j) = log10(T(i,j))
            call read_vals(.true.) ! flxrad
            call read_vals(.true.) ! terad
            call read_vals(.true.) ! bkmean
            call read_vals(.true.) ! jkmean
            call read_vals(.true.) ! fkmean
            call read_vals(.true.) ! rkmean
            call read_vals(.true.) ! pgas
            Pgas(i,j) = vals(n)
            logPgas(i,j) = log10(Pgas(i,j))
            call read_vals(.true.) ! pe
            Pe(i,j) = vals(n)
            logPe(i,j) = log10(Pe(i,j))
            Psum(i,j) = Pgas(i,j) + Pe(i,j)
            logPsum(i,j) = log10(Psum(i,j))
         end subroutine read_nextgen
         
         
         subroutine read_vals(skip)
            logical, intent(in) :: skip
            include 'tau100_T2s.dek'
            if (skip) read(io,*) ! skip line
            do k=1,16
               read(io,*) vals(1+(k-1)*3:k*3)
            end do
            read(io,*) vals(49:50)
         end subroutine read_vals
         
                  
      end subroutine make_nextgen
         
      
      end module nextgen_support


      program nextgen
      use nextgen_support
      call make_nextgen
      end program nextgen
