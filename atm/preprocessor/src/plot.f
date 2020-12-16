
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
      module plot_support
      use const_def, only: dp

      implicit none

      contains
      
      subroutine make_plot_files

         character (len=256) :: filename
         integer, parameter :: nZ=8, ng=13, nT=82
         integer :: table_atm_version = 1
         character(len=40) :: table_atm_files(nZ)
         character(len=8):: atm_mix(nZ) !not used currently -- nor is alphaFe
      
         real(dp) :: Teff_array(nT), logg_array(ng), Teff_bound(ng)
         real(dp) :: logZ(nZ), alphaFe(nZ), Pgas(ng,nT,nZ)
         real(dp) :: Teff_tmp(nT), logg_tmp(ng)
         integer :: ierr, iounit, i, j, ii, jj, ibound_tmp(ng), nZ_tmp, nT_tmp, ng_tmp, &
            tmp_version(nZ), ibound(ng,nZ), iZ, text_file_version

         real(dp), parameter :: Pfill = 1d2 !used for interpolation in filling missing values
         real(dp) :: d0, d1, Pinterp_max, Pinterp
         integer :: i_max, j_max
         
         logical, parameter :: dbg = .false.

         include 'formats.dek'
         
         iounit = 33

         filename = 'atm_data/table_summary.txt'
         open(iounit,file=trim(filename),action='read',status='old',iostat=ierr)
         if(ierr/=0) then
            write(*,*) 'table_atm_init: missing atm data'
            write(*,*) filename
            error stop 1
         endif

         !read first line and (nZ, nT, ng)
         read(iounit,*)            !first line is text, skip it
         read(iounit,*) nZ_tmp, nT_tmp, ng_tmp
         if(nZ_tmp /= nZ .or. nT_tmp /= nT .or. ng_tmp /= ng) then
            write(*,*) 'table_atm_init: problem with table dimensions'
            error stop 2
         endif

         !read filenames and headers
         read(iounit,*)            !text
         do i=1,nZ
            read(iounit,'(a)') table_atm_files(i)
            read(iounit,'(14x,i4)') tmp_version(i)
            read(iounit,'(13x,f5.2,8x,f4.1,1x,a8,1x,15x,13i4)') &
               logZ(i), alphaFe(i), atm_mix(i), ibound(:,i)
         enddo

         !read Teff_array
         read(iounit,*)            !text
         read(iounit,'(13f7.0)') Teff_array(:)

         !read logg_array
         read(iounit,*)            !text
         read(iounit,'(13f7.2)') logg_array(:)

         close(iounit)

         !determine table boundaries
         do i=1,ng                 ! -- for each logg, smallest Teff at which Pgas->0
            Teff_bound(i) = Teff_array(ibound(i,1))
            do j=2,nZ
               Teff_bound(i) = min( Teff_bound(i) , Teff_array(ibound(i,j)) )
            enddo
         enddo
         
         do iZ=1,nZ
            filename = 'atm_data/' // trim(table_atm_files(iZ))
            open(iounit,file=trim(filename),action='read',status='old',iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'load_atm_table: missing atm data:'
               write(*,*) filename
               error stop 1
            end if
            read(iounit,'(14x,i4)') text_file_version
            if(text_file_version /= table_atm_version) then
               write(*,*) 'load_atm_table: mismatch in table versions'
               write(*,*)
               error stop 1
            endif
            
            read(iounit,'(13x,f5.2,8x,f4.1,1x,a8,1x,15x,13i4)') &
               logZ(iZ), alphaFe(iZ), atm_mix(iZ), ibound_tmp(:)
            read(iounit,'(15x,13(9x,f5.2,1x))') logg_tmp(:)
            do j=1,nT
               read(iounit,'(14e15.7)') Teff_tmp(j), Pgas(:,j,iZ)
            enddo
            close(iounit)
         end do
         
         do iZ=1, nZ
            filename = 'plot_data/' // trim(table_atm_files(iZ)) // '.data'
            write(*,*) 'write ' // trim(filename)
            open(iounit,file=trim(filename),action='write',iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open:'
               write(*,*) filename
               error stop 1
            end if
            write(iounit,'(e20.10)') log10(max(1d-99,Pgas(:,:,iZ)))
            close(iounit)
         end do
         
         filename = 'plot_data/logg.data'
         write(*,*) 'write ' // trim(filename)
         open(iounit,file=trim(filename),action='write',iostat=ierr)
         write(iounit,'(e20.10)') logg_array(:)
         close(iounit)
         
         filename = 'plot_data/Teff.data'
         write(*,*) 'write ' // trim(filename)
         open(iounit,file=trim(filename),action='write',iostat=ierr)
         write(iounit,'(e20.10)') Teff_array(:)
         close(iounit)
                  
      end subroutine make_plot_files
         
      
      end module plot_support


      program plot
      use plot_support
      call make_plot_files
      end program plot
