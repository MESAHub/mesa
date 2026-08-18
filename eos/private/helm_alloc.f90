! ***********************************************************************
!
!   Copyright (C) 2006, 2007-2019  Frank Timmes & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module helm_alloc

         use const_def, only: dp, use_mesa_temp_cache
      use math_lib

      implicit none

      contains

      subroutine alloc_helm_table(h, imax, jmax)
         ! This routine allocates a Helm_Table and places pointer to it in h.
         ! It also allocates the arrays in the Helm_Table record.

         use eos_def

         type (Helm_Table), allocatable, intent(out) :: h
         integer, intent(in) :: imax, jmax

         allocate(h)

         h% with_coulomb_corrections = .true.
         h% imax = imax
         h% jmax = jmax

         allocate(h% d(imax))
         allocate(h% t(jmax))

         !..for the helmholtz free energy tables
         allocate(h% f(imax, jmax))
         allocate(h% fd(imax, jmax))
         allocate(h% ft(imax, jmax))
         allocate(h% fdd(imax, jmax))
         allocate(h% ftt(imax, jmax))
         allocate(h% fdt(imax, jmax))
         allocate(h% fddt(imax, jmax))
         allocate(h% fdtt(imax, jmax))
         allocate(h% fddtt(imax, jmax))

         !..for the pressure derivative with density tables
         allocate(h% dpdf(imax, jmax))
         allocate(h% dpdfd(imax, jmax))
         allocate(h% dpdft(imax, jmax))
         allocate(h% dpdfdt(imax, jmax))

         !..for chemical potential tables
         allocate(h% ef(imax, jmax))
         allocate(h% efd(imax, jmax))
         allocate(h% eft(imax, jmax))
         allocate(h% efdt(imax, jmax))

         !..for the number density tables
         allocate(h% xf(imax, jmax))
         allocate(h% xfd(imax, jmax))
         allocate(h% xft(imax, jmax))
         allocate(h% xfdt(imax, jmax))

         !..for storing the differences
         allocate(h% dt_sav(jmax - 1))
         allocate(h% dt2_sav(jmax - 1))
         allocate(h% dti_sav(jmax - 1))
         allocate(h% dt2i_sav(jmax - 1))
         allocate(h% dt3i_sav(jmax - 1))

         allocate(h% dd_sav(imax - 1))
         allocate(h% dd2_sav(imax - 1))
         allocate(h% ddi_sav(imax - 1))
         allocate(h% dd2i_sav(imax - 1))
         allocate(h% dd3i_sav(imax - 1))

      end subroutine alloc_helm_table


      subroutine read_helm_table(h, data_dir, ierr)
      use eos_def
      use forum_m, only: hdf5io_t, OPEN_FILE_RO
      use utils_lib, only: mv, switch_str


      type (Helm_Table), intent(inout) :: h
      character(*), intent(IN) :: data_dir
      integer, intent(out) :: ierr

!..this routine reads the helmholtz eos file, and
!..must be called once before the helmeos routine is invoked.

!..declare local variables
      character (len=256) :: filename
      integer          :: i
      type(hdf5io_t) :: hi

       ierr = 0

!..read the normal helmholtz free energy table
       h% logtlo   = 3.0d0
       h% logthi   = 13.0d0
       h% logdlo   = -12.0d0
       h% logdhi   = 15.0d0

       h% templo = exp10(h% logtlo)
       h% temphi = exp10(h% logthi)
       h% denlo = exp10(h% logdlo)
       h% denhi = exp10(h% logdhi)

       h% logtstp  = (h% logthi - h% logtlo)/real(h% jmax-1,kind=dp)
       h% logtstpi = 1.0d0/h% logtstp
       h% logdstp  = (h% logdhi - h% logdlo)/real(h% imax-1,kind=dp)
       h% logdstpi = 1.0d0/h% logdstp

       write(filename,'(2a)') trim(data_dir), '/helm-table.hdf5'

       hi = hdf5io_t(filename, OPEN_FILE_RO)

       call hi% read_dset('f', h% f)
       call hi% read_dset('fd', h% fd)
       call hi% read_dset('ft', h% ft)
       call hi% read_dset('fdd', h% fdd)
       call hi% read_dset('ftt', h% ftt)
       call hi% read_dset('fdt', h% fdt)
       call hi% read_dset('fddt', h% fddt)
       call hi% read_dset('fdtt', h% fdtt)
       call hi% read_dset('fddtt', h% fddtt)

       call hi% read_dset('dpdf', h% dpdf)
       call hi% read_dset('dpdfd', h% dpdfd)
       call hi% read_dset('dpdft', h% dpdft)
       call hi% read_dset('dpdfdt', h% dpdfdt)

       call hi% read_dset('ef', h% ef)
       call hi% read_dset('efd', h% efd)
       call hi% read_dset('eft', h% eft)
       call hi% read_dset('efdt', h% efdt)

       call hi% read_dset('xf', h% xf)
       call hi% read_dset('xfd', h% xfd)
       call hi% read_dset('xft', h% xft)
       call hi% read_dset('xfdt', h% xfdt)

       call hi% final()

       do i=1,h% jmax
          h% t(i) = exp10(h% logtlo + (i-1)*h% logtstp)
       end do

       do i=1,h% imax
          h% d(i) = exp10(h% logdlo + (i-1)*h% logdstp)
       end do

       h% dt_sav(:) = h% t(2:h% jmax) - h% t(1:h% jmax-1)
       h% dt2_sav(:) = h% dt_sav * h% dt_sav
       h% dti_sav(:) = 1d0 / h% dt_sav
       h% dt2i_sav(:) = 1d0 / h% dt2_sav
       h% dt3i_sav(:) = h% dt2i_sav * h% dti_sav

       h% dd_sav(:) = h% d(2:h% imax) - h% d(1:h% imax-1)
       h% dd2_sav(:) = h% dd_sav * h% dd_sav
       h% ddi_sav(:) = 1d0 / h% dd_sav
       h% dd2i_sav(:) = 1d0 / h% dd2_sav
       h% dd3i_sav(:) = h% dd2i_sav * h% ddi_sav

      end subroutine read_helm_table

      end module helm_alloc
