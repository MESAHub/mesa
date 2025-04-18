! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

      module pgstar_astero_plots
      use star_lib
      use star_def
      use star_pgstar

      implicit none

      contains

      subroutine astero_pgstar_plots_info(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine astero_pgstar_plots_info


      subroutine write_plot_to_file(s, p, file_prefix, number, ierr)
         type (star_info), pointer :: s
         type (pgstar_win_file_data), pointer :: p
         character (len=*), intent(in) :: file_prefix
         integer, intent(in) :: number
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine write_plot_to_file

      end module pgstar_astero_plots
