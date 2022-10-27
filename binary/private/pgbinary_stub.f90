! ***********************************************************************
!
!   Copyright (C) 2010-2022  Bill Paxton, Matthias Fabry & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module pgbinary

   use binary_def
   use const_def
   use chem_def, only : category_name
   use rates_def, only : i_rate

   implicit none


contains

   ! pgbinary interface
   subroutine start_new_run_for_pgbinary(b, ierr) ! reset logs
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine start_new_run_for_pgbinary


   subroutine restart_run_for_pgbinary(b, ierr)
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine restart_run_for_pgbinary


   subroutine read_pgbinary_controls(b, ierr)
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine read_pgbinary_controls


   subroutine read_pgbinary_inlist(b, inlist_fname, ierr)
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character(*), intent(in) :: inlist_fname
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine read_pgbinary_inlist

   subroutine update_pgbinary_plots(b, must_write_files, ierr)
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      logical, intent(in) :: must_write_files
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine update_pgbinary_plots


   subroutine do_create_file_name(b, dir, prefix, name)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: dir, prefix
      character (len = *), intent(out) :: name
      name = ''
   end subroutine do_create_file_name
   


   subroutine do_write_plot_to_file(b, p, filename, ierr)
      use binary_def, only : binary_info, pgbinary_win_file_data
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      character (len = *), intent(in) :: filename
      integer, intent(out) :: ierr
      
   end subroutine do_write_plot_to_file


   subroutine do_show_pgbinary_annotations(&
      b, show_annotation1, show_annotation2, show_annotation3)
      type (binary_info), pointer :: b
      logical, intent(in) :: &
         show_annotation1, show_annotation2, show_annotation3

   end subroutine do_show_pgbinary_annotations


   subroutine do_start_new_run_for_pgbinary(b, ierr) ! reset logs
      use utils_lib
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      integer :: iounit
      character (len = strlen) :: fname
      logical :: fexist
      ierr = 0
   end subroutine do_start_new_run_for_pgbinary


   subroutine do_restart_run_for_pgbinary(b, ierr)
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      logical :: fexists
      ierr = 0
      
   end subroutine do_restart_run_for_pgbinary


   subroutine do_read_pgbinary_controls(b, inlist_fname, ierr)
      type (binary_info), pointer :: b
      character(*), intent(in) :: inlist_fname
      integer, intent(out) :: ierr
      ierr = 0
     
   end subroutine do_read_pgbinary_controls


   subroutine set_win_file_data(b, ierr)
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr

      
   end subroutine set_win_file_data


   subroutine do_pgbinary_plots(b, must_write_files, ierr)
      type (binary_info), pointer :: b
      logical, intent(in) :: must_write_files
      integer, intent(out) :: ierr

      integer :: i
      integer(8) :: time0, time1, clock_rate
      logical :: pause

      include 'formats'

      ierr = 0

   end subroutine do_pgbinary_plots


   ! pgbinary driver, called after each timestep
   subroutine onScreen_Plots(b, must_write_files_in, ierr)
      use utils_lib

      type (binary_info), pointer :: b
      logical :: must_write_files_in
      integer, intent(out) :: ierr

      integer :: i
      logical, parameter :: dbg = .false.
      real(dp) :: dlgL, dlgTeff, dHR
      logical :: must_write_files, show_plot_now, save_plot_now

      include 'formats'
      ierr = 0

   end subroutine onScreen_Plots


   subroutine update_pgbinary_history_file(b, ierr)
      use utils_lib
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr

      integer :: iounit, i, n
      character (len = 1024) :: fname

      logical, parameter :: dbg = .false.

      include 'formats'

      ierr = 0
     
   end subroutine update_pgbinary_history_file


   subroutine read_pgbinary_data(b, ierr)
      use utils_lib
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr

      logical :: fexist
      integer :: iounit, i, n
      character (len = 1024) :: fname

      logical, parameter :: dbg = .false.

      include 'formats'
      ierr = 0

   end subroutine read_pgbinary_data


   subroutine update_pgbinary_data(b, ierr)

      type (binary_info), pointer :: b
      integer, intent(out) :: ierr

      integer :: num, i

      include 'formats'

      ierr = 0

   end subroutine update_pgbinary_data

   subroutine shutdown_pgbinary(b)
      type (binary_info), pointer :: b

   end subroutine shutdown_pgbinary

end module pgbinary
