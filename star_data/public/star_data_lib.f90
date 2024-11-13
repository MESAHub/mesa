! ***********************************************************************
!
!   Copyright (C) 2019  The MESA Team
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
!
! ***********************************************************************

      module star_data_lib

      implicit none

      contains

      subroutine star_data_init(mesa_dir_init,ierr)
         use star_data_def, only: do_star_def_init
         character (len=*), intent(in) :: mesa_dir_init
         integer, intent(out) :: ierr
         ierr = 0
         call do_star_def_init(mesa_dir_init,ierr)
      end subroutine star_data_init

      end module star_data_lib

