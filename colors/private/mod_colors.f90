! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
module mod_colors
  use colors_def
  use math_lib
  use const_def
  use utils_lib
  use colors_ctrls_io

  implicit none
  private

  public :: init_colors, shutdown_colors

  contains

  subroutine init_colors(ierr)
    integer, intent(out) :: ierr
    
    ! Initialize colors by reading controls
    call read_colors_controls(colors_controls, ierr)
    if (ierr /= 0) then
      write(*,*) 'Error initializing colors module'
      return
    end if
    
    ! Output the configuration
    call write_colors_controls_info(colors_controls, 6) ! 6 is stdout
  end subroutine init_colors
  
  subroutine shutdown_colors()
    ! Any cleanup needed
  end subroutine shutdown_colors

end module mod_colors

