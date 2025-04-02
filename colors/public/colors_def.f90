! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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



module colors_def
  use const_def, only : strlen, dp
  implicit none
  
  ! Define the controls type
  type :: colors_controls_type
    ! Paths and filenames
    character(len=256) :: instrument
    character(len=256) :: vega_sed
    character(len=256) :: stellar_atm
    
    ! Numeric parameters
    real(dp) :: metallicity
    real(dp) :: distance
    
    ! Boolean controls
    logical :: make_csv
    logical :: use_colors
    
    ! Bookkeeping (similar to kap)
    integer :: handle
    logical :: in_use
  end type colors_controls_type

  ! Handle management (similar to kap)
  integer, parameter :: max_colors_handles = 10
  logical :: colors_is_initialized = .false.
  ! Global array of handles (like kap_handles)
  type (colors_controls_type), target :: colors_handles(max_colors_handles)
  
  ! Global instance accessible everywhere (needed for compatibility)
  type(colors_controls_type), target :: colors_controls
  public :: alloc_colors_handle, free_colors_handle, colors_ptr



contains
  
  integer function alloc_colors_handle(ierr) result(handle)
    integer, intent(out) :: ierr
    integer :: i
    ierr = 0
    handle = -1
    do i = 1, max_colors_handles
       if (.not. colors_handles(i)% in_use) then
          colors_handles(i)% in_use = .true.
          colors_handles(i)% handle = i
          handle = i
          exit
       end if
    end do
    if (handle == -1) then
       ierr = -1
       return
    end if
  end function alloc_colors_handle

  subroutine free_colors_handle(handle)
    integer, intent(in) :: handle
    if (handle >= 1 .and. handle <= max_colors_handles) then
       colors_handles(handle)% in_use = .false.
    end if
  end subroutine free_colors_handle

  subroutine colors_ptr(handle, ctrl, ierr)
    integer, intent(in) :: handle
    type (colors_controls_type), pointer :: ctrl
    integer, intent(out) :: ierr
    
    ierr = 0
    if (handle < 1 .or. handle > max_colors_handles) then
       ierr = -1
       return
    end if
    
    ctrl => colors_handles(handle)
  end subroutine colors_ptr

end module colors_def