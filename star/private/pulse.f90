! ***********************************************************************
!
!   Copyright (C) 2010-2019 Bill Paxton, Rich Townsend
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

module pulse

  ! Uses

  use star_def
  use utils_lib

  use pulse_cafein
  use pulse_fgong
  use pulse_osc
  use pulse_gyre
  use pulse_saio
  use pulse_gr1d

  ! No implicit typing
  
  implicit none

  ! Access specifiers

  private

  public :: export_pulse_data
  public :: get_pulse_data
  public :: write_pulse_data

contains

  subroutine export_pulse_data (id, data_format, filename, &
       add_center_point, keep_surface_point, add_atmosphere, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: data_format
    character(*), intent(in) :: filename
    logical, intent(in)      :: add_center_point
    logical, intent(in)      :: keep_surface_point
    logical, intent(in)      :: add_atmosphere
    integer, intent(out)     :: ierr

    type(star_info), pointer :: s
    real(dp), allocatable    :: global_data(:)
    real(dp), allocatable    :: point_data(:,:)

    ! Export pulsation data to a file with the specified format

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for export_pulse_data'
       return
    end if

    ! If necessary, hand off to the user hook

    if (s%use_other_export_pulse_data .AND. ASSOCIATED(s%other_export_pulse_data)) then
       call s%other_export_pulse_data(id, data_format, filename, &
            add_center_point, keep_surface_point, add_atmosphere, ierr)
       return
    end if

    ! Get the pulsation data

    call get_pulse_data(id, data_format, &
         add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    if (ierr /= 0) return

    ! Write the pulsation data

    call write_pulse_data(id, data_format, filename, global_data, point_data, ierr)
    if (ierr /= 0) return

    ! Finish

    return

  end subroutine export_pulse_data

  !****

  subroutine get_pulse_data (id, data_format, &
       add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

    integer, intent(in)                :: id
    character(*), intent(in)           :: data_format
    logical, intent(in)                :: add_center_point
    logical, intent(in)                :: keep_surface_point
    logical, intent(in)                :: add_atmosphere
    real(dp), allocatable, intent(out) :: global_data(:)
    real(dp), allocatable, intent(out) :: point_data(:,:)
    integer, intent(out)               :: ierr

    type(star_info), pointer :: s

    ! Get pulsation data

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for export_pulse_data'
       return
    end if

    ! If necessary, hand off to the user hook

    if (s%use_other_get_pulse_data .AND. ASSOCIATED(s%other_get_pulse_data)) then
       call s%other_get_pulse_data(id, data_format, &
            add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
       return
    end if

    select case (StrLowCase(data_format))
    case ('cafein')
       call get_cafein_data(id, add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    case ('fgong')
       call get_fgong_data(id, add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    case ('osc')
       call get_osc_data(id, add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    case ('gyre')
       call get_gyre_data(id, add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    case ('saio')
       call get_saio_data(id, keep_surface_point, add_atmosphere, global_data, point_data, ierr)
    case ('gr1d')
       call get_gr1d_data(id, global_data, point_data, ierr)
    case default
       write(*,*) 'unknown format in get_pulse_data: '//TRIM(data_format)
       ierr = -1
    end select

    ! Edit the data

    if (s%use_other_edit_pulse_data .AND. ASSOCIATED(s%other_edit_pulse_data)) then
       call s%other_edit_pulse_data(s%id, data_format, global_data, point_data, ierr)
    end if

    ! Finish

    return

  end subroutine get_pulse_data

  !****

  subroutine write_pulse_data (id, data_format, filename, global_data, point_data, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: data_format
    character(*), intent(in) :: filename
    real(dp), intent(in)     :: global_data(:)
    real(dp), intent(in)     :: point_data(:,:)
    integer, intent(out)     :: ierr
    
    ! Write pulsation data

    select case (StrLowCase(data_format))
    case ('cafein')
       call write_cafein_data(id, filename, global_data, point_data, ierr)
    case ('fgong')
       call write_fgong_data(id, filename, global_data, point_data, ierr)
    case ('osc')
       call write_osc_data(id, filename, global_data, point_data, ierr)
    case ('gyre')
       call write_gyre_data(id, filename, global_data, point_data, ierr)
    case ('saio')
       call write_saio_data(id, filename, global_data, point_data, ierr)
    case ('gr1d')
       call write_saio_data(id, filename, global_data, point_data, ierr)
    case default
       write(*,*) 'unknown format in write_pulse_data: '//TRIM(data_format)
       ierr = -1
    end select

    ! Finish

    return

  end subroutine write_pulse_data

end module pulse
