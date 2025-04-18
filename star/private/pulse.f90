! ***********************************************************************
!
!   Copyright (C) 2010-2019  Rich Townsend & The MESA Team
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

module pulse

  use star_def
  use utils_lib

  use pulse_cafein
  use pulse_fgong
  use pulse_osc
  use pulse_gyre
  use pulse_gsm
  use pulse_saio
  use pulse_gr1d

  implicit none

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

    return

  end subroutine export_pulse_data


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
    case ('gsm')
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

    return

  end subroutine get_pulse_data


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
    case ('gsm')
       call write_gsm_data(id, filename, global_data, point_data, ierr)
    case ('saio')
       call write_saio_data(id, filename, global_data, point_data, ierr)
    case ('gr1d')
       call write_saio_data(id, filename, global_data, point_data, ierr)
    case default
       write(*,*) 'unknown format in write_pulse_data: '//TRIM(data_format)
       ierr = -1
    end select

    return

  end subroutine write_pulse_data

end module pulse
