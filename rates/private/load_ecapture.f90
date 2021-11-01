! ***********************************************************************
!
!   Copyright (C) 2013  Josiah Schwab
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


module load_ecapture

  use rates_def
  use const_def, only: dp
  use chem_def, only: iso_name_length

  implicit none

contains

  subroutine load_ecapture_data(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    if (do_ecapture) then
       call load_ecapture_states_list(ierr)
       if (ierr /= 0) return
       call load_ecapture_transitions_list(ierr)
    endif
  end subroutine load_ecapture_data


  subroutine load_ecapture_states_list(ierr)
    use utils_lib
    use chem_lib, only: chem_get_iso_id
    use chem_def, only: iso_name_length

    integer, intent(out) :: ierr
    integer :: iounit, i, j, id
    character (len=256) :: filename, string
    integer :: nstates
    character(len=iso_name_length) :: iso

    real(dp) :: Ei, Ji

    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0

    filename = trim(ecapture_states_file)

    if (dbg) then
       write(*,'(A)')
       write(*,*) 'ecapture states filename <' // trim(filename) // '>'
       write(*,'(A)')
    end if

    ierr = 0
    open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open special_weak_states_file ' // trim(filename)
       return
    end if

    allocate(ecapture_nuclide_name(max_num_ecapture_nuclei), &
             ecapture_nuclide_id(max_num_ecapture_nuclei))

    allocate(ecapture_states_data(max_num_ecapture_nuclei * &
         max_num_ecapture_states, num_states_data))

    nullify(ecapture_states_number_dict)
    nullify(ecapture_states_offset_dict)

    num_ecapture_nuclei = 0
    num_ecapture_states = 0
    do i = 1, max_num_ecapture_nuclei ! keep reading until end of file or max # nuclei

       do ! skip any blank or comment lines
          read(iounit,'(a)',iostat=ierr) string
          if (ierr /= 0) exit
          if ((index(string,"!") /= 0) .or. (len_trim(string) == 0)) then
             cycle ! comment or blank line
          else
             exit ! good line
          end if
       end do
       if (ierr /= 0) then
          ierr = 0; exit
       end if

       ! string now holds the first good line
       read(string,fmt=*,iostat=ierr) iso, nstates
       if (ierr /= 0) then
          ierr = 0; exit
       end if
       num_ecapture_nuclei = i

       if (nstates .gt. max_num_ecapture_states) stop "ecapture: too many states"
       call integer_dict_define(ecapture_states_number_dict, iso, nstates, ierr)

       id = chem_get_iso_id(iso)
       if (id <= 0) then
          write(*,*) 'ecapture FATAL ERROR: unknown nuclide ' // iso
          call mesa_error(__FILE__,__LINE__)
       end if

       ecapture_nuclide_id(i) = id
       ecapture_nuclide_name(i) = iso

       if (dbg) write(*,'(a)') 'ecapture nucleus ' // trim(iso)

       ! store where this list of states starts
       call integer_dict_define(ecapture_states_offset_dict, iso, num_ecapture_states, ierr)
       if (failed('integer_dict_define')) return

       do j = 1, nstates
          read(iounit,fmt=*,iostat=ierr) Ei, Ji
          if (ierr /= 0) then
             ierr = 0; exit
          end if
          num_ecapture_states = num_ecapture_states + 1

          ! pack the data into the array
          ecapture_states_data(num_ecapture_states, i_E) = Ei
          ecapture_states_data(num_ecapture_states, i_J) = Ji

       end do

    end do

    close(iounit)

    if (num_ecapture_nuclei == 0) then
       ierr = -1
       write(*,*) 'failed trying to read special_weak_states_file -- no nuclei?'
       return
    end if

    if (num_ecapture_nuclei == max_num_ecapture_nuclei) then
       ierr = -1
       write(*,*) 'failed trying to read special_weak_states_file -- too many nuclei?'
       return
    end if

    call integer_dict_create_hash(ecapture_states_number_dict, ierr)
    if (ierr /= 0) return

    call integer_dict_create_hash(ecapture_states_offset_dict, ierr)
    if (ierr /= 0) return

    call realloc_double2(ecapture_states_data, num_ecapture_states, num_states_data, ierr)
    if (ierr /= 0) return

    call realloc_integer(ecapture_nuclide_id, num_ecapture_nuclei, ierr)
    if (ierr /= 0) return

  contains

    logical function failed(str)
      character (len=*) :: str
      failed = (ierr /= 0)
      if (failed) then
         write(*,*) 'failed: ' // trim(str)
      end if
    end function failed


  end subroutine load_ecapture_states_list


  subroutine load_ecapture_transitions_list(ierr)
    use utils_lib
    use chem_lib, only: chem_get_iso_id
    use chem_def, only: iso_name_length

    integer, intent(out) :: ierr
    integer :: iounit, i, j, id
    character (len=256) :: filename, string
    character(len=iso_name_length) :: lhs, rhs
    integer :: ntrans
    character(len=2*iso_name_length+1) :: key

    integer :: Si, Sf
    real(dp) :: logft

    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0

    filename = ecapture_transitions_file

    if (dbg) then
       write(*,'(A)')
       write(*,*) 'ecapture transitions filename <' // trim(filename) // '>'
       write(*,'(A)')
    end if

    ierr = 0
    open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open special_weak_transitions_file ' // trim(filename)
       write(*,*) trim(string)
       return
    end if

    allocate(ecapture_lhs_nuclide_name(max_num_ecapture_reactions), &
         ecapture_rhs_nuclide_name(max_num_ecapture_reactions), &
         ecapture_lhs_nuclide_id(max_num_ecapture_reactions), &
         ecapture_rhs_nuclide_id(max_num_ecapture_reactions))

    allocate(ecapture_transitions_data(max_num_ecapture_reactions * &
         max_num_ecapture_transitions, &
         num_transitions_data))
    allocate(ecapture_logft_data(max_num_ecapture_reactions * &
         max_num_ecapture_transitions))

    nullify(ecapture_reactions_dict)
    nullify(ecapture_transitions_number_dict)
    nullify(ecapture_transitions_offset_dict)
    num_ecapture_reactions = 0
    num_ecapture_transitions = 0

    do i = 1, max_num_ecapture_reactions ! keep reading until end of file

       do ! skip any blank or comment lines
          read(iounit,'(a)',iostat=ierr) string
          if (ierr /= 0) exit
          if ((index(string,"!") /= 0) .or. (len_trim(string) == 0)) then
             cycle ! comment or blank line
          else
             exit ! good line
          end if
       end do
       if (ierr /= 0) then
          ierr = 0; exit
       end if

       ! string now holds the first good line
       read(string,fmt=*,iostat=ierr) lhs, rhs, ntrans
       if (ierr /= 0) then
          ierr = 0; exit
       end if

       id = chem_get_iso_id(lhs)
       if (id <= 0) then
          write(*,*) 'ecapture FATAL ERROR: unknown nuclide ' // lhs
          call mesa_error(__FILE__,__LINE__)
       end if
       ecapture_lhs_nuclide_id(i) = id
       id = chem_get_iso_id(rhs)
       if (id <= 0) then
          write(*,*) 'ecapture FATAL ERROR: unknown nuclide ' // rhs
          call mesa_error(__FILE__,__LINE__)
       end if
       ecapture_rhs_nuclide_id(i) = id
       ecapture_lhs_nuclide_name(i) = lhs
       ecapture_rhs_nuclide_name(i) = rhs

       call create_ecapture_dict_key(lhs, rhs, key)

       if (dbg) write(*,'(a)') 'ecapture transitions ' // trim(key)

       call integer_dict_define(ecapture_reactions_dict, key, i, ierr)
       if (failed('integer_dict_define')) return
       num_ecapture_reactions = i

       call integer_dict_define(ecapture_transitions_number_dict, key, ntrans, ierr)
       if (failed('integer_dict_define')) return

       call integer_dict_define(ecapture_transitions_offset_dict, key, num_ecapture_transitions, ierr)
       if (failed('integer_dict_define')) return

       do j = 1, ntrans
          read(iounit,fmt=*,iostat=ierr) Si, Sf, logft
          if (ierr /= 0) then
             ierr = 0; exit
          end if
          num_ecapture_transitions = num_ecapture_transitions + 1

          ! pack the data into the array
          ecapture_transitions_data(num_ecapture_transitions, i_Si) = Si
          ecapture_transitions_data(num_ecapture_transitions, i_Sf) = Sf

          ecapture_logft_data(num_ecapture_transitions) = logft

       end do

    end do

    close(iounit)

    if (num_ecapture_reactions == 0) then
       ierr = -1
       write(*,*) 'failed trying to read special_weak_transitions_file -- no reactions?'
       return
    end if

    if (num_ecapture_reactions == max_num_ecapture_reactions) then
       ierr = -1
       write(*,*) 'failed trying to read special_weak_transitions_file -- too many reactions?'
       return
    end if

    call integer_dict_create_hash(ecapture_reactions_dict, ierr)
    if (ierr /= 0) return

    call integer_dict_create_hash(ecapture_transitions_number_dict, ierr)
    if (ierr /= 0) return

    call integer_dict_create_hash(ecapture_transitions_offset_dict, ierr)
    if (ierr /= 0) return

    call realloc_integer2(ecapture_transitions_data, num_ecapture_transitions, num_transitions_data, ierr)
    if (ierr /= 0) return

    call realloc_double(ecapture_logft_data, num_ecapture_transitions, ierr)
    if (ierr /= 0) return

    call realloc_integer(ecapture_lhs_nuclide_id, num_ecapture_reactions, ierr)
    if (ierr /= 0) return

    call realloc_integer(ecapture_rhs_nuclide_id, num_ecapture_reactions, ierr)
    if (ierr /= 0) return

  contains

    logical function failed(str)
      character (len=*) :: str
      failed = (ierr /= 0)
      if (failed) then
         write(*,*) 'failed: ' // trim(str)
      end if
    end function failed


  end subroutine load_ecapture_transitions_list



end module load_ecapture
