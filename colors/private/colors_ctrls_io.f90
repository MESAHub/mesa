! ***********************************************************************
!
!   Copyright (C) 2010-2025  The MESA Team
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
module colors_ctrls_io
   use const_def
   use colors_def

   implicit none

   ! Module variables
   integer, parameter :: colors_namelist_id = 17  ! Choose a unique ID

contains

   !************************************
   !******* Read Colors Controls *******
   !************************************

   subroutine read_colors_controls(col, ierr)
      type(colors_controls_type), intent(inout) :: col
      integer, intent(out) :: ierr

      ! Local variables for namelist (MUST be declared at the beginning)
      character(len=256) :: color_instrument
      character(len=256) :: color_vega_sed
      character(len=256) :: color_atm
      real(dp) :: color_z
      real(dp) :: color_d
      logical :: color_make_csv
      logical :: do_colors

      ! Define the namelist
      namelist /colors_controls/ color_instrument, color_vega_sed, color_make_csv, &
                                color_atm, color_z, color_d, do_colors

      ! Set default ierr
      ierr = 0

      ! First try to read from the defaults file
      call read_colors_defaults('colors.defaults', col, ierr)
      if (ierr /= 0) then
        ! If fails, fall back to hardcoded defaults
        call set_default_colors_controls(col)
      end if

      ! Set namelist variables from col values
      color_instrument = col% instrument
      color_vega_sed = col% vega_sed
      color_atm = col% stellar_atm
      color_z = col% metallicity
      color_d = col% distance
      color_make_csv = col% make_csv
      do_colors = col% use_colors

      ! Read controls from input file
      read(colors_namelist_id, nml=colors_controls, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'error reading colors_controls namelist'
        return
      end if

      ! Update col from local variables
      col% instrument = color_instrument
      col% vega_sed = color_vega_sed
      col% stellar_atm = color_atm
      col% metallicity = color_z
      col% distance = color_d
      col% make_csv = color_make_csv
      col% use_colors = do_colors

      ! Additional validation could go here

   end subroutine read_colors_controls

   !************************************
   !****** Set Default Controls ********
   !************************************

   subroutine set_default_colors_controls(col)
      type (colors_controls_type), intent(out) :: col

      ! Set default values
      col% instrument = 'data/filters/GAIA/GAIA'
      col% vega_sed = 'data/stellar_models/vega_flam.csv'
      col% stellar_atm = 'data/stellar_models/Kurucz2003all/'
      col% metallicity = 0.0d0
      col% distance = 3.0857d17  ! 10pc for abs mag
      col% make_csv = .false.
      col% use_colors = .true.

   end subroutine set_default_colors_controls

   !*************************************************
   !********* Write Colors Controls Info ************
   !*************************************************

   subroutine write_colors_controls_info(col, io)
      type (colors_controls_type), intent(in) :: col
      integer, intent(in) :: io

      write(io,'(A)') '================ Colors Controls ================'
      write(io,'(A,A)') ' instrument: ', trim(col% instrument)
      write(io,'(A,A)') ' vega_sed: ', trim(col% vega_sed)
      write(io,'(A,A)') ' stellar_atm: ', trim(col% stellar_atm)
      write(io,'(A,1PE26.16)') ' metallicity: ', col% metallicity
      write(io,'(A,1PE26.16)') ' distance: ', col% distance
      write(io,'(A,L1)') ' make_csv: ', col% make_csv
      write(io,'(A,L1)') ' use_colors: ', col% use_colors
      write(io,'(A)') '================================================='

   end subroutine write_colors_controls_info

   !*************************************************
   !********* Read Colors Defaults ************
   !*************************************************

   subroutine read_colors_defaults(filename, col, ierr)
     character(len=*), intent(in) :: filename
     type(colors_controls_type), intent(out) :: col
     integer, intent(out) :: ierr

     integer :: unit, status
     character(len=256) :: line, key, value

     ! Set defaults first
     call set_default_colors_controls(col)

     ! Open the defaults file
     unit = 20
     open(unit, file=filename, status='old', action='read', iostat=status)
     if (status /= 0) then
       write(*,*) 'Error opening defaults file: ', trim(filename)
       ierr = 1
       return
     end if

     ! Read key-value pairs
     do
       read(unit, '(A)', iostat=status) line
       if (status /= 0) exit

       ! Skip comments and empty lines
       line = adjustl(line)
       if (line == '' .or. line(1:1) == '!' .or. line(1:1) == '#') cycle

       ! Parse key = value
       if (index(line, '=') > 0) then
         key = line(1:index(line, '=')-1)
         value = line(index(line, '=')+1:)
         key = adjustl(trim(key))
         value = adjustl(trim(value))

         ! Set the appropriate control value
         if (key == 'color_instrument') then
           col% instrument = value
         else if (key == 'color_vega_sed') then
           col% vega_sed = value
         else if (key == 'color_atm') then
           col% stellar_atm = value
         else if (key == 'color_z') then
           read(value, *) col% metallicity
         else if (key == 'color_d') then
           read(value, *) col% distance
         else if (key == 'color_make_csv') then
           col% make_csv = (value == 'true' .or. value == 'TRUE')
         end if
       end if
     end do

     close(unit)
     ierr = 0
   end subroutine read_colors_defaults

end module colors_ctrls_io