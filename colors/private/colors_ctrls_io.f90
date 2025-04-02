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

module colors_controls_io

   use const_def
   use colors_def
   
   implicit none
   
   ! Define the colors_controls structure
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
   end type colors_controls_type
   
   ! Module variables
   integer, parameter :: colors_namelist_id = 17  ! Choose a unique ID
   
contains

   !************************************
   !******* Read Colors Controls *******
   !************************************
   
   subroutine read_colors_controls(ctrl, ierr)
      type (colors_controls_type), intent(inout) :: ctrl
      integer, intent(out) :: ierr
      
      ! Local variables for namelist
      character(len=256) :: color_instrument
      character(len=256) :: color_vega_sed
      character(len=256) :: color_atm
      real(dp) :: color_z
      real(dp) :: color_d
      logical :: color_make_csv
      logical :: do_colors
      
      namelist /colors_controls/ &
         color_instrument, color_vega_sed, color_make_csv, &
         color_atm, color_z, color_d, do_colors
      
      ierr = 0
      
      ! Set defaults
      call set_default_colors_controls(ctrl)
      
      ! Set local variables from ctrl values
      color_instrument = ctrl% instrument
      color_vega_sed = ctrl% vega_sed
      color_atm = ctrl% stellar_atm
      color_z = ctrl% metallicity
      color_d = ctrl% distance
      color_make_csv = ctrl% make_csv
      do_colors = ctrl% use_colors
      
      ! Read controls from input file
      read(colors_namelist_id, nml=colors_controls, iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'error reading colors_controls namelist'
         return
      end if
      
      ! Update ctrl from local variables
      ctrl% instrument = color_instrument
      ctrl% vega_sed = color_vega_sed
      ctrl% stellar_atm = color_atm
      ctrl% metallicity = color_z
      ctrl% distance = color_d
      ctrl% make_csv = color_make_csv
      ctrl% use_colors = do_colors
      
      ! Additional validation could go here
      
   end subroutine read_colors_controls
   
   !************************************
   !****** Set Default Controls ********
   !************************************
   
   subroutine set_default_colors_controls(ctrl)
      type (colors_controls_type), intent(out) :: ctrl
      
      ! Set default values
      ctrl% instrument = 'data/filters/GAIA/GAIA'
      ctrl% vega_sed = 'data/stellar_models/vega_flam.csv'
      ctrl% stellar_atm = 'data/stellar_models/Kurucz2003all/'
      ctrl% metallicity = 0.0d0
      ctrl% distance = 3.0857d17  ! 10pc for abs mag
      ctrl% make_csv = .false.
      ctrl% use_colors = .true.
      
   end subroutine set_default_colors_controls
   
   !*************************************************
   !********* Write Colors Controls Info ************
   !*************************************************
   
   subroutine write_colors_controls_info(ctrl, io)
      type (colors_controls_type), intent(in) :: ctrl
      integer, intent(in) :: io
      
      write(io,'(A)') '================ Colors Controls ================'
      write(io,'(A,A)') ' instrument: ', trim(ctrl% instrument)
      write(io,'(A,A)') ' vega_sed: ', trim(ctrl% vega_sed)
      write(io,'(A,A)') ' stellar_atm: ', trim(ctrl% stellar_atm)
      write(io,'(A,1PE26.16)') ' metallicity: ', ctrl% metallicity
      write(io,'(A,1PE26.16)') ' distance: ', ctrl% distance
      write(io,'(A,L1)') ' make_csv: ', ctrl% make_csv
      write(io,'(A,L1)') ' use_colors: ', ctrl% use_colors
      write(io,'(A)') '================================================='
      
   end subroutine write_colors_controls_info

end module colors_controls_io