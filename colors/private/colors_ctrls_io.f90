module colors_ctrls_io
   use const_def
   ! Remove the import of star_def
   
   implicit none
   
   ! Module variables
   integer, parameter :: colors_namelist_id = 17  ! Choose a unique ID
   
   ! Define a simple structure for color controls
   type :: colors_controls_type
      character(len=256) :: instrument
      character(len=256) :: vega_sed
      character(len=256) :: atm
      real(dp) :: z
      real(dp) :: d
      logical :: make_csv
      logical :: use_colors
   end type colors_controls_type
   
   ! Global instance - renamed to avoid conflict
   type(colors_controls_type), save :: color_settings
   
contains
   !************************************
   !******* Read Colors Controls *******
   !************************************
   subroutine read_colors_controls(col, ierr)
      type(colors_controls_type), intent(inout) :: col
      integer, intent(out) :: ierr
      
      ! Local variables for namelist
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
      
      ! Set defaults
      call set_default_colors_controls(col)
      
      ! Set namelist variables from col values
      color_instrument = col% instrument
      color_vega_sed = col% vega_sed
      color_atm = col% atm
      color_z = col% z
      color_d = col% d
      color_make_csv = col% make_csv
      do_colors = col% use_colors
      
      ! Read controls from input file
      read(colors_namelist_id, nml=colors_controls, iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'error reading colors_controls namelist'
        return
      end if
      
      ! Update col directly from local variables
      col% instrument = color_instrument
      col% vega_sed = color_vega_sed
      col% atm = color_atm
      col% z = color_z
      col% d = color_d
      col% make_csv = color_make_csv
      col% use_colors = do_colors
   end subroutine read_colors_controls
   
   !************************************
   !****** Set Default Controls ********
   !************************************
   subroutine set_default_colors_controls(col)
      type(colors_controls_type), intent(out) :: col
      
      ! Set default values
      col% instrument = 'data/filters/GAIA/GAIA'
      col% vega_sed = 'data/stellar_models/vega_flam.csv'
      col% atm = 'data/stellar_models/Kurucz2003all/'
      col% z = 0.0d0
      col% d = 3.0857d17  ! 10pc for abs mag
      col% make_csv = .false.
      col% use_colors = .true.
   end subroutine set_default_colors_controls
   
   !*************************************************
   !********* Write Colors Controls Info ************
   !*************************************************
   subroutine write_colors_controls_info(col, io)
      type(colors_controls_type), intent(in) :: col
      integer, intent(in) :: io
      
      write(io,'(A)') '================ Colors Controls ================'
      write(io,'(A,A)') ' instrument: ', trim(col% instrument)
      write(io,'(A,A)') ' vega_sed: ', trim(col% vega_sed)
      write(io,'(A,A)') ' stellar_atm: ', trim(col% atm)
      write(io,'(A,1PE26.16)') ' metallicity: ', col% z
      write(io,'(A,1PE26.16)') ' distance: ', col% d
      write(io,'(A,L1)') ' make_csv: ', col% make_csv
      write(io,'(A)') '================================================='
   end subroutine write_colors_controls_info
end module colors_ctrls_io