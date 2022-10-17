program create_table
   use chem_support
   use iso_fortran_env, only : error_unit
   character(len = 64) :: inlist_fname
   integer :: ios
   
   ios = 0
   call get_command_argument(1, inlist_fname, status = ios)
   if (ios /= 0) then
      write (error_unit, *) 'unable to get inlist filename from command line'
      stop
   end if
   
   call read_input_parameters(inlist_fname)
   call init_preprocessor
   call read_mass_table
   call process_winvn_table
   call cleanup
end program create_table
