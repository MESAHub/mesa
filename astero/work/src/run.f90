program run

   use astero_lib, only: run_star_astero
   use run_star_extras, only: extras_controls

   implicit none

   call run_star_astero(extras_controls)

end program run
