module run_binary
   implicit none

contains
   
   subroutine do_run_binary(tst)
      use binary_lib, only : run1_binary
      use run_star_extras
      use run_binary_extras
      
      logical, intent(in) :: tst
      
      integer :: ierr
      
      call run1_binary(tst, &
         ! star extras
         extras_controls, &
         ! binary extras
         extras_binary_controls, &
         ierr)
   
   end subroutine do_run_binary

end module run_binary
      
