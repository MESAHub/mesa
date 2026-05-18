      program binary_run
      use binary_lib, only: run1_binary
      use run_star_extras
      use run_binary_extras

      integer :: ierr

      call run1_binary(.true., extras_controls, extras_binary_controls, &
         ierr, 'inlist_binary')
      
      end program
