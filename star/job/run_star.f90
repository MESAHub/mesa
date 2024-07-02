 
      module run_star
         use caliper_mod
         implicit none
      
      contains
      
      subroutine do_run_star(inlist_fname_arg)
         use run_star_support, only: run1_star
         use run_star_extras, only: extras_controls
         use star_lib, only: starlib_shutdown
         character (len=*) :: inlist_fname_arg
         optional inlist_fname_arg
         logical, parameter :: &
            do_alloc_star = .true., &
            do_free_star = .true., &
            restart_okay = .true.
         integer :: id, ierr
         logical :: restart

         ! A scope annotation. Start region 'main'
         call cali_begin_region('do_run_star')

         call run1_star( &
            do_alloc_star, &
            do_free_star, &
            restart_okay, &
            id, restart, &
            extras_controls, &
            ierr, &
            inlist_fname_arg)
         call starlib_shutdown

         call cali_end_region('do_run_star')
      end subroutine do_run_star

      end module run_star
      
