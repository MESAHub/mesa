      program run
      use run_star_support, only: do_read_star_job_and_return_id
      use run_star, only: do_run_star
      
      implicit none
      
      integer :: id, ierr
      
      ierr = 0
      call do_read_star_job_and_return_id('inlist', id, ierr)
      if (ierr /= 0) stop 1
      
      call do_run_star
      
      end program
