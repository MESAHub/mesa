      program run
      use run_star_support, only: do_read_star_job
      use run_star, only: do_run_star
      
      implicit none
      
      integer :: ierr
      character (len=32) :: inlist_fname
      
      ierr = 0
      inlist_fname = 'inlist'
      
      call do_read_star_job(inlist_fname, ierr)
      if (ierr /= 0) stop 1
      
      call do_run_star(inlist_fname)
      
      end program
