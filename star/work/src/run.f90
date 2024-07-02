      program run
      use caliper_mod
      use run_star_support, only: do_read_star_job
      use run_star, only: do_run_star
      
      implicit none
      
      integer :: ierr
      character (len=32) :: inlist_fname

    
      type(ConfigManager)   :: mgr
      !integer               :: argc
      !logical               :: ret
      !character(len=:), allocatable :: errmsg
      !character(len=256)    :: arg
    
      ! (Optional) create a ConfigManager object to control profiling.
      ! Users can provide a configuration string (e.g., 'runtime-report')
      ! on the command line.
      print *, 'XXXXXXX'
      mgr = ConfigManager_new()
      call mgr%set_default_parameter('aggregate_across_ranks', 'false')
      call mgr%add('runtime-report')
      !argc = command_argument_count()
      !if (argc .ge. 1) then
      !    call get_command_argument(1, arg)
      !    call mgr%add(arg)
      !    ret = mgr%error()
      !    if (ret) then
      !        errmsg = mgr%error_msg()
      !        write(*,*) 'ConfigManager: ', errmsg
      !    endif
      !endif
    
      ! Start configured profiling channels
      call mgr%start
    
      ! A scope annotation. Start region 'main'
      call cali_begin_region('run')

      
      ierr = 0
      inlist_fname = 'inlist'
      
      call do_read_star_job(inlist_fname, ierr)
      if (ierr /= 0) stop 1
      
      call do_run_star(inlist_fname)

      call cali_end_region('run')
      ! Compute and flush output for the ConfigManager profiles.
      call mgr%flush
      call ConfigManager_delete(mgr)     
      end program
