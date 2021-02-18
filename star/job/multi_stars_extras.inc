
      
      subroutine do_run_multi_stars
         integer :: id, ierr, i, i_prev, result, result_reason, model_number
         type (star_info), pointer :: s
         character (len=64) :: inlist_fname
         logical :: first_try, continue_evolve_loop
         real(dp) :: sum_times
         real(dp) :: dt
         logical :: restart, pgstar_ok
         logical, parameter :: &
            do_alloc_star = .true., &
            do_free_star = .true., &
            okay_to_restart = .true., &
            dbg = .false.
         
         include 'formats'

         ierr = 0
         call test_suite_startup(s, .false., ierr)

         inlist_fname = 'inlist_multi_stars_job'    
         call read_controls(inlist_fname,ierr)
         if (ierr /= 0) return
         
         if (num_stars < 1) then
            write(*,*) 'need to set num_stars >= 1'
            return
         end if
         
         if (num_stars > max_num_stars) then
            write(*,*) 'need to set num_stars <= max_num_stars or rebuild with larger max'
            return
         end if
         
         do i = 1, num_stars
            
            id_from_read_star_job = 0
            call do_read_star_job_and_return_id(inlist_names(i), id, ierr)
            if (failed('do_read_star_job',ierr)) return
            star_ids(i) = id
            star_done(i) = .false.
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            
            pgstar_ok = (i == which_for_pgstar .or. which_for_pgstar < 0)

            call start_run1_star( &
               do_alloc_star, do_free_star, okay_to_restart, &
               id, restart, pgstar_ok, dbg, &
               extras_controls, &
               ierr, inlist_names(i))
            if (failed('before_evolve_loop',ierr)) return

            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            
            s% doing_timing = .false.

         end do

         continue_evolve_loop = .true.
         i_prev = 0

         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            
            i = select_youngest_star()
            if (i == 0) exit ! all done
            
            call star_ptr(star_ids(i), s, ierr)
            if (failed('star_ptr',ierr)) return
            
            if (s% star_age >= stopping_age) then
               write(*,*) 'stars have reached stopping age'
               exit
            end if
         
            continue_evolve_loop = do_evolve_one_step(s, dbg, ierr)
            if (failed('do_evolve_one_step',ierr)) return
            
            if (.not. continue_evolve_loop) then
               if (s% result_reason == result_reason_normal) then
                  star_done(i) = .true.
                  if (.not. stop_when_1st_finishes) continue_evolve_loop = .true.
               end if
            end if
            
         end do evolve_loop

         do i = 1, num_stars

            call after_evolve_loop(star_ids(i), do_free_star, ierr)
            if (failed('after_evolve_loop',ierr)) return
            
         end do
         
         call starlib_shutdown
         
         call test_suite_after_evolve(s, ierr)


         contains
         

         integer function select_youngest_star()
            integer :: i
            real(dp) :: age_min
            type (star_info), pointer :: s
            select_youngest_star = 0
            age_min = 1d99
            do i = 1, num_stars
               if (star_done(i)) cycle
               call star_ptr(star_ids(i), s, ierr)
               if (failed('star_ptr',ierr)) return
               if (s% star_age < age_min) then
                  age_min = s% star_age
                  select_youngest_star = i
               end if
            end do
         end function select_youngest_star

         subroutine read_controls(filename,ierr)
            use utils_lib
            character (len=*) :: filename
            integer, intent(out) :: ierr         
            character (len=256) :: message
            integer :: unit
            ! set defaults
            num_stars = 0
            restart_flag = .false.
            restart_names(:) = 'undefined'
            inlist_names(:) = 'undefined'
            which_for_pgstar = -1
            stopping_age = 1d99
            stop_when_1st_finishes = .true.
            open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'Failed to open control namelist file ', trim(filename)
            else
               read(unit, nml=multi_stars_job, iostat=ierr)  
               close(unit)
               if (ierr /= 0) then
                  write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
                  write(*, '(a)') &
                     'The following runtime error message might help you find the problem'
                  write(*, *) 
                  open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
                  read(unit, nml=multi_stars_job)
                  close(unit)
               end if  
            end if
         end subroutine read_controls

      end subroutine do_run_multi_stars   
