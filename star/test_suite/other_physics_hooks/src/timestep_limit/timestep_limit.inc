
      
      subroutine read_inlist_other_timestep_limit(ierr)
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit         
         filename = 'inlist_other_timestep_limit'
         
         write(*,*) 'read_inlist_other_timestep_limit'
         
         ! set defaults
         delta_lgT_limit = 0.5d0
         delta_lgT_hard_limit = -1
         delta_lgT_limit_min_lgT = 1d99
         rotation_steps_limit = 500
         rotation_steps_hard_limit = 700
         delta_mdot_atol = 1d-3
         delta_mdot_rtol = 0.5d0
         delta_mdot_limit = -1
         delta_mdot_hard_limit = -1

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=other_timestep_limit, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=other_timestep_limit)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_other_timestep_limit


      integer function set_other_timestep_limit( &
            id, skip_hard_limit, dt, dt_limit_ratio)
         ! check max change in log10(temperature)
         integer, intent(in) :: id
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         type (star_info), pointer :: s
         real(dp) :: relative_excess, max_dlnT
         integer :: i, ierr
         include 'formats'

         set_other_timestep_limit = keep_going

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            set_other_timestep_limit = terminate
            return
         end if
         if (s% model_number < 29 .or. s% model_number > 42) return

         !write(*,*) 'check_rotation_steps_limit'
         set_other_timestep_limit = check_rotation_steps_limit( &
            s, skip_hard_limit, s% dt, dt_limit_ratio)
         if (return_now('check_rotation_steps_limit')) return

         !write(*,*) 'check_delta_mdot'
         set_other_timestep_limit = check_delta_mdot( &
            s, skip_hard_limit, s% dt, dt_limit_ratio)
         if (return_now('check_delta_mdot')) return
         
         !write(*,*) 'check_dlgT_change'
         set_other_timestep_limit = check_dlgT_change( &
            s, skip_hard_limit, dt_limit_ratio)
         if (return_now('check_dlgT_change')) return

         contains
         
         logical function return_now(str)
            character (len=*), intent(in) :: str            
            if (set_other_timestep_limit == keep_going) then
               return_now = .false.
               return
            end if            
            write(*, '(a, i8)') &
               'set_other_timestep retry because of timestep limit ' // &
               trim(str), s% model_number            
            return_now = .true.            
         end function return_now
         
      end function set_other_timestep_limit


      integer function check_dlgT_change(s, skip_hard_limit, dt_limit_ratio)
         ! check max change in log10(temperature)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnT
         integer :: i
         include 'formats'
         check_dlgT_change = keep_going
         call get_dlgT_info(s, i, max_dlnT)
         if (i == 0) return
         check_dlgT_change = check_change_timestep_limit(s% id, max_dlnT/ln10, &
            delta_lgT_limit, delta_lgT_hard_limit, &
            i, 'check_dlgT_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         
         contains
         
         subroutine get_dlgT_info(s, i, max_dlnT)
            use const_def, only:ln10
            type (star_info), pointer :: s
            integer, intent(out) :: i
            real(dp), intent(out) :: max_dlnT
            real(dp) :: lim, dlnT
            integer :: k
            include 'formats'
            lim = ln10*delta_lgT_limit_min_lgT
            i = 0
            max_dlnT = 0
            do k=1,s% nz
               if (s% lnT(k) < lim) cycle
               dlnT = abs(s% lnT(k) - s% lnT_start(k))
               if (dlnT > max_dlnT) then
                  max_dlnT = dlnT
                  i = k
               end if
            end do
         end subroutine get_dlgT_info

      end function check_dlgT_change


      integer function check_rotation_steps_limit(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         include 'formats'
         check_rotation_steps_limit = keep_going
         if (s% num_retries == 0) then  ! JUST CHECKING FORCED RETRY
            check_rotation_steps_limit = retry
            return
         end if
         if (.not. s% rotation_flag) return
         check_rotation_steps_limit = check_change_integer_timestep_limit( &
           s% id, rotation_steps_limit, rotation_steps_hard_limit, &
           s% num_rotation_solver_steps,  &
           'num_rotation_solver_steps', skip_hard_limit, dt, dt_limit_ratio)
      end function check_rotation_steps_limit


      integer function check_delta_mdot(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: mdot, mdot_old, delta_mdot
         check_delta_mdot = keep_going
         mdot = s% mstar_dot
         mdot_old = s% mstar_dot_old
         delta_mdot = abs(mdot - mdot_old)/ &
            (delta_mdot_atol*Msun/secyer + &
               delta_mdot_rtol*max(abs(mdot),abs(mdot_old)))
         if (delta_mdot == 0) return
         if (delta_mdot_hard_limit > 0 .and. (.not. skip_hard_limit)) then
            if (delta_mdot > delta_mdot_hard_limit) then
               if (s% report_all_dt_limits) &
                  write(*, '(a30, f20.10, 99e20.10)') 'delta_mdot_hard_limit', &
                     delta_mdot, delta_mdot_hard_limit
                  check_delta_mdot= retry
               return
            end if
         end if
         if (delta_mdot_limit <= 0) return
         dt_limit_ratio = delta_mdot/delta_mdot_limit
         if (dt_limit_ratio <= 1d0) then
            dt_limit_ratio = 0
         else if (s% report_all_dt_limits) then
            write(*, '(a30, f20.10, 99e20.10)') 'delta_mdot_limit', &
               dt_limit_ratio, delta_mdot, delta_mdot_limit
         end if
      end function check_delta_mdot
      
      
