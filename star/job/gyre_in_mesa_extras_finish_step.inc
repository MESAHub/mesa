!      example lines from inlist

!      x_integer_ctrl(1) = 10 ! output GYRE info at this step interval
!      x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile

!      x_integer_ctrl(2) = 2 ! max number of modes to output per call
!      x_logical_ctrl(2) = .false. ! output eigenfunction files

!      x_integer_ctrl(3) = 1 ! mode l (e.g. 0 for p modes, 1 for g modes)
         ! should match gyre.in mode l
!      x_ctrl(1) = 0.6d-4 ! freq < this (Hz)
!      x_ctrl(2) = 0.3d4 ! growth < this (days)

  integer function gyre_in_mesa_extras_finish_step(id)
    integer, intent(in) :: id

    integer                   :: ierr
    logical                   :: call_gyre
    type (star_info), pointer :: s
    real(dp), allocatable     :: global_data(:)
    real(dp), allocatable     :: point_data(:,:)
    integer                   :: ipar(2), step_interval, mode_l
    real(dp)                  :: rpar(1), growth_lim
    
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    gyre_in_mesa_extras_finish_step = terminate

    ! Extract parameters

    step_interval = s% x_integer_ctrl(1)
    mode_l = s% x_integer_ctrl(3)
    growth_lim = s% x_ctrl(2)

    if (s% need_to_save_profiles_now .and. s% x_logical_ctrl(1)) then
       call_gyre = .TRUE.
    elseif (step_interval > 0) then
       call_gyre = MOD(s% model_number, step_interval) == 0
    else
       call_gyre = .FALSE.
    endif
    
    s% ixtra3_array(1) = 0 ! number of results

    ! If necessary, call GYRE

    if (call_gyre) then

       call star_get_pulse_data(s%id, 'GYRE', .FALSE., .FALSE., .FALSE., global_data, point_data, ierr)
       if (ierr /= 0) then
          print *,'Failed when calling star_get_pulse_data'
          return
       end if

       call gyre_set_model(global_data, point_data, 101)
       
       if (growth_lim > 0d0) then
          write(*, 100) 'model', 'order', 'freq (Hz)', &
             'P (sec)', 'P (day)', 'growth (day)', 'growth', 'cycles to double'
100       format(2A8,99A20)
       else
          write(*, 101) 'model', 'order', 'freq (Hz)', &
             'P (sec)', 'P (day)'
101       format(2A8,5A20)
       end if

       ipar(1) = s% id
       ipar(2) = 0 ! num_written

       call gyre_get_modes(mode_l, process_mode_, ipar, rpar)
       
       s% ixtra3_array(1) = ipar(2)

    end if

    gyre_in_mesa_extras_finish_step = keep_going

  contains

    subroutine process_mode_ (md, ipar, rpar, retcode)

      type(mode_t), intent(in) :: md
      integer, intent(inout)   :: ipar(:)
      real(dp), intent(inout)  :: rpar(:)
      integer, intent(out)     :: retcode

      character(LEN=strlen) :: filename
      integer               :: ierr, unit, k, model_number, num_written, max_to_write !, order_target
      complex(dp)           :: cfreq
      real(dp)              :: freq, growth, period, freq_target, growth_lim
      logical               :: write_flag, freq_okay
      type(grid_t)          :: gr
      type (star_info), pointer :: s
      

      ierr = 0
      call star_ptr(ipar(1), s, ierr)
      if (ierr /= 0) return

      num_written = ipar(2)
      max_to_write = s% x_integer_ctrl(2)
      if (num_written >= max_to_write) return
      num_written = num_written + 1
      ipar(2) = num_written

      model_number = s% model_number
      freq_target = s% x_ctrl(1)
      growth_lim = s% x_ctrl(2)
      write_flag = s% x_logical_ctrl(2)

      cfreq = md% freq('HZ')
      growth = AIMAG(cfreq) ! in seconds
      freq = REAL(cfreq) ! in seconds
      freq_okay = (abs(freq - freq_target) < freq_target*3d-2)
      period = 0 ! days
      if (growth_lim > 0d0) then ! report growth
         if (growth > 0d0) then ! unstable
            period = 1d0/freq ! in seconds
            write(*, 100) model_number, md%n_pg, &
               freq, period, period/(24*3600), 1d0/(2*pi*24*3600*AIMAG(cfreq)), &
               (2d0*pi*growth)/freq, freq/(2d0*pi*growth)
100         format(2I8,E20.4,5F20.4)
            if (freq_okay .and. growth < growth_lim) & ! .and. md%n_pg == order_target) &
               write(*,*) 'matched target frequency'
         else ! stable
            write(*, 110) model_number, md%n_pg, &
               freq, period, period/(24*3600), 'stable'
110         format(2I8,E20.4,2F20.4,A20)
         end if
      else ! growth_lim <= 0 means ignore it
         write(*, 111) model_number, md%n_pg, freq, period, period/(24*3600)
111      format(2I8,E20.4,2F20.4)
         if (freq_okay) & ! .and. md%n_pg == order_target) &
            write(*,*) 'matched target frequency'
      endif

      s% ixtra1_array(num_written) = md%n_pg
      s% xtra1_array(num_written) = growth
      s% xtra2_array(num_written) = freq
      if (growth > 0d0) then
         s% xtra3_array(num_written) = freq/(2d0*pi*growth)
      else
         s% xtra3_array(num_written) = -1d0
      end if

      if (write_flag) then

      ! Write the mode radial & horizontal eigenfunctions, together with the differential work
       
         write(filename, 120) 'eigfunc.', model_number, '.', md%n_pg, '.dat'
120      format(A,I0,A,I0,A)

         print *,'Writing eigenfunction to file:', TRIM(filename)
         write(*,*)

         open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')

         write(unit, 130) 'x=r/R', 'Real(xi_r/R)', 'Imag(xi_r/R)', 'Real(xi_h/R)', 'Imag(xi_h/R)', 'dW/dx'
130      format(6(1X,A24))

         gr = md%grid()

         do k = 1, md%n
            write(unit, 140) gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k)
140         format(6(1X,E24.16))
         end do

         close(unit)

      end if

      retcode = 0

    end subroutine process_mode_

  end function gyre_in_mesa_extras_finish_step
