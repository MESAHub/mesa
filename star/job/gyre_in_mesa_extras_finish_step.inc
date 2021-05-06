!      example lines from inlist

!      x_integer_ctrl(1) = 10 ! output GYRE info at this step interval
!      x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile

!      x_integer_ctrl(2) = 2 ! max number of modes to output per call
!      x_logical_ctrl(2) = .false. ! output eigenfunction files

!      x_integer_ctrl(3) = 1 ! mode l (e.g. 0 for p modes, 1 for g modes)
         ! should match gyre.in mode l
!      x_integer_ctrl(4) = 1 ! order
!      x_ctrl(1) = 0.6d-4 ! freq < this (Hz)
!      x_ctrl(2) = 0.3d4 ! growth < this (days)

  integer function gyre_in_mesa_extras_finish_step(id)
    integer, intent(in) :: id

    integer                   :: ierr
    logical                   :: call_gyre
    type (star_info), pointer :: s
    real(dp), allocatable     :: global_data(:)
    real(dp), allocatable     :: point_data(:,:)
    integer                   :: ipar(5), step_interval, mode_l
    real(dp)                  :: rpar(2), growth_lim
    
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
             'P (sec)', 'P (min)', 'P (day)', 'growth (day)', 'growth/P', '(4pi*im/re)'
100       format(2A8,7A20)
       else
          write(*, 101) 'model', 'order', 'freq (Hz)', &
             'P (sec)', 'P (min)', 'P (day)', '(4pi*im/re)'
101       format(2A8,5A20)
       end if

       ipar(1) = s% model_number
       ipar(2) = s% x_integer_ctrl(4) ! order_target
       rpar(1) = s% x_ctrl(1)         ! freq_target
       rpar(2) = growth_lim
       ! write_flag
       if (s% x_logical_ctrl(2)) then
          ipar(3) = 1
       else
          ipar(3) = 0
       endif
       ipar(4) = s% x_integer_ctrl(2) ! max_to_write
       ipar(5) = 0 ! num_written

       call gyre_get_modes(mode_l, process_mode_, ipar, rpar)

    end if

    gyre_in_mesa_extras_finish_step = keep_going

  contains

    subroutine process_mode_ (md, ipar, rpar, retcode)

      type(mode_t), intent(in) :: md
      integer, intent(inout)   :: ipar(:)
      real(dp), intent(inout)  :: rpar(:)
      integer, intent(out)     :: retcode

      character(LEN=strlen) :: filename
      integer               :: unit, k, model_number, num_written, max_to_write !, order_target
      complex(dp)           :: cfreq
      real(dp)              :: freq, growth, freq_target, growth_lim
      logical               :: write_flag, freq_okay
      type(grid_t)          :: gr
      
      max_to_write = ipar(4)
      num_written = ipar(5)
      if (num_written >= max_to_write) return
      ipar(5) = num_written + 1
      
      model_number = ipar(1)
      !order_target = ipar(2)
      freq_target = rpar(1)
      growth_lim = rpar(2)
      write_flag = (ipar(3) == 1)

      cfreq = md% freq('HZ')
      freq = REAL(cfreq)
      freq_okay = (abs(freq - freq_target) < freq_target*3d-2)

      if (growth_lim > 0d0) then ! report growth
         if (AIMAG(cfreq) > 0._dp) then ! unstable
            growth = 1d0/(2*pi*24*3600*AIMAG(cfreq))
            write(*, 100) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), growth, &
               growth*freq*24*3600, 4*pi*AIMAG(cfreq)/freq
100         format(2I8,E20.4,3F20.4,3E20.4)
            if (freq_okay .and. growth < growth_lim) & ! .and. md%n_pg == order_target) &
               write(*,*) 'matched target frequency'
         else ! stable
            write(*, 110) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), 'stable'
110         format(2I8,E20.4,3F20.4,A20)
         end if
      else ! growth_lim <= 0 means ignore it
         write(*, 111) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600)
111      format(2I8,E20.4,3F20.4)
         if (freq_okay) & ! .and. md%n_pg == order_target) &
            write(*,*) 'matched target frequency'
      endif

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

         do k = 1, md%n_k
            write(unit, 140) gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k)
140         format(6(1X,E24.16))
         end do

         close(unit)

      end if

      retcode = 0

    end subroutine process_mode_

  end function gyre_in_mesa_extras_finish_step
