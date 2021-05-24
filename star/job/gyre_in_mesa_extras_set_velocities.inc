!      example lines from inlist

! GYRE set starting velocities
!   x_logical_ctrl(5) = .true. 
!   x_ctrl(4) = 0d0 ! fraction_1st_overtone (order 2)
!   x_ctrl(5) = 0d0 ! fraction_2nd_overtone (order 3)
!   x_ctrl(6) = 0.1d0 ! initial vsurf (kms)


      subroutine gyre_in_mesa_extras_set_velocities(s,restart,ierr)
         use const_def
         use math_lib
         use gyre_lib
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector_pm
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         
         real(dp), allocatable     :: global_data(:)
         real(dp), allocatable     :: point_data(:,:)
         integer                   :: ipar(5), mode_l
         real(dp)                  :: rpar(1)
         
         integer, parameter :: modes = 3
         integer :: npts(modes), nz, nz_max, i, k, i_v, j, nz_gyre, mode
         real(dp), pointer :: vel(:)
         real(dp), allocatable, dimension(:,:) :: r, v
         real(dp) :: v_surf, v1, period(modes)
         real(dp), allocatable, dimension(:) :: &
            r_old, r_new, v_old, v_new
         real(dp), pointer :: work1(:) ! =(nz_max, pm_work_size)
         
         include 'formats'
         
         if (restart) return
         
         write(*,*) 'set gyre starting velocities'
         
         nz = s% nz
         nz_max = nz+100
         allocate( &
            r(modes,nz_max), v(modes,nz_max), &
            r_old(nz_max), r_new(nz_max), v_old(nz_max), v_new(nz_max))
         npts(:) = 0
         
         call gyre_init('gyre.in')

         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)

         call gyre_set_constant('M_SUN', Msun)
         call gyre_set_constant('R_SUN', Rsun)
         call gyre_set_constant('L_SUN', Lsun)

         call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')
         
         mode_l = 0 ! mode l (e.g. 0 for p modes, 1 for g modes)
                        ! should match gyre.in mode l
         
         !write(*,*) 'call star_get_pulse_data'
         call star_get_pulse_data(s%id, 'GYRE', &
            .FALSE., .FALSE., .FALSE., global_data, point_data, ierr)
         if (ierr /= 0) then
            write(*,*) 'Failed when calling star_get_pulse_data'
            return
         end if
         
         !write(*,*) 'call star_write_pulse_data'
         call star_write_pulse_data(s%id, &
            'GYRE', 'gyre.data', global_data, point_data, ierr)
         if (ierr /= 0) then
            write(*,*) 'Failed when calling star_write_pulse_data'
            return
         end if

         !write(*,*) 'call gyre_set_model'
         call gyre_set_model(global_data, point_data, 101)

          write(*, 100) 'order', 'freq (Hz)', &
             'P (sec)', 'P (day)', 'cycles to double'
100       format(A8,8A20)

         rpar(1) = 0.5d-6 ! freq < this (Hz)
         ipar(1) = s% model_number
         ipar(2) = 1 ! order_target
         ipar(3) = 1 ! 1 means output eigenfunction files
         ipar(4) = 3 ! max number of modes to output per call
         ipar(5) = 0 ! num_written

         call gyre_get_modes(mode_l, process_mode_, ipar, rpar)
         
         if (s% x_ctrl(4) == 1d0) then
            mode = 2
         else if (s% x_ctrl(5) == 1d0) then
            mode = 3
         else
            mode = 1
         end if
         nz_gyre = npts(mode)
         
         if (nz_gyre == 0) then
            ierr = -1
            write(*,*) 'Failed to find requested mode for setting velocities'
            return
         end if
         
         v_surf = v(mode,nz_gyre)
         v1 = 1d5/v_surf ! rescale v's
         if (s% x_ctrl(6) > 0d0) v1 = v1*s% x_ctrl(6)

         do i=1,nz_gyre ! change to star index order
            k = nz_gyre-i+1
            v_old(k) = v1*v(mode,i)
            r_old(k) = r(mode,i)*s% r(1)
         end do
         
         do k=1,nz
            r_new(k) = s% r(k)
         end do
         
         allocate(work1(nz_max*pm_work_size))
         call interpolate_vector_pm( &
            nz_gyre, r_old, nz, r_new, v_old, v_new, work1, 'get_pulse_data', ierr)
         deallocate(work1)
         if (ierr /= 0) then
            write(*,*) 'Failed when calling interpolate_vector_pm for pulse data'
            return
         end if
                  
         if (s% v_flag) then
            do k=1,nz
               s% xh(s% i_v,k) = v_new(k)
            end do
         else if (s% u_flag) then
            do k=1,nz
               s% xh(s% i_u,k) = v_new(k)
            end do
         else
            stop 'set_gyre_linear_analysis failed to set velocities'
         end if
         
         return
         
         do k=1,nz_gyre
            write(*,2) 'gyre r v', k, r_old(k)/Rsun, v_old(k)
            if (k <= nz) write(*,2) 'star r v', k, r_new(k)/Rsun, v_new(k)
            write(*,*)
         end do
         
         do k=nz_gyre+1,nz
            write(*,2) 'star r v', k, r_new(k), v_new(k)
         end do
         
         stop 'gyre_in_mesa_extras_set_velocities'

         contains

         subroutine process_mode_ (md, ipar, rpar, retcode)

            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            character(LEN=128) :: filename
            integer               :: unit, k, model_number, order_target, num_written, max_to_write
            complex(dp)           :: cfreq
            real(dp)              :: freq, growth, freq_lim, per
            logical               :: write_flag
            type(grid_t)          :: gr

            max_to_write = ipar(4)
            num_written = ipar(5)
            if (num_written >= max_to_write) return
            ipar(5) = num_written + 1

            model_number = ipar(1)
            order_target = ipar(2)
            freq_lim = rpar(1)
            write_flag = (ipar(3) == 1)

            cfreq = md% freq('HZ')
            freq = REAL(cfreq)
            growth = AIMAG(cfreq)
            if (freq /= 0d0) then
               per = 1d0/freq
            else
               per = 0d0
            end if

            if (growth > 0._dp) then ! unstable
               write(*, 100) md%n_pg, &
                  freq, per, per/(24*3600), freq/(2d0*pi*growth)
100            format(I8,E20.4,3F20.4)
            else ! stable
               write(*, 110) md%n_pg, &
                  freq, per, per/(24*3600), 'stable'
110            format(I8,E20.4,2F20.4,A20)
            end if
            
            if (md%n_pg > modes) return

            gr = md%grid()

            period(md%n_pg) = 1d0/freq
            npts(md%n_pg) = md%n_k
            do k = 1, md%n_k
               r(md%n_pg,k) = gr%pt(k)%x
               v(md%n_pg,k) = md%xi_r(k)            
            end do

            if (write_flag) then
               ! Write the mode radial & horizontal eigenfunctions, together with the differential work
               write(filename, 120) 'eigfunc.', md%n_pg, '.dat'
120            format(A,I0,A)
               !print *,'Writing eigenfunction to ', TRIM(filename)
               !write(*,*)
               open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')
               write(unit, 130) 'x=r/R', 'Real(xi_r/R)', 'Imag(xi_r/R)', 'Real(xi_h/R)', 'Imag(xi_h/R)', 'dW/dx'
130               format(6(1X,A24))
               do k = 1, md%n_k
                  write(unit, 140) gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k)
140               format(6(1X,E24.16))
               end do
               close(unit)
            end if

            retcode = 0

         end subroutine process_mode_
         
      end subroutine gyre_in_mesa_extras_set_velocities
