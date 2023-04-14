! ***********************************************************************
!
!   Copyright (C) 2018-2019  Rich Townsend, The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      use gyre_lib

      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      
      include 'run_star_extras.inc'
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         if (.not. s% use_other_RSP_linear_analysis) return

         s% other_rsp_linear_analysis => rsp_set_gyre_linear_analysis

      end subroutine extras_controls
      

      subroutine rsp_set_gyre_linear_analysis(id,restart,ierr)
         use const_def
         use math_lib
         use gyre_lib
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         real(dp), allocatable     :: global_data(:)
         real(dp), allocatable     :: point_data(:,:)
         integer                   :: ipar(5), mode_l
         real(dp)                  :: rpar(1)
         
         integer, parameter :: modes = 3
         integer :: npts(modes), nz, i, k
         real(dp), allocatable, dimension(:,:) :: r, v
         real(dp) :: velkm, v_surf, amix1, amix2, amixF, &
            period(modes)
         
         include 'formats'
         
         if (restart) return
         
         write(*,*) 'set gyre starting velocities'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         allocate(r(modes,nz+10), v(modes,nz+10))
         npts = 0

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
         
         call star_get_pulse_data(s%id, 'GYRE', &
            .FALSE., .FALSE., .FALSE., global_data, point_data, ierr)
         if (ierr /= 0) then
            print *,'Failed when calling get_pulse_data'
            return
         end if
         
         call star_write_pulse_data(s%id, &
            'GYRE', 'gyre.data', global_data, point_data, ierr)
         if (ierr /= 0) return

         call gyre_set_model(global_data, point_data, 101)

         write(*, 100) 'order', 'freq (Hz)', 'P (sec)', &
           'P (min)', 'P (day)', 'growth (day)', '(4pi*im/re)'
100      format(A8,A16,A16,A14,A12,A16,A14)

         rpar(1) = 0.2d-3 ! freq < this (Hz)
         ipar(1) = s% model_number
         ipar(2) = 1 ! order_target
         ipar(3) = 1 ! 1 means output eigenfunction files
         ipar(4) = 3 ! max number of modes to output per call
         ipar(5) = 0 ! num_written

         call gyre_get_modes(mode_l, process_mode_, ipar, rpar)

         call gyre_final()
         
         amix1 = s% x_ctrl(4) ! s% RSP_fraction_1st_overtone
         amix2 = s% x_ctrl(5) ! s% RSP_fraction_2nd_overtone
         if((amix1+amix2) > 1d0) then
            write(*,*) 'AMIX DO NOT ADD UP RIGHT' 
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
         end if
         velkm = s% x_ctrl(6) ! s% RSP_kick_vsurf_km_per_sec
         amixF = 1d0 - (amix1 + amix2)
         
         if (amixF > 0d0 .and. npts(1) /= nz-1) then
            write(*,3) 'amixF > 0d0 .and. npts(1) /= nz-1', npts(1)
            write(*,*) 'cannot use fundamental for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         if (AMIX1 > 0d0 .and. npts(2) /= nz-1) then
            write(*,3) 'AMIX1 > 0d0 .and. npts(2) /= nz-1', npts(2)
            write(*,*) 'cannot use 1st overtone for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         if (AMIX2 > 0d0 .and. npts(2) /= nz-1) then
            write(*,3) 'AMIX2 > 0d0 .and. npts(3) /= nz-1', npts(3)
            write(*,*) 'cannot use 2nd overtone for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         v_surf = amixF*v(1,nz-1) + AMIX1*v(2,nz-1) + AMIX2*v(3,nz-1)

         do i=1,nz-1
            k = nz+1-i ! v(1) from gyre => s% v(nz) in star
            s% v(k)=1.0d5*VELKM/v_surf* &
                (amixF*v(1,i) + AMIX1*v(2,i) + AMIX2*v(3,i))
         end do
         s% v(1) = s% v(2)
         s% v_center = 0d0
         
         do k=1,nz
            s% xh(s% i_v,k) = s% v(k)
         end do
         
         s% rsp_period = period(s% RSP_mode_for_setting_PERIODLIN + 1)
         
         !write(*,*) 'amix1 amix2 amixF velkm v_surf period', amix1, amix2, amixF, velkm, v_surf, s% rsp_period

         contains

         subroutine process_mode_ (md, ipar, rpar, retcode)

            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            character(LEN=128) :: filename
            integer               :: unit, k, model_number, order_target, num_written, max_to_write
            complex(dp)           :: cfreq
            real(dp)              :: freq, growth, freq_lim
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

            if (AIMAG(cfreq) > 0._dp) then ! unstable
               growth = 1d0/(2*pi*24*3600*AIMAG(cfreq))
               write(*, 100)  md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), &
                  growth, 4*pi*AIMAG(cfreq)/freq
100               format(I8,E16.4,F16.4,F14.4,F12.4,E16.4,E14.4)
            else ! stable
               write(*, 110) md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), 'stable'
110         format(I8,E16.4,F16.4,F14.4,F12.4,A16)
            end if
            
            if (md%n_pg > modes) return

            gr = md%grid()

            period(md%n_pg) = 1d0/freq
            npts(md%n_pg) = md%n
            do k = 1, md%n
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
               do k = 1, md%n
                  write(unit, 140) gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k)
140               format(6(1X,E24.16))
               end do
               close(unit)
            end if

            retcode = 0

         end subroutine process_mode_
         
      end subroutine rsp_set_gyre_linear_analysis
      

      end module run_star_extras
      
