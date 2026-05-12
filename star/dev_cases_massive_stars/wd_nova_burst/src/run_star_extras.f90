! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
      use interp_1d_lib

      implicit none

      include "test_suite_extras_def.inc"

      integer :: num_bursts
      logical :: waiting_for_burst
      real(dp) :: L_burst = 1d4, L_between = 1d3 ! Lsun units


      contains

      include "test_suite_extras.inc"

      subroutine nova_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         real(dp) :: b
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         w = 0
         ierr = 0
         b = nova_wind_b(Msurf / Msun)
         if (Lsurf / Lsun > 1d3) then
            w = exp10(-1.49 * safe_log10(Tsurf / 1d5) + b) / (Msun / secyer)
            if (w < 1d-6) then
               w = 0
            else
                write(*,'(A,ES10.3,A)') ' nova wind ', w, ' Msun/year'
            end if
         end if
      end subroutine nova_wind

      function nova_wind_b(m) result(b_out)
         use interp_1d_def, only: pm_work_size
          real(dp), intent(in) :: m
          real(dp) :: b_out
         integer, parameter :: n = 9
         real(dp) :: m_interp(n), b_interp(n)
         real(dp) :: v_new(1), x_new(1)
         real(dp), pointer :: work1(:) => null()
         integer :: ierr

         m_interp = [0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.33d0]
         b_interp = [10.01d0, 20.24d0, 20.38d0, 20.49d0, 20.55d0, 20.62d0, 20.67d0, 20.72d0, 20.78d0]

         x_new(1) = m

         allocate(work1(n * pm_work_size))
         call interpolate_vector_pm( &
            n, m_interp, 1, x_new, b_interp, v_new, work1, 'nova_wind_b', ierr)
         if (ierr /= 0) then
            write(*,*) 'ERROR in nova_wind_b interpolation', ierr
            b_out = 0d0
            return
         end if

         b_out = v_new(1)
         deallocate(work1)
      end function nova_wind_b

      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit,iostat=ierr) num_bursts, waiting_for_burst
      end subroutine extras_photo_read


      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_bursts, waiting_for_burst
      end subroutine extras_photo_write


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write
         s% other_wind => nova_wind

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (.not. restart) then
            num_bursts = 0
            waiting_for_burst = .true.
         end if
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         include 'formats'
         extras_check_model = keep_going
!         if (s% L_phot > L_burst) then
!            if (waiting_for_burst) then
!               num_bursts = num_bursts + 1
!               write(*,2) 'num_bursts', num_bursts
!               waiting_for_burst = .false.
!            end if
!         else if (s% L_phot < L_between) then
!            if (num_bursts >= 1) then
!               write(*,*) 'have finished burst'
!               extras_check_model = terminate
!               s% termination_code = t_extras_check_model
!            end if
!            waiting_for_burst = .true.
!         end if
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'zbar_div_abar'
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

      end function extras_finish_step


      end module run_star_extras
