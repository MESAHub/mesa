! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
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

      implicit none

      include "test_suite_extras_def.inc"
      real(dp) :: t_spindown


      contains

      include "test_suite_extras.inc"


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_torque => magnetic_braking
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
      end subroutine extras_controls



      subroutine magnetic_braking(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         real(dp) :: bfield, vinf, eta, factor
         real(dp) :: j_dot, check_delta_j, i_tot, j_tot, delta_j
         real(dp) :: residual_jdot, torque, j_average


         type (star_info), pointer :: s
         integer :: k, j
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !Parameters

         bfield = s% x_ctrl(1) ! Surface Magnetic Fields in Gauss

         !Initialize variables

         s% extra_jdot(:) = 0d0
         s% extra_omegadot(:) = 0d0
         i_tot = 0d0
         j_dot = 0d0
         t_spindown = 0d0
         residual_jdot = 0d0
         check_delta_j = 0d0
         torque = 0d0
         factor = 0d0
         j_average = 0d0

         ! Calculate total specific moment of inertia and angular momentum

         i_tot = sum(s% i_rot(1:s% nz))
         j_tot = dot_product(s% j_rot(1:s% nz),s% dm_bar(1:s% nz)) ! g cm^2/s Total Stellar Angular Momentum Content

         if ((s% mstar_dot /= 0) .and. (j_tot .gt. 1d49)) then ! Only 'brake' when mass is lost and star has non-negligible amount of angular momentum
          !Calculate V_inf of stellar wind (e.g. Vinf = 1.92 Vesc, see Lamers & Cassinelli 2000)
          !N.B. This is good for line-driven winds in hot stars. For different types of Vinf = Vesc might be a better choice?
          vinf = 1.92d0 * sqrt(2.0d0 * standard_cgrav * s% mstar / (s% photosphere_r * Rsun))
          !Calculate Wind Confinement parameter 'Eta' (See ud-doula & Owocki 2002)
          eta = pow2( s% photosphere_r * Rsun * bfield ) / (abs( s% mstar_dot ) * vinf)
          !Calculate Jdot !Eq. 7 in  Ud-Doula et al. 2008
          j_dot = (2d0/3d0) * s% mstar_dot * s% omega_avg_surf*pow2( s% photosphere_r *Rsun)
          j_dot = j_dot * eta ! This is g*cm^2/s^2
          ! This could be improved using
          ! Eq. 20 in Ud-Doula et al. (2008) MNRAS (Dipole Weber-Davis)
          ! delta_j = delta_j * (0.29 *(eta+0.25)**0.25)**2.0
          ! But note that this formula doesnt work in the B -> 0 limit

          ! Note that other_torque uses extra_jdot, which is the rate at which specific angular momentum is changed
          ! i.e. units of extra_jdot are cm^2/s^2

          ! Fraction of angular momentum lost during timestep
          delta_j = j_dot * s% dt / j_tot

          ! Check if spindown timescale is shorter than timestep. Print a warning in case.
          ! In extras_finish_step check that dt < t_spindown. If not, decrease timestep
          t_spindown = abs(j_tot / j_dot) ! Estimate spindown timescale
          if (s% x_logical_ctrl(1)) then
             write(*,*) 'Spindown Timescale (Myr)', t_spindown / (1d6*secyer)
             write(*,*) 'Spindown Timescale / dt: ', t_spindown / s% dt
          end if



          ! Let's assume the magnetic field applies a ~uniform torque through the star
          ! Strategy:
          ! 1) Start from the surface. Apply weighted torque (jdot) to each gridpoint
          ! 2) If j_dot * dt < j_rot(k), all good. If j_dot * dt > j_rot(k) then set j_dot = -j_rot(k)/dt. This is to avoid flip/flops in omega
          ! 3) To conserve angular momentum one has to redistribute the residual torque
          j = 0
          do k = 1, s% nz
            s% extra_jdot(k) =  j_dot / (s% nz * s% dm_bar(k))  ! Specific torque cm^2/s^2. Divide by shell mass and total number of shells
            if (abs(s% extra_jdot(k)) .gt. abs(s% j_rot(k)/ s% dt)) then
              residual_jdot = residual_jdot - (abs( s% extra_jdot(k)) - abs( s% j_rot(k) / s% dt )) * s% dm_bar(k) ! Residual J_dot cm^2/s^2 * g
              s% extra_jdot(k) = - s% j_rot(k)/ s% dt ! Set torque = - s% j_rot(k)/ s% dt. Note this way we're not conserving angular momentum, need to distribute residual torque
              j=j+1
            endif
          end do

          ! Redistribute residual J_dot (only to gridpoints that can accomodate more torque)
          ! j number of cells that can not take anymore torque ()
          do k = 1, s% nz
            if (abs(s% extra_jdot(k)) .lt. abs(s% j_rot(k)/ s% dt)) then
              s% extra_jdot(k) = s% extra_jdot(k) + (residual_jdot / ((s% nz - j) * s% dm_bar(k)))
            endif
          end do

          torque = dot_product(s% extra_jdot(1:s% nz),s% dm_bar(1:s% nz)) ! Total applied Torque

          ! Wind Diagnostics
          !write(*,*) 'Rotational Velocity: ',s% omega_avg_surf * s% photosphere_r *Rsun / 1d5, 'km/s' !Rot Vel.
          !write(*,*) 'Wind terminal Velocity: ',vinf/1d5, 'km/s'
          !write(*,*) 'Wind confinement parameter: ', eta
          ! Torque Diagnostics
          !write(*,*) 'Torque we want to apply: ', j_dot
          !write(*,*) 'Angular Momentum we want to remove', j_dot* s% dt
          !write(*,*) 'Stellar Angular Momentum: ', j_tot
          !write(*,*) 'Specific Torque not allocated: ', residual_jdot

          ! Angular Momentum Conservation check
          if (s% x_logical_ctrl(1)) then
             write(*,*) 'Fraction of total angular momentum to remove', (j_dot * s% dt) / j_tot
             write(*,*) 'Torque / J_dot (if = 1.0 angular momentum is conserved): ', ((torque) / j_dot)
          end if

        else
          t_spindown = 100 * s% dt ! To avoid decreasing the timestep in extras_finish_step when j_tot < 1d49
        endif


      end subroutine magnetic_braking




      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% v_rot_avg_surf < 1d5) then
            write(*,*) 'Test is ok: surf_avg_v_rot < 1'
         else
            write(*,*) 'Test failed: surf_avg_v_rot > 1'
         end if
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
         extras_check_model = keep_going
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
         !write(*,*) t_spindown/s% dt
         if (s% dt > 0) then
            if ((t_spindown / s% dt) .lt. 10) then
               s% dt_next = s% dt * 0.5d0
               write(*,*) "Warning: Torque too large. Decreasing timestep to ",s% dt_next
            end if
         end if


      end function extras_finish_step



      end module run_star_extras
